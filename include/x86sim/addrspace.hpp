#ifndef X86SIM_ADDRSPACE_HPP
#define X86SIM_ADDRSPACE_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <format>
#include <memory>
#include <span>
#include <unordered_map>

namespace x86sim {

using address_t = std::uint64_t;

enum class Protection : std::uint8_t { none = 0, read = 1, write = 2, execute = 4 };
enum class MemoryError { unaligned_address, zero_size, unmapped_address, out_of_memory, mapping_failed };

[[nodiscard]] constexpr Protection operator|(Protection lhs, Protection rhs) noexcept {
  return static_cast<Protection>(static_cast<std::uint8_t>(lhs) | static_cast<std::uint8_t>(rhs));
}

[[nodiscard]] constexpr Protection operator&(Protection lhs, Protection rhs) noexcept {
  return static_cast<Protection>(static_cast<std::uint8_t>(lhs) & static_cast<std::uint8_t>(rhs));
}

[[nodiscard]] constexpr bool has_protection(Protection value, Protection flag) noexcept {
  return (static_cast<std::uint8_t>(value) & static_cast<std::uint8_t>(flag)) == static_cast<std::uint8_t>(flag);
}

// A guest virtual address space: a sparse set of mapped 4 KiB pages with
// per-page protection bits. The library user constructs and mutates these
// directly (map/unmap/read/write) and hands one to Machine::run alongside a
// CpuState; multiple address spaces can be kept and swapped independently.
class AddressSpace {
public:
  static constexpr std::uint64_t kPageSize = 4096;
  using Page = std::array<std::byte, kPageSize>;

  AddressSpace();
  AddressSpace(const AddressSpace&);
  AddressSpace& operator=(const AddressSpace&);
  AddressSpace(AddressSpace&&) noexcept;
  AddressSpace& operator=(AddressSpace&&) noexcept;
  ~AddressSpace();

  [[nodiscard]] std::expected<void, MemoryError> map(address_t start, std::uint64_t size, Protection) noexcept;
  void unmap(address_t start, std::uint64_t size) noexcept;
  [[nodiscard]] std::expected<void, MemoryError> read(address_t start, std::span<std::byte> out) const noexcept;
  [[nodiscard]] std::expected<void, MemoryError> write(address_t start, std::span<const std::byte>) noexcept;
  [[nodiscard]] std::unique_ptr<AddressSpace> clone_deep() const;

  [[nodiscard]] void* page_virt_to_mapped(address_t addr) noexcept;
  [[nodiscard]] const void* page_virt_to_mapped(address_t addr) const noexcept;

  // Accessibility test for the simulator core: returns true only if every
  // protection bit requested in `prot` is granted for the page at `address`.
  [[nodiscard]] bool check(address_t address, Protection prot) const noexcept;

  // Self-modifying-code dirty tracking, keyed by machine frame number.
  [[nodiscard]] bool isdirty(address_t mfn) const noexcept;
  void setdirty(address_t mfn);
  void cleardirty(address_t mfn);

private:
  // Shadow page-attribute tables (read/write/exec/dirty). These are an
  // implementation detail of the accessibility checks above and never escape
  // the class.
  struct ShadowMap;
  using spat_t = ShadowMap*;
  [[nodiscard]] bool fastcheck(address_t addr, spat_t top) const noexcept;
  [[nodiscard]] spat_t read_map() const noexcept { return readmap_.get(); }
  [[nodiscard]] spat_t write_map() const noexcept { return writemap_.get(); }
  [[nodiscard]] spat_t exec_map() const noexcept { return execmap_.get(); }
  [[nodiscard]] Protection protection_at(address_t address) const noexcept;

  using PageMap = std::unordered_map<address_t, std::unique_ptr<Page>>;

  void map_or_throw(address_t start, std::uint64_t size, Protection);
  void setattr(address_t start, std::uint64_t size, Protection);
  void make_accessible(address_t address, std::uint64_t size, ShadowMap&);
  void make_inaccessible(address_t address, std::uint64_t size, ShadowMap&);
  void make_page_accessible(address_t address, ShadowMap&);
  void make_page_inaccessible(address_t address, ShadowMap&);
  [[nodiscard]] static address_t pageid(address_t address) noexcept;

  PageMap mapped_mem_;
  std::unique_ptr<ShadowMap> readmap_;
  std::unique_ptr<ShadowMap> writemap_;
  std::unique_ptr<ShadowMap> execmap_;
  std::unique_ptr<ShadowMap> dirtymap_;
};

} // namespace x86sim

template<>
struct std::formatter<x86sim::Protection> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(x86sim::Protection protection, FormatContext& ctx) const {
    const auto bits = static_cast<unsigned>(protection);
    if ((bits & ~0x7u) != 0)
      return std::format_to(ctx.out(), "unknown(0x{:02x})", bits);
    if (protection == x86sim::Protection::none)
      return std::format_to(ctx.out(), "none");

    auto out = ctx.out();
    if (x86sim::has_protection(protection, x86sim::Protection::read))
      out = std::format_to(out, "r");
    if (x86sim::has_protection(protection, x86sim::Protection::write))
      out = std::format_to(out, "w");
    if (x86sim::has_protection(protection, x86sim::Protection::execute))
      out = std::format_to(out, "x");
    return out;
  }
};

template<>
struct std::formatter<x86sim::MemoryError> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(x86sim::MemoryError error, FormatContext& ctx) const {
    using enum x86sim::MemoryError;
    switch (error) {
    case unaligned_address:
      return std::format_to(ctx.out(), "address is not page-aligned");
    case zero_size:
      return std::format_to(ctx.out(), "range has zero size");
    case unmapped_address:
      return std::format_to(ctx.out(), "address is not mapped");
    case out_of_memory:
      return std::format_to(ctx.out(), "out of memory");
    case mapping_failed:
      return std::format_to(ctx.out(), "mapping failed");
    }
    return std::format_to(ctx.out(), "unknown({})", static_cast<int>(error));
  }
};

#endif
