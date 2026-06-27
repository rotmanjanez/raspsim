// Portable half of the Linux syscall support library.
//
// This translation unit contains ONLY portable logic: the shared decode/error
// helpers, the brk/mmap heap emulation, and the register-only syscalls. It
// touches solely CpuState / AddressSpace / synthetic data and therefore must not
// include any host header (<unistd.h>, <fcntl.h>, <sys/socket.h>,
// <netinet/in.h>, <sys/stat.h>, <sys/time.h>, <dirent.h>, <time.h>, <utime.h>,
// <sys/un.h>). The host-POSIX handlers live in syscall-linux-posix.cpp.
#include "syscall-linux-detail.hpp"
#include "x86sim-support/syscall-linux.hpp"
#include "x86sim/addrspace.hpp"
#include "x86sim/registerfile.hpp"
#include "x86sim/x86sim.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

namespace x86sim::linux_syscalls {

std::shared_ptr<ProcessTable> make_process_table() {
  return std::make_shared<ProcessTable>();
}

namespace detail {

[[nodiscard]] word_t syscall_number(const RegisterFile& context) noexcept {
  return context[Register::rax];
}

[[nodiscard]] word_t syscall_arg(const CpuState& context, std::size_t index) noexcept {
  static constexpr Register registers[] = {Register::rdi, Register::rsi, Register::rdx,
                                           Register::r10, Register::r8,  Register::r9};
  return context[registers[index]];
}

[[nodiscard]] bool handles(const CpuState& context, SyscallKind kind, word_t number) noexcept {
  return kind == SyscallKind::syscall64 && syscall_number(context) == number;
}

namespace {

[[nodiscard]] SyscallResult continue_result() {
  return {.reason = StopReason::host_request, .continue_execution = true, .message = {}};
}

} // namespace

[[nodiscard]] SyscallResult return_value(CpuState& context, std::int64_t value) {
  context[Register::rax] = static_cast<word_t>(value);
  return continue_result();
}

[[nodiscard]] SyscallResult return_error(CpuState& context, int error) {
  return return_value(context, -static_cast<std::int64_t>(error));
}

SyscallResult unsupported_syscall(const CpuState& context, SyscallKind kind) {
  return {.reason = StopReason::unsupported_syscall,
          .continue_execution = false,
          .message = "unsupported Linux syscall " + std::to_string(syscall_number(context)) + " via " +
                     (kind == SyscallKind::syscall64 ? "syscall64"
                      : kind == SyscallKind::int80   ? "int80"
                                                     : "sysenter")};
}

[[nodiscard]] int memory_error_to_linux(MemoryError error) noexcept {
  switch (error) {
  case MemoryError::unaligned_address:
  case MemoryError::zero_size:
    return linux_einval;
  case MemoryError::unmapped_address:
    return linux_efault;
  case MemoryError::out_of_memory:
  case MemoryError::mapping_failed:
    return linux_enomem;
  }
  return linux_eio;
}

[[nodiscard]] std::optional<int> checked_fd(word_t raw_fd) noexcept {
  if (raw_fd > static_cast<word_t>(std::numeric_limits<int>::max()))
    return std::nullopt;
  return static_cast<int>(raw_fd);
}

[[nodiscard]] bool fits_host_transfer(word_t count) noexcept {
  return count <= static_cast<word_t>(std::numeric_limits<std::intptr_t>::max());
}

[[nodiscard]] bool range_overflows(address_t start, word_t length) noexcept {
  return length != 0 && start > std::numeric_limits<address_t>::max() - (length - 1);
}

[[nodiscard]] bool is_linux_at_fdcwd(word_t raw_fd) noexcept {
  return static_cast<std::int32_t>(raw_fd & 0xffffffffu) == -100;
}

[[nodiscard]] std::uint64_t read_le(std::span<const std::byte> bytes, std::size_t offset, std::size_t width) noexcept {
  std::uint64_t value = 0;
  for (std::size_t i = 0; i < width; ++i)
    value |= static_cast<std::uint64_t>(std::to_integer<unsigned char>(bytes[offset + i])) << (i * 8);
  return value;
}

[[nodiscard]] CStringResult read_c_string(AddressSpace& space, address_t address, std::size_t max_length) noexcept {
  if (address == 0)
    return {.ok = false, .error = linux_efault, .value = {}};

  try {
    std::string result;
    result.reserve(128);

    std::array<std::byte, AddressSpace::kPageSize> buffer{};
    std::size_t total = 0;
    while (total < max_length) {
      if (range_overflows(address, total + 1))
        return {.ok = false, .error = linux_efault, .value = {}};

      const address_t current = address + total;
      const std::size_t page_remaining = AddressSpace::kPageSize - (current & (AddressSpace::kPageSize - 1));
      const std::size_t chunk_size = std::min(max_length - total, page_remaining);
      if (range_overflows(address, total + chunk_size))
        return {.ok = false, .error = linux_efault, .value = {}};

      auto read = space.read(current, std::span<std::byte>(buffer.data(), chunk_size));
      if (!read)
        return {.ok = false, .error = memory_error_to_linux(read.error()), .value = {}};

      for (std::size_t i = 0; i < chunk_size; ++i) {
        const char ch = static_cast<char>(std::to_integer<unsigned char>(buffer[i]));
        if (ch == '\0')
          return {.ok = true, .error = 0, .value = std::move(result)};

        result.push_back(ch);
      }
      total += chunk_size;
    }
  } catch (const std::bad_alloc&) {
    return {.ok = false, .error = linux_enomem, .value = {}};
  } catch (...) {
    return {.ok = false, .error = linux_eio, .value = {}};
  }

  return {.ok = false, .error = linux_enametoolong, .value = {}};
}

// --- portable-only helpers (not part of the shared cross-TU surface) -------

[[nodiscard]] std::optional<std::uint64_t> page_align_up(std::uint64_t value) noexcept {
  constexpr std::uint64_t mask = AddressSpace::kPageSize - 1;
  if (value > std::numeric_limits<std::uint64_t>::max() - mask)
    return std::nullopt;
  return (value + mask) & ~mask;
}

[[nodiscard]] bool page_aligned(std::uint64_t value) noexcept {
  return (value & (AddressSpace::kPageSize - 1)) == 0;
}

[[nodiscard]] bool range_is_unmapped(AddressSpace& space, address_t start, word_t length) noexcept {
  if (length == 0)
    return true;
  if (range_overflows(start, length))
    return false;

  std::array<std::byte, 1> byte{};
  const address_t end = start + length;
  for (address_t page = start; page < end;) {
    if (space.read(page, byte))
      return false;

    const address_t next_page = (page + AddressSpace::kPageSize) & ~(AddressSpace::kPageSize - 1);
    if (next_page <= page)
      return false;
    page = next_page;
  }
  return true;
}

[[nodiscard]] std::optional<address_t> find_unmapped_range(AddressSpace& space, address_t& cursor,
                                                           word_t length) noexcept {
  auto candidate = page_align_up(cursor);
  if (!candidate)
    return std::nullopt;

  while (*candidate <= std::numeric_limits<address_t>::max() - length) {
    if (range_is_unmapped(space, *candidate, length)) {
      cursor = *candidate + length;
      return *candidate;
    }
    if (*candidate > std::numeric_limits<address_t>::max() - AddressSpace::kPageSize)
      return std::nullopt;
    *candidate += AddressSpace::kPageSize;
  }
  return std::nullopt;
}

[[nodiscard]] std::optional<Protection> protection_from_linux(word_t prot) noexcept {
  if ((prot & ~(linux_prot_read | linux_prot_write | linux_prot_exec)) != 0)
    return std::nullopt;

  Protection result = Protection::none;
  if ((prot & linux_prot_read) != 0)
    result = result | Protection::read;
  if ((prot & linux_prot_write) != 0)
    result = result | Protection::write;
  if ((prot & linux_prot_exec) != 0)
    result = result | Protection::execute;
  return result;
}

[[nodiscard]] std::optional<int> write_word(AddressSpace& space, address_t address, word_t value) noexcept {
  std::array<std::byte, sizeof(word_t)> bytes{};
  for (std::size_t i = 0; i < bytes.size(); ++i)
    bytes[i] = static_cast<std::byte>((value >> (i * 8)) & 0xff);

  auto written = space.write(address, bytes);
  if (!written)
    return memory_error_to_linux(written.error());
  return std::nullopt;
}

template<std::size_t Size>
void copy_c_string(std::array<std::byte, Size>& bytes, std::size_t offset, const char* value) noexcept {
  const std::size_t length = std::min<std::size_t>(std::strlen(value), 64);
  for (std::size_t i = 0; i < length; ++i)
    bytes[offset + i] = static_cast<std::byte>(value[i]);
}

[[nodiscard]] SysPrlimit64::Limit read_rlimit64(std::span<const std::byte> bytes) noexcept {
  return {.current = read_le(bytes, 0, 8), .maximum = read_le(bytes, 8, 8)};
}

[[nodiscard]] std::optional<int> write_rlimit64(AddressSpace& space, address_t address,
                                                SysPrlimit64::Limit limit) noexcept {
  std::array<std::byte, 16> bytes{};
  write_le(bytes, 0, limit.current, 8);
  write_le(bytes, 8, limit.maximum, 8);

  auto written = space.write(address, bytes);
  if (!written)
    return memory_error_to_linux(written.error());
  return std::nullopt;
}

// --- brk bookkeeping ------------------------------------------------------

namespace {

std::unordered_map<ProcessId, BrkState> brk_states;

} // namespace

[[nodiscard]] BrkState& brk_state(ProcessId context_id, address_t initial_break) {
  auto [it, inserted] = brk_states.try_emplace(context_id);
  if (inserted)
    it->second = BrkState{.minimum_break = initial_break, .current_break = initial_break, .mapped_end = initial_break};
  return it->second;
}

void clone_brk_state(ProcessId parent, ProcessId child) {
  brk_states[child] = brk_state(parent);
}

void reset_brk_state(ProcessId context_id, address_t initial_break) {
  brk_states[context_id] =
      BrkState{.minimum_break = initial_break, .current_break = initial_break, .mapped_end = initial_break};
}

void erase_brk_state(ProcessId context_id) {
  brk_states.erase(context_id);
}

// --- file-backed mmap seam ------------------------------------------------

namespace {

MmapFileLoader g_mmap_file_loader = nullptr;

} // namespace

[[nodiscard]] MmapFileLoader mmap_file_loader() noexcept {
  return g_mmap_file_loader;
}

void set_mmap_file_loader(MmapFileLoader loader) noexcept {
  g_mmap_file_loader = loader;
}

} // namespace detail

namespace abi {

[[nodiscard]] word_t syscall_number(const RegisterFile& context) noexcept {
  return detail::syscall_number(context);
}

[[nodiscard]] word_t syscall_arg(const CpuState& context, std::size_t index) noexcept {
  return detail::syscall_arg(context, index);
}

[[nodiscard]] SyscallResult return_value(CpuState& context, std::int64_t value) noexcept {
  return detail::return_value(context, value);
}

} // namespace abi

std::optional<SyscallResult> SysBrk::try_syscall(Machine& machine, ProcessId context_id, CpuState& context,
                                                 AddressSpace& space, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_brk))
    return std::nullopt;

  detail::BrkState& state = detail::brk_state(context_id, initial_break);
  const address_t requested_break = detail::syscall_arg(context, 0);
  if (requested_break == 0)
    return detail::return_value(context, static_cast<std::int64_t>(state.current_break));
  if (requested_break < state.minimum_break)
    return detail::return_value(context, static_cast<std::int64_t>(state.current_break));

  auto requested_end = detail::page_align_up(requested_break);
  if (!requested_end)
    return detail::return_value(context, static_cast<std::int64_t>(state.current_break));

  if (*requested_end > state.mapped_end) {
    auto mapped = space.map(state.mapped_end, *requested_end - state.mapped_end, Protection::read | Protection::write);
    if (!mapped)
      return detail::return_value(context, static_cast<std::int64_t>(state.current_break));
    state.mapped_end = *requested_end;
  }

  state.current_break = requested_break;
  return detail::return_value(context, static_cast<std::int64_t>(state.current_break));
}

std::optional<SyscallResult> SysMmap::try_syscall(Machine& machine, ProcessId context_id, CpuState& context,
                                                  AddressSpace& space, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_mmap))
    return std::nullopt;

  const address_t requested_address = detail::syscall_arg(context, 0);
  const word_t length = detail::syscall_arg(context, 1);
  const word_t prot = detail::syscall_arg(context, 2);
  const word_t flags = detail::syscall_arg(context, 3);
  const word_t offset = detail::syscall_arg(context, 5);

  if (length == 0 || !detail::page_aligned(offset))
    return detail::return_error(context, detail::linux_einval);

  const bool private_mapping = (flags & detail::linux_map_private) != 0;
  const bool shared_mapping = (flags & detail::linux_map_shared) != 0;
  if (private_mapping == shared_mapping)
    return detail::return_error(context, detail::linux_einval);

  constexpr word_t supported_flags = detail::linux_map_shared | detail::linux_map_private | detail::linux_map_fixed |
                                     detail::linux_map_anonymous | detail::linux_map_denied_write |
                                     detail::linux_map_executable | detail::linux_map_noreserve |
                                     detail::linux_map_populate | detail::linux_map_stack;
  if ((flags & ~supported_flags) != 0)
    return detail::return_error(context, detail::linux_einval);

  auto protection = detail::protection_from_linux(prot);
  if (!protection)
    return detail::return_error(context, detail::linux_einval);

  auto aligned_length = detail::page_align_up(length);
  if (!aligned_length)
    return detail::return_error(context, detail::linux_enomem);

  address_t mapped_address = 0;
  if ((flags & detail::linux_map_fixed) != 0) {
    if (!detail::page_aligned(requested_address))
      return detail::return_error(context, detail::linux_einval);
    mapped_address = requested_address;
  } else {
    auto allocated = detail::find_unmapped_range(space, next_mapping_address, *aligned_length);
    if (!allocated)
      return detail::return_error(context, detail::linux_enomem);
    mapped_address = *allocated;
  }

  if (mapped_address > std::numeric_limits<address_t>::max() - *aligned_length)
    return detail::return_error(context, detail::linux_einval);

  auto mapped = space.map(mapped_address, *aligned_length, *protection);
  if (!mapped)
    return detail::return_error(context, detail::memory_error_to_linux(mapped.error()));

  if ((flags & detail::linux_map_anonymous) == 0) {
    // File-backed mappings copy bytes from a host fd, which only the POSIX
    // library can do. The portable library exposes a loader hook the POSIX TU
    // registers at startup; if it is absent (portable-only build), file-backed
    // mmap is unsupported.
    auto loader = detail::mmap_file_loader();
    if (!loader) {
      space.unmap(mapped_address, *aligned_length);
      return detail::return_error(context, detail::linux_ebadf);
    }
    if (auto error = loader(context_id, space, mapped_address, length, detail::syscall_arg(context, 4), offset)) {
      space.unmap(mapped_address, *aligned_length);
      return detail::return_error(context, *error);
    }
  }

  return detail::return_value(context, static_cast<std::int64_t>(mapped_address));
}

std::optional<SyscallResult> SysMunmap::try_syscall(Machine& machine, ProcessId context_id, CpuState& context,
                                                    AddressSpace& space, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_munmap))
    return std::nullopt;

  const address_t address = detail::syscall_arg(context, 0);
  const word_t length = detail::syscall_arg(context, 1);
  if (length == 0 || !detail::page_aligned(address))
    return detail::return_error(context, detail::linux_einval);

  auto aligned_length = detail::page_align_up(length);
  if (!aligned_length)
    return detail::return_error(context, detail::linux_einval);

  space.unmap(address, *aligned_length);
  return detail::return_value(context, 0);
}

std::optional<SyscallResult> SysMremap::try_syscall(Machine& machine, ProcessId context_id, CpuState& context,
                                                    AddressSpace& space, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_mremap))
    return std::nullopt;

  const address_t old_address = detail::syscall_arg(context, 0);
  const word_t old_size = detail::syscall_arg(context, 1);
  const word_t new_size = detail::syscall_arg(context, 2);
  const word_t flags = detail::syscall_arg(context, 3);
  const address_t fixed_address = detail::syscall_arg(context, 4);

  constexpr word_t supported_flags =
      detail::linux_mremap_maymove | detail::linux_mremap_fixed | detail::linux_mremap_dontunmap;
  if ((flags & ~supported_flags) != 0)
    return detail::return_error(context, detail::linux_einval);
  if ((flags & detail::linux_mremap_fixed) != 0 && (flags & detail::linux_mremap_maymove) == 0)
    return detail::return_error(context, detail::linux_einval);
  if ((flags & detail::linux_mremap_dontunmap) != 0)
    return detail::return_error(context, detail::linux_einval);
  if (old_size == 0 || new_size == 0 || !detail::page_aligned(old_address))
    return detail::return_error(context, detail::linux_einval);

  auto old_length = detail::page_align_up(old_size);
  auto new_length = detail::page_align_up(new_size);
  if (!old_length || !new_length)
    return detail::return_error(context, detail::linux_enomem);
  if (detail::range_overflows(old_address, *old_length))
    return detail::return_error(context, detail::linux_efault);

  try {
    const word_t copy_length = std::min(*old_length, *new_length);
    std::vector<std::byte> copied(static_cast<std::size_t>(copy_length));
    if (copy_length != 0) {
      auto read = space.read(old_address, copied);
      if (!read)
        return detail::return_error(context, detail::memory_error_to_linux(read.error()));
    }

    if (*new_length <= *old_length) {
      space.unmap(old_address + *new_length, *old_length - *new_length);
      return detail::return_value(context, static_cast<std::int64_t>(old_address));
    }

    const address_t old_end = old_address + *old_length;
    const word_t growth = *new_length - *old_length;
    if (old_end <= std::numeric_limits<address_t>::max() - growth &&
        detail::range_is_unmapped(space, old_end, growth)) {
      auto mapped = space.map(old_end, growth, Protection::read | Protection::write);
      if (mapped) {
        next_mapping_address = std::max(next_mapping_address, old_end + growth);
        return detail::return_value(context, static_cast<std::int64_t>(old_address));
      }
    }

    if ((flags & detail::linux_mremap_maymove) == 0)
      return detail::return_error(context, detail::linux_enomem);

    address_t new_address = 0;
    if ((flags & detail::linux_mremap_fixed) != 0) {
      if (!detail::page_aligned(fixed_address))
        return detail::return_error(context, detail::linux_einval);
      if (!detail::range_is_unmapped(space, fixed_address, *new_length))
        space.unmap(fixed_address, *new_length);
      new_address = fixed_address;
    } else {
      auto allocated = detail::find_unmapped_range(space, next_mapping_address, *new_length);
      if (!allocated)
        return detail::return_error(context, detail::linux_enomem);
      new_address = *allocated;
    }

    auto mapped = space.map(new_address, *new_length, Protection::read | Protection::write);
    if (!mapped)
      return detail::return_error(context, detail::memory_error_to_linux(mapped.error()));

    if (!copied.empty()) {
      auto written = space.write(new_address, copied);
      if (!written) {
        space.unmap(new_address, *new_length);
        return detail::return_error(context, detail::memory_error_to_linux(written.error()));
      }
    }

    space.unmap(old_address, *old_length);
    return detail::return_value(context, static_cast<std::int64_t>(new_address));
  } catch (const std::bad_alloc&) {
    return detail::return_error(context, detail::linux_enomem);
  } catch (...) {
    return detail::return_error(context, detail::linux_eio);
  }
}

std::optional<SyscallResult> SysArchPrctl::try_syscall(Machine& machine, ProcessId context_id, CpuState& context,
                                                       AddressSpace& space, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_arch_prctl))
    return std::nullopt;

  const word_t code = detail::syscall_arg(context, 0);
  const address_t address = detail::syscall_arg(context, 1);

  // Segment bases live directly in the CpuState the caller owns; arch_prctl just
  // reads/writes them.
  constexpr auto fs = segment_register_index(SegmentRegister::fs);
  constexpr auto gs = segment_register_index(SegmentRegister::gs);
  switch (code) {
  case detail::linux_arch_set_fs:
    context.segment_bases[fs] = address;
    return detail::return_value(context, 0);
  case detail::linux_arch_set_gs:
    context.segment_bases[gs] = address;
    return detail::return_value(context, 0);
  case detail::linux_arch_get_fs:
    if (auto error = detail::write_word(space, address, context.segment_bases[fs]))
      return detail::return_error(context, *error);
    return detail::return_value(context, 0);
  case detail::linux_arch_get_gs:
    if (auto error = detail::write_word(space, address, context.segment_bases[gs]))
      return detail::return_error(context, *error);
    return detail::return_value(context, 0);
  default:
    return detail::return_error(context, detail::linux_einval);
  }
}

std::optional<SyscallResult> SysGetIdentity::try_syscall(Machine&, ProcessId context_id, CpuState& context,
                                                         AddressSpace& space, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_getpid) &&
      !detail::handles(context, kind, detail::syscall_gettid) &&
      !detail::handles(context, kind, detail::syscall_getuid) &&
      !detail::handles(context, kind, detail::syscall_geteuid) &&
      !detail::handles(context, kind, detail::syscall_getgid) &&
      !detail::handles(context, kind, detail::syscall_getegid)) {
    return std::nullopt;
  }

  switch (detail::syscall_number(context)) {
  case detail::syscall_getpid:
  case detail::syscall_gettid:
    return detail::return_value(context, static_cast<std::int64_t>(context_id));
  case detail::syscall_getuid:
  case detail::syscall_geteuid:
    return detail::return_value(context, detail::synthetic_uid);
  case detail::syscall_getgid:
  case detail::syscall_getegid:
    return detail::return_value(context, detail::synthetic_gid);
  default:
    return std::nullopt;
  }
}

std::optional<SyscallResult> SysUname::try_syscall(Machine& machine, ProcessId context_id, CpuState& context,
                                                   AddressSpace& space, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_uname))
    return std::nullopt;

  constexpr std::size_t field_size = 65;
  std::array<std::byte, field_size * 6> utsname{};
  detail::copy_c_string(utsname, field_size * 0, "Linux");
  detail::copy_c_string(utsname, field_size * 1, "x86sim");
  detail::copy_c_string(utsname, field_size * 2, "6.0.0");
  detail::copy_c_string(utsname, field_size * 3, "#1");
  detail::copy_c_string(utsname, field_size * 4, "x86_64");
  detail::copy_c_string(utsname, field_size * 5, "localdomain");

  auto written = space.write(detail::syscall_arg(context, 0), utsname);
  if (!written)
    return detail::return_error(context, detail::memory_error_to_linux(written.error()));
  return detail::return_value(context, 0);
}

std::optional<SyscallResult> SysGetcwd::try_syscall(Machine& machine, ProcessId context_id, CpuState& context,
                                                    AddressSpace& space, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_getcwd))
    return std::nullopt;

  const address_t buffer_address = detail::syscall_arg(context, 0);
  const word_t size = detail::syscall_arg(context, 1);
  if (buffer_address == 0)
    return detail::return_error(context, detail::linux_efault);
  if (size == 0)
    return detail::return_error(context, detail::linux_einval);
  if (detail::range_overflows(buffer_address, size))
    return detail::return_error(context, detail::linux_efault);

  try {
    std::string cwd = std::filesystem::current_path().string();
    if (cwd.size() + 1 > size)
      return detail::return_error(context, ERANGE);

    std::vector<std::byte> bytes(cwd.size() + 1);
    std::memcpy(bytes.data(), cwd.data(), cwd.size());
    auto written = space.write(buffer_address, bytes);
    if (!written)
      return detail::return_error(context, detail::memory_error_to_linux(written.error()));
    return detail::return_value(context, static_cast<std::int64_t>(cwd.size() + 1));
  } catch (const std::bad_alloc&) {
    return detail::return_error(context, detail::linux_enomem);
  } catch (...) {
    return detail::return_error(context, detail::linux_eio);
  }
}

std::optional<SyscallResult> SysSetTidAddress::try_syscall(Machine&, ProcessId, CpuState& context, AddressSpace& space,
                                                           SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_set_tid_address))
    return std::nullopt;

  return detail::return_value(context, 1);
}

std::optional<SyscallResult> SysSetRobustList::try_syscall(Machine&, ProcessId, CpuState& context, AddressSpace& space,
                                                           SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_set_robust_list))
    return std::nullopt;

  return detail::return_value(context, 0);
}

std::optional<SyscallResult> SysRseq::try_syscall(Machine&, ProcessId, CpuState& context, AddressSpace& space,
                                                  SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_rseq))
    return std::nullopt;

  return detail::return_error(context, detail::linux_enosys);
}

std::optional<SyscallResult> SysPrlimit64::try_syscall(Machine& machine, ProcessId context_id, CpuState& context,
                                                       AddressSpace& space, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_prlimit64))
    return std::nullopt;

  const word_t pid = detail::syscall_arg(context, 0);
  if (pid != 0 && pid != context_id)
    return detail::return_error(context, detail::linux_esrch);

  const word_t resource = detail::syscall_arg(context, 1);
  if (resource >= limits.size())
    return detail::return_error(context, detail::linux_einval);

  const address_t new_limit_address = detail::syscall_arg(context, 2);
  const address_t old_limit_address = detail::syscall_arg(context, 3);
  const Limit old_limit = limits[resource];
  Limit new_limit = old_limit;

  if (new_limit_address != 0) {
    if (detail::range_overflows(new_limit_address, 16))
      return detail::return_error(context, detail::linux_efault);

    std::array<std::byte, 16> bytes{};
    if (auto read_error = detail::read_guest_memory(space, new_limit_address, bytes))
      return detail::return_error(context, *read_error);

    new_limit = detail::read_rlimit64(bytes);
    if (new_limit.current > new_limit.maximum)
      return detail::return_error(context, detail::linux_einval);
  }

  if (old_limit_address != 0) {
    if (detail::range_overflows(old_limit_address, 16))
      return detail::return_error(context, detail::linux_efault);
    if (auto error = detail::write_rlimit64(space, old_limit_address, old_limit))
      return detail::return_error(context, *error);
  }

  if (new_limit_address != 0)
    limits[resource] = new_limit;

  return detail::return_value(context, 0);
}

std::optional<SyscallResult> SysSignals::try_syscall(Machine&, ProcessId, CpuState& context, AddressSpace& space,
                                                     SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_rt_sigaction) &&
      !detail::handles(context, kind, detail::syscall_rt_sigprocmask)) {
    return std::nullopt;
  }

  return detail::return_value(context, 0);
}

std::optional<SyscallResult> SysFutex::try_syscall(Machine& machine, ProcessId context_id, CpuState& context,
                                                   AddressSpace& space, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_futex))
    return std::nullopt;

  const address_t uaddr = detail::syscall_arg(context, 0);
  const word_t op = detail::syscall_arg(context, 1);
  const std::uint32_t expected = static_cast<std::uint32_t>(detail::syscall_arg(context, 2));
  const address_t timeout = detail::syscall_arg(context, 3);
  const word_t command = op & detail::linux_futex_cmd_mask;

  switch (command) {
  case detail::linux_futex_wait:
  case detail::linux_futex_wait_bitset: {
    if (uaddr == 0 || detail::range_overflows(uaddr, 4))
      return detail::return_error(context, detail::linux_efault);

    std::array<std::byte, 4> bytes{};
    if (auto read_error = detail::read_guest_memory(space, uaddr, bytes))
      return detail::return_error(context, *read_error);

    const auto current = static_cast<std::uint32_t>(detail::read_le(bytes, 0, 4));
    if (current != expected)
      return detail::return_error(context, detail::linux_eagain);

    if (timeout != 0)
      return detail::return_error(context, detail::linux_etimedout);

    return detail::return_value(context, 0);
  }
  case detail::linux_futex_wake:
  case detail::linux_futex_wake_bitset:
    return detail::return_value(context, 0);
  default:
    return detail::return_error(context, detail::linux_enosys);
  }
}

std::optional<SyscallResult> SysProcess::try_syscall(Machine& machine, ProcessId context_id, CpuState& context,
                                                     AddressSpace& space, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_clone) && !detail::handles(context, kind, detail::syscall_fork) &&
      !detail::handles(context, kind, detail::syscall_vfork) &&
      !detail::handles(context, kind, detail::syscall_execve) &&
      !detail::handles(context, kind, detail::syscall_wait4)) {
    return std::nullopt;
  }

  const word_t number = detail::syscall_number(context);

  // Phase 1 runs a single process. fork/clone/vfork/execve/wait4 all require the
  // multi-process scheduler that Phase 2 rebuilds on top of run(CpuState&,
  // AddressSpace&); until then they degrade to single-process behavior.
  // TODO(phase2): restore real fork/exec/wait once the support library owns the
  // process table + scheduler.
  (void)machine;
  (void)context_id;
  (void)space;
  if (number == detail::syscall_wait4)
    return detail::return_error(context, detail::linux_echild);
  return detail::return_error(context, detail::linux_enosys);
}

namespace {

// Phase 1 never has a parent/child relationship (fork is disabled), so this is
// never reached; it is kept as the seam where Phase 2's scheduler will deliver a
// child's exit status to a waiting parent process.
// TODO(phase2): reinstate parent notification / wait wakeup on the scheduler.
void publish_child_exit(Machine&, ProcessTable&, ProcessId, int) {}

} // namespace

std::optional<SyscallResult> SysExit::try_syscall(Machine& machine, ProcessId context_id, CpuState& context,
                                                  AddressSpace& space, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_exit))
    return std::nullopt;

  const int status = static_cast<int>(detail::syscall_arg(context, 0) & 0xff);
  if (processes->parent.contains(context_id)) {
    publish_child_exit(machine, *processes, context_id, status);
    return detail::return_value(context, 0);
  }
  if (exit_status)
    *exit_status = status;
  return SyscallResult{.reason = StopReason::guest_exit, .continue_execution = false, .message = {}};
}

std::optional<SyscallResult> SysExitGroup::try_syscall(Machine& machine, ProcessId context_id, CpuState& context,
                                                       AddressSpace& space, SyscallKind kind) noexcept {
  if (!detail::handles(context, kind, detail::syscall_exit_group))
    return std::nullopt;

  const int status = static_cast<int>(detail::syscall_arg(context, 0) & 0xff);
  if (processes->parent.contains(context_id)) {
    publish_child_exit(machine, *processes, context_id, status);
    return detail::return_value(context, 0);
  }
  if (exit_status)
    *exit_status = status;
  return SyscallResult{.reason = StopReason::guest_exit, .continue_execution = false, .message = {}};
}

} // namespace x86sim::linux_syscalls
