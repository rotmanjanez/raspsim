#include "x86sim/addrspace.hpp"

#include "globals.h"
#include "x86sim/logging.hpp"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <new>

namespace x86sim {
namespace {

constexpr address_t address_space_bits = 48;
constexpr address_t address_space_size = address_t{1} << address_space_bits;
constexpr address_t address_space_mask = address_space_size - 1;
constexpr address_t page_bits = 12;
constexpr address_t page_size = AddressSpace::kPageSize;
constexpr address_t page_mask = page_size - 1;
constexpr std::size_t spat_toplevel_chunk_bits = 17;
constexpr std::size_t spat_pages_per_chunk_bits = 19;
constexpr std::size_t spat_toplevel_chunks = std::size_t{1} << spat_toplevel_chunk_bits;
constexpr std::size_t spat_pages_per_chunk = std::size_t{1} << spat_pages_per_chunk_bits;
constexpr std::size_t spat_bytes_per_chunk = spat_pages_per_chunk / 8;

[[nodiscard]] address_t page_floor(address_t value) noexcept {
  return value & ~page_mask;
}

[[nodiscard]] std::optional<address_t> page_ceil(address_t value) noexcept {
  if (value > std::numeric_limits<address_t>::max() - page_mask)
    return std::nullopt;
  return (value + page_mask) & ~page_mask;
}

[[nodiscard]] bool range_overflows(address_t start, std::uint64_t size) noexcept {
  return size != 0 && start > std::numeric_limits<address_t>::max() - (size - 1);
}

} // namespace

struct AddressSpace::ShadowMap {
  using Chunk = std::array<std::byte, spat_bytes_per_chunk>;

  std::array<std::unique_ptr<Chunk>, spat_toplevel_chunks> chunks;

  [[nodiscard]] std::byte& byte_for_page(address_t pageid) {
    const address_t chunkid = pageid >> spat_pages_per_chunk_bits;
    if (!chunks[chunkid])
      chunks[chunkid] = std::make_unique<Chunk>();
    const address_t byteid = (pageid >> 3) & (spat_bytes_per_chunk - 1);
    return (*chunks[chunkid])[byteid];
  }

  [[nodiscard]] bool bit_for_page(address_t pageid) const noexcept {
    const address_t chunkid = pageid >> spat_pages_per_chunk_bits;
    if (!chunks[chunkid])
      return false;
    const address_t byteid = (pageid >> 3) & (spat_bytes_per_chunk - 1);
    const auto value = std::to_integer<unsigned>((*chunks[chunkid])[byteid]);
    return (value & (1u << (pageid & 7u))) != 0;
  }

  void set_page(address_t pageid) {
    std::byte& value = byte_for_page(pageid);
    value = static_cast<std::byte>(std::to_integer<unsigned>(value) | (1u << (pageid & 7u)));
  }

  void clear_page(address_t pageid) {
    std::byte& value = byte_for_page(pageid);
    value = static_cast<std::byte>(std::to_integer<unsigned>(value) & ~(1u << (pageid & 7u)));
  }
};

AddressSpace::AddressSpace()
    : readmap_(std::make_unique<ShadowMap>()), writemap_(std::make_unique<ShadowMap>()),
      execmap_(std::make_unique<ShadowMap>()), dirtymap_(std::make_unique<ShadowMap>()) {}

AddressSpace::AddressSpace(const AddressSpace& other) : AddressSpace() {
  for (const auto& [address, page] : other.mapped_mem_) {
    map_or_throw(address, kPageSize, other.protection_at(address));
    std::memcpy(mapped_mem_.at(address)->data(), page->data(), kPageSize);
  }
}

AddressSpace& AddressSpace::operator=(const AddressSpace& other) {
  if (this == &other)
    return *this;
  AddressSpace copy(other);
  *this = std::move(copy);
  return *this;
}

AddressSpace::AddressSpace(AddressSpace&&) noexcept = default;

AddressSpace& AddressSpace::operator=(AddressSpace&&) noexcept = default;

AddressSpace::~AddressSpace() = default;

std::expected<void, MemoryError> AddressSpace::map(address_t start, std::uint64_t size, Protection prot) noexcept {
  if (size == 0)
    return std::unexpected(MemoryError::zero_size);
  if (range_overflows(start, size))
    return std::unexpected(MemoryError::mapping_failed);

  try {
    map_or_throw(start, size, prot);
  } catch (const std::bad_alloc&) {
    return std::unexpected(MemoryError::out_of_memory);
  } catch (...) {
    return std::unexpected(MemoryError::mapping_failed);
  }
  return {};
}

void AddressSpace::unmap(address_t start, std::uint64_t size) noexcept {
  if (size == 0 || range_overflows(start, size))
    return;

  start = page_floor(start);
  const auto end = page_ceil(start + size);
  if (!end)
    return;

  for (address_t page = start; page < *end; page += kPageSize)
    mapped_mem_.erase(page);
  setattr(start, *end - start, Protection::none);
}

std::expected<void, MemoryError> AddressSpace::read(address_t start, std::span<std::byte> out) const noexcept {
  if (out.empty())
    return std::unexpected(MemoryError::zero_size);
  if (range_overflows(start, out.size()))
    return std::unexpected(MemoryError::unmapped_address);

  std::uint64_t processed = 0;
  while (processed < out.size()) {
    const address_t addr = start + processed;
    const auto* src = static_cast<const std::byte*>(page_virt_to_mapped(addr));
    if (!src)
      return std::unexpected(MemoryError::unmapped_address);

    const auto chunk = std::min<std::uint64_t>(out.size() - processed, kPageSize - (addr & page_mask));
    std::memcpy(out.data() + processed, src, chunk);
    processed += chunk;
  }
  return {};
}

std::expected<void, MemoryError> AddressSpace::write(address_t start, std::span<const std::byte> bytes) noexcept {
  if (bytes.empty())
    return std::unexpected(MemoryError::zero_size);
  if (range_overflows(start, bytes.size()))
    return std::unexpected(MemoryError::unmapped_address);

  std::uint64_t processed = 0;
  while (processed < bytes.size()) {
    const address_t addr = start + processed;
    auto* dst = static_cast<std::byte*>(page_virt_to_mapped(addr));
    if (!dst)
      return std::unexpected(MemoryError::unmapped_address);

    const auto chunk = std::min<std::uint64_t>(bytes.size() - processed, kPageSize - (addr & page_mask));
    std::memcpy(dst, bytes.data() + processed, chunk);
    processed += chunk;
  }
  return {};
}

Protection AddressSpace::protection_at(address_t address) const noexcept {
  address &= address_space_mask;
  const address_t page = pageid(address);

  Protection prot = Protection::none;
  if (readmap_->bit_for_page(page))
    prot = prot | Protection::read;
  if (writemap_->bit_for_page(page))
    prot = prot | Protection::write;
  if (execmap_->bit_for_page(page))
    prot = prot | Protection::execute;
  return prot;
}

std::unique_ptr<AddressSpace> AddressSpace::clone_deep() const {
  return std::make_unique<AddressSpace>(*this);
}

void* AddressSpace::page_virt_to_mapped(address_t addr) noexcept {
  auto it = mapped_mem_.find(page_floor(addr));
  if (it == mapped_mem_.end())
    return nullptr;

  return it->second->data() + (addr & page_mask);
}

const void* AddressSpace::page_virt_to_mapped(address_t addr) const noexcept {
  auto it = mapped_mem_.find(page_floor(addr));
  if (it == mapped_mem_.end())
    return nullptr;

  return it->second->data() + (addr & page_mask);
}

bool AddressSpace::fastcheck(address_t addr, spat_t top) const noexcept {
  if ((addr >> address_space_bits) != 0 || top == nullptr)
    return false;
  return top->bit_for_page(pageid(addr));
}

bool AddressSpace::check(address_t address, Protection prot) const noexcept {
  if (has_protection(prot, Protection::read) && !fastcheck(address, read_map()))
    return false;
  if (has_protection(prot, Protection::write) && !fastcheck(address, write_map()))
    return false;
  if (has_protection(prot, Protection::execute) && !fastcheck(address, exec_map()))
    return false;
  return true;
}

bool AddressSpace::isdirty(address_t mfn) const noexcept {
  return fastcheck(mfn << page_bits, dirtymap_.get());
}

void AddressSpace::setdirty(address_t mfn) {
  make_page_accessible(mfn << page_bits, *dirtymap_);
}

void AddressSpace::cleardirty(address_t mfn) {
  make_page_inaccessible(mfn << page_bits, *dirtymap_);
}

void AddressSpace::map_or_throw(address_t start, std::uint64_t size, Protection prot) {
  start = page_floor(start);
  const auto aligned_size = page_ceil(size);
  if (!aligned_size)
    throw std::bad_alloc();

  const address_t num_pages = *aligned_size / kPageSize;
  for (address_t i = 0; i < num_pages; ++i)
    mapped_mem_.insert_or_assign(start + i * kPageSize, std::make_unique<Page>());
  setattr(start, *aligned_size, prot);
}

void AddressSpace::setattr(address_t start, std::uint64_t size, Protection prot) {
  logging::println(
      "setattr: region {} to {} ({} KB) has user-visible attributes {}{}{}", reinterpret_cast<void*>(start),
      reinterpret_cast<void*>(start + size), size >> 10, (has_protection(prot, Protection::read) ? 'r' : '-'),
      (has_protection(prot, Protection::write) ? 'w' : '-'), (has_protection(prot, Protection::execute) ? 'x' : '-'));

  if (has_protection(prot, Protection::read))
    make_accessible(start, size, *readmap_);
  else
    make_inaccessible(start, size, *readmap_);

  if (has_protection(prot, Protection::write))
    make_accessible(start, size, *writemap_);
  else
    make_inaccessible(start, size, *writemap_);

  if (has_protection(prot, Protection::execute))
    make_accessible(start, size, *execmap_);
  else
    make_inaccessible(start, size, *execmap_);
}

void AddressSpace::make_accessible(address_t address, std::uint64_t size, ShadowMap& top) {
  address &= address_space_mask;
  const address_t firstpage = address >> page_bits;
  const address_t lastpage = (address + size - 1) >> page_bits;
  assert(((address + size + page_mask) & ~page_mask) <= address_space_size);
  for (address_t page = firstpage; page <= lastpage; ++page)
    top.set_page(page);
}

void AddressSpace::make_inaccessible(address_t address, std::uint64_t size, ShadowMap& top) {
  address &= address_space_mask;
  const address_t firstpage = address >> page_bits;
  const address_t lastpage = (address + size - 1) >> page_bits;
  assert(((address + size + page_mask) & ~page_mask) <= address_space_size);
  for (address_t page = firstpage; page <= lastpage; ++page)
    top.clear_page(page);
}

void AddressSpace::make_page_accessible(address_t address, ShadowMap& top) {
  top.set_page(pageid(address));
}

void AddressSpace::make_page_inaccessible(address_t address, ShadowMap& top) {
  top.clear_page(pageid(address));
}

address_t AddressSpace::pageid(address_t address) noexcept {
  return (address & address_space_mask) >> page_bits;
}

} // namespace x86sim
