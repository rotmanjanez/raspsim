#include "x86sim/x86sim.hpp"

#include "x86sim/addrspace.h"
#include "ooocore.h"
#include "ptlhwdef.h"
#include "ptlsim.h"
#include "ptlsim-api.h"
#include "seqcore.h"
#include "stats.h"
#include "x86sim/logging.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <format>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

namespace x86sim {
namespace {

[[nodiscard]] SyscallKind to_syscall_kind(int semantics) {
  switch (semantics) {
  case SYSCALL_SEMANTICS_INT80:
    return SyscallKind::int80;
  case SYSCALL_SEMANTICS_SYSENTER:
    return SyscallKind::sysenter;
  default:
    throw std::runtime_error("unknown syscall semantics");
  }
}

void advance_default_syscall_rip_if_unchanged(Context& context, SyscallKind kind, W64 old_rip) noexcept {
  if (context[Register::rip] != old_rip)
    return;

  switch (kind) {
  case SyscallKind::int80:
  case SyscallKind::sysenter:
    context[Register::rip] = context.commitarf[REG_nextrip];
    break;
  case SyscallKind::syscall64:
    context[Register::rip] = context.commitarf[REG_rcx];
    break;
  }
}

} // namespace

bool requested_switch_to_native = false;

W64 loadphys(Waddr addr) {
  return *reinterpret_cast<W64*>(addr);
}

W64 storemask(Waddr addr, W64 data, byte bytemask) {
  W64& mem = *reinterpret_cast<W64*>(addr);
  mem = mux64(expand_8bit_to_64bit_lut[bytemask], mem, data);
  return data;
}

int Context::copy_from_user(void* target, Waddr addr, int bytes, PageFaultErrorCode& pfec, Waddr& faultaddr,
                            bool forexec, Level1PTE& ptelo, Level1PTE& ptehi) {
  assert(machine);
  AddressSpace& asp = machine->address_space();
  logging::println("VMEM: Read from user {} ({})", reinterpret_cast<void*>(addr), bytes);
  logging::flush();

  int n = 0;
  pfec = 0;
  ptelo = 0;
  ptehi = 0;

  bool readable = asp.fastcheck(reinterpret_cast<byte*>(addr), asp.readmap);
  bool executable = true;
  if likely (forexec)
    executable = asp.fastcheck(reinterpret_cast<byte*>(addr), asp.execmap);

  if unlikely ((!readable) | (forexec & !executable)) {
    faultaddr = addr;
    pfec.p = readable;
    pfec.nx = (forexec & (!executable));
    pfec.us = 1;
    return n;
  }

  n = std::min(static_cast<Waddr>(4096 - lowbits(addr, 12)), static_cast<Waddr>(bytes));

  void* mapped_addr = asp.page_virt_to_mapped(addr);
  assert(mapped_addr);
  logging::println("VMEM: Read {} = {}[{:016x}]", mapped_addr, *reinterpret_cast<W8*>(mapped_addr),
                   *reinterpret_cast<W8*>(mapped_addr));
  logging::flush();
  std::memcpy(target, mapped_addr, n);

  if likely (n == bytes)
    return n;

  readable = asp.fastcheck(reinterpret_cast<byte*>(addr + n), asp.readmap);
  executable = true;
  if likely (forexec)
    executable = asp.fastcheck(reinterpret_cast<byte*>(addr + n), asp.execmap);

  if unlikely ((!readable) | (forexec & !executable)) {
    faultaddr = addr + n;
    pfec.p = readable;
    pfec.nx = (forexec & (!executable));
    pfec.us = 1;
    return n;
  }

  std::memcpy(reinterpret_cast<byte*>(target) + n, asp.page_virt_to_mapped(addr + n), bytes - n);
  return bytes;
}

int Context::copy_to_user(Waddr target, void* source, int bytes, PageFaultErrorCode& pfec, Waddr& faultaddr) {
  assert(machine);
  AddressSpace& asp = machine->address_space();
  logging::println("VMEM: Write to user {} ({})", reinterpret_cast<void*>(target), bytes);
  logging::flush();

  pfec = 0;
  bool writable = asp.fastcheck(reinterpret_cast<byte*>(target), asp.writemap);
  if unlikely (!writable) {
    faultaddr = target;
    pfec.p = asp.fastcheck(reinterpret_cast<byte*>(target), asp.readmap);
    pfec.rw = 1;
    return 0;
  }

  byte* targetlo = reinterpret_cast<byte*>(asp.page_virt_to_mapped(target));
  int nlo = std::min(static_cast<Waddr>(4096 - lowbits(target, 12)), static_cast<Waddr>(bytes));

  asp.setdirty(target >> 12);

  if likely (nlo == bytes) {
    std::memcpy(targetlo, source, nlo);
    return bytes;
  }

  writable = asp.fastcheck(reinterpret_cast<byte*>(target + nlo), asp.writemap);
  if unlikely (!writable) {
    faultaddr = target + nlo;
    pfec.p = asp.fastcheck(reinterpret_cast<byte*>(target + nlo), asp.readmap);
    pfec.rw = 1;
    pfec.us = 1;
    return nlo;
  }

  std::memcpy(asp.page_virt_to_mapped(target + nlo), reinterpret_cast<byte*>(source) + nlo, bytes - nlo);
  std::memcpy(targetlo, source, nlo);

  asp.setdirty((target + nlo) >> 12);
  return bytes;
}

Waddr Context::check_and_translate(Waddr virtaddr, int sizeshift, bool store, bool internal, int& exception,
                                   PageFaultErrorCode& pfec, PTEUpdate& pteupdate, Level1PTE& pteused) {
  assert(machine);
  AddressSpace& asp = machine->address_space();
  exception = 0;
  pteupdate = 0;
  pteused = 0;
  pfec = 0;

  if unlikely (lowbits(virtaddr, sizeshift)) {
    exception = EXCEPTION_UnalignedAccess;
    return INVALID_PHYSADDR;
  }

  if unlikely (internal)
    return virtaddr;

  AddressSpace::spat_t top = store ? asp.writemap : asp.readmap;
  if unlikely (!asp.fastcheck(virtaddr, top)) {
    exception = store ? EXCEPTION_PageFaultOnWrite : EXCEPTION_PageFaultOnRead;
    pfec.p = asp.fastcheck(virtaddr, asp.readmap);
    pfec.rw = store;
    pfec.us = 1;
    return 0;
  }

  return reinterpret_cast<Waddr>(asp.page_virt_to_mapped(floor(signext64(virtaddr, 48), 8)));
}

int Context::write_segreg(unsigned int, W16) {
  return EXCEPTION_x86_gp_fault;
}

void Context::update_shadow_segment_descriptors() {
  W64 limit = use64 ? 0xffffffffffffffffULL : 0xffffffffULL;

  for (SegmentDescriptorCache* seg_cache : {&seg[SEGID_CS], &seg[SEGID_SS], &seg[SEGID_DS], &seg[SEGID_ES],
                                            &seg[SEGID_FS], &seg[SEGID_GS]}) {
    seg_cache->present = 1;
    seg_cache->base = 0;
    seg_cache->limit = limit;
  }

  virt_addr_mask = limit;
}

RIPVirtPhys& RIPVirtPhys::update(Context& ctx, int bytes) {
  use64 = ctx.use64;
  kernel = 0;
  df = ((ctx.internal_eflags & FLAG_DF) != 0);
  padlo = 0;
  padhi = 0;
  mfnlo = rip >> 12;
  mfnhi = (rip + (bytes - 1)) >> 12;
  return *this;
}

void assist_ptlcall(Context& ctx) {
  requested_switch_to_native = 1;
  ctx.commitarf[REG_rip] = ctx.commitarf[REG_nextrip];
}

Machine::Machine(HostCallbacks& callbacks, Options options)
    : callbacks_(callbacks), options_(std::move(options)), address_space_(std::make_unique<AddressSpace>()) {
  init_uops();

  handle_config_change(options_, 0, nullptr);

  switch (options_.core) {
  case CoreModel::sequential:
    machine_ = std::make_unique<SequentialMachine>(*this, options_);
    break;
  case CoreModel::out_of_order:
    machine_ = std::make_unique<OutOfOrderModel::OutOfOrderMachine>(*this, options_);
    break;
  }
}

Machine::~Machine() = default;

std::expected<void, MemoryError> Machine::map(address_t start, std::uint64_t size, Protection prot) noexcept {
  if (size == 0)
    return std::unexpected(MemoryError::zero_size);

  try {
    address_space_->map(start, size, prot);
  } catch (const std::bad_alloc&) {
    return std::unexpected(MemoryError::out_of_memory);
  } catch (...) {
    return std::unexpected(MemoryError::mapping_failed);
  }
  return {};
}

void Machine::unmap(address_t start, std::uint64_t size) noexcept {
  if (size == 0)
    return;
  address_space_->unmap(start, size);
}

std::expected<std::span<const std::byte>, MemoryError> Machine::read_memory(address_t start,
                                                                           std::size_t num_bytes) const noexcept {
  if (num_bytes == 0)
    return std::span<const std::byte>{};

  const auto page_remaining = kPageSize - (start & (kPageSize - 1));
  if (num_bytes > page_remaining)
    return std::unexpected(MemoryError::unmapped_address);

  auto* src = static_cast<const std::byte*>(address_space_->page_virt_to_mapped(start));
  if (!src)
    return std::unexpected(MemoryError::unmapped_address);

  return std::span<const std::byte>(src, num_bytes);
}

std::expected<void, MemoryError> Machine::write_memory(address_t start, std::span<const std::byte> bytes) noexcept {
  std::uint64_t processed = 0;
  while (processed < bytes.size()) {
    address_t addr = start + processed;
    auto* dst = static_cast<std::byte*>(address_space_->page_virt_to_mapped(addr));
    if (!dst)
      return std::unexpected(MemoryError::unmapped_address);

    auto chunk = std::min<std::uint64_t>(bytes.size() - processed, kPageSize - (addr & (kPageSize - 1)));
    std::memcpy(dst, bytes.data() + processed, chunk);
    processed += chunk;
  }
  return {};
}

RunResult Machine::run(RunOptions options) {
  pending_stop_.reset();

  const auto stop_at = options.instruction_limit ? machine_->total_user_insns_committed + *options.instruction_limit
                                                 : std::numeric_limits<W64>::max();
  machine_->config.stop_at_user_insns = stop_at;

  try {
    machine_->run();
  } catch (const X86Exception& e) {
    return {
        .reason = StopReason::x86_exception,
        .stats = stats(),
        .x86_exception = e,
        .message = e.message,
    };
  }

  if (options.instruction_limit && stats().instructions >= stop_at)
    return {.reason = StopReason::instruction_limit, .stats = stats(), .x86_exception = std::nullopt, .message = {}};

  if (pending_stop_) {
    pending_stop_->stats = stats();
    return *pending_stop_;
  }

  return {.reason = StopReason::guest_exit, .stats = stats(), .x86_exception = std::nullopt, .message = {}};
}

Stats Machine::stats() const noexcept {
  return machine_ ? Stats{.cycles = machine_->sim_cycle, .instructions = machine_->total_user_insns_committed} : Stats{};
}

const Options& Machine::options() const noexcept {
  return options_;
}

RegisterFile& Machine::register_file(std::size_t core_index) noexcept {
  return machine_->register_file(core_index);
}

const RegisterFile& Machine::register_file(std::size_t core_index) const noexcept {
  return machine_->register_file(core_index);
}

AddressSpace& Machine::address_space() noexcept {
  return *address_space_;
}

const AddressSpace& Machine::address_space() const noexcept {
  return *address_space_;
}

SyscallResult Machine::dispatch_syscall(Context& context, SyscallKind kind) noexcept {
  return callbacks_.syscall(*this, context, kind);
}

CpuidResult Machine::dispatch_cpuid(Context& context, CpuidRequest request) noexcept {
  return callbacks_.cpuid(*this, context, request);
}

void Machine::set_pending_stop(RunResult result) {
  pending_stop_ = std::move(result);
}

[[noreturn]] void throw_x86_exception(Context& context, byte exception, W32 errorcode, Waddr virtaddr) {
  X86Exception e{
      .vector = exception,
      .error_code = errorcode,
      .virtual_address = virtaddr,
      .rip = context.commitarf[REG_rip],
      .message = std::format("Exception {} code={} addr={} @ rip {}", exception, errorcode,
                             reinterpret_cast<void*>(static_cast<uintptr_t>(virtaddr)),
                             reinterpret_cast<void*>(static_cast<uintptr_t>(context.commitarf[REG_rip]))),
      .context = std::format("{}", context),
  };

  throw e;
}

void dispatch_syscall_64bit(Context& context) {
  assert(context.machine);
  Machine& machine = *context.machine;
  const W64 old_rip = context.commitarf[REG_rip];
  SyscallResult result = machine.dispatch_syscall(context, SyscallKind::syscall64);
  advance_default_syscall_rip_if_unchanged(context, SyscallKind::syscall64, old_rip);

  if (!result.continue_execution) {
    machine.set_pending_stop(
        RunResult{.reason = result.reason, .stats = machine.stats(), .x86_exception = std::nullopt, .message = result.message});
    requested_switch_to_native = 1;
  }
}

void dispatch_syscall_32bit(Context& context, int semantics) {
  assert(context.machine);
  Machine& machine = *context.machine;
  const SyscallKind kind = to_syscall_kind(semantics);
  const W64 old_rip = context.commitarf[REG_rip];
  SyscallResult result = machine.dispatch_syscall(context, kind);
  advance_default_syscall_rip_if_unchanged(context, kind, old_rip);

  if (!result.continue_execution) {
    machine.set_pending_stop(
        RunResult{.reason = result.reason, .stats = machine.stats(), .x86_exception = std::nullopt, .message = result.message});
    requested_switch_to_native = 1;
  }
}

CpuidResult dispatch_cpuid(Context& context, W32 func, W32 subfunc) {
  assert(context.machine);
  x86sim::CpuidResult result = context.machine->dispatch_cpuid(context, {.function = func, .subfunction = subfunc});
  return {.eax = result.eax, .ebx = result.ebx, .ecx = result.ecx, .edx = result.edx};
}

} // namespace x86sim

namespace x86sim {

void Context::propagate_x86_exception(byte exception, W32 errorcode, Waddr virtaddr) {
  throw_x86_exception(*this, exception, errorcode, virtaddr);
}

void handle_syscall_64bit(Context& context) {
  dispatch_syscall_64bit(context);
}

void handle_syscall_32bit(Context& context, int semantics) {
  dispatch_syscall_32bit(context, semantics);
}

CpuidResult handle_cpuid(Context& context, W32 func, W32 subfunc) {
  return dispatch_cpuid(context, func, subfunc);
}

} // namespace x86sim
