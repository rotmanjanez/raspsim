#include "x86sim/x86sim.hpp"

#include "x86sim/addrspace.hpp"
#include "ooocore.h"
#include "ptlhwdef.h"
#include "ptlsim.h"
#include "ptlsim-api.h"
#include "seqcore.h"
#include "stats.h"
#include "x86sim/logging.hpp"

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
  assert(address_space);
  AddressSpace& asp = *address_space;
  logging::println("VMEM: Read from user {} ({})", reinterpret_cast<void*>(addr), bytes);
  logging::flush();

  int n = 0;
  pfec = 0;
  ptelo = 0;
  ptehi = 0;

  bool readable = asp.check(addr, Protection::read);
  bool executable = true;
  if likely (forexec)
    executable = asp.check(addr, Protection::execute);

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

  readable = asp.check(addr + n, Protection::read);
  executable = true;
  if likely (forexec)
    executable = asp.check(addr + n, Protection::execute);

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
  assert(address_space);
  AddressSpace& asp = *address_space;
  logging::println("VMEM: Write to user {} ({})", reinterpret_cast<void*>(target), bytes);
  logging::flush();

  pfec = 0;
  bool writable = asp.check(target, Protection::write);
  if unlikely (!writable) {
    faultaddr = target;
    pfec.p = asp.check(target, Protection::read);
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

  writable = asp.check(target + nlo, Protection::write);
  if unlikely (!writable) {
    faultaddr = target + nlo;
    pfec.p = asp.check(target + nlo, Protection::read);
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
  assert(address_space);
  AddressSpace& asp = *address_space;
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

  if unlikely (!asp.check(virtaddr, store ? Protection::write : Protection::read)) {
    exception = store ? EXCEPTION_PageFaultOnWrite : EXCEPTION_PageFaultOnRead;
    pfec.p = asp.check(virtaddr, Protection::read);
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
  // Userspace always runs with flat segments; the only derived quantity is the
  // virtual address mask, which depends on the active code size.
  virt_addr_mask = use64 ? 0xffffffffffffffffULL : 0xffffffffULL;
}

void Context::load_cpu_state(const CpuState& state) {
  static_cast<CpuState&>(*this) = state;
  // REG_ctx and REG_fpstack are part of the register file and were copied
  // verbatim from the snapshot; repoint them at this context's own storage.
  commitarf[REG_ctx] = reinterpret_cast<Waddr>(this);
  commitarf[REG_fpstack] = reinterpret_cast<Waddr>(fpstack.data());
}

CpuState Context::syscall_return_state(SyscallKind kind) const {
  CpuState state = to_cpu_state();
  switch (kind) {
  case SyscallKind::int80:
  case SyscallKind::sysenter:
    state[Register::rip] = commitarf[REG_nextrip];
    break;
  case SyscallKind::syscall64:
    state[Register::rip] = commitarf[REG_rcx];
    break;
  }
  return state;
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

Machine::Machine(HostCallbacks& callbacks, Options options) : callbacks_(callbacks), options_(std::move(options)) {
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

RunResult Machine::run(CpuState& state, AddressSpace& space, RunOptions options) {
  machine_->pending_stop.reset();

  machine_->state = &state;
  machine_->address_space = &space;

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

  // A pending stop (guest exit, unsupported syscall, host request) takes
  // precedence over the instruction limit: if the guest exits on the very
  // instruction that reaches the limit, report the exit rather than masking it
  // as a plain limit stop. This matters for single-stepping, where every step
  // sets a one-instruction limit.
  if (machine_->pending_stop) {
    machine_->pending_stop->stats = stats();
    return *machine_->pending_stop;
  }

  if (options.instruction_limit && stats().instructions >= stop_at)
    return {.reason = StopReason::instruction_limit, .stats = stats(), .x86_exception = std::nullopt, .message = {}};

  return {.reason = StopReason::guest_exit, .stats = stats(), .x86_exception = std::nullopt, .message = {}};
}

Stats Machine::stats() const noexcept {
  return machine_ ? Stats{.cycles = machine_->sim_cycle, .instructions = machine_->total_user_insns_committed}
                  : Stats{};
}

const Options& Machine::options() const noexcept {
  return options_;
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

// Deliver a guest syscall to the host callback and, if the host asks to stop,
// record the reason on the core so Machine::run can report it.
void dispatch_syscall(Context& context, SyscallKind kind) {
  assert(context.machine_impl);
  assert(context.address_space);
  MachineImpl& impl = *context.machine_impl;
  const W64 old_rip = context.commitarf[REG_rip];
  SyscallResult result = impl.callbacks.syscall(impl.owner, context, *context.address_space, kind);
  advance_default_syscall_rip_if_unchanged(context, kind, old_rip);

  if (!result.continue_execution) {
    impl.pending_stop = RunResult{
        .reason = result.reason, .stats = impl.owner.stats(), .x86_exception = std::nullopt, .message = result.message};
    requested_switch_to_native = 1;
  }
}

} // namespace x86sim

namespace x86sim {

void Context::propagate_x86_exception(byte exception, W32 errorcode, Waddr virtaddr) {
  throw_x86_exception(*this, exception, errorcode, virtaddr);
}

void handle_syscall_64bit(Context& context) {
  dispatch_syscall(context, SyscallKind::syscall64);
}

void handle_syscall_32bit(Context& context, int semantics) {
  dispatch_syscall(context, to_syscall_kind(semantics));
}

CpuidResult handle_cpuid(Context& context, W32 func, W32 subfunc) {
  assert(context.machine_impl);
  assert(context.address_space);
  MachineImpl& impl = *context.machine_impl;
  return impl.callbacks.cpuid(impl.owner, context, *context.address_space, {.function = func, .subfunction = subfunc});
}

} // namespace x86sim
