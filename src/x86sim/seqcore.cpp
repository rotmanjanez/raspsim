//
// PTLsim: Cycle Accurate x86-64 Simulator
// Sequential Core Simulator
//
// Copyright 2003-2008 Matt T. Yourst <yourst@yourst.com>
//

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <format>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <utility>

#include "decode.h"
#include "globals.h"
#include "ptlsim-api.h"
#include "ptlsim.h"
#include "seqcore.h"
#include "branchpred.h"
#include "stats.h"
#include "ptlhwdef.h"
#include "superstl.h"
#include "typedefs.h"
#include "x86sim/addrspace.hpp"
#include "x86sim/logging.hpp"

namespace x86sim {

W64 suppress_total_user_insn_count_updates_in_seqcore;

static const byte archreg_remap_table[TRANSREG_COUNT] = {
    REG_rax,
    REG_rcx,
    REG_rdx,
    REG_rbx,
    REG_rsp,
    REG_rbp,
    REG_rsi,
    REG_rdi,
    REG_r8,
    REG_r9,
    REG_r10,
    REG_r11,
    REG_r12,
    REG_r13,
    REG_r14,
    REG_r15,

    REG_xmml0,
    REG_xmmh0,
    REG_xmml1,
    REG_xmmh1,
    REG_xmml2,
    REG_xmmh2,
    REG_xmml3,
    REG_xmmh3,
    REG_xmml4,
    REG_xmmh4,
    REG_xmml5,
    REG_xmmh5,
    REG_xmml6,
    REG_xmmh6,
    REG_xmml7,
    REG_xmmh7,

    REG_xmml8,
    REG_xmmh8,
    REG_xmml9,
    REG_xmmh9,
    REG_xmml10,
    REG_xmmh10,
    REG_xmml11,
    REG_xmmh11,
    REG_xmml12,
    REG_xmmh12,
    REG_xmml13,
    REG_xmmh13,
    REG_xmml14,
    REG_xmmh14,
    REG_xmml15,
    REG_xmmh15,

    REG_fptos,
    REG_fpsw,
    REG_fptags,
    REG_fpstack,
    REG_msr,
    REG_dlptr,
    REG_trace,
    REG_ctx,

    REG_rip,
    REG_flags,
    REG_dlend,
    REG_selfrip,
    REG_nextrip,
    REG_ar1,
    REG_ar2,
    REG_zero,

    REG_temp0,
    REG_temp1,
    REG_temp2,
    REG_temp3,
    REG_temp4,
    REG_temp5,
    REG_temp6,
    REG_temp7,

    // Notice how these (REG_zf, REG_cf, REG_of) are all mapped to REG_flags in an in-order processor:
    REG_flags,
    REG_flags,
    REG_flags,
    REG_imm,
    REG_mem,
    REG_temp8,
    REG_temp9,
    REG_temp10,
};

const char* seqexec_result_names[SEQEXEC_RESULT_COUNT] = {
    "ok", "early-exit", "smc", "check", "unaligned", "exception", "invalidrip", "skipblock", "barrier", "interrupt",
};

template<int N, int setcount>
struct TransactionalMemory {
  W64 addr_list[N];
  W64 data_list[N];
  W16s next_list[N];
  W8 bytemask_list[N];

  W16s sets[setcount];
  int count;

  // Set by the owning SequentialCore so the (currently unused) commit path can
  // log against that machine's logger rather than a global.
  logging::Logger* logger = nullptr;

  TransactionalMemory() { reset(); }

  void reset() {
    /*
    foreach (i, N) {
      addr_list[i] = 0;
      data_list[i] = 0;
      next_list[i] = -1;
      bytemask_list[i] = 0;
    }
    */

    count = 0;
    memset(sets, 0xff, sizeof(sets));
  }

  static int setof(W64 addr) {
    W64 set = 0;
    addr >>= 3; // cut off subword bits

    return foldbits<log2(setcount)>(addr);
  }

  int lookup(W64 addr, int& set) const {
    set = setof(addr);
    W16s slot = sets[set];
    while (slot >= 0) {
      if (addr_list[slot] == addr)
        return slot;
      slot = next_list[slot];
    }
    return -1;
  }

  int lookup(W64 addr) const {
    int dummy;
    return lookup(addr, dummy);
  }

  bool store(W64 addr, W64 data, byte bytemask) {
    int set;
    int slot = lookup(addr, set);
    if likely (slot < 0) {
      slot = count++;
      assert(slot < N);
      addr_list[slot] = addr;
      data_list[slot] = data;
      bytemask_list[slot] = bytemask;
      next_list[slot] = sets[set];
      sets[set] = slot;
      return true;
    }

    W64& d = data_list[slot];
    d = mux64(expand_8bit_to_64bit_lut[bytemask], d, data);
    bytemask_list[slot] |= bytemask;
    return false;
  }

  bool load(W64 addr, W64& data, byte& bytemask) const {
    bytemask = 0;
    int slot = lookup(addr);
    if likely (slot < 0)
      return false;
    data = data_list[slot];
    bytemask = bytemask_list[slot];
    return true;
  }

  W64 load(W64 addr) const {
    W64 data;
    byte bytemask;

    W64 memdata = loadimpl(addr);
    if likely (!load(addr, data, bytemask))
      return memdata;

    if likely (bytemask == 0xff)
      return data;
    data = mux64(expand_8bit_to_64bit_lut[bytemask], memdata, data);
    return data;
  }

  void rollback() {
    memset(sets, 0xff, sizeof(sets));
    count = 0;
  }

  void commit() {
    foreach (i, count) {
      storeimpl(addr_list[i], data_list[i], bytemask_list[i]);
    }
    memset(sets, 0xff, sizeof(sets));
    count = 0;
  }

  static W64 loadimpl(W64 addr);
  W64 storeimpl(W64 addr, W64 data, byte bytemask);

  /*

  Sample implementations for straight virtual addresses:

  static W64 loadimpl(W64 addr) {
  return *(const W64*)addr;
  }

  static W64 storeimpl(W64 addr, W64 data, byte bytemask) {
  addr = signext64(addr, 48);
  W64& mem = *(W64*)(Waddr)addr;
  mem = mux64(expand_8bit_to_64bit_lut[bytemask], mem, data);
  return mem;
  }

  */

  int update(CommitRecord& cmtrec) {
    foreach (i, count) {
      SFR& sfr = cmtrec.stores[i];
      setzero(sfr);
      sfr.data = data_list[i];
      sfr.physaddr = addr_list[i] >> 3;
      sfr.bytemask = bytemask_list[i];
    }
    cmtrec.store_count = count;
    return count;
  }
};

} // namespace x86sim

template<int N, int setcount>
struct std::formatter<x86sim::TransactionalMemory<N, setcount>> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }
  auto format(const x86sim::TransactionalMemory<N, setcount>& tm, std::format_context& ctx) const;
};


template<int N, int setcount>
auto std::formatter<x86sim::TransactionalMemory<N, setcount>>::format(
    const x86sim::TransactionalMemory<N, setcount>& tm, std::format_context& ctx) const {
  using namespace x86sim;
  auto out = ctx.out();
  out = std::format_to(out, "x86sim::TransactionalMemory containing {} stores:\n", tm.count);
  foreach (i, tm.count) {
    W64 data = tm.data_list[i];
    out = std::format_to(out, "  {:>4}: 0x{:x} <= {}\n", i, W64(tm.addr_list[i]),
                         bytemaskstring((const byte*)&data, tm.bytemask_list[i], 8));
  }
  out = std::format_to(out, "  Hash chains:\n");
  foreach (set, tm.setcount) {
    if (tm.sets[set] < 0)
      continue;
    out = std::format_to(out, "  Set {:>2}:", set);
    W16s slot = tm.sets[set];
    while (slot >= 0) {
      out = std::format_to(out, " {}", slot);
      slot = tm.next_list[slot];
    }
    out = std::format_to(out, "\n");
  }
  return out;
}

namespace x86sim {

template<int N, int setcount>
W64 TransactionalMemory<N, setcount>::loadimpl(W64 physaddr) {
  return loadphys(physaddr);
}

template<int N, int setcount>
W64 TransactionalMemory<N, setcount>::storeimpl(W64 physaddr, W64 data, byte bytemask) {
  if (logger)
    logger->println(logging::DEBUG, "Before storeimpl: physaddr {} => mfn {}, data {}", (void*)physaddr,
                    (physaddr >> 12), bytemaskstring((const byte*)&data, bytemask, 8));
  W64 rc = storemask(physaddr, data, bytemask);
  if (logger)
    logger->println(logging::DEBUG, "After storeimpl: physaddr {} => mfn {}, data {} => rc {}", (void*)physaddr,
                    (physaddr >> 12), bytemaskstring((const byte*)&data, bytemask, 8), (void*)rc);
  return rc;
}

enum {
  EVENT_INVALID = 0,
  EVENT_TRANSLATE,
  EVENT_EXECUTE_BB,
  EVENT_ISSUE,
  EVENT_BRANCH,
  EVENT_LOAD,
  EVENT_STORE,
  EVENT_LOAD_STORE_UNALIGNED,
  EVENT_LOAD_ANNUL,
  EVENT_SKIPBLOCK,
  EVENT_ALIGNMENT_FIXUP,
  EVENT_PTE_UPDATE,
  EVENT_ASSIST,
  EVENT_SMC,
  EVENT_COUNT,
};

//
// Event that gets written to the trace buffer
//
// In the interest of minimizing space, the cycle counters
// and uuids are only 32-bits; in practice wraparound is
// not likely to be a problem.
//
struct SequentialCoreEvent {
  W32 cycle;
  W32 uuid;
  W64 rip;
  W32 eomid;
  byte type;
  byte coreid;
  byte threadid;
  byte uopid;
  TransOp uop;

  union {
    struct {
      IssueState state;
    } issue;
    struct {
      SFR sfr;
      W64 virtaddr;
      W64 origaddr;
      W64 pteused;
      W32 pfec;
    } loadstore;
    struct {
      int uopindex;
    } alignfixup;
    struct {
      W64 chk_recovery_rip;
      byte bytes_in_current_insn;
    } skipblock;
    struct {
      W64 virtaddr;
      W64 pteupdate;
    } pteupdate;
    struct {
      RIPVirtPhysBase rvp;
      void* bb;
      byte bbcount;
    } bb;
    struct {
      W64 rip;
      W64 ptl_pip;
      W64 next_rip;
      W64 real_target_rip;
      W16 id;
    } assist;
  };
};

} // namespace x86sim

template<>
struct std::formatter<x86sim::SequentialCoreEvent> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }
  auto format(const x86sim::SequentialCoreEvent& ev, std::format_context& ctx) const {
    using namespace x86sim;
    auto out = ctx.out();

    if (ev.uuid > 0)
      out = std::format_to(out, "{:>20}", ev.uuid);
    else
      out = std::format_to(out, "{:>20}", "-");
    out = std::format_to(out, " c{} t{} ", ev.coreid, ev.threadid);

    bool st = isstore(ev.uop.opcode);
    bool ld = isload(ev.uop.opcode);
    bool br = isbranch(ev.uop.opcode);
    out = std::format_to(out, "[type {}]", ev.type);

    switch (ev.type) {
    case EVENT_ISSUE:
    case EVENT_BRANCH: {
      std::string rdstr = format_value_and_flags(ev.issue.state.reg.rddata, ev.issue.state.reg.rdflags);
      out = std::format_to(
          out, "{} rip {}:{:<2}  {} {}",
          ((ev.issue.state.reg.rdflags & FLAG_INV) ? ((br) ? "brxcpt" : "except") : ((br) ? "branch" : "issue ")),
          (void*)ev.rip, ev.uopid, ev.uop, rdstr.c_str());
      break;
    }
    case EVENT_LOAD:
    case EVENT_STORE: {
      out = std::format_to(out, "{} rip {}:{:<2}  {} {} (virt 0x{:012x})",
                           ((ev.loadstore.sfr.invalid) ? ((st) ? "stxcpt" : "ldxcpt") : ((st) ? "store " : "load ")),
                           (void*)ev.rip, ev.uopid, ev.uop, ev.loadstore.sfr, ev.loadstore.virtaddr);
      if (ev.loadstore.origaddr != ev.loadstore.virtaddr)
        out = std::format_to(out, " (orig 0x{:012x})", ev.loadstore.origaddr);
      if (ev.loadstore.sfr.invalid)
        out = std::format_to(out, " (PFEC {}, PTE {})", x86sim::PageFaultErrorCode(ev.loadstore.pfec),
                             Level1PTE(ev.loadstore.pteused));
      break;
    }
    case EVENT_LOAD_ANNUL: {
      out = std::format_to(out, "ldanul rip {}:{:<2}  {} {} (virt 0x{:012x})", (void*)ev.rip, ev.uopid, ev.uop,
                           ev.loadstore.sfr, ev.loadstore.virtaddr);
      if (ev.loadstore.origaddr != ev.loadstore.virtaddr)
        out = std::format_to(out, " (orig 0x{:012x})", ev.loadstore.origaddr);
      out = std::format_to(out, " was annulled (high unaligned load)");
      break;
    }
    case EVENT_LOAD_STORE_UNALIGNED: {
      out = std::format_to(out, "{} rip {}:{:<2}  {} virt 0x{:012x} (size {})", ((st) ? "stalgn" : "ldalgn"),
                           (void*)ev.rip, ev.uopid, ev.uop, ev.loadstore.virtaddr, (1 << ev.uop.size));
      break;
    }
    case EVENT_SKIPBLOCK: {
      out = std::format_to(out, "skip   rip {}:{:<2}  {}: advance by {} bytes to {}", (void*)ev.rip, ev.uopid, ev.uop,
                           ev.skipblock.bytes_in_current_insn, (void*)ev.skipblock.chk_recovery_rip);
      break;
    }
    case EVENT_ALIGNMENT_FIXUP: {
      out = std::format_to(out, "algnfx rip {}:{:<2} {} set unaligned bit for uop index {}", (void*)ev.rip, ev.uopid,
                           ev.uop, ev.alignfixup.uopindex);
      break;
    }
    case EVENT_TRANSLATE: {
      out = std::format_to(out, "xlate  rip {} (rvp {}): BB of {} uops", (void*)ev.rip, ev.bb.rvp, ev.bb.bbcount);
      break;
    }
    case EVENT_EXECUTE_BB: {
      out = std::format_to(out, "execbb rip {} (rvp {}): BB of {} uops", (void*)ev.rip, ev.bb.rvp, ev.bb.bbcount);
      break;
    }
    case EVENT_PTE_UPDATE: {
      out = std::format_to(out, "pteupd 0x{:012x}: 0x{:02x}", ev.pteupdate.virtaddr, ev.pteupdate.pteupdate);
      break;
    }
    default: {
      out = std::format_to(out, "Unknown event type {}", ev.type);
      break;
    }
    }

    if (ev.uop.eom)
      out = std::format_to(out, " [EOM #{}]", ev.eomid);
    out = std::format_to(out, "\n");

    return out;
  }
};

namespace x86sim {

struct SequentialCoreEventLog {
  MachineImpl* machine = nullptr;
  SequentialCoreEvent* start;
  SequentialCoreEvent* end;
  SequentialCoreEvent* tail;

  SequentialCoreEventLog() = default;

  explicit SequentialCoreEventLog(MachineImpl& machine_) : machine(&machine_) {
    start = null;
    end = null;
    tail = null;
  }

  bool init(size_t bufsize);
  void reset();

  SequentialCoreEvent* add() {
    if unlikely (tail >= end) {
      assert(machine);
      Options& config = machine->config;
      tail = start;
      if likely ((config.log.loglevel >= 6) || config.log.flush_event_log_every_cycle)
        flush();
    }
    SequentialCoreEvent* event = tail;
    tail++;
    return event;
  }

  void flush(bool only_to_tail = false);

  void clear() { tail = start; }

  SequentialCoreEvent* add(int type, int coreid, const TransOp& uop, W64 rip, int uopid, W64 uuid, W64 eomid) {
    assert(machine);

    SequentialCoreEvent* event = add();
    event->type = type;
    event->cycle = machine->sim_cycle;
    event->uuid = uuid;
    event->rip = rip;
    event->eomid = eomid;
    event->uopid = uopid;
    event->coreid = coreid;
    event->threadid = 0;
    event->uop = uop;
    return event;
  }

  void print(bool only_to_tail = false);
};

bool SequentialCoreEventLog::init(size_t bufsize) {
  reset();
  size_t bytes = bufsize * sizeof(SequentialCoreEvent);
  start = (SequentialCoreEvent*)std::aligned_alloc(PAGE_SIZE, ceil(bytes, PAGE_SIZE));
  if unlikely (!start)
    return false;
  end = start + bufsize;
  tail = start;

  foreach (i, bufsize)
    start[i].type = EVENT_INVALID;
  return true;
}

void SequentialCoreEventLog::reset() {
  if (!start)
    return;

  std::free(start);
  start = null;
  end = null;
  tail = null;
}

void SequentialCoreEventLog::flush(bool only_to_tail) {
  print(only_to_tail);
  tail = start;
}

void SequentialCoreEventLog::print(bool only_to_tail) {
  if (!start)
    return;

  assert(machine);
  Options& config = machine->config;

  if (tail > end)
    tail = start;
  if (tail < start)
    tail = end;

  SequentialCoreEvent* p = (only_to_tail) ? start : tail;
  size_t bufsize = end - start;
  size_t limit = (only_to_tail) ? (tail - start) : bufsize;

  if (!config.log.flush_event_log_every_cycle)
    machine->logger.println(logging::INFO, "#-------- Start of event log --------");

  W64 cycle = std::numeric_limits<W64>::max();
  foreach (i, limit) {
    if unlikely (p >= end)
      p = start;
    if unlikely (p < start)
      p = end - 1;

    if unlikely (p->type == EVENT_INVALID) {
      p++;
      continue;
    }

    if (p->type == EVENT_EXECUTE_BB) {
      foreach (i, 24)
        machine->logger.print(logging::INFO, "--------");
      machine->logger.println(logging::INFO, "");
    }

    if unlikely (p->cycle != cycle) {
      cycle = p->cycle;
      machine->logger.println(logging::INFO, "Cycle {}:", cycle);
    }

    machine->logger.print(logging::INFO, "{}", *p);
    p++;
  }

  if (!config.log.flush_event_log_every_cycle)
    machine->logger.println(logging::INFO, "#-------- End of event log --------");
}

static SequentialCoreEventLog eventlog;

struct SequentialCore {
  SequentialMachine& machine;
  Options& config;
  Context ctx;
  CommitRecord* cmtrec;
  bool runnable = true;

  SequentialCore(SequentialMachine& machine_, int core_id, CommitRecord* cmtrec_ = null)
      : machine(machine_), config(machine_.config), ctx(config, machine_, core_id), cmtrec(cmtrec_) {
    transactmem.logger = &machine.logger;
  }

  BasicBlock* current_basic_block;
  int bytes_in_current_insn;
  int current_uop_in_macro_op;
  W64 current_uuid;

  // (n/a):
  W64 fetch_blocks_fetched;
  W64 fetch_uops_fetched;
  W64 fetch_user_insns_fetched;

  W64 bbcache_inserts;
  W64 bbcache_removes;

  W64 seq_total_basic_blocks;
  W64 seq_total_uops_committed;
  W64 seq_total_user_insns_committed;
  W64 seq_total_cycles;

  //
  // Shadow flags are maintained for each archreg to simulate renaming,
  // since the x86 decoder assumes renaming will be done and hence may
  // specify some uops as "don't update user flags".
  //
  W64 arf[TRANSREG_COUNT];
  W16 arflags[TRANSREG_COUNT];

  TransactionalMemory<MAX_STORES_IN_COMMIT_RECORD, 16> transactmem;

  void print_state() {
    machine.logger.println("General state:");
    machine.logger.println("  RIP:                {}", (void*)(Waddr)arf[REG_rip]);
    machine.logger.println("  Flags:              {:04x} {}", arf[REG_flags], flagstring(arf[REG_flags]));
    machine.logger.println("  UUID:               {}", current_uuid);
    machine.logger.println("  Bytes in macro-op:  {}", bytes_in_current_insn);
    machine.logger.println("  Uop in macro-op:    {}", current_uop_in_macro_op);
    machine.logger.println("Basic block state:");
    machine.logger.println("  BBcache block:      {}", (void*)current_basic_block);
    machine.logger.println("  uop count in block: {}", ((current_basic_block) ? current_basic_block->count : 0));
    machine.logger.println("Register state:");

    static const int width = 4;
    foreach (i, TRANSREG_COUNT) {
      std::string flagsb = std::format("{}", flagstring(arflags[i]));
      machine.logger.print("    {:<6} 0x{:016x}|{:<6}  ", arch_reg_names[i], arf[i], flagsb);
      if ((i % width) == (width - 1))
        machine.logger.println("");
    }
  }

  void reset_fetch(W64 realrip) {
    arf[REG_rip] = realrip;
    current_basic_block = null;
  }

  enum {
    ISSUE_COMPLETED = 1,
    ISSUE_REFETCH = 0,
    ISSUE_EXCEPTION = -1,
  };

  //
  // Address generation common to both loads and stores
  //
  template<int STORE>
  Waddr addrgen(const TransOp& uop, SFR& state, Waddr& origaddr, W64 ra, W64 rb, W64 rc, PTEUpdate& pteupdate,
                Waddr& addr, int& exception, PageFaultErrorCode& pfec, Level1PTE& pteused, bool& annul) {
    int sizeshift = uop.size;
    int aligntype = uop.cond;
    bool internal = uop.internal;
    bool signext = (uop.opcode == OP_ldx);

    Waddr rip = arf[REG_rip];

    addr = (STORE) ? (ra + rb) : ((aligntype == LDST_ALIGN_NORMAL) ? (ra + rb) : ra);
    //
    // x86-64 requires virtual addresses to be canonical: if bit 47 is set,
    // all upper 16 bits must be set. If this is not true, we need to signal
    // a general protection fault.
    //
    addr = (W64)signext64(addr, 48);
    addr &= ctx.virt_addr_mask;
    origaddr = addr;
    annul = 0;

    switch (aligntype) {
    case LDST_ALIGN_NORMAL:
      break;
    case LDST_ALIGN_LO:
      addr = floor(addr, 8);
      break;
    case LDST_ALIGN_HI:
      //
      // Is the high load ever even used? If not, don't check for exceptions;
      // otherwise we may erroneously flag page boundary conditions as invalid
      //
      addr = floor(addr, 8);
      annul = (floor(origaddr + ((1 << sizeshift) - 1), 8) == addr);
      addr += 8;
      break;
    }

    state.physaddr = addr >> 3;
    state.invalid = 0;
    state.addrvalid = 1;
    state.datavalid = 0;

    //
    // Special case: if no part of the actual user load/store falls inside
    // of the high 64 bits, do not perform the access and do not signal
    // any exceptions if that page was invalid.
    //
    // However, we must be extremely careful if we're inheriting an SFR
    // from an earlier store: the earlier store may have updated some
    // bytes in the high 64-bit chunk even though we're not updating
    // any bytes. In this case we still must do the write since it
    // could very well be the final commit to that address. In any
    // case, the SFR mismatch and LSAT must still be checked.
    //
    // The store commit code checks if the bytemask is zero and does
    // not attempt the actual store if so. This will always be correct
    // for high stores as described in this scenario.
    //

    exception = 0;

    W64 physaddr =
        (annul) ? INVALID_PHYSADDR
                : ctx.check_and_translate(addr, uop.size, STORE, uop.internal, exception, pfec, pteupdate, pteused);
    return physaddr;
  }

  //
  // Handle exceptions common to both loads and stores
  //
  template<bool STORE>
  int handle_common_exceptions(const TransOp& uop, SFR& state, Waddr& origaddr, Waddr& addr, int& exception,
                               PageFaultErrorCode& pfec, Level1PTE& pteused) {
    if likely (!exception)
      return ISSUE_COMPLETED;

    int aligntype = uop.cond;

    state.invalid = 1;
    state.data = exception | ((W64)pfec << 32);
    state.datavalid = 1;

    if likely (exception == EXCEPTION_UnalignedAccess) {
      //
      // If we have an unaligned access, retry this dynamic load/store as
      // a split access. Alignment depends on the current register values,
      // so the sequential core must not cache this decision in the basic block.
      //
      // The real hardware would detect unaligned uops in the fetch stage
      // and split them up on the fly.
      //
      if unlikely (config.log.event_log_enabled) {
        SequentialCoreEvent* event =
            eventlog.add(EVENT_LOAD_STORE_UNALIGNED, ctx.vcpuid, uop, arf[REG_rip], current_uop_in_macro_op,
                         current_uuid, machine.total_user_insns_committed);
        event->loadstore.virtaddr = origaddr;
      }

      return ISSUE_REFETCH;
    }

    if unlikely (((exception == EXCEPTION_PageFaultOnRead) | (exception == EXCEPTION_PageFaultOnWrite)) &&
                 (aligntype == LDST_ALIGN_HI)) {
      //
      // If we have a page fault on an unaligned access, and this is the high
      // half (ld.hi / st.hi) of that access, the page fault address recorded
      // in CR2 must be at the very first byte of the second page the access
      // overlapped onto (otherwise the kernel will repeatedly fault in the
      // first page, even though that one is already present.
      //
      origaddr = addr;
    }

    if unlikely (config.log.event_log_enabled) {
      SequentialCoreEvent* event =
          eventlog.add((STORE) ? EVENT_STORE : EVENT_LOAD, ctx.vcpuid, uop, arf[REG_rip], current_uop_in_macro_op,
                       current_uuid, machine.total_user_insns_committed);
      event->loadstore.sfr = state;
      event->loadstore.virtaddr = addr;
      event->loadstore.origaddr = origaddr;
      event->loadstore.pfec = pfec;
      event->loadstore.pteused = pteused;
    }

    return ISSUE_EXCEPTION;
  }

  int issuestore(const TransOp& uop, SFR& state, Waddr& origaddr, W64 ra, W64 rb, W64 rc, PTEUpdate& pteupdate) {
    int status;
    Waddr rip = arf[REG_rip];
    int sizeshift = uop.size;
    int aligntype = uop.cond;

    Waddr addr;
    int exception = 0;
    PageFaultErrorCode pfec;
    Level1PTE pteused;
    bool annul;

    W64 physaddr = addrgen<1>(uop, state, origaddr, ra, rb, rc, pteupdate, addr, exception, pfec, pteused, annul);

    if unlikely ((status = handle_common_exceptions<1>(uop, state, origaddr, addr, exception, pfec, pteused)) !=
                 ISSUE_COMPLETED)
      return status;

    //
    // At this point all operands are valid, so merge the data and mark the store as valid.
    //
    state.physaddr = (annul) ? INVALID_PHYSADDR : (physaddr >> 3);

    bool ready;
    byte bytemask;

    switch (aligntype) {
    case LDST_ALIGN_NORMAL:
    case LDST_ALIGN_LO:
      bytemask = ((1 << (1 << sizeshift)) - 1) << (lowbits(origaddr, 3));
      rc <<= 8 * lowbits(origaddr, 3);
      break;
    case LDST_ALIGN_HI:
      bytemask = ((1 << (1 << sizeshift)) - 1) >> (8 - lowbits(origaddr, 3));
      rc >>= 8 * (8 - lowbits(origaddr, 3));
    }

    state.invalid = 0;
    state.data = rc;
    state.bytemask = bytemask;
    state.datavalid = !annul;

    if unlikely (config.log.event_log_enabled) {
      SequentialCoreEvent* event = eventlog.add(EVENT_STORE, ctx.vcpuid, uop, rip, current_uop_in_macro_op,
                                                current_uuid, machine.total_user_insns_committed);
      event->loadstore.sfr = state;
      event->loadstore.virtaddr = addr;
      event->loadstore.origaddr = origaddr;
      event->loadstore.pfec = pfec;
      event->loadstore.pteused = pteused;
    }

    return ISSUE_COMPLETED;
  }

  int issueload(const TransOp& uop, SFR& state, Waddr& origaddr, W64 ra, W64 rb, W64 rc, PTEUpdate& pteupdate) {
    int status;
    Waddr rip = arf[REG_rip];
    int sizeshift = uop.size;
    int aligntype = uop.cond;
    bool signext = (uop.opcode == OP_ldx);
    Waddr addr;
    int exception = 0;
    PageFaultErrorCode pfec;
    Level1PTE pteused;
    bool annul;

    W64 physaddr = addrgen<0>(uop, state, origaddr, ra, rb, rc, pteupdate, addr, exception, pfec, pteused, annul);

    if unlikely ((status = handle_common_exceptions<0>(uop, state, origaddr, addr, exception, pfec, pteused)) !=
                 ISSUE_COMPLETED)
      return status;

    state.physaddr = (annul) ? 0xffffffffffffffffULL : (physaddr >> 3);

    W64 data = 0;
    if likely (!annul) {
      if unlikely (cmtrec) {
        data = transactmem.load(state.physaddr << 3);
      } else {
        machine.logger.println("[cycle {}] load from physaddr {} for virtaddr {}", machine.sim_cycle, (void*)physaddr,
                               (void*)origaddr);
        data = loadphys(physaddr);
      }
    }

    if unlikely (aligntype == LDST_ALIGN_HI) {
      //
      // Concatenate the aligned data from a previous ld.lo uop provided in rb
      // with the currently loaded data D as follows:
      //
      // rb | D
      //
      // Example:
      //
      // floor(a) floor(a)+8
      // ---rb--  --DD---
      // 0123456701234567
      //    XXXXXXXX
      //    ^ origaddr
      //
      if likely (!annul) {
        struct {
          W64 lo;
          W64 hi;
        } aligner;

        aligner.lo = rb;
        aligner.hi = data;

        W64 offset = lowbits(origaddr - floor(origaddr, 8), 4);

        data = extract_bytes(((byte*)&aligner) + offset, sizeshift, signext);
      } else {
        //
        // annulled: we need no data from the high load anyway; only use the low data
        // that was already checked for exceptions and forwarding:
        //
        W64 offset = lowbits(origaddr, 3);
        state.data = extract_bytes(((byte*)&rb) + offset, sizeshift, signext);
        state.invalid = 0;
        state.datavalid = 1;

        if unlikely (config.log.event_log_enabled) {
          SequentialCoreEvent* event = eventlog.add(EVENT_LOAD_ANNUL, ctx.vcpuid, uop, rip, current_uop_in_macro_op,
                                                    current_uuid, machine.total_user_insns_committed);
          event->loadstore.sfr = state;
          event->loadstore.virtaddr = addr;
          event->loadstore.origaddr = origaddr;
          event->loadstore.pfec = 0;
        }

        return ISSUE_COMPLETED;
      }
    } else {
      data = extract_bytes(((byte*)&data) + lowbits(addr, 3), sizeshift, signext);
    }

    //
    // NOTE: Technically the data is valid right now for simulation purposes
    // only; in reality it may still be arriving from the cache.
    //
    state.data = data;
    state.invalid = 0;
    state.datavalid = 1;
    state.bytemask = 0xff;

    if unlikely (config.log.event_log_enabled) {
      SequentialCoreEvent* event = eventlog.add(EVENT_LOAD, ctx.vcpuid, uop, rip, current_uop_in_macro_op, current_uuid,
                                                machine.total_user_insns_committed);
      event->loadstore.sfr = state;
      event->loadstore.virtaddr = addr;
      event->loadstore.origaddr = origaddr;
      event->loadstore.pfec = pfec;
      event->loadstore.pteused = pteused;
    }

    return ISSUE_COMPLETED;
  }

  void external_to_core_state(const Context& ctx) {
    foreach (i, ARCHREG_COUNT) {
      arf[i] = ctx.commitarf[i];
      arflags[i] = 0;
    }
    for (int i = ARCHREG_COUNT; i < TRANSREG_COUNT; i++) {
      arf[i] = 0;
      arflags[i] = 0;
    }

    arflags[REG_flags] = ctx.commitarf[REG_flags];
  }

  void core_to_external_state(Context& ctx) {
    foreach (i, ARCHREG_COUNT) {
      ctx.commitarf[i] = arf[i];
    }
  }

  bool handle_barrier() {
    core_to_external_state(ctx);

    int assistid = ctx.commitarf[REG_rip];
    assist_func_t assist = (assist_func_t)(Waddr)assistid_to_func[assistid];

    machine.logger.println(
        logging::DEBUG, "[vcpu {}] Barrier (#{} -> {} {}) called from {}; return to {}) at {} cycles, {} commits",
        ctx.vcpuid, assistid, (void*)assist, assist_name(assist), (RIPVirtPhys(ctx.commitarf[REG_selfrip]).update(ctx)),
        (void*)(Waddr)ctx.commitarf[REG_nextrip], machine.sim_cycle, machine.total_user_insns_committed);
    machine.logger.flush();

    machine.logger.println(logging::TRACE, "Calling assist function at {}...", (void*)assist);
    machine.logger.flush();

    update_assist_stats(assist);
    machine.logger.println(logging::TRACE, "Before assist:");
    machine.logger.println(logging::TRACE, "{}", ctx);

    requested_switch_to_native = false;
    assist(ctx);
    const bool switch_to_native = requested_switch_to_native;
    requested_switch_to_native = false;

    machine.logger.println(logging::TRACE, "Done with assist");
    machine.logger.println(logging::TRACE, "New state:");
    machine.logger.println(logging::TRACE, "{}", ctx);

    reset_fetch(ctx.commitarf[REG_rip]);
    external_to_core_state(ctx);
    if (switch_to_native) {
      machine.logger.println("PTL call requested switch to native mode at rip {}",
                             (void*)(Waddr)ctx.commitarf[REG_rip]);
      return false;
    }
    return true;
  }

  bool handle_exception() {
    core_to_external_state(ctx);

    machine.logger.println(logging::INFO, "PTL Exception {} called from rip {} at {} cycles, {} commits", ctx.exception,
                           (void*)(Waddr)ctx.commitarf[REG_rip], machine.sim_cycle, machine.total_user_insns_committed);
    machine.logger.flush();

    //
    // Map PTL internal hardware exceptions to their x86 equivalents,
    // depending on the context. The error_code field should already
    // be filled out.
    //
    switch (ctx.exception) {
    case EXCEPTION_PageFaultOnRead:
    case EXCEPTION_PageFaultOnWrite:
    case EXCEPTION_PageFaultOnExec:
      ctx.x86_exception = EXCEPTION_x86_page_fault;
      break;
    case EXCEPTION_FloatingPointNotAvailable:
      ctx.x86_exception = EXCEPTION_x86_fpu_not_avail;
      break;
    case EXCEPTION_FloatingPoint:
      ctx.x86_exception = EXCEPTION_x86_fpu;
      break;
    default:
      machine.logger.println("Unsupported internal exception type {}", ctx.exception);
      machine.logger.flush();
      assert(false);
    }

    if unlikely ((ctx.x86_exception == EXCEPTION_x86_page_fault) && (ctx.cr2 == 0xffffffff00000018)) {
      eventlog.print();
    }
    machine.logger.println(logging::DEBUG, "{}", ctx);

    ctx.propagate_x86_exception(ctx.x86_exception, ctx.error_code, ctx.cr2);

    external_to_core_state(ctx);

    return true;
  }

  BasicBlock* fetch_or_translate_basic_block(Waddr rip) {
    RIPVirtPhys rvp(rip);

    rvp.update(ctx);

    BasicBlockCache& bbcache = *machine.bbcache;
    BasicBlock* bb = bbcache(rvp);

    if likely (bb) {
      current_basic_block = bb;
    } else {
      current_basic_block = bbcache.translate(ctx, rvp);
      assert(current_basic_block);
      synth_uops_for_bb(*current_basic_block);

      if unlikely (config.log.event_log_enabled) {
        TransOp dummyuop;
        setzero(dummyuop);
        SequentialCoreEvent* event =
            eventlog.add(EVENT_TRANSLATE, ctx.vcpuid, dummyuop, rip, 0, 0, machine.total_user_insns_committed);
        event->bb.rvp = rvp;
        event->bb.bb = current_basic_block;
        event->bb.bbcount = current_basic_block->count;
      }

      bbcache_inserts++;
    }

    current_basic_block->use(machine.sim_cycle);
    return current_basic_block;
  }

  //
  // Execute one basic block sequentially
  //

  int execute(BasicBlock* bb, W64 insnlimit) {
    arf[REG_rip] = bb->rip;

    //
    // Fetch
    //

    bool barrier = 0;

    machine.logger.println(logging::DEBUG, "[vcpu {}] Sequentially executing basic block {} ({} uops), insn limit {}",
                           ctx.vcpuid, W64(bb->rip), bb->count, insnlimit);

    if unlikely (config.log.event_log_enabled) {
      TransOp dummyuop;
      setzero(dummyuop);
      SequentialCoreEvent* event =
          eventlog.add(EVENT_EXECUTE_BB, ctx.vcpuid, dummyuop, bb->rip, 0, 0, machine.total_user_insns_committed);
      event->bb.rvp = bb->rip;
      event->bb.bb = bb;
      event->bb.bbcount = bb->count;
    }

    if unlikely (!bb->synthops)
      synth_uops_for_bb(*bb);
    bb->hitcount++;

    TransOpBuffer unaligned_ldst_buf;
    unaligned_ldst_buf.index = -1;

    int uopindex = 0;
    int current_uop_in_macro_op = 0;

    int user_insns = 0;

    seq_total_basic_blocks++;
    machine.total_basic_blocks_committed++;

    RIPVirtPhys rvp(arf[REG_rip]);

    assert(bb->rip.rip == arf[REG_rip]);

    // See comment below about idempotent updates
    W64 saved_flags = 0;

    while ((uopindex < bb->count) & (user_insns < insnlimit)) {
      TransOp uop;
      UopImpl synthop;

      if unlikely (arf[REG_rip] == config.stop_at_rip) {
        return SEQEXEC_EARLY_EXIT;
      }

      if unlikely (arf[REG_rip] == config.log.start_log_at_rip) {
        config.log.start_log_at_iteration = 0;
        machine.logger.set_enabled(true);
      }

      if likely (!unaligned_ldst_buf.get(uop, synthop)) {
        uop = bb->transops[uopindex];
        synthop = bb->synthops[uopindex];
      }

      assert(uopindex < bb->count);

      if likely (uop.som) {
        current_uop_in_macro_op = 0;
        bytes_in_current_insn = uop.bytes;
        fetch_user_insns_fetched++;
        // Update the span of bytes to watch for SMC:
        rvp.update(ctx, uop.bytes);
        //
        // Save the flags at the start of this x86 insn in
        // case an ALU uop inside the macro-op updates the
        // flags before all exceptions (i.e. from stores)
        // can be detected. All other registers are updated
        // idempotently.
        //
        saved_flags = arf[REG_flags];
      }

      //
      // Check for self modifying code (SMC) by checking if any previous
      // instruction has dirtied the page(s) on which the current instruction
      // resides. The SMC check is done first since it's perfectly legal for a
      // store to overwrite its own instruction bytes, but this update only
      // becomes visible after the store has committed.
      //
      AddressSpace& asp = *ctx.address_space;
      if unlikely (asp.isdirty(rvp.mfnlo) | (asp.isdirty(rvp.mfnhi))) {
        machine.logger.println("Self-modifying code at rip {} detected: mfn was dirty (invalidate and retry)", rvp);
        machine.bbcache->invalidate_page(asp, rvp.mfnlo, INVALIDATE_REASON_SMC);
        if (rvp.mfnlo != rvp.mfnhi)
          machine.bbcache->invalidate_page(asp, rvp.mfnhi, INVALIDATE_REASON_SMC);
        return SEQEXEC_SMC;
      }

      fetch_uops_fetched++;

      //
      // Issue
      //
      IssueState state;
      state.reg.rddata = 0;
      state.reg.addr = 0;
      state.reg.rdflags = 0;
      ctx.exception = 0;

      W64 radata = arf[archreg_remap_table[uop.ra]];
      W64 rbdata = (uop.rb == REG_imm) ? uop.rbimm : arf[archreg_remap_table[uop.rb]];
      W64 rcdata = (uop.rc == REG_imm) ? uop.rcimm : arf[archreg_remap_table[uop.rc]];

      W16 raflags = arflags[archreg_remap_table[uop.ra]];
      W16 rbflags = arflags[archreg_remap_table[uop.rb]];
      W16 rcflags = arflags[archreg_remap_table[uop.rc]];

      const UopInputs inputs{.ra = radata,
                             .rb = rbdata,
                             .rc = rcdata,
                             .raflags = raflags,
                             .rbflags = rbflags,
                             .rcflags = rcflags,
                             .riptaken = uop.riptaken,
                             .ripseq = uop.ripseq,
                             .logger = &machine.logger};

      bool ld = isload(uop.opcode);
      bool st = isstore(uop.opcode);
      bool br = isbranch(uop.opcode);

      SFR sfr;

      bool refetch = 0;

      PTEUpdate pteupdate = 0;
      Waddr origvirt = 0;
      PageFaultErrorCode pfec = 0;

      bool force_fpu_not_avail_fault = 0;

      Waddr rip = arf[REG_rip];

      force_fpu_not_avail_fault = (uop.is_sse & ctx.no_sse) | (uop.is_x87 & ctx.no_x87);
      if unlikely (force_fpu_not_avail_fault) {
        if unlikely (config.log.event_log_enabled) {
          SequentialCoreEvent* event = eventlog.add(EVENT_ISSUE, ctx.vcpuid, uop, rip, current_uop_in_macro_op,
                                                    current_uuid, machine.total_user_insns_committed);
          IssueState state;
          state.reg.rdflags = FLAG_INV;
          state.reg.rddata = EXCEPTION_FloatingPointNotAvailable;
        }

        ctx.exception = EXCEPTION_FloatingPointNotAvailable;
        ctx.error_code = 0;
        arf[REG_flags] = saved_flags;
        return SEQEXEC_EXCEPTION;
      } else if unlikely (ld | st) {
        int status;
        if likely (ld) {
          status = issueload(uop, sfr, origvirt, radata, rbdata, rcdata, pteupdate);
        } else if unlikely (uop.opcode == OP_mf) {
          // Memory fences are NOPs on the in-order core:
          status = ISSUE_COMPLETED;
          sfr.data = 0;
        } else {
          assert(st);
          status = issuestore(uop, sfr, origvirt, radata, rbdata, rcdata, pteupdate);
        }

        state.reg.rddata = sfr.data;
        state.reg.rdflags = 0;

        if (status == ISSUE_EXCEPTION) {
          ctx.exception = LO32(state.reg.rddata);
          ctx.error_code = HI32(state.reg.rddata); // page fault error code
          ctx.cr2 = origvirt;
          arf[REG_flags] = saved_flags;
          return SEQEXEC_EXCEPTION;
        } else if (status == ISSUE_REFETCH) {
          if unlikely (config.log.event_log_enabled) {
            SequentialCoreEvent* event =
                eventlog.add(EVENT_ALIGNMENT_FIXUP, ctx.vcpuid, uop, rip, current_uop_in_macro_op, current_uuid,
                             machine.total_user_insns_committed);
            event->alignfixup.uopindex = uopindex;
          }
          uop.unaligned = 1;
          machine.logger.println(logging::INFO, "{:20} fetch  rip {}: split unaligned load or store {}", "",
                                 (void*)(Waddr)arf[REG_rip], uop);
          split_unaligned(uop, unaligned_ldst_buf);
          continue;
        }
      } else if unlikely (br) {
        assert(synthop);
        state = issue_state_from(synthop(inputs));

        if unlikely (config.log.event_log_enabled) {
          SequentialCoreEvent* event = eventlog.add(EVENT_BRANCH, ctx.vcpuid, uop, rip, current_uop_in_macro_op,
                                                    current_uuid, machine.total_user_insns_committed);
          event->issue.state = state;
        }

        bb->predcount +=
            (uop.opcode == OP_jmp) ? (state.reg.rddata == bb->lasttarget) : (state.reg.rddata == uop.riptaken);
        bb->lasttarget = state.reg.rddata;
      } else {
        assert(synthop);
        state = issue_state_from(synthop(inputs));
        if unlikely (state.reg.rdflags & FLAG_INV)
          ctx.exception = LO32(state.reg.rddata);

        if unlikely (config.log.event_log_enabled) {
          SequentialCoreEvent* event = eventlog.add(EVENT_ISSUE, ctx.vcpuid, uop, rip, current_uop_in_macro_op,
                                                    current_uuid, machine.total_user_insns_committed);
          event->issue.state = state;
        }

        if unlikely (ctx.exception) {
          if (isclass(uop.opcode, OPCLASS_CHECK) & (ctx.exception == EXCEPTION_SkipBlock)) {
            W64 chk_recovery_rip = arf[REG_rip] + bytes_in_current_insn;
            if unlikely (config.log.event_log_enabled) {
              SequentialCoreEvent* event = eventlog.add(EVENT_SKIPBLOCK, ctx.vcpuid, uop, rip, current_uop_in_macro_op,
                                                        current_uuid, machine.total_user_insns_committed);
              event->skipblock.bytes_in_current_insn = bytes_in_current_insn;
              event->skipblock.chk_recovery_rip = chk_recovery_rip;
            }
            current_uuid++;
            arf[REG_rip] = chk_recovery_rip;

            seq_total_user_insns_committed++;
            machine.total_user_insns_committed += (!suppress_total_user_insn_count_updates_in_seqcore);
            user_insns++;
            return SEQEXEC_OK;
          } else {
            arf[REG_flags] = saved_flags;
            return SEQEXEC_EXCEPTION;
          }
        }
      }

      //
      // Commit
      //

      machine.total_uops_committed++;
      seq_total_uops_committed++;

      assert(!ctx.exception);

      if unlikely (uop.opcode == OP_st || uop.opcode == OP_st_a16) {
        if (sfr.bytemask) {
          if unlikely (cmtrec) {
            transactmem.store(sfr.physaddr << 3, sfr.data, sfr.bytemask);
          } else {
            storemask(sfr.physaddr << 3, sfr.data, sfr.bytemask);
          }

          /* see ooopipe.cpp for a short comment why we can't use sfr->physaddr */
          asp.setdirty(sfr.smc_mfn);
        }
      } else if likely (uop.rd != REG_zero) {
        arf[uop.rd] = state.reg.rddata;
        arflags[uop.rd] = state.reg.rdflags;

        if (!uop.nouserflags) {
          W64 flagmask = setflags_to_x86_flags[uop.setflags];
          arf[REG_flags] = (arf[REG_flags] & ~flagmask) | (state.reg.rdflags & flagmask);
          arflags[REG_flags] = arf[REG_flags];
        }
      }

      if unlikely (pteupdate) {
        if unlikely (config.log.event_log_enabled) {
          SequentialCoreEvent* event = eventlog.add(EVENT_PTE_UPDATE, ctx.vcpuid, uop, rip, current_uop_in_macro_op,
                                                    current_uuid, machine.total_user_insns_committed);
          event->pteupdate.virtaddr = origvirt;
          event->pteupdate.pteupdate = pteupdate;
        }

        if unlikely (cmtrec) {
          /*
          assert(cmtrec->pte_update_count < lengthof(cmtrec->pte_update_list));
          cmtrec->pte_update_list[cmtrec->pte_update_count] = pteupdate;
          cmtrec->pte_update_virt[cmtrec->pte_update_count] = origvirt;
          cmtrec->pte_update_count++;
          */
        } else {
          ctx.update_pte_acc_dirty(origvirt, pteupdate);
        }
      }

      barrier = isclass(uop.opcode, OPCLASS_BARRIER);

      if unlikely ((arf[REG_rip] == config.log.log_backwards_from_trigger_rip) && (config.log.event_log_enabled)) {
        machine.logger.println("Hit trigger rip {}; printing event ring buffer:",
                               (void*)(Waddr)config.log.log_backwards_from_trigger_rip);
        machine.logger.flush();
        eventlog.print();
        machine.logger.println("End of triggered event dump");
      }

      if likely (uop.eom) {
        arf[REG_rip] = (uop.rd == REG_rip) ? state.reg.rddata : (arf[REG_rip] + bytes_in_current_insn);
        // Do not commit transactional memory: that's up to the caller:
        // if unlikely (cmtrec) transactmem.commit();
      }

      if unlikely (barrier) {
        if unlikely (config.log.event_log_enabled) {
          int assistid = arf[REG_rip];
          assist_func_t assist = (assist_func_t)(Waddr)assistid_to_func[assistid];
          TransOp dummyuop;
          setzero(dummyuop);
          SequentialCoreEvent* event = eventlog.add(EVENT_ASSIST, ctx.vcpuid, dummyuop, rip, current_uop_in_macro_op,
                                                    current_uuid, machine.total_user_insns_committed);
          event->assist.id = assistid;
          event->assist.rip = arf[REG_selfrip];
          event->assist.ptl_pip = (W64)assist;
          event->assist.next_rip = arf[REG_nextrip];
          event->assist.real_target_rip = arf[REG_rip];
        }
      }

      seq_total_user_insns_committed += uop.eom;
      machine.total_user_insns_committed += uop.eom && (!suppress_total_user_insn_count_updates_in_seqcore);
      user_insns += uop.eom;
      stats.summary.insns += uop.eom;
      stats.summary.uops++;

      current_uuid++;
      // Don't advance on cracked loads/stores:
      uopindex += unaligned_ldst_buf.empty();
      current_uop_in_macro_op++;

      if unlikely (br) {
        core_to_external_state(ctx);
        machine.logger.println(logging::VERBOSE, "Core State after branch:");
        machine.logger.println(logging::VERBOSE, "{}", ctx);
      }
    }

    if (barrier)
      return SEQEXEC_BARRIER;

    return (insnlimit < bb->user_insn_count) ? SEQEXEC_EARLY_EXIT : SEQEXEC_OK;
  }

  int execute() {
    Waddr rip = arf[REG_rip];

    current_basic_block = fetch_or_translate_basic_block(rip);

    bool exiting = 0;

    int result = execute(current_basic_block, (config.stop_at_user_insns - machine.total_user_insns_committed));

    switch (result) {
    case SEQEXEC_OK:
    case SEQEXEC_SMC:
    case SEQEXEC_SKIPBLOCK:
      // no action required
      break;
    case SEQEXEC_EARLY_EXIT:
      exiting = 1;
      break;
    case SEQEXEC_EXCEPTION:
    case SEQEXEC_INVALIDRIP:
      exiting = (!handle_exception());
      break;
    case SEQEXEC_BARRIER:
      exiting = (!handle_barrier());
      break;
    default:
      assert(false);
    }

    return exiting;
  }
};

std::string_view SequentialMachine::name() const {
  return "seq";
}

SequentialMachine::~SequentialMachine() = default;

SequentialMachine::SequentialMachine(Machine& machine, Options config)
    : MachineImpl(machine, std::move(config)), core(std::make_unique<SequentialCore>(*this, 0)) {}

Context& SequentialMachine::cpu_context() {
  return core->ctx;
}

const Context& SequentialMachine::cpu_context() const {
  return core->ctx;
}

void SequentialMachine::flush_tlb(Context& ctx) {}

void SequentialMachine::flush_tlb_virt(Context& ctx, Waddr virtaddr) {}

int SequentialMachine::run() {
  eventlog.machine = this;

  logger.println("Starting sequential core toplevel loop at {} cycles and {} commits", sim_cycle,
                 total_user_insns_committed);
  logger.flush();

  if unlikely (config.log.event_log_enabled && (!eventlog.start)) {
    eventlog.init(config.log.event_log_ring_buffer_size);
  }

  SequentialCore& seq = *core;
  Context& ctx = seq.ctx;
  ctx.machine = &owner;
  ctx.address_space = address_space;
  ctx.load_cpu_state(*state);
  seq.runnable = true;
  seq.external_to_core_state(ctx);

  logger.println("Initial state:");
  logger.println("{}", ctx);

  bool exiting = false;

  logger.println(logging::INFO, "Current logenable = {}, start_log_at_iteration = {}, loglevel {}", logger.enabled(),
                 config.log.start_log_at_iteration, config.log.loglevel);

  for (;;) {
    if unlikely (iterations >= config.log.start_log_at_iteration)
      logger.set_enabled(true);

    int running_thread_count = 0;
    if (seq.runnable) {
      running_thread_count++;
      exiting |= seq.execute();
    }

    exiting |= running_thread_count == 0;
    exiting |= iterations >= config.stop_at_iteration || total_user_insns_committed >= config.stop_at_user_insns;

    if unlikely (config.log.event_log_enabled) {
      if unlikely (config.log.flush_event_log_every_cycle)
        eventlog.flush(true);
    }

    iterations++;
    sim_cycle++;
    unhalted_cycle_count += (running_thread_count > 0);
    stats.summary.cycles++;

    if unlikely (exiting)
      break;
  }

  logger.println(logging::INFO, "Exiting sequential mode at {} commits, {} uops and {} iterations (cycles)",
                 total_user_insns_committed, total_uops_committed, iterations);

  seq.core_to_external_state(ctx);
  *state = ctx.to_cpu_state();

  logger.println(logging::VERBOSE, "Core State at end:");
  logger.println(logging::VERBOSE, "{}", ctx);

  return exiting;
}

void SequentialMachine::dump_state() {
  logger.println("Dumping event log for sequential core:");
  eventlog.print();

  core->print_state();
}

void SequentialMachine::update_stats(PTLsimStats& stats) {
  // (nop)
}


} // namespace x86sim

auto std::formatter<x86sim::CommitRecord>::format(const x86sim::CommitRecord& cr, std::format_context& ctx) const {
  using namespace x86sim;
  auto out = ctx.out();
  out = std::format_to(out, "x86sim::CommitRecord: {} stores, {} PTE updates\n", cr.store_count, cr.pte_update_count);
  foreach (i, cr.store_count) {
    out = std::format_to(out, "  Store {:>3}: {}\n", i, cr.stores[i]);
  }
  return out;
}
