//
// PTLsim: Cycle Accurate x86-64 Simulator
// Out-of-Order Core Simulator
// Core Structures
//
// Copyright 2003-2008 Matt T. Yourst <yourst@yourst.com>
// Copyright 2006-2008 Hui Zeng <hzeng@cs.binghamton.edu>
//

#include <cstdlib>
#include <format>
#include <print>
#include <cstdio>
#include <utility>
#include "globals.h"
#include "ptlsim.h"
#include "branchpred.h"
#include "logic.h"
#include "dcache.h"
#include "x86sim/logging.h"

#define INSIDE_OOOCORE
#define DECLARE_STRUCTURES
#include "ooocore.h"
#include "stats.h"

namespace x86sim {

#ifndef ENABLE_CHECKS
#undef assert
#define assert(x) (x)
#endif

#ifndef ENABLE_LOGGING
#endif

using namespace OutOfOrderModel;

namespace OutOfOrderModel {
byte uop_executable_on_cluster[OP_MAX_OPCODE];
W32 forward_at_cycle_lut[MAX_CLUSTERS][MAX_FORWARDING_LATENCY + 1];
}; // namespace OutOfOrderModel

void StateList::init(const char* name, ListOfStateLists& lol, W32 flags) {
  this->name = strdup(name);
  this->flags = flags;
  listid = lol.add(this);
  reset();
}

void StateList::reset() {
  selfqueuelink::reset();
  count = 0;
  dispatch_source_counter = 0;
  issue_source_counter = 0;
}

int ListOfStateLists::add(StateList* list) {
  assert(count < size());
  (*this)[count] = list;
  return count++;
}

void ListOfStateLists::reset() {
  foreach (i, count) {
    (*this)[i]->reset();
  }
}

//
// Initialize lookup tables used by the simulation
//
static void init_luts() {
  // Initialize opcode maps
  foreach (i, OP_MAX_OPCODE) {
    W32 allowedfu = fuinfo[i].fu;
    W32 allowedcl = 0;
    foreach (cl, MAX_CLUSTERS) {
      if (clusters[cl].fu_mask & allowedfu)
        setbit(allowedcl, cl);
    }
    uop_executable_on_cluster[i] = allowedcl;
  }

  // Initialize forward-at-cycle LUTs
  foreach (srcc, MAX_CLUSTERS) {
    foreach (destc, MAX_CLUSTERS) {
      foreach (lat, MAX_FORWARDING_LATENCY + 1) {
        if (lat == intercluster_latency_map[srcc][destc]) {
          setbit(forward_at_cycle_lut[srcc][lat], destc);
        }
      }
    }
  }
}

void ThreadContext::reset() {
  setzero(specrrt);
  setzero(commitrrt);

  setzero(fetchrip);
  current_basic_block = null;
  current_basic_block_transop_index = -1;
  stall_frontend = false;
  waiting_for_icache_fill = false;
  waiting_for_icache_fill_physaddr = 0;
  fetch_uuid = 0;
  current_icache_block = 0;
  loads_in_flight = 0;
  stores_in_flight = 0;
  prev_interrupts_pending = false;
  handle_interrupt_at_next_eom = false;
  stop_at_next_eom = false;

  last_commit_at_cycle = 0;
  smc_invalidate_pending = 0;
  setzero(smc_invalidate_rvp);

  chk_recovery_rip = 0;
  unaligned_ldst_buf.reset();
  consecutive_commits_inside_spinlock = 0;

  total_uops_committed = 0;
  total_insns_committed = 0;
  dispatch_deadlock_countdown = 0;
  issueq_count = 0;
  queued_mem_lock_release_count = 0;
  branchpred.init();
}

void ThreadContext::init() {
  rob_states.reset();
  //
  // ROB states
  //
  rob_free_list("free", rob_states, 0);
  rob_frontend_list("frontend", rob_states, ROB_STATE_PRE_READY_TO_DISPATCH);
  rob_ready_to_dispatch_list("ready-to-dispatch", rob_states, 0);
  InitClusteredROBList(rob_dispatched_list, "dispatched", ROB_STATE_IN_ISSUE_QUEUE);
  InitClusteredROBList(rob_ready_to_issue_list, "ready-to-issue", ROB_STATE_IN_ISSUE_QUEUE);
  InitClusteredROBList(rob_ready_to_store_list, "ready-to-store", ROB_STATE_IN_ISSUE_QUEUE);
  InitClusteredROBList(rob_ready_to_load_list, "ready-to-load", ROB_STATE_IN_ISSUE_QUEUE);
  InitClusteredROBList(rob_issued_list, "issued", 0);
  InitClusteredROBList(rob_completed_list, "completed", ROB_STATE_READY);
  InitClusteredROBList(rob_ready_to_writeback_list, "ready-to-write", ROB_STATE_READY);
  rob_cache_miss_list("cache-miss", rob_states, 0);
  rob_tlb_miss_list("tlb-miss", rob_states, 0);
  rob_memory_fence_list("memory-fence", rob_states, 0);
  rob_ready_to_commit_queue("ready-to-commit", rob_states, ROB_STATE_READY);

  reset();
}

void OutOfOrderCore::reset() {
  round_robin_tid = 0;
  round_robin_reg_file_offset = 0;
  caches.reset();
  caches.callback = &cache_callbacks;
  setzero(robs_on_fu);
  foreach_issueq(reset(*this));

  reserved_iq_entries = (int)std::sqrt((double)ISSUE_QUEUE_SIZE / MAX_THREADS_PER_CORE);
  assert(reserved_iq_entries && reserved_iq_entries < ISSUE_QUEUE_SIZE);

  foreach_issueq(set_reserved_entries(reserved_iq_entries * MAX_THREADS_PER_CORE));
  foreach_issueq(reset_shared_entries());

  unaligned_predictor.reset();

  foreach (i, threadcount)
    threads[i]->reset();
}

void OutOfOrderCore::init_generic() {
  reset();
}

template<typename T>
static void print_list_of_state_lists(const ListOfStateLists& lol, const char* title) {
  logging::println("{}:", title);
  foreach (i, lol.count) {
    StateList& list = *lol[i];
    logging::println("{} ({} entries):", list.name, list.count);
    int n = 0;
    T* obj;
    foreach_list_mutable(list, obj, entry, nextentry) {
      if ((n % 16) == 0)
        logging::print(" ");
      logging::print(" {:<3}", obj->index());
      if (((n % 16) == 15) || (n == list.count - 1))
        logging::println("");
      n++;
    }
    assert(n == list.count);
    logging::println("");
    // list.validate();
  }
}

void StateList::checkvalid() {
#if 0
  int realcount = 0;
  selfqueuelink* obj;
  foreach_list_mutable(*this, obj, entry, nextentry) {
    realcount++;
  }
  assert(count == realcount);
#endif
}

void PhysicalRegister::init(OutOfOrderCore& core, int rfid, int idx) {
  this->core = &core;
  this->coreid = core.coreid;
  this->rfid = rfid;
  this->idx = idx;
  reset();
}

void PhysicalRegisterFile::init(const char* name, OutOfOrderCore& core, int rfid, int size) {
  assert(rfid < PHYS_REG_FILE_COUNT);
  assert(size <= MAX_PHYS_REG_FILE_SIZE);
  this->size = size;
  this->core = &core;
  this->coreid = core.coreid;
  this->rfid = rfid;
  this->name = name;
  this->allocations = 0;
  this->frees = 0;

  foreach (i, MAX_PHYSREG_STATE) {
    std::string name_str = std::format("{}-{}", name, physreg_state_names[i]);
    states[i].init(name_str.c_str(), getcore().physreg_states);
  }

  foreach (i, size) {
    (*this)[i].init(core, rfid, i);
  }
}

PhysicalRegister* PhysicalRegisterFile::alloc(W8 threadid, int r) {
  PhysicalRegister* physreg = (PhysicalRegister*)((r == 0) ? &(*this)[r] : states[PHYSREG_FREE].peek());
  if unlikely (!physreg)
    return null;
  physreg->changestate(PHYSREG_WAITING);
  physreg->flags = FLAG_WAIT;
  physreg->threadid = threadid;
  allocations++;

  assert(states[PHYSREG_FREE].count >= 0);
  return physreg;
}

void PhysicalRegisterFile::reset(W8 threadid) {
  foreach (i, size) {
    if ((*this)[i].threadid == threadid) {
      (*this)[i].reset(threadid);
    }
  }
}

void PhysicalRegisterFile::reset() {
  foreach (i, MAX_PHYSREG_STATE) {
    states[i].reset();
  }

  foreach (i, size) {
    (*this)[i].reset(0, false);
  }
}

StateList& PhysicalRegister::get_state_list(int s) const {
  return getcore().physregfiles[rfid].states[s];
}

//
// Get the thread priority, with lower numbers receiving higher priority.
// This is used to regulate the order in which fetch, rename, frontend
// and dispatch slots are filled in each cycle.
//
// The well known ICOUNT algorithm adds up the number of uops in
// the frontend pipeline stages and gives highest priority to
// the thread with the lowest number, since this thread is moving
// uops through very quickly and can make more progress.
//
int ThreadContext::get_priority() const {
  int priority = fetchq.count + rob_frontend_list.count + rob_ready_to_dispatch_list.count;

  for_each_cluster(cluster) {
    priority += rob_dispatched_list[cluster].count + rob_ready_to_issue_list[cluster].count +
                rob_ready_to_store_list[cluster].count + rob_ready_to_load_list[cluster].count;
  }

  return priority;
}

//
// Execute one cycle of the entire core state machine
//
bool OutOfOrderCore::runcycle() {
  bool exiting = 0;

  //
  // Compute reserved issue queue entries to avoid starvation:
  //
#ifdef ENABLE_CHECKS
  int total_issueq_count = 0;
  int total_issueq_reserved_free = 0;

  foreach (i, MAX_THREADS_PER_CORE) {
    auto& thread = threads[i];

    if unlikely (!thread) {
      total_issueq_reserved_free += reserved_iq_entries;
    } else {
      total_issueq_count += thread->issueq_count;
      if (thread->issueq_count < reserved_iq_entries) {
        total_issueq_reserved_free += reserved_iq_entries - thread->issueq_count;
      }
    }
  }

  // assert (total_issueq_count == issueq_all.count);
  // assert((ISSUE_QUEUE_SIZE - issueq_all.count) == (issueq_all.shared_entries + total_issueq_reserved_free));
#endif

  foreach (i, threadcount)
    threads[i]->loads_in_this_cycle = 0;

  fu_avail = bitmask(FU_COUNT);
  caches.clock();

  //
  // Backend and issue pipe stages run with round robin priority
  //
  int commitrc[MAX_THREADS_PER_CORE];
  commitcount = 0;
  writecount = 0;

  foreach (permute, threadcount) {
    int tid = add_index_modulo(round_robin_tid, +permute, threadcount);
    auto& thread = threads[tid];
    if unlikely (!thread->ctx.running)
      continue;

    commitrc[tid] = thread->commit();
    for_each_cluster(j) thread->writeback(j);
    for_each_cluster(j) thread->transfer(j);
  }


  //
  // Always clock the issue queues: they're independent of all threads
  // NOTE(AE): do it here so that forwarding can happen without one cycle delay
  //
  foreach_issueq(clock());

  //
  // Issue whatever is ready
  //
  for_each_cluster(i) {
    issue(i);
  }

  //
  // Most of the frontend (except fetch!) also works with round robin priority
  //
  int dispatchrc[MAX_THREADS_PER_CORE];
  dispatchcount = 0;
  foreach (permute, threadcount) {
    int tid = add_index_modulo(round_robin_tid, +permute, threadcount);
    auto& thread = threads[tid];
    if unlikely (!thread->ctx.running)
      continue;

    dispatchrc[tid] = thread->dispatch();

    for_each_cluster(j) {
      thread->complete(j);
    }

    if likely (dispatchrc[tid] >= 0) {
      thread->frontend();
      thread->rename();
    }
  }

  //
  // Compute fetch priorities (default is ICOUNT algorithm)
  //
  // This means we sort in ascending order, with any unused threads
  // (if any) given the lowest priority.
  //

  int priority_value[MAX_THREADS_PER_CORE];
  int priority_index[MAX_THREADS_PER_CORE];

  if likely (threadcount == 1) {
    priority_value[0] = 0;
    priority_index[0] = 0;
  } else {
    foreach (i, threadcount) {
      priority_index[i] = i;
      auto& thread = threads[i];
      priority_value[i] = thread->get_priority();
      if unlikely (!thread->ctx.running)
        priority_value[i] = std::numeric_limits<int>::max();
    }

    sort(priority_index, threadcount, SortPrecomputedIndexListComparator<int, false>(priority_value));
  }

  //
  // Fetch in thread priority order
  //
  // NOTE: True ICOUNT only fetches the highest priority
  // thread per cycle, since there is usually only one
  // instruction cache port. In a banked i-cache, we can
  // fetch from multiple threads every cycle.
  //
  foreach (j, threadcount) {
    int i = priority_index[j];
    auto& thread = threads[i];
    assert(thread);
    if unlikely (!thread->ctx.running) {
      continue;
    }

    if likely (dispatchrc[i] >= 0) {
      thread->fetch();
    }
  }

  //
  // Advance the round robin priority index
  //
  round_robin_tid = add_index_modulo(round_robin_tid, +1, threadcount);

  //
  // Flush event log ring buffer
  //
  if unlikely (config.log.event_log_enabled) {
    logging::println(logging::INFO, "[cycle {}] Miss buffer contents:", machine.sim_cycle);
    logging::println(logging::INFO, "{}", caches.missbuf);
    if unlikely (config.log.flush_event_log_every_cycle) {
      eventlog.flush(true);
    }
  }

#ifdef ENABLE_CHECKS
  // This significantly slows down simulation; only enable it if absolutely needed:
  // check_refcounts();
#endif

  foreach (i, threadcount) {
    auto& thread = threads[i];
    if unlikely (!thread->ctx.running)
      continue;
    int rc = commitrc[i];
    if likely ((rc == COMMIT_RESULT_OK) | (rc == COMMIT_RESULT_NONE))
      continue;

    switch (rc) {
    case COMMIT_RESULT_SMC: {
      logging::println(logging::INFO,
                       "Potentially cross-modifying SMC detected: global flush required (cycle {}, {} commits)",
                       machine.sim_cycle, machine.total_user_insns_committed);
      logging::flush();
      //
      // DO NOT GLOBALLY FLUSH! It will cut off the other thread(s) in the
      // middle of their currently committing x86 instruction, causing massive
      // internal corruption on any VCPUs that happen to be straddling the
      // instruction boundary.
      //
      // BAD: machine.flush_all_pipelines();
      //
      // This is a temporary fix: in the *extremely* rare case where both
      // threads have the same basic block in their pipelines and that
      // BB is being invalidated, the BB cache will forbid us from
      // freeing it (and will print a warning to that effect).
      //
      // I'm working on a solution to this, to put some BBs on an
      // "invisible" list, where they cannot be looked up anymore,
      // but their memory is not freed until the lock is released.
      //
      foreach (i, threadcount) {
        auto& t = threads[i];
        if unlikely (!t)
          continue;

        logging::println(logging::DEBUG, "[vcpu {}] current_basic_block = {}: ", i, (void*)t->current_basic_block);
        if (t->current_basic_block)
          logging::println(logging::DEBUG, "  [vcpu {}] current_basic_block = {}: {}", i, (void*)t->current_basic_block,
                           (void*)(W64)t->current_basic_block->rip);
        else
          logging::println(logging::DEBUG, "  [vcpu {}] current_basic_block = {}", i, (void*)t->current_basic_block);
      }

      thread->flush_pipeline();
      thread->invalidate_smc();
      break;
    }
    case COMMIT_RESULT_EXCEPTION: {
      exiting = !thread->handle_exception();
      break;
    }
    case COMMIT_RESULT_BARRIER: {
      exiting = !thread->handle_barrier();
      break;
    }
    case COMMIT_RESULT_INTERRUPT: {
      thread->handle_interrupt();
      break;
    }
    case COMMIT_RESULT_STOP: {
      thread->flush_pipeline();
      thread->stall_frontend = 1;
      machine.stopped[thread->ctx.vcpuid] = 1;
      // Wait for other cores to sync up, so don't exit right away
      break;
    }
    }
  }

  foreach (i, threadcount) {
    auto& thread = threads[i];
    if unlikely (!thread->ctx.running)
      break;

    if unlikely ((machine.sim_cycle - thread->last_commit_at_cycle) > 4096) {
      logging::println(logging::ERROR,
                       "[vcpu {}] thread {}: WARNING: At cycle {}, {} user commits: no instructions have "
                       "committed for {} cycles; the pipeline could be deadlocked",
                       thread->ctx.vcpuid, thread->threadid, machine.sim_cycle, machine.total_user_insns_committed,
                       (machine.sim_cycle - thread->last_commit_at_cycle));
      logging::flush();
      exiting = 1;
    }
  }

  return exiting;
}

//
// ReorderBufferEntry
//
void ReorderBufferEntry::init(int idx) {
  this->idx = idx;
  entry_valid = 0;
  selfqueuelink::reset();
  current_state_list = null;
  reset();
}

//
// Clean out various fields from the ROB entry that are
// expected to be zero when allocating a new ROB entry.
//
void ReorderBufferEntry::reset() {
  int latency, operand;
  // Deallocate ROB entry
  entry_valid = false;
  cycles_left = 0;
  physreg = (PhysicalRegister*)null;
  lfrqslot = -1;
  lsq = 0;
  load_store_second_phase = 0;
  lock_acquired = 0;
  consumer_count = 0;
  executable_on_cluster_mask = 0;
  pteupdate = 0;
  cluster = -1;
#ifdef ENABLE_TRANSIENT_VALUE_TRACKING
  dest_renamed_before_writeback = 0;
  no_branches_between_renamings = 0;
#endif
  issued = 0;
}

bool ReorderBufferEntry::ready_to_issue() const {
  bool raready = operands[0]->ready();
  bool rbready = operands[1]->ready();
  bool rcready = operands[2]->ready();
  bool rsready = operands[3]->ready();

  if (isstore(uop.opcode)) {
    return (load_store_second_phase) ? (raready & rbready & rcready & rsready) : (raready & rbready);
  } else if (isload(uop.opcode)) {
    return (load_store_second_phase) ? (raready & rbready & rcready & rsready) : (raready & rbready & rcready);
  } else {
    return (raready & rbready & rcready & rsready);
  }
}

bool ReorderBufferEntry::ready_to_commit() const {
  return (current_state_list == &getthread().rob_ready_to_commit_queue);
}

StateList& ReorderBufferEntry::get_ready_to_issue_list() const {
  OutOfOrderCore& core = getcore();
  ThreadContext& thread = getthread();
  return isload(uop.opcode)    ? thread.rob_ready_to_load_list[cluster]
         : isstore(uop.opcode) ? thread.rob_ready_to_store_list[cluster]
                               : thread.rob_ready_to_issue_list[cluster];
}

//
// Reorder Buffer
//
std::string ReorderBufferEntry::get_operand_info(int operand) const {
  PhysicalRegister& physreg = *operands[operand];
  ReorderBufferEntry& sourcerob = *physreg.rob;

  std::string result = std::format("r{}", physreg.index());
  if (PHYS_REG_FILE_COUNT > 1)
    result += std::format("@{}", getcore().physregfiles[physreg.rfid].name);

  switch (physreg.state) {
  case PHYSREG_WRITTEN:
    result += " (written)";
    break;
  case PHYSREG_BYPASS:
    result += " (ready)";
    break;
  case PHYSREG_WAITING:
    result += std::format(" (wait rob {} uuid {})", sourcerob.index(), sourcerob.uop.uuid);
    break;
  case PHYSREG_ARCH:
    break;
    if (physreg.index() == PHYS_REG_NULL)
      result += " (zero)";
    else
      result += std::format(" (arch {})", arch_reg_names[physreg.archreg]);
    break;
  case PHYSREG_PENDINGFREE:
    result += std::format(" (pending free for {})", arch_reg_names[physreg.archreg]);
    break;
  default:
    // Cannot be in free state!
    result += " (FREE)";
    break;
  }

  return result;
}

ThreadContext& ReorderBufferEntry::getthread() const {
  return *getcore().threads[threadid];
}

issueq_tag_t ReorderBufferEntry::get_tag() {
  int mask = ((1 << MAX_THREADS_BIT) - 1) << MAX_ROB_IDX_BIT;
  logging::println(logging::VERBOSE, " get_tag() thread {} rob idx {} mask 0x{:x}", threadid, idx, mask);

  assert(!(idx & mask));
  assert(!(threadid >> MAX_THREADS_BIT));
  //  int threadid = 1;
  issueq_tag_t rc = (idx | (threadid << MAX_ROB_IDX_BIT));
  logging::println(logging::VERBOSE, " tag 0x{:x}", rc);
  return rc;
}


void OutOfOrderCore::print_smt_state() {
  logging::println("Print SMT statistics:");

  foreach (i, threadcount) {
    auto& thread = threads[i];
    logging::println("Thread {}:", i);
    logging::println("  total_uops_committed {}", thread->total_uops_committed);
    logging::println("  uipc {}", double(thread->total_uops_committed) / double(machine.iterations));
    logging::println("  total_insns_committed {}", thread->total_insns_committed);
    logging::println("  ipc {}", double(thread->total_insns_committed) / double(machine.iterations));
  }
}

} // namespace x86sim

auto std::formatter<x86sim::OutOfOrderModel::ThreadContext>::format(
    const x86sim::OutOfOrderModel::ThreadContext& thread, std::format_context& fctx) const {
  using namespace x86sim;
  auto out = fctx.out();
  out = std::format_to(out, "SMT per-thread state for t{}:\n", thread.threadid);
  out = std::format_to(out, "SpecRRT:\n");
  out = std::format_to(out, "{}", thread.specrrt);
  out = std::format_to(out, "CommitRRT:\n");
  out = std::format_to(out, "{}\n", thread.commitrrt);
  out = std::format_to(out, "{}\n", thread.ROB);
  out = std::format_to(out, "{}\n", thread.LSQ);
  return out;
}


auto std::formatter<x86sim::OutOfOrderModel::OutOfOrderCore>::format(
    const x86sim::OutOfOrderModel::OutOfOrderCore& core, std::format_context& ctx) const {
  using namespace x86sim;
  auto out = ctx.out();

  out = std::format_to(out, "SMT common structures:\n");
  foreach (i, PHYS_REG_FILE_COUNT) {
    out = std::format_to(out, "{}\n", core.physregfiles[i]);
  }
  out = std::format_to(out, "Issue Queues: int0={} int1={} ld={} fp={}\n", core.issueq_int0.count,
                       core.issueq_int1.count, core.issueq_ld.count, core.issueq_fp.count);

  out = std::format_to(out, "Caches: lfrq_count={} missbuf_free={}\n", core.caches.lfrq.count,
                       core.caches.missbuf.freemap.count());

  out = std::format_to(out, "Unaligned predictor:\n");
  out = std::format_to(out, "  {} unaligned bits out of {} bits\n", core.unaligned_predictor.count(),
                       UNALIGNED_PREDICTOR_SIZE);
  out = std::format_to(out, "  Raw data: {}\n", core.unaligned_predictor.to_string());

  foreach (i, core.threadcount) {
    out = std::format_to(out, "{}\n", *core.threads[i]);
  }
  return out;
}

namespace x86sim {

//
// Validate the physical register reference counters against what
// is really accessible from the various tables and operand fields.
//
// This is for debugging only.
//
void OutOfOrderCore::check_refcounts() {
  // this should be for each thread instead of whole core:
  // for now, we just work on thread[0];
  ThreadContext& thread = *threads[0];
  Queue<ReorderBufferEntry, ROB_SIZE>& ROB = thread.ROB;
  RegisterRenameTable& specrrt = thread.specrrt;
  RegisterRenameTable& commitrrt = thread.commitrrt;

  int refcounts[PHYS_REG_FILE_COUNT][MAX_PHYS_REG_FILE_SIZE];
  memset(refcounts, 0, sizeof(refcounts));

  foreach (rfid, PHYS_REG_FILE_COUNT) {
    // Null physreg in each register file is special and can never be freed:
    refcounts[rfid][PHYS_REG_NULL]++;
  }

  foreach_forward(ROB, i) {
    ReorderBufferEntry& rob = ROB[i];
    foreach (j, MAX_OPERANDS) {
      refcounts[rob.operands[j]->rfid][rob.operands[j]->index()]++;
    }
  }

  foreach (i, TRANSREG_COUNT) {
    refcounts[commitrrt[i]->rfid][commitrrt[i]->index()]++;
    refcounts[specrrt[i]->rfid][specrrt[i]->index()]++;
  }

  bool errors = 0;

  foreach (rfid, PHYS_REG_FILE_COUNT) {
    PhysicalRegisterFile& physregs = physregfiles[rfid];
    foreach (i, physregs.size) {
      if unlikely (physregs[i].refcount != refcounts[rfid][i]) {
        logging::println(logging::ERROR, "ERROR: r{} refcount is {} but should be {}", i, physregs[i].refcount,
                         refcounts[rfid][i]);

        foreach_forward(ROB, r) {
          ReorderBufferEntry& rob = ROB[r];
          foreach (j, MAX_OPERANDS) {
            if ((rob.operands[j]->index() == i) & (rob.operands[j]->rfid == rfid))
              logging::println(logging::ERROR, "  ROB {} operand {}", r, j);
          }
        }

        foreach (j, TRANSREG_COUNT) {
          if ((commitrrt[j]->index() == i) & (commitrrt[j]->rfid == rfid))
            logging::println(logging::ERROR, "  CommitRRT {}", arch_reg_names[j]);
          if ((specrrt[j]->index() == i) & (specrrt[j]->rfid == rfid))
            logging::println(logging::ERROR, "  SpecRRT {}", arch_reg_names[j]);
        }

        errors = 1;
      }
    }
  }

  if (errors)
    assert(false);
}

void OutOfOrderCore::check_rob() {
  // this should be for each thread instead of whole core:
  // for now, we just work on thread[0];
  ThreadContext& thread = *threads[0];
  Queue<ReorderBufferEntry, ROB_SIZE>& ROB = thread.ROB;

  foreach (i, ROB_SIZE) {
    ReorderBufferEntry& rob = ROB[i];
    if (!rob.entry_valid)
      continue;
    assert(inrange((int)rob.forward_cycle, 0, (MAX_FORWARDING_LATENCY + 1) - 1));
  }

  foreach (i, threadcount) {
    auto& thread = threads[i];
    foreach (i, rob_states.count) {
      StateList& list = *(thread->rob_states[i]);
      ReorderBufferEntry* rob;
      foreach_list_mutable(list, rob, entry, nextentry) {
        assert(inrange(rob->index(), 0, ROB_SIZE - 1));
        assert(rob->current_state_list == &list);
        if (!((rob->current_state_list != &thread->rob_free_list) ? rob->entry_valid : (!rob->entry_valid))) {
          logging::println(logging::ERROR, "ROB {} list = {} entry_valid {}", rob->index(),
                           rob->current_state_list->name, static_cast<int>(rob->entry_valid));
          logging::flush();
          dump_smt_state();
          assert(false);
        }
      }
    }
  }
}

//
// Barriers must flush the fetchq and stall the frontend until
// after the barrier is consumed. Execution resumes at the address
// in internal register nextrip (rip after the instruction) after
// handling the barrier in microcode.
//
bool ThreadContext::handle_barrier() {
  // Release resources of everything in the pipeline:

  core_to_external_state();
  flush_pipeline();

  int assistid = ctx.commitarf[REG_rip];
  assist_func_t assist = (assist_func_t)(Waddr)assistid_to_func[assistid];

  logging::println(
      logging::INFO, "[vcpu {}] Barrier (#{} -> {} {} called from {}; return to {}) at {} cycles, {} commits",
      ctx.vcpuid, assistid, (void*)assist, assist_name(assist), RIPVirtPhys(ctx.commitarf[REG_selfrip]).update(ctx),
      (void*)(Waddr)ctx.commitarf[REG_nextrip], core.machine.sim_cycle, core.machine.total_user_insns_committed);
  logging::flush();

  logging::println(logging::DEBUG, "Calling assist function at {}...", (void*)assist);
  logging::flush();

  update_assist_stats(assist);
  logging::println(logging::DEBUG, "Before assist:");
  logging::println(logging::DEBUG, "{}", ctx);

  assist(ctx);

  logging::println(logging::DEBUG, "Done with assist");
  logging::println(logging::DEBUG, "New state:");
  logging::println(logging::DEBUG, "{}", ctx);

  // Flush again, but restart at possibly modified rip
  flush_pipeline();

  if (requested_switch_to_native) {
    logging::println(logging::INFO, "PTL call requested switch to native mode at rip {}",
                     (void*)(Waddr)ctx.commitarf[REG_rip]);
    return false;
  }
  return true;
}

bool ThreadContext::handle_exception() {
  // Release resources of everything in the pipeline:
  core_to_external_state();
  flush_pipeline();

  logging::println(logging::INFO, "[vcpu {}] Exception {} called from rip {} at {} cycles, {} commits", ctx.vcpuid,
                   ctx.exception, (void*)(Waddr)ctx.commitarf[REG_rip], core.machine.sim_cycle,
                   core.machine.total_user_insns_committed);
  logging::flush();

  //
  // CheckFailed and SkipBlock exceptions are raised by the chk uop.
  // This uop is used at the start of microcoded instructions to assert
  // that certain conditions are true so complex corrective actions can
  // be taken if the check fails.
  //
  // SkipBlock is a special case used for checks at the top of REP loops.
  // Specifically, if the %rcx register is zero on entry to the REP, no
  // action at all is to be taken; the rip should simply advance to
  // whatever is in chk_recovery_rip and execution should resume.
  //
  // CheckFailed exceptions usually indicate the processor needs to take
  // evasive action to avoid a user visible exception. For instance,
  // CheckFailed is raised when an inlined floating point operand is
  // denormal or otherwise cannot be handled by inlined fastpath uops,
  // or when some unexpected segmentation or page table conditions
  // arise.
  //
  if (ctx.exception == EXCEPTION_SkipBlock) {
    ctx.commitarf[REG_rip] = chk_recovery_rip;
    logging::println(logging::DEBUG, "SkipBlock pseudo-exception: skipping to {}",
                     (void*)(Waddr)ctx.commitarf[REG_rip]);
    logging::flush();
    flush_pipeline();
    return true;
  }

  //
  // Map PTL internal hardware exceptions to their x86 equivalents,
  // depending on the context. The error_code field should already
  // be filled out.
  //
  // Exceptions not listed here are propagated by microcode
  // rather than the processor itself.
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
  case EXCEPTION_DivideOverflow:
    ctx.x86_exception = EXCEPTION_x86_divide;
    break;
  case EXCEPTION_InvalidAddr:
    ctx.x86_exception = EXCEPTION_x86_gp_fault;
    break;
  default:
    logging::println(logging::ERROR, "Unsupported internal exception type {}", static_cast<int>(ctx.exception));
    logging::flush();
    assert(false);
  }

  logging::println(logging::INFO, "{}", ctx);

  ctx.propagate_x86_exception(ctx.x86_exception, ctx.error_code, ctx.cr2);

  // Flush again, but restart at modified rip
  flush_pipeline();

  return true;
}

bool ThreadContext::handle_interrupt() {
  return true;
}

//
// Event Formatting
//
void PhysicalRegister::fill_operand_info(PhysicalRegisterOperandInfo& opinfo) {
  opinfo.physreg = index();
  opinfo.state = state;
  opinfo.rfid = rfid;
  opinfo.archreg = archreg;
  if (rob) {
    opinfo.rob = rob->index();
    opinfo.uuid = rob->uop.uuid;
  }
}

bool EventLog::init(size_t bufsize) {
  reset();
  size_t bytes = bufsize * sizeof(OutOfOrderCoreEvent);
  start = (OutOfOrderCoreEvent*)std::aligned_alloc(PAGE_SIZE, ceil(bytes, PAGE_SIZE));
  if unlikely (!start)
    return false;
  end = start + bufsize;
  tail = start;

  foreach (i, bufsize)
    start[i].type = EVENT_INVALID;
  return true;
}

void EventLog::reset() {
  if (!start)
    return;

  std::free(start);
  start = null;
  end = null;
  tail = null;
}

void EventLog::flush(bool only_to_tail) {
  print(only_to_tail);
  tail = start;
}

void EventLog::print(bool only_to_tail) {
  Options& config = machine.config;

  if (tail >= end)
    tail = start;
  if (tail < start)
    tail = end;

  OutOfOrderCoreEvent* p = (only_to_tail) ? start : tail;

  W64 cycle = std::numeric_limits<W64>::max();
  size_t bufsize = end - start;

  if (!config.log.flush_event_log_every_cycle)
    logging::println(logging::INFO, "#-------- Start of event log --------");

  foreach (i, (only_to_tail ? (tail - start) : bufsize)) {
    if unlikely (p >= end)
      p = start;
    if unlikely (p < start)
      p = end - 1;
    if unlikely (p->type == EVENT_INVALID) {
      p++;
      continue;
    }

    if unlikely (p->cycle != cycle) {
      cycle = p->cycle;
      logging::println(logging::INFO, "Cycle {}:", cycle);
    }

    logging::print(logging::INFO, "{}", *p);
    p++;
  }

  if (!config.log.flush_event_log_every_cycle)
    logging::println(logging::INFO, "#-------- End of event log --------");
}

} // namespace x86sim

auto std::formatter<x86sim::OutOfOrderModel::OutOfOrderCoreEvent>::format(
    const x86sim::OutOfOrderModel::OutOfOrderCoreEvent& ev, std::format_context& ctx) const {
  using namespace x86sim;
  auto out = ctx.out();
  bool ld = isload(ev.uop.opcode);
  bool st = isstore(ev.uop.opcode);
  bool br = isbranch(ev.uop.opcode);

  std::string uopname = nameof(ev.uop);

  out = std::format_to(out, "{:>20} t{} ", ev.uuid, ev.threadid);
  switch (ev.type) {
    //
    // Fetch Events
    //
  case EVENT_FETCH_STALLED:
    out = std::format_to(out, "fetch  frontend stalled");
    break;
  case EVENT_FETCH_ICACHE_WAIT:
    out = std::format_to(out, "fetch  rip {}: wait for icache fill", ev.rip);
    break;
  case EVENT_FETCH_FETCHQ_FULL:
    out = std::format_to(out, "fetch  rip {}: fetchq full", ev.rip);
    break;
  case EVENT_FETCH_IQ_QUOTA_FULL:
    out = std::format_to(out, "fetch  rip {}: issue queue quota full = {} ", ev.rip, ev.issueq_count);
    break;
  case EVENT_FETCH_BOGUS_RIP:
    out = std::format_to(out, "fetch  rip {}: bogus RIP or decode failed", ev.rip);
    break;
  case EVENT_FETCH_ICACHE_MISS:
    out = std::format_to(out, "fetch  rip {}: wait for icache fill of phys {} on missbuf {}", ev.rip,
                         (void*)(Waddr)((ev.rip.mfnlo << 12) + lowbits(ev.rip.rip, 12)), ev.fetch.missbuf);
    break;
  case EVENT_FETCH_SPLIT:
    out = std::format_to(out, "fetch  rip {}: split unaligned load or store {}", ev.rip, ev.uop);
    break;
  case EVENT_FETCH_ASSIST:
    out = std::format_to(out, "fetch  rip {}: branch into assist microcode: {}", ev.rip, ev.uop);
    break;
  case EVENT_FETCH_TRANSLATE:
    out = std::format_to(out, "xlate  rip {}: {} uops", ev.rip, ev.fetch.bb_uop_count);
    break;
  case EVENT_FETCH_OK: {
    out = std::format_to(out, "fetch  rip {}: {} (uopid {}", ev.rip, ev.uop, ev.uop.bbindex);
    if (ev.uop.som)
      out = std::format_to(out, "; SOM");
    if (ev.uop.eom)
      out = std::format_to(out, "; EOM {} bytes", ev.uop.bytes);
    out = std::format_to(out, ")");
    if (ev.uop.eom && ev.fetch.predrip)
      out = std::format_to(out, " -> pred {}", (void*)ev.fetch.predrip);
    if (isload(ev.uop.opcode) | isstore(ev.uop.opcode)) {
      out = std::format_to(out, "; unaligned pred slot {} -> {}", OutOfOrderCore::hash_unaligned_predictor_slot(ev.rip),
                           ev.uop.unaligned);
    }
    break;
  }
    //
    // Rename Events
    //
  case EVENT_RENAME_FETCHQ_EMPTY:
    out = std::format_to(out, "rename fetchq empty");
    break;
  case EVENT_RENAME_ROB_FULL:
    out = std::format_to(out, "rename ROB full");
    break;
  case EVENT_RENAME_PHYSREGS_FULL:
    out = std::format_to(out, "rename physical register file full");
    break;
  case EVENT_RENAME_LDQ_FULL:
    out = std::format_to(out, "rename load queue full");
    break;
  case EVENT_RENAME_STQ_FULL:
    out = std::format_to(out, "rename store queue full");
    break;
  case EVENT_RENAME_MEMQ_FULL:
    out = std::format_to(out, "rename memory queue full");
    break;
  case EVENT_RENAME_OK: {
    out = std::format_to(out, "rename rob {:<3}({:<5}) r{:<3}@{}", ev.rob, uopname, ev.physreg,
                         phys_reg_file_names[ev.rfid]);
    if (ld | st)
      out = std::format_to(out, " lsq{}", ev.lsq);
    out = std::format_to(out, " = ");
    foreach (i, MAX_OPERANDS)
      out = std::format_to(out, "{}{}", ev.rename.opinfo[i], ((i < MAX_OPERANDS - 1) ? " " : ""));
    out = std::format_to(out, "; renamed");
    out = std::format_to(out, " {} (old r{})", arch_reg_names[ev.uop.rd], ev.rename.oldphys);
    if unlikely (!ev.uop.nouserflags) {
      if likely (ev.uop.setflags & SETFLAG_ZF)
        out = std::format_to(out, " zf (old r{})", ev.rename.oldzf);
      if likely (ev.uop.setflags & SETFLAG_CF)
        out = std::format_to(out, " cf (old r{})", ev.rename.oldcf);
      if likely (ev.uop.setflags & SETFLAG_OF)
        out = std::format_to(out, " of (old r{})", ev.rename.oldof);
    }
    break;
  }
  case EVENT_FRONTEND:
    out = std::format_to(out, "front  rob {:<3}({:<5}) frontend stage {} of {}", ev.rob, uopname,
                         (FRONTEND_STAGES - ev.frontend.cycles_left), FRONTEND_STAGES);
    break;
  case EVENT_CLUSTER_NO_CLUSTER:
  case EVENT_CLUSTER_OK: {
    out = std::format_to(out, "{} rob {:<3}({:<5}) allowed FUs = {} -> clusters {} avail",
                         ((ev.type == EVENT_CLUSTER_OK) ? "clustr" : "noclus"), ev.rob, uopname,
                         bitstring(fuinfo[ev.uop.opcode].fu, FU_COUNT, true),
                         bitstring(ev.select_cluster.allowed_clusters, MAX_CLUSTERS, true));
    foreach (i, MAX_CLUSTERS)
      out = std::format_to(out, " {}", ev.select_cluster.iq_avail[i]);
    out = std::format_to(out, "-> ");
    if (ev.type == EVENT_CLUSTER_OK)
      out = std::format_to(out, "cluster {}", clusters[ev.cluster].name);
    else
      out = std::format_to(out, "-> none");
    break;
  }
  case EVENT_DISPATCH_NO_CLUSTER:
  case EVENT_DISPATCH_OK: {
    out = std::format_to(out, "{} rob {:<3}({:<5}) operands ", ((ev.type == EVENT_DISPATCH_OK) ? "disptc" : "nodisp"),
                         ev.rob, uopname);
    foreach (i, MAX_OPERANDS)
      out = std::format_to(out, "{}{}", ev.dispatch.opinfo[i], ((i < MAX_OPERANDS - 1) ? " " : ""));
    if (ev.type == EVENT_DISPATCH_OK)
      out = std::format_to(out, " -> cluster {}", clusters[ev.cluster].name);
    else
      out = std::format_to(out, " -> none");
    break;
  }
  case EVENT_ISSUE_NO_FU: {
    out = std::format_to(out, "issue  rob {:<3}({:<5})", ev.rob, uopname);
    out = std::format_to(out, "no FUs available in cluster {}: fu_avail = {}, op_fu = {}, fu_cl_mask = {}",
                         clusters[ev.cluster].name, bitstring(ev.issue.fu_avail, FU_COUNT, true),
                         bitstring(fuinfo[ev.uop.opcode].fu, FU_COUNT, true),
                         bitstring(clusters[ev.cluster].fu_mask, FU_COUNT, true));
    break;
  }
  case EVENT_ISSUE_OK: {
    out = std::format_to(out, "issue  rob {:<3}({:<5})", ev.rob, uopname);
    out = std::format_to(out, " on {:<4} in {:<4}: r{:<3}@{}", fu_names[ev.fu], cluster_names[ev.cluster], ev.physreg,
                         phys_reg_file_names[ev.rfid]);
    out = std::format_to(out, " {}", format_value_and_flags(ev.issue.state.reg.rddata, ev.issue.state.reg.rdflags));
    out = std::format_to(out, " = {}", format_value_and_flags(ev.issue.operand_data[RA], ev.issue.operand_flags[RA]));
    out = std::format_to(out, ",  {}", format_value_and_flags(ev.issue.operand_data[RB], ev.issue.operand_flags[RB]));
    out = std::format_to(out, ",  {}", format_value_and_flags(ev.issue.operand_data[RC], ev.issue.operand_flags[RC]));
    out = std::format_to(out, " ({} cycles left)", ev.issue.cycles_left);
    if (ev.issue.mispredicted)
      out = std::format_to(out, "; mispredicted (real {} vs expected {})", (void*)(Waddr)ev.issue.state.reg.rddata,
                           (void*)(Waddr)ev.issue.predrip);
    break;
  }
  case EVENT_REPLAY: {
    out = std::format_to(out, "replay rob {:<3}({:<5}) r{:<3}@{} on cluster {}: waiting on", ev.rob, uopname,
                         ev.physreg, phys_reg_file_names[ev.rfid], clusters[ev.cluster].name);
    foreach (i, MAX_OPERANDS) {
      if (!bit(ev.replay.ready, i))
        out = std::format_to(out, " {}", ev.replay.opinfo[i]);
    }
    break;
  }
  case EVENT_STORE_WAIT: {
    out = std::format_to(out, "store{} rob {:<3}({:<5}) stq {} r{:<3} on {:<4} @ {} (phys {}): ",
                         (ev.loadstore.load_store_second_phase ? "2" : " "), ev.rob, uopname, ev.lsq, ev.physreg,
                         fu_names[ev.fu], (void*)(Waddr)ev.loadstore.virtaddr,
                         (void*)(Waddr)(ev.loadstore.sfr.physaddr << 3));
    out = std::format_to(out, "wait on ");
    if (!ev.loadstore.rcready)
      out = std::format_to(out, " rc");
    if (ev.loadstore.inherit_sfr_used) {
      out = std::format_to(out, "{}{} (uuid {}, stq {}, rob {}, r{})", ((ev.loadstore.rcready) ? "" : " and "),
                           ev.loadstore.inherit_sfr, ev.loadstore.inherit_sfr_uuid, ev.loadstore.inherit_sfr_lsq,
                           ev.loadstore.inherit_sfr_rob, ev.loadstore.inherit_sfr_physreg);
    }
    break;
  }
  case EVENT_STORE_PARALLEL_FORWARDING_MATCH: {
    out = std::format_to(out, "store{} rob {:<3}({:<5}) stq {} r{:<3} on {:<4} @ {} (phys {}): ",
                         (ev.loadstore.load_store_second_phase ? "2" : " "), ev.rob, uopname, ev.lsq, ev.physreg,
                         fu_names[ev.fu], (void*)(Waddr)ev.loadstore.virtaddr,
                         (void*)(Waddr)(ev.loadstore.sfr.physaddr << 3));
    out = std::format_to(out, "ignored parallel forwarding match with ldq {} (uuid {} rob{} r{})",
                         ev.loadstore.inherit_sfr_lsq, ev.loadstore.inherit_sfr_uuid, ev.loadstore.inherit_sfr_rob,
                         ev.loadstore.inherit_sfr_physreg);
    break;
  }
  case EVENT_STORE_ALIASED_LOAD: {
    out = std::format_to(out, "store{} rob {:<3}({:<5}) stq {} r{:<3} on {:<4} @ {} (phys {}): ",
                         (ev.loadstore.load_store_second_phase ? "2" : " "), ev.rob, uopname, ev.lsq, ev.physreg,
                         fu_names[ev.fu], (void*)(Waddr)ev.loadstore.virtaddr,
                         (void*)(Waddr)(ev.loadstore.sfr.physaddr << 3));
    out = std::format_to(out,
                         "aliased with ldbuf {} (uuid {} rob{} r{}); (add colliding load rip {}; replay from rip {})",
                         ev.loadstore.inherit_sfr_lsq, ev.loadstore.inherit_sfr_uuid, ev.loadstore.inherit_sfr_rob,
                         ev.loadstore.inherit_sfr_physreg, (void*)(Waddr)ev.loadstore.inherit_sfr_rip, ev.rip);
    break;
  }
  case EVENT_STORE_ISSUED: {
    out = std::format_to(out, "store{} rob {:<3}({:<5}) stq {} r{:<3} on {:<4} @ {} (phys {}): ",
                         (ev.loadstore.load_store_second_phase ? "2" : " "), ev.rob, uopname, ev.lsq, ev.physreg,
                         fu_names[ev.fu], (void*)(Waddr)ev.loadstore.virtaddr,
                         (void*)(Waddr)(ev.loadstore.sfr.physaddr << 3));
    if (ev.loadstore.inherit_sfr_used) {
      out = std::format_to(out, "inherit from {} (uuid {}, rob{}, lsq {}, r{});", ev.loadstore.inherit_sfr,
                           ev.loadstore.inherit_sfr_uuid, ev.loadstore.inherit_sfr_rob, ev.loadstore.inherit_sfr_lsq,
                           ev.loadstore.inherit_sfr_physreg);
    }
    out = std::format_to(out, " <= {} = {}", hexstring(ev.loadstore.data_to_store, 8 * (1 << ev.uop.size)),
                         ev.loadstore.sfr);
    break;
  }
  case EVENT_STORE_LOCK_RELEASED: {
    out = std::format_to(out, "lk-rel rob {:<3}({:<5}) stq {} r{:<3} on {:<4} @ {} (phys {}): ", ev.rob, uopname,
                         ev.lsq, ev.physreg, fu_names[ev.fu], (void*)(Waddr)ev.loadstore.virtaddr,
                         (void*)(Waddr)(ev.loadstore.sfr.physaddr << 3));
    out = std::format_to(out, "lock released (original ld.acq uuid {} rob {} on vcpu {})", ev.loadstore.locking_uuid,
                         ev.loadstore.locking_rob, ev.loadstore.locking_vcpuid);
    break;
  }
  case EVENT_STORE_LOCK_ANNULLED: {
    out = std::format_to(out, "lk-anl rob {:<3}({:<5}) stq {} r{:<3} on {:<4} @ {} (phys {}): ", ev.rob, uopname,
                         ev.lsq, ev.physreg, fu_names[ev.fu], (void*)(Waddr)ev.loadstore.virtaddr,
                         (void*)(Waddr)(ev.loadstore.sfr.physaddr << 3));
    out = std::format_to(out, "lock annulled (original ld.acq uuid {} rob {} on vcpu {})", ev.loadstore.locking_uuid,
                         ev.loadstore.locking_rob, ev.loadstore.locking_vcpuid);
    break;
  }
  case EVENT_STORE_LOCK_REPLAY: {
    out = std::format_to(out, "store{} rob {:<3}({:<5}) stq {} r{:<3} on {:<4} @ {} (phys {}): ",
                         (ev.loadstore.load_store_second_phase ? "2" : " "), ev.rob, uopname, ev.lsq, ev.physreg,
                         fu_names[ev.fu], (void*)(Waddr)ev.loadstore.virtaddr,
                         (void*)(Waddr)(ev.loadstore.sfr.physaddr << 3));
    out = std::format_to(out, "replay because vcpuid {} uop uuid {} has lock", ev.loadstore.locking_vcpuid,
                         ev.loadstore.locking_uuid);
    break;
  }

  case EVENT_LOAD_WAIT: {
    out = std::format_to(out, "{} rob {:<3}({:<5}) ldq {} r{:<3} on {:<4} @ {} (phys {}): ",
                         (ev.loadstore.load_store_second_phase ? "load2 " : "load "), ev.rob, uopname, ev.lsq,
                         ev.physreg, fu_names[ev.fu], (void*)(Waddr)ev.loadstore.virtaddr,
                         (void*)(Waddr)(ev.loadstore.sfr.physaddr << 3));
    out = std::format_to(out, "wait on sfr {} (uuid {}, stq {}, rob {}, r{})", ev.loadstore.inherit_sfr,
                         ev.loadstore.inherit_sfr_uuid, ev.loadstore.inherit_sfr_lsq, ev.loadstore.inherit_sfr_rob,
                         ev.loadstore.inherit_sfr_physreg);
    if (ev.loadstore.predicted_alias)
      out = std::format_to(out, "; stalled by predicted aliasing");
    break;
  }
  case EVENT_LOAD_HIT:
  case EVENT_LOAD_MISS: {
    if (ev.type == EVENT_LOAD_HIT)
      out = std::format_to(out, "{}", (ev.loadstore.load_store_second_phase ? "load2 " : "load  "));
    else
      out = std::format_to(out, "{}", (ev.loadstore.load_store_second_phase ? "ldmis2" : "ldmiss"));

    out = std::format_to(out, " rob {:<3}({:<5}) ldq {} r{:<3} on {:<4} @ {} (phys {}): ", ev.rob, uopname, ev.lsq,
                         ev.physreg, fu_names[ev.fu], (void*)(Waddr)ev.loadstore.virtaddr,
                         (void*)(Waddr)(ev.loadstore.sfr.physaddr << 3));
    if (ev.loadstore.inherit_sfr_used) {
      out = std::format_to(out, "inherit from {} (uuid {}, rob{}, lsq {}, r{}); ", ev.loadstore.inherit_sfr,
                           ev.loadstore.inherit_sfr_uuid, ev.loadstore.inherit_sfr_rob, ev.loadstore.inherit_sfr_lsq,
                           ev.loadstore.inherit_sfr_physreg);
    }
    if (ev.type == EVENT_LOAD_HIT)
      out = std::format_to(out, "hit L1: value 0x{:016x}", ev.loadstore.sfr.data);
    else
      out =
          std::format_to(out, "missed L1 (lfrqslot {}) [value would be 0x{:016x}]", ev.lfrqslot, ev.loadstore.sfr.data);
    break;
  }
  case EVENT_LOAD_BANK_CONFLICT: {
    out = std::format_to(out, "ldbank rob {:<3}({:<5}) ldq {} r{:<3} on {:<4} @ {} (phys {}): ", ev.rob, uopname,
                         ev.lsq, ev.physreg, fu_names[ev.fu], (void*)(Waddr)ev.loadstore.virtaddr,
                         (void*)(Waddr)(ev.loadstore.sfr.physaddr << 3));
    out = std::format_to(out, "L1 bank conflict over bank {}",
                         lowbits(ev.loadstore.sfr.physaddr, log2(CacheSubsystem::L1_DCACHE_BANKS)));
    break;
  }
  case EVENT_LOAD_TLB_MISS: {
    out = std::format_to(out, "{}", (ev.loadstore.load_store_second_phase ? "ldtlb2" : "ldtlb "));
    out = std::format_to(out, " rob {:<3}({:<5}) ldq {} r{:<3} on {:<4} @ {} (phys {}): ", ev.rob, uopname, ev.lsq,
                         ev.physreg, fu_names[ev.fu], (void*)(Waddr)ev.loadstore.virtaddr,
                         (void*)(Waddr)(ev.loadstore.sfr.physaddr << 3));
    if (ev.loadstore.inherit_sfr_used) {
      out = std::format_to(out, "inherit from {} (uuid {}, rob{}, lsq {}, r{}); ", ev.loadstore.inherit_sfr,
                           ev.loadstore.inherit_sfr_uuid, ev.loadstore.inherit_sfr_rob, ev.loadstore.inherit_sfr_lsq,
                           ev.loadstore.inherit_sfr_physreg);
    } else
      out = std::format_to(out, "DTLB miss [value would be 0x{:016x}]", ev.loadstore.sfr.data);
    break;
  }
  case EVENT_LOAD_LOCK_REPLAY: {
    out = std::format_to(out, "{} rob {:<3}({:<5}) ldq {} r{:<3} on {:<4} @ {} (phys {}): ",
                         (ev.loadstore.load_store_second_phase ? "load2 " : "load "), ev.rob, uopname, ev.lsq,
                         ev.physreg, fu_names[ev.fu], (void*)(Waddr)ev.loadstore.virtaddr,
                         (void*)(Waddr)(ev.loadstore.sfr.physaddr << 3));
    out = std::format_to(out, "replay because vcpuid {} uop uuid {} has lock", ev.loadstore.locking_vcpuid,
                         ev.loadstore.locking_uuid);
    break;
  }
  case EVENT_LOAD_LOCK_OVERFLOW: {
    out = std::format_to(out, "{} rob {:<3}({:<5}) ldq {} r{:<3} on {:<4} @ {} (phys {}): ",
                         (ev.loadstore.load_store_second_phase ? "load2 " : "load "), ev.rob, uopname, ev.lsq,
                         ev.physreg, fu_names[ev.fu], (void*)(Waddr)ev.loadstore.virtaddr,
                         (void*)(Waddr)(ev.loadstore.sfr.physaddr << 3));
    out = std::format_to(out, "replay because locking required but no free interlock buffers\n");
    break;
  }
  case EVENT_LOAD_LOCK_ACQUIRED: {
    out = std::format_to(out, "lk-acq rob {:<3}({:<5}) ldq {} r{:<3} on {:<4} @ {} (phys {}): ", ev.rob, uopname,
                         ev.lsq, ev.physreg, fu_names[ev.fu], (void*)(Waddr)ev.loadstore.virtaddr,
                         (void*)(Waddr)(ev.loadstore.sfr.physaddr << 3));
    out = std::format_to(out, "lock acquired");
    break;
  }
  case EVENT_LOAD_LFRQ_FULL:
    out = std::format_to(out, "load   rob {:<3}({:<5}) ldq {} r{:<3}: LFRQ or miss buffer full; replaying", ev.rob,
                         uopname, ev.lsq, ev.physreg);
    break;
  case EVENT_LOAD_HIGH_ANNULLED: {
    out = std::format_to(out, "{} rob {:<3}({:<5}) ldq {} r{:<3} on {:<4} @ {} (phys {}): ",
                         (ev.loadstore.load_store_second_phase ? "load2 " : "load "), ev.rob, uopname, ev.lsq,
                         ev.physreg, fu_names[ev.fu], (void*)(Waddr)ev.loadstore.virtaddr,
                         (void*)(Waddr)(ev.loadstore.sfr.physaddr << 3));
    out = std::format_to(out, "load was annulled (high unaligned load)");
    break;
  }
  case EVENT_LOAD_WAKEUP:
    out = std::format_to(out, "ldwake rob {:<3}({:<5}) ldq {} r{:<3} wakeup load via lfrq slot {}", ev.rob, uopname,
                         ev.lsq, ev.physreg, ev.lfrqslot);
    break;
  case EVENT_TLBWALK_HIT: {
    out =
        std::format_to(out, "wlkhit rob {:<3}({:<5}) ldq {} r{:<3} page table walk (level {}): hit for PTE at phys {}",
                       ev.rob, uopname, ev.lsq, ev.physreg, ev.loadstore.tlb_walk_level, (void*)ev.loadstore.virtaddr);
    break;
    break;
  }
  case EVENT_TLBWALK_MISS: {
    out = std::format_to(
        out, "wlkmis rob {:<3}({:<5}) ldq {} r{:<3} page table walk (level {}): miss for PTE at phys {}: lfrq {}",
        ev.rob, uopname, ev.lsq, ev.physreg, ev.loadstore.tlb_walk_level, (void*)ev.loadstore.virtaddr, ev.lfrqslot);
    break;
    break;
  }
  case EVENT_TLBWALK_WAKEUP: {
    out = std::format_to(
        out,
        "wlkwak rob {:<3}({:<5}) ldq {} r{:<3} page table walk (level {}): wakeup from cache miss for phys {}: lfrq {}",
        ev.rob, uopname, ev.lsq, ev.physreg, ev.loadstore.tlb_walk_level, (void*)ev.loadstore.virtaddr, ev.lfrqslot);
    break;
    break;
  }
  case EVENT_TLBWALK_NO_LFRQ_MB: {
    out = std::format_to(
        out,
        "wlknml rob {:<3}({:<5}) ldq {} r{:<3} page table walk (level {}): no LFRQ or MB for PTE at phys {}: lfrq {}",
        ev.rob, uopname, ev.lsq, ev.physreg, ev.loadstore.tlb_walk_level, (void*)ev.loadstore.virtaddr, ev.lfrqslot);
    break;
    break;
  }
  case EVENT_TLBWALK_COMPLETE: {
    out = std::format_to(out, "wlkhit rob {:<3}({:<5}) ldq {} r{:<3} page table walk (level {}): complete!", ev.rob,
                         uopname, ev.lsq, ev.physreg, ev.loadstore.tlb_walk_level);
    break;
    break;
  }
  case EVENT_LOAD_EXCEPTION: {
    out = std::format_to(out, "{} rob {:<3}({:<5}) stq {} r{:<3} on {:<4} @ {}: exception {}, pfec {}",
                         (ev.loadstore.load_store_second_phase ? "load2 " : "load "), ev.rob, uopname, ev.lsq,
                         ev.physreg, fu_names[ev.fu], (void*)(Waddr)ev.loadstore.virtaddr,
                         exception_name(LO32(ev.loadstore.sfr.data)), PageFaultErrorCode(HI32(ev.loadstore.sfr.data)));
    break;
  }
  case EVENT_STORE_EXCEPTION: {
    out = std::format_to(out, "store{} rob {:<3}({:<5}) stq {} r{:<3} on {:<4} @ {}: exception {}, pfec {}",
                         (ev.loadstore.load_store_second_phase ? "2" : " "), ev.rob, uopname, ev.lsq, ev.physreg,
                         fu_names[ev.fu], (void*)(Waddr)ev.loadstore.virtaddr,
                         exception_name(LO32(ev.loadstore.sfr.data)), PageFaultErrorCode(HI32(ev.loadstore.sfr.data)));
    break;
  }
  case EVENT_ALIGNMENT_FIXUP:
    out = std::format_to(out, "algnfx rip {}: set unaligned bit for uop {} (unaligned predictor slot {}) and refetch",
                         ev.rip, ev.uop.bbindex, OutOfOrderCore::hash_unaligned_predictor_slot(ev.rip));
    break;
  case EVENT_FENCE_ISSUED:
    out = std::format_to(out, "mfence rob {:<3}({:<5}) lsq {} r{:<3}: memory fence ({})", ev.rob, uopname, ev.lsq,
                         ev.physreg, ev.uop);
    break;
  case EVENT_ANNUL_NO_FUTURE_UOPS:
    out = std::format_to(out, "misspc rob {:<3}({:<5}): SOM rob {}, EOM rob {}: no future uops to annul", ev.rob,
                         uopname, ev.annul.somidx, ev.annul.eomidx);
    break;
  case EVENT_ANNUL_MISSPECULATION: {
    out = std::format_to(out, "misspc rob {:<3}({:<5}): SOM rob {}, EOM rob {}: annul from rob {} to rob {}", ev.rob,
                         uopname, ev.annul.somidx, ev.annul.eomidx, ev.annul.startidx, ev.annul.endidx);
    break;
  }
  case EVENT_ANNUL_EACH_ROB: {
    out = std::format_to(out, "annul  rob {:<3}({:<5}): annul rip {}", ev.rob, uopname, ev.rip);
    out = std::format_to(out, "{}", (ev.uop.som ? " SOM" : "    "));
    out = std::format_to(out, "{}", (ev.uop.eom ? " EOM" : "    "));
    out = std::format_to(out, ": free");
    out = std::format_to(out, " r{}", ev.physreg);
    if (ld | st)
      out = std::format_to(out, " lsq{}", ev.lsq);
    if (ev.lfrqslot >= 0)
      out = std::format_to(out, " lfrq{}", ev.lfrqslot);
    if (ev.annul.annulras)
      out = std::format_to(out, " ras");
    break;
  }
  case EVENT_ANNUL_PSEUDOCOMMIT: {
    out = std::format_to(out, "pseucm rob {:<3}({:<5}): r{} rebuild rrt:", ev.rob, uopname, ev.physreg);
    out = std::format_to(out, " arch {}", arch_reg_names[ev.uop.rd]);
    if likely (!ev.uop.nouserflags) {
      if (ev.uop.setflags & SETFLAG_ZF)
        out = std::format_to(out, " zf");
      if (ev.uop.setflags & SETFLAG_CF)
        out = std::format_to(out, " cf");
      if (ev.uop.setflags & SETFLAG_OF)
        out = std::format_to(out, " of");
    }
    out = std::format_to(out, " = r{}", ev.physreg);
    break;
  }
  case EVENT_ANNUL_FETCHQ_RAS:
    out = std::format_to(out, "anlras rip {}: annul RAS update still in fetchq", ev.rip);
    break;
  case EVENT_ANNUL_FLUSH:
    out = std::format_to(out, "flush  rob {:<3}({:<5}) rip {}", ev.rob, uopname, ev.rip);
    break;
  case EVENT_REDISPATCH_DEPENDENTS:
    out = std::format_to(out, "redisp rob {:<3}({:<5}) find all dependents", ev.rob, uopname);
    break;
  case EVENT_REDISPATCH_DEPENDENTS_DONE:
    out = std::format_to(out, "redisp rob {:<3}({:<5}) redispatched {} dependent uops", ev.rob, uopname,
                         (ev.redispatch.count - 1));
    break;
  case EVENT_REDISPATCH_EACH_ROB: {
    out = std::format_to(out, "redisp rob {:<3}({:<5}) from state {}: dep on ", ev.rob, uopname,
                         ev.redispatch.current_state_list->name);
    if (ev.redispatch.dependent_operands.none()) {
      out = std::format_to(out, " [self]");
    } else {
      foreach (i, MAX_OPERANDS) {
        if (ev.redispatch.dependent_operands[i])
          out = std::format_to(out, " {}", ev.redispatch.opinfo[i]);
      }
    }

    out = std::format_to(out, "; redispatch ");
    out = std::format_to(out, " [rob {}]", ev.rob);
    out = std::format_to(out, " [physreg {}]", ev.physreg);
    if (ld | st)
      out = std::format_to(out, " [lsq {}]", ev.lsq);
    if (ev.redispatch.iqslot)
      out = std::format_to(out, " [iqslot]");
    if (ev.lfrqslot >= 0)
      out = std::format_to(out, " [lfrqslot {}]", ev.lfrqslot);
    if (ev.redispatch.opinfo[RS].physreg != PHYS_REG_NULL)
      out = std::format_to(out, " [inheritsfr {}]", ev.redispatch.opinfo[RS]);

    break;
  }
  case EVENT_COMPLETE:
    out = std::format_to(out, "complt rob {:<3}({:<5}) on {:<4}: r{:<3}", ev.rob, uopname, fu_names[ev.fu], ev.physreg);
    break;
  case EVENT_FORWARD: {
    out = std::format_to(
        out, "forwd{} rob {:<3}({:<5}) ({}) r{:<3} => uuid {} rob {} ({}) r{} operand {}", ev.forwarding.forward_cycle,
        ev.rob, uopname, clusters[ev.cluster].name, ev.physreg, ev.forwarding.target_uuid, ev.forwarding.target_rob,
        clusters[ev.forwarding.target_cluster].name, ev.forwarding.target_physreg, ev.forwarding.operand);
    if (ev.forwarding.target_st)
      out = std::format_to(out, " => st{}", ev.forwarding.target_lsq);
    out = std::format_to(out, " [still waiting?");
    foreach (i, MAX_OPERANDS) {
      if (!bit(ev.forwarding.target_operands_ready, i))
        out = std::format_to(out, " r{}", (char)('a' + i));
    }
    if (ev.forwarding.target_all_operands_ready)
      out = std::format_to(out, " READY");
    out = std::format_to(out, "]");
    break;
  }
  case EVENT_BROADCAST: {
    out = std::format_to(out, "brcst{} rob {:<3}({:<5}) from cluster {} to cluster {} on forwarding cycle {}",
                         ev.forwarding.forward_cycle, ev.rob, uopname, clusters[ev.cluster].name,
                         clusters[ev.forwarding.target_cluster].name, ev.forwarding.forward_cycle);
    break;
  }
  case EVENT_WRITEBACK: {
    out = std::format_to(out, "write  rob {:<3}({:<5}) (cluster {}) r{:<3}@{} = 0x{:016x} {}", ev.rob, uopname,
                         clusters[ev.cluster].name, ev.physreg, phys_reg_file_names[ev.rfid], ev.writeback.data,
                         flagstring(ev.writeback.flags));
    if (ev.writeback.transient)
      out = std::format_to(out, " (transient)");
    out = std::format_to(out, " ({} consumers", ev.writeback.consumer_count);
    if (ev.writeback.all_consumers_sourced_from_bypass)
      out = std::format_to(out, ", all from bypass");
    if (ev.writeback.no_branches_between_renamings)
      out = std::format_to(out, ", no intervening branches");
    if (ev.writeback.dest_renamed_before_writeback)
      out = std::format_to(out, ", dest renamed before writeback");
    out = std::format_to(out, ")");
    break;
  }
  case EVENT_COMMIT_FENCE_COMPLETED:
    out = std::format_to(out, "mfcmit rob {:<3}({:<5}) fence committed: wake up waiting memory uops", ev.rob, uopname);
    break;
  case EVENT_COMMIT_EXCEPTION_DETECTED:
    out = std::format_to(out, "detect rob {:<3}({:<5}) exception {} ({}), error code {:04x}, origvirt {}", ev.rob,
                         uopname, exception_name(LO32(ev.commit.state.reg.rddata)), LO32(ev.commit.state.reg.rddata),
                         HI32(ev.commit.state.reg.rddata), (void*)(Waddr)ev.commit.origvirt);
    break;
  case EVENT_COMMIT_EXCEPTION_ACKNOWLEDGED:
    out = std::format_to(out, "except rob {:<3}({:<5}) exception {} [EOM #{}]", ev.rob, uopname,
                         exception_name(LO32(ev.commit.state.reg.rddata)), ev.commit.total_user_insns_committed);
    break;
  case EVENT_COMMIT_SKIPBLOCK:
    out = std::format_to(out, "skipbk rob {:<3}({:<5}) skip block: advance rip by {} to {} [EOM #{}]", ev.rob, uopname,
                         ev.uop.bytes, (void*)(Waddr)(ev.rip.rip + ev.uop.bytes), ev.commit.total_user_insns_committed);
    break;
  case EVENT_COMMIT_SMC_DETECTED:
    out = std::format_to(out,
                         "smcdet rob {:<3}({:<5}) self-modifying code at rip {} detected (mfn was dirty); invalidate "
                         "and retry [EOM #{}]",
                         ev.rob, uopname, ev.rip, ev.commit.total_user_insns_committed);
    break;
  case EVENT_COMMIT_MEM_LOCKED:
    out = std::format_to(out, "waitlk rob {:<3}({:<5}) wait for lock on physaddr {} to be released", ev.rob, uopname,
                         (void*)(ev.commit.state.st.physaddr << 3));
    break;
  case EVENT_COMMIT_OK: {
    out = std::format_to(out, "commit rob {:<3}({:<5})", ev.rob, uopname);
    if likely (archdest_can_commit[ev.uop.rd])
      out = std::format_to(out, " [rrt {} = r{} 0x{:016x}]", arch_reg_names[ev.uop.rd], ev.physreg,
                           ev.commit.state.reg.rddata);

    if ((!ev.uop.nouserflags) && ev.uop.setflags) {
      out = std::format_to(out, " [flags {}{}{} -> {}]", ((ev.uop.setflags & SETFLAG_ZF) ? "z" : ""),
                           ((ev.uop.setflags & SETFLAG_CF) ? "c" : ""), ((ev.uop.setflags & SETFLAG_OF) ? "o" : ""),
                           flagstring(ev.commit.state.reg.rdflags));
    }

    if (ev.uop.eom)
      out = std::format_to(out, " [rip = {}]", (void*)(Waddr)ev.commit.target_rip);

    if unlikely (st && (ev.commit.state.st.bytemask != 0))
      out =
          std::format_to(out, " [mem {} = {} mask {}]", (void*)(Waddr)(ev.commit.state.st.physaddr << 3),
                         bytemaskstring((const unsigned char*)&ev.commit.state.st.data, ev.commit.state.st.bytemask, 8),
                         bitstring(ev.commit.state.st.bytemask, 8, true));

    if unlikely (ev.commit.pteupdate.a | ev.commit.pteupdate.d | ev.commit.pteupdate.ptwrite) {
      out = std::format_to(out, " [pte:");
      if (ev.commit.pteupdate.a)
        out = std::format_to(out, " a");
      if (ev.commit.pteupdate.d)
        out = std::format_to(out, " d");
      if (ev.commit.pteupdate.ptwrite)
        out = std::format_to(out, " w");
      out = std::format_to(out, "]");
    }

    if unlikely (ld | st) {
      out = std::format_to(out, " [lsq {}]", ev.lsq);
      out = std::format_to(out, " [upslot {} = {}]", OutOfOrderCore::hash_unaligned_predictor_slot(ev.rip),
                           ev.commit.ld_st_truly_unaligned);
    }

    if likely (ev.commit.oldphysreg > 0) {
      if unlikely (ev.commit.oldphysreg_refcount) {
        out = std::format_to(out, " [pending free old r{} ref by", ev.commit.oldphysreg);
        out = std::format_to(out, " refcount {}", ev.commit.oldphysreg_refcount);
        out = std::format_to(out, "]");
      } else {
        out = std::format_to(out, " [free old r{}]", ev.commit.oldphysreg);
      }
    }

    out = std::format_to(out, " [commit r{}]", ev.physreg);

    foreach (i, MAX_OPERANDS) {
      if unlikely (ev.commit.operand_physregs[i] != PHYS_REG_NULL)
        out = std::format_to(out, " [unref r{}]", ev.commit.operand_physregs[i]);
    }

    if unlikely (br) {
      out = std::format_to(out, " [brupdate{}{} {}]", (ev.commit.taken ? " tk" : " nt"),
                           (ev.commit.predtaken ? " pt" : " np"),
                           ((ev.commit.taken == ev.commit.predtaken) ? " ok" : " MP"));
    }

    if (ev.uop.eom)
      out = std::format_to(out, " [EOM #{}]", ev.commit.total_user_insns_committed);
    break;
  }
  case EVENT_COMMIT_ASSIST: {
    out = std::format_to(out, "assist rob {:<3} calling assist {} (#{}: {})", ev.rob, (void*)ev.rip.rip,
                         assist_index((assist_func_t)ev.rip.rip), assist_name((assist_func_t)ev.rip.rip));
    break;
  }
  case EVENT_RECLAIM_PHYSREG:
    out = std::format_to(out, "free   r{} no longer referenced; moving to free state", ev.physreg);
    break;
  case EVENT_RELEASE_MEM_LOCK: {
    out = std::format_to(out, "unlkcm phys {}: lock release committed", (void*)(ev.loadstore.sfr.physaddr << 3));
    break;
  }
  default:
    out = std::format_to(out, "?????? unknown event type {}", ev.type);
    break;
  }

  out = std::format_to(out, "\n");
  return out;
}

namespace x86sim {

//
// Construct all the structures necessary to configure
// the cores. This function is only called once, after
// all other PTLsim subsystems are brought up.
//

std::string_view OutOfOrderMachine::name() const {
  return "ooo";
}

OutOfOrderMachine::~OutOfOrderMachine() = default;

OutOfOrderMachine::OutOfOrderMachine(Machine& machine, Options config) : MachineImpl(machine, std::move(config)) {
  // Note: we only create a single core for all contexts for now.
  cores[0] = std::make_unique<OutOfOrderCore>(0, *this);

  foreach (i, this->config.core_count) {
    OutOfOrderCore& core = *cores[0];
    core.threadcount++;
    core.threads[i] = std::make_unique<ThreadContext>(core, i);
    core.threads[i]->init();

    //
    // Note: in a multi-processor model, config may
    // specify various ways of slicing contextcount up
    // into threads, cores and sockets; the appropriate
    // interconnect and cache hierarchy parameters may
    // be specified here.
    //
  }

  cores[0]->init();
  init_luts();
}

RegisterFile& OutOfOrderMachine::register_file(std::size_t core_index) noexcept {
  return cores[0]->threads[core_index]->ctx;
}

const RegisterFile& OutOfOrderMachine::register_file(std::size_t core_index) const noexcept {
  return cores[0]->threads[core_index]->ctx;
}

//
// Run the processor model, until a stopping point
// is hit (as configured elsewhere in config).
//
int OutOfOrderMachine::run() {
  logging::println("Starting out-of-order core toplevel loop");
  logging::flush();

  // All VCPUs are running:
  stopped = 0;

  if unlikely (this->iterations >= config.log.start_log_at_iteration) {
    logging::println("Start logging at level {} in cycle {}", config.log.loglevel, this->iterations);
    logging::flush();

    logenable = 1;
  }

  cores[0]->reset();
  cores[0]->flush_pipeline_all();

  logging::println("IssueQueue states:");

  if unlikely (config.log.event_log_enabled && (!cores[0]->eventlog.start)) {
    cores[0]->eventlog.init(config.log.event_log_ring_buffer_size);
  }

  bool exiting = false;
  bool stopping = false;

  for (;;) {
    if unlikely (this->iterations >= config.log.start_log_at_iteration) {
      if unlikely (!logenable) {
        logging::println("Start logging at level {} in cycle {}", config.log.loglevel, this->iterations);
        logging::flush();
      }
      logenable = 1;
    }

    OutOfOrderCore& core = *cores[0]; // only one core for now
    int running_thread_count = 0;
    foreach (i, core.threadcount) {
      auto& thread = core.threads[i];
    }

    exiting |= core.runcycle();

    if unlikely (((this->iterations >= config.stop_at_iteration ||
                   this->total_user_insns_committed >= config.stop_at_user_insns) &&
                  (!stopping))) {
      logging::println("Waiting for all VCPUs to reach stopping point, starting at cycle {}", this->sim_cycle);
      // force_logging_enabled();
      OutOfOrderCore& core = *cores[0];
      foreach (i, core.threadcount)
        core.threads[i]->stop_at_next_eom = 1;
      if (config.debug.abort_at_end) {
        config.debug.abort_at_end = 0;
        logging::println("Abort immediately: do not wait for next x86 boundary nor flush pipelines");
        stopped = 1;
        exiting = 1;
      }
      stopping = 1;
    }

    stats.summary.cycles++;
    stats.ooocore.cycles++;
    this->sim_cycle++;
    unhalted_cycle_count += (running_thread_count > 0);
    this->iterations++;

    if unlikely (stopping) {
      logging::println(logging::TRACE, "Waiting for all VCPUs to stop at {}: mask = {} (need {} VCPUs)",
                       this->sim_cycle, stopped.to_string(), contextcount);
      exiting |= (stopped.count() == contextcount);
    }

    if unlikely (exiting)
      break;
  }

  logging::println("Exiting out-of-order core at {} commits, {} uops and {} iterations (cycles)",
                   this->total_user_insns_committed, this->total_uops_committed, this->iterations);

  OutOfOrderCore& core = *cores[0]; /// only one core for now.

  foreach (i, core.threadcount) {
    auto& thread = core.threads[i];

    thread->core_to_external_state();

    if (((this->sim_cycle - thread->last_commit_at_cycle) > 1024) | config.debug.dump_state_now) {
      logging::println(logging::TRACE, "Core State at end for thread {}:", thread->threadid);
      logging::println(logging::TRACE, "{}", thread->ctx);
    }
  }

  config.debug.dump_state_now = 0;

  dump_state();

  // Flush everything to remove any remaining refs to basic blocks
  flush_all_pipelines();

  return exiting;
}

void OutOfOrderCore::flush_tlb(Context& ctx, int threadid, bool selective, Waddr virtaddr) {
  ThreadContext& thread = *threads[threadid];

  logging::print(logging::DEBUG, "[vcpu {}] core {}, thread {}: Flush TLBs", ctx.vcpuid, coreid, threadid);
  if (selective)
    logging::println(logging::DEBUG, " for virtaddr {}", (void*)virtaddr);
  logging::println(logging::DEBUG, "");
  logging::println(logging::TRACE, "DTLB before:");
  logging::println(logging::TRACE, "{}", caches.dtlb);
  logging::println(logging::TRACE, "ITLB before:");
  logging::println(logging::TRACE, "{}", caches.itlb);

  int dn;
  int in;

  if unlikely (selective) {
    dn = caches.dtlb.flush_virt(virtaddr, threadid);
    in = caches.itlb.flush_virt(virtaddr, threadid);
  } else {
    dn = caches.dtlb.flush_thread(threadid);
    in = caches.itlb.flush_thread(threadid);
  }

  logging::println(logging::DEBUG, "Flushed {} DTLB slots and {} ITLB slots", dn, in);
  logging::println(logging::TRACE, "DTLB after:");
  logging::println(logging::TRACE, "{}", caches.dtlb);
  logging::println(logging::TRACE, "ITLB after:");
  logging::println(logging::TRACE, "{}", caches.itlb);
}

void OutOfOrderMachine::flush_tlb(Context& ctx) {
  // This assumes all VCPUs are mapped as threads in a single SMT core
  int coreid = 0;
  int threadid = ctx.vcpuid;
  cores[coreid]->flush_tlb(ctx, threadid);
}

void OutOfOrderMachine::flush_tlb_virt(Context& ctx, Waddr virtaddr) {
  // This assumes all VCPUs are mapped as threads in a single SMT core
  int coreid = 0;
  int threadid = ctx.vcpuid;
  cores[coreid]->flush_tlb(ctx, threadid, true, virtaddr);
}

void OutOfOrderMachine::dump_state() {
  logging::println(logging::INFO, "dump_state include event if -ringbuf enabled:");
  foreach (i, MAX_SMT_CORES) {
    logging::println(logging::INFO, "dump_state for core {}", i);
    logging::flush();
    if (!cores[i])
      continue;
    OutOfOrderCore& core = *cores[i];
    if unlikely (config.log.event_log_enabled)
      core.eventlog.print();
    else
      logging::println(logging::INFO, "config.log.event_log_enabled is not enabled: {}", config.log.event_log_enabled);

    core.dump_smt_state();
    core.print_smt_state();
  }
  logging::println(logging::INFO, "Memory interlock buffer:");
  logging::flush();
  logging::println(logging::INFO, "{}", interlocks);
  //
  // For debugging only:
  //
  foreach (i, cores[0]->threadcount) {
    auto& thread = cores[0]->threads[i];
    logging::println(logging::TRACE, "Thread {}:", i);
    logging::println(logging::TRACE, "  rip:                                 {}",
                     (void*)thread->ctx.commitarf[REG_rip]);
    logging::println(logging::TRACE, "  consecutive_commits_inside_spinlock: {}",
                     thread->consecutive_commits_inside_spinlock);
    logging::println(logging::TRACE, "  State:");
    logging::println(logging::TRACE, "{}", thread->ctx);
  }
}

// Stub implementation for dump_smt_state - prints core state using logging
void OutOfOrderCore::dump_smt_state() {
  logging::println(logging::INFO, "Core {} SMT state:", coreid);
  foreach (i, threadcount) {
    if (threads[i]) {
      logging::println(logging::INFO, "  Thread {}: ROB count {}, LSQ count {}", i, threads[i]->ROB.count,
                       threads[i]->LSQ.count);
    }
  }
}

void OutOfOrderMachine::update_stats(PTLsimStats& stats) {
  foreach (vcpuid, contextcount) {
    PerContextOutOfOrderCoreStats& s = per_context_ooocore_stats_ref(vcpuid);
    s.issue.uipc = s.issue.uops / (double)stats.ooocore.cycles;
    s.commit.uipc = (double)s.commit.uops / (double)stats.ooocore.cycles;
    s.commit.ipc = (double)s.commit.insns / (double)stats.ooocore.cycles;
  }

  PerContextOutOfOrderCoreStats& s = stats.ooocore.total;
  s.issue.uipc = s.issue.uops / (double)stats.ooocore.cycles;
  s.commit.uipc = (double)s.commit.uops / (double)stats.ooocore.cycles;
  s.commit.ipc = (double)s.commit.insns / (double)stats.ooocore.cycles;
}

//
// Flush all pipelines in every core, and process any
// pending BB cache invalidates.
//
// Typically this is in response to some infrequent event
// like cross-modifying SMC or cache coherence deadlocks.
//
void OutOfOrderMachine::flush_all_pipelines() {
  assert(cores[0]);
  OutOfOrderCore* core = cores[0].get();

  //
  // Make sure all pipelines are flushed BEFORE
  // we try to invalidate the dirty page!
  // Otherwise there will still be some remaining
  // references to to the basic block
  //
  core->flush_pipeline_all();

  foreach (i, core->threadcount) {
    auto& thread = core->threads[i];
    thread->invalidate_smc();
  }
}

// Formatter implementations
} // namespace x86sim

auto std::formatter<x86sim::OutOfOrderModel::PhysicalRegister>::format(
    const x86sim::OutOfOrderModel::PhysicalRegister& physreg, std::format_context& ctx) const {
  using namespace x86sim;
  using namespace OutOfOrderModel;
  std::string vf = format_value_and_flags(physreg.data, physreg.flags);
  auto out = std::format_to(ctx.out(), "TH {} rfid {}", physreg.threadid, physreg.rfid);
  out = std::format_to(out, "  r{:<3} state {:<12} {}", physreg.index(), physreg.get_state_list().name, vf);
  if (physreg.rob)
    out = std::format_to(out, " rob {} (uuid {})", physreg.rob->index(), physreg.rob->uop.uuid);
  out = std::format_to(out, " refcount {}", physreg.refcount);
  return out;
}


auto std::formatter<x86sim::OutOfOrderModel::PhysicalRegisterFile>::format(
    const x86sim::OutOfOrderModel::PhysicalRegisterFile& prf, std::format_context& ctx) const {
  using namespace x86sim;
  auto out = std::format_to(ctx.out(), "PhysicalRegisterFile<{}, rfid {}, size {}>:\n", prf.name, prf.rfid, prf.size);
  for (int i = 0; i < prf.size; i++) {
    out = std::format_to(out, "{}\n", prf[i]);
  }
  return out;
}


auto std::formatter<x86sim::OutOfOrderModel::ReorderBufferEntry>::format(
    const x86sim::OutOfOrderModel::ReorderBufferEntry& rob, std::format_context& ctx) const {
  using namespace x86sim;
  using namespace OutOfOrderModel;
  std::string name = nameof(rob.uop);
  std::string rainfo = rob.get_operand_info(0);
  std::string rbinfo = rob.get_operand_info(1);
  std::string rcinfo = rob.get_operand_info(2);

  auto out = std::format_to(
      ctx.out(), "rob {:<3} uuid {:>16} rip 0x{:x} {:<24} {} {} @ {:<4} {:<12} r{:<3} {:<6}", rob.index(), rob.uop.uuid,
      (W64)rob.uop.rip, rob.current_state_list->name, (rob.uop.som ? "SOM" : "   "), (rob.uop.eom ? "EOM" : "   "),
      (rob.cluster >= 0) ? clusters[rob.cluster].name : "???", name, rob.physreg->index(), arch_reg_names[rob.uop.rd]);

  if (isload(rob.uop.opcode))
    out = std::format_to(out, " ld{:<3}", rob.lsq->index());
  else if (isstore(rob.uop.opcode))
    out = std::format_to(out, " st{:<3}", rob.lsq->index());
  else
    out = std::format_to(out, "      ");

  out = std::format_to(out, " = ");
  out = std::format_to(out, "{:<30}{:<30}{:<30}", rainfo, rbinfo, rcinfo);
  return out;
}


auto std::formatter<x86sim::OutOfOrderModel::LoadStoreQueueEntry>::format(
    const x86sim::OutOfOrderModel::LoadStoreQueueEntry& lsq, std::format_context& ctx) const {
  using namespace x86sim;
  using namespace OutOfOrderModel;
  auto out = std::format_to(ctx.out(), "{}{:<3} uuid {:>10} rob {:<3} r{:<3}", (lsq.store ? "st" : "ld"), lsq.index(),
                            lsq.rob->uop.uuid, lsq.rob->index(), lsq.rob->physreg->index());

  if (PHYS_REG_FILE_COUNT > 1)
    out = std::format_to(out, "@{}", lsq.getcore().physregfiles[lsq.rob->physreg->rfid].name);

  out = std::format_to(out, " ");

  if (lsq.invalid) {
    out = std::format_to(out, "< Invalid: fault 0x{:02x} > ", lsq.data);
  } else {
    if (lsq.datavalid)
      out = std::format_to(out, "{}", bytemaskstring((const unsigned char*)&lsq.data, lsq.bytemask, 8));
    else
      out = std::format_to(out, "<    Data Invalid     >");
    out = std::format_to(out, " @ ");
    if (lsq.addrvalid)
      out = std::format_to(out, "0x{:012x}", lsq.physaddr << 3);
    else
      out = std::format_to(out, "< Addr Inval >");
  }
  return out;
}
