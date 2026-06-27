//
// PTLsim: Cycle Accurate x86-64 Simulator
// L1 and L2 Data Caches
//
// Copyright 2005-2008 Matt T. Yourst <yourst@yourst.com>
//

#include "dcache.h"
#include "globals.h"
#include "ptlhwdef.h"
#include "ptlsim-api.h"
#include "stats.h"
#include "typedefs.h"
#include "x86sim/logging.hpp"
#include <bitset>
#include <cassert>
#include <format>
#include <utility>

namespace x86sim {

using namespace CacheSubsystem;

#if 0
#define starttimer(timer) timer.start()
#define stoptimer(timer) timer.stop()
#else
#define starttimer(timer) (1)
#define stoptimer(timer) (1)
#endif

#ifdef TRACK_LINE_USAGE
// Lifetime
template<>
W64 L1StatsCollectorBase::line_lifetime_histogram[DCACHE_L1_LINE_LIFETIME_SLOTS] = {};
template<>
W64 L1IStatsCollectorBase::line_lifetime_histogram[DCACHE_L1I_LINE_LIFETIME_SLOTS] = {};
template<>
W64 L2StatsCollectorBase::line_lifetime_histogram[DCACHE_L2_LINE_LIFETIME_SLOTS] = {};
template<>
W64 L3StatsCollectorBase::line_lifetime_histogram[DCACHE_L3_LINE_LIFETIME_SLOTS] = {};

// Deadtime
template<>
W64 L1StatsCollectorBase::line_deadtime_histogram[DCACHE_L1_LINE_DEADTIME_SLOTS] = {};
template<>
W64 L1IStatsCollectorBase::line_deadtime_histogram[DCACHE_L1I_LINE_DEADTIME_SLOTS] = {};
template<>
W64 L2StatsCollectorBase::line_deadtime_histogram[DCACHE_L2_LINE_DEADTIME_SLOTS] = {};
template<>
W64 L3StatsCollectorBase::line_deadtime_histogram[DCACHE_L3_LINE_DEADTIME_SLOTS] = {};

// Hit count
template<>
W64 L1StatsCollectorBase::line_hitcount_histogram[DCACHE_L1_LINE_HITCOUNT_SLOTS] = {};
template<>
W64 L1IStatsCollectorBase::line_hitcount_histogram[DCACHE_L1I_LINE_HITCOUNT_SLOTS] = {};
template<>
W64 L2StatsCollectorBase::line_hitcount_histogram[DCACHE_L2_LINE_HITCOUNT_SLOTS] = {};
template<>
W64 L3StatsCollectorBase::line_hitcount_histogram[DCACHE_L3_LINE_HITCOUNT_SLOTS] = {};
#endif

//
// Load Fill Request Queue
//

template<int size>
void LoadFillReqQueue<size>::restart() {
  while (!freemap.all()) {
    int idx = lsb(~freemap);
    LoadFillReq& req = reqs[idx];
    logging::println(logging::DEBUG, "iter {}: force final wakeup/reset of LFRQ slot {}", machine.iterations, idx);
    annul(idx);
  }
  reset();
  stats.dcache.lfrq.resets++;
}

template<int size>
void LoadFillReqQueue<size>::reset(int threadid) {
  foreach (i, SIZE) {
    LoadFillReq& req = reqs[i];
    if likely ((!freemap[i]) && (req.lsi.threadid == threadid)) {
      logging::println(logging::DEBUG, "[vcpu {}] reset lfrq slot {}", threadid, i);
      waiting[i] = 0;
      ready[i] = 0;
      freemap[i] = 1;
      count--;
      assert(count >= 0);
    }
  }

  stats.dcache.lfrq.resets++;
}

template<int size>
void LoadFillReqQueue<size>::annul(int lfrqslot) {
  LoadFillReq& req = reqs[lfrqslot];
  logging::println(logging::TRACE, "  Annul LFRQ slot {}", lfrqslot);
  stats.dcache.lfrq.annuls++;
  hierarchy->missbuf.annul_lfrq(lfrqslot);
  reqs[lfrqslot].mbidx = -1;
  assert(!freemap[lfrqslot]);
  changestate(lfrqslot, ready, freemap);
  count--;
  assert(count >= 0);
}

//
// Add an entry to the LFRQ in the waiting state.
//
template<int size>
int LoadFillReqQueue<size>::add(const LoadFillReq& req) {
  if unlikely (full())
    return -1;
#ifndef NDEBUG
  // Sanity check: make sure (tid, rob) is unique:
  foreach (i, size) {
    if likely (freemap[i])
      continue;
    const LoadFillReq& old = reqs[i];
    if ((old.lsi.threadid == req.lsi.threadid) && (old.lsi.rob == req.lsi.rob)) {
      logging::println(logging::ERROR, "ERROR: during add LFRQ req {}, entry {} ({}) already matches at cycle {}", req,
                       i, old, machine.sim_cycle);
      logging::println(logging::ERROR, "{}", *this);
      logging::println(logging::ERROR, "{}", hierarchy->missbuf);
      // assert(false);
    }
  }
#endif
  int idx = lsb(freemap);
  changestate(idx, freemap, waiting);
  reqs[idx] = req;
  assert(count < size);
  count++;
  stats.dcache.lfrq.inserts++;
  return idx;
}

//
// Move any LFRQ entries in <mask> to the ready state
// in response to the arrival of the corresponding
// line at the L1 level. Once a line is delivered,
// it is copied into the L1 cache and the corresponding
// miss buffer can be freed.
//
template<int size>
void LoadFillReqQueue<size>::wakeup(W64 address, const std::bitset<LFRQ_SIZE>& lfrqmask) {
  logging::println(logging::DEBUG, "LFRQ.wakeup({})", (void*)(Waddr)address);
  //assert(L2.probe(address));
  waiting &= ~lfrqmask;
  ready |= lfrqmask;
}

//
// Find the first N requests (N = 2) in the READY state,
// and extract, sign extend and write into their target
// register, then mark that register as ready.
//
// Also mark the entire cache line containing each load
// as fully valid.
//
// Loads will always be allocated a physical register
// since if the load misses the L1, it will have fallen
// off the end of the pipeline and into the register file
// by the earliest time we can receive the data from the
// L2 cache and/or lower levels.
//
template<int size>
void LoadFillReqQueue<size>::clock() {
  //
  // Process up to MAX_WAKEUPS_PER_CYCLE missed loads per cycle:
  //
  int wakeupcount = 0;
  foreach (i, MAX_WAKEUPS_PER_CYCLE) {
    if unlikely (ready.none())
      break;

    int idx = lsb(ready);
    LoadFillReq& req = reqs[idx];

    logging::println(logging::DEBUG, "[vcpu {}] at cycle {}: wakeup LFRQ slot {}", req.lsi.threadid, machine.sim_cycle,
                     idx);

    W64 delta = LO32(machine.sim_cycle) - LO32(req.initcycle);
    if unlikely (delta >= 65536) {
      // avoid overflow induced erroneous values:
      logging::println(
          logging::WARNING,
          "LFRQ: warning: cycle counter wraparound in initcycle latency (current {} vs init {} = delta {})",
          machine.sim_cycle, req.initcycle, delta);
    } else {
      stats.dcache.lfrq.total_latency += delta;
    }

    stats.dcache.lfrq.wakeups++;
    wakeupcount++;
    if likely (hierarchy->callback)
      hierarchy->callback->dcache_wakeup(req.lsi, req.addr);

    assert(!freemap[idx]);
    changestate(idx, ready, freemap);
    count--;
    assert(count >= 0);
  }

  stats.dcache.lfrq.width[wakeupcount]++;
}

LoadFillReq::LoadFillReq(W64 addr, W64 data, byte mask, LoadStoreInfo lsi, W64 initcycle) {
  this->addr = addr;
  this->data = data;
  this->mask = mask;
  this->lsi = lsi;
  this->lsi.threadid = lsi.threadid;
  this->fillL1 = 1;
  this->fillL2 = 1;
  this->initcycle = initcycle;
  this->mbidx = -1;
}

} // namespace x86sim

auto std::formatter<x86sim::CacheSubsystem::LoadFillReq>::format(const x86sim::CacheSubsystem::LoadFillReq& req,
                                                                 std::format_context& ctx) const {
  using namespace x86sim;
  auto out = ctx.out();
  out = std::format_to(out, "0x{:016x} @ {} -> rob {} @ t{}", req.data, (void*)(Waddr)req.addr, req.lsi.rob,
                       req.lsi.threadid);
  out = std::format_to(out, ": shift {}, signext {}, mask {}", req.lsi.sizeshift, req.lsi.signext,
                       bitstring(req.mask, 8, true));
  return out;
}

namespace x86sim {

//
// Miss Buffer
//

template<int SIZE>
void MissBuffer<SIZE>::reset() {
  foreach (i, SIZE) {
    missbufs[i].reset();
  }
  freemap.set();
  count = 0;
}


template<int SIZE>
void MissBuffer<SIZE>::reset(int threadid) {
  foreach (i, SIZE) {
    Entry& mb = missbufs[i];
    // NOTE SD: This check is broken. A MBE may be shared by LFRQs from different threads.
#if (0)
    if likely (mb.threadid == threadid) {
      logging::println(logging::DEBUG, "[vcpu {}] reset missbuf slot {}: for rob{}", threadid, i, mb.rob);
      assert(!freemap[i]);
      mb.reset();
      freemap[i] = 1;
      count--;
      assert(count >= 0);

      //
      // If multiple threads depend on the same missbuf but one thread is
      // flushed, we'll wake up a stale LFRQ. We have to make sure after
      // a missbuf reset, all the entries point to a valid lfrqmap.
      //
      // NOTE SD: mb.lfrqmap should be zero by now (mb.reset() above)
      if (*mb.lfrqmap) {
        std::bitset<LFRQ_SIZE> tmp_lfrqmap = mb.lfrqmap ^ hierarchy->lfrq.waiting;
        if (*tmp_lfrqmap) {
          logging::println(logging::DEBUG, "Multithread share same missbufs[{}] : lfrqmap conflict detected", i);
          mb.lfrqmap &= ~tmp_lfrqmap;
          logging::println(logging::DEBUG, "after remove stale lfrq entries for missbuf[{}]", i);
        }
        // NB SD: This is the same as mb.lfrqmap &= hierarchy->lfrq.waiting; which actually makes sense :)
      }
    }
#endif
    if (freemap[i])
      continue;

    logging::println(logging::DEBUG, "Adjusting LFR wakeups for missbuf[{}] : before", i);
    if (mb.lfrqmap.any())
      for (size_t l = 0; l < mb.lfrqmap.size(); ++l) {
        if (!mb.lfrqmap[l])
          continue;

        // Go through all associated LFRs and check their thread ID
        // NOTE: This could also be done in LFRQ::reset or through several
        //       assumptions about the waiting bitmap, as above.
        const LoadFillReq& lfr = hierarchy->lfrq.reqs[l];
        if (lfr.lsi.threadid == threadid)
          mb.lfrqmap[l] = 0;
      }
    logging::println(logging::DEBUG, "Adjusting LFR wakeups for missbuf[{}] : after", i);

    if likely (mb.lfrqmap.none() && (mb.threadid == threadid)) {
      // Drop empty MBEs that had only wakeups for the flushed thread
      logging::println(logging::DEBUG, "[vcpu {}] reset missbuf slot {}: for rob{}", threadid, i, mb.rob);
      assert(!freemap[i]);
      mb.reset();
      freemap[i] = 1;
      count--;
      assert(count >= 0);
    }
  }
}

template<int SIZE>
void MissBuffer<SIZE>::restart() {
  if likely (!(freemap.all())) {
    foreach (i, SIZE) {
      missbufs[i].lfrqmap = 0;
    }
  }
}

template<int SIZE>
int MissBuffer<SIZE>::find(W64 addr) {
  W64 match = 0;
  foreach (i, SIZE) {
    if ((missbufs[i].addr == addr) && !freemap[i])
      return i;
  }
  return -1;
}

//
// Request fully or partially missed both the L2 and L1
// caches and needs service from below.
//
template<int SIZE>
int MissBuffer<SIZE>::initiate_miss(W64 addr, bool hit_in_L2, bool icache, int rob, int threadid) {
  addr = floor(addr, L1_LINE_SIZE);

  int idx = find(addr);

  // NOTE SD: This is a fundamental question here: can MBEs be shared by two threads?
  //          Care has to betaken with respect to reset / annullment if so.
  //          Further note that prefetches and store commits have their own bogus thread IDs
  // if unlikely (idx >= 0 && threadid == missbufs[idx].threadid) {
  if unlikely (idx >= 0) {
    // Handle case where dcache miss is already in progress but some
    // code needed in icache is also stored in that line:
    Entry& mb = missbufs[idx];
    mb.icache |= icache;
    mb.dcache |= (!icache);
    // Handle case where icache miss is already in progress but some
    // data needed in dcache is also stored in that line:
    logging::println(logging::DEBUG, "[vcpu {}] miss buffer hit for address {}: returning old slot {}", threadid,
                     (void*)(Waddr)addr, idx);
    return idx;
  }

  if unlikely (full()) {
    logging::println(logging::DEBUG, "[vcpu {}] miss buffer full while allocating slot for address {}", threadid,
                     (void*)(Waddr)addr);
    return -1;
  }

  idx = lsb(freemap);
  freemap[idx] = 0;
  assert(count < SIZE);
  count++;

  stats.dcache.missbuf.inserts++;
  Entry& mb = missbufs[idx];
  mb.addr = addr;
  mb.lfrqmap = 0;
  mb.icache = icache;
  mb.dcache = (!icache);
  mb.rob = rob;
  mb.threadid = threadid;

  logging::println(logging::DEBUG, "[vcpu {}] mb{}: allocated for address {} (iter {})", mb.threadid, idx,
                   (void*)(Waddr)addr, machine.iterations);

  if likely (hit_in_L2) {
    logging::println(logging::DEBUG, "[vcpu {}] mb{}: enter state deliver to L1 on {} (iter {})", mb.threadid, idx,
                     (void*)(Waddr)addr, machine.iterations);
    mb.state = STATE_DELIVER_TO_L1;
    mb.cycles = L2_LATENCY;

    if (mb.threadid <= 31) {
      if unlikely (icache)
        per_context_dcache_stats_update(mb.threadid, fetch.hit.L2++);
      else
        per_context_dcache_stats_update(mb.threadid, load.hit.L2++);
    }
    return idx;
  }
#ifdef ENABLE_L3_CACHE
  bool L3hit = hierarchy->L3.probe(addr);
  if likely (L3hit) {
    logging::println(logging::DEBUG, "[vcpu {}] mb{}: enter state deliver to L2 on {} (iter {})", mb.threadid, idx,
                     (void*)(Waddr)addr, machine.iterations);
    mb.state = STATE_DELIVER_TO_L2;
    mb.cycles = L3_LATENCY;
    if (mb.threadid <= 31) {
      if (icache)
        per_context_dcache_stats_update(mb.threadid, fetch.hit.L3++);
      else
        per_context_dcache_stats_update(mb.threadid, load.hit.L3++);
    };
    return idx;
  }

  logging::println(logging::DEBUG, "[vcpu {}] mb{}: enter state deliver to L3 on {} (iter {})", mb.threadid, idx,
                   (void*)(Waddr)addr, machine.iterations);
  mb.state = STATE_DELIVER_TO_L3;
  mb.cycles = MAIN_MEM_LATENCY;
#else
  // L3 cache disabled
  logging::println(logging::DEBUG, "[vcpu {}] mb{}: enter state deliver to L2 on {} (iter {})", mb.threadid, idx,
                   (void*)(Waddr)addr, machine.iterations);
  mb.state = STATE_DELIVER_TO_L2;
  mb.cycles = MAIN_MEM_LATENCY;
#endif
  if (mb.threadid <= 31) {
    if unlikely (icache)
      per_context_dcache_stats_update(mb.threadid, fetch.hit.mem++);
    else
      per_context_dcache_stats_update(mb.threadid, load.hit.mem++);
  };

  return idx;
}

template<int SIZE>
int MissBuffer<SIZE>::initiate_miss(LoadFillReq& req, bool hit_in_L2, int rob) {
  int lfrqslot = hierarchy->lfrq.add(req);

  logging::println(logging::DEBUG, "[vcpu {}] missbuf.initiate_miss(L2hit? {}) -> lfrqslot {}", req.lsi.threadid,
                   hit_in_L2, lfrqslot);

  if unlikely (lfrqslot < 0)
    return -1;

  int mbidx = initiate_miss(req.addr, hit_in_L2, 0, rob, req.lsi.threadid);
  if unlikely (mbidx < 0) {
    hierarchy->lfrq.free(lfrqslot);
    return -1;
  }

  Entry& missbuf = missbufs[mbidx];
  missbuf.lfrqmap[lfrqslot] = 1;
  hierarchy->lfrq[lfrqslot].mbidx = mbidx;
  // missbuf.threadid = req.lsi.threadid;       //NOTE SD: Multiple threads can share the same MBE

  return lfrqslot;
}

template<int SIZE>
void MissBuffer<SIZE>::clock() {
  if likely (freemap.all())
    return;

  foreach (i, SIZE) {
    Entry& mb = missbufs[i];
    switch (mb.state) {
    case STATE_IDLE:
      break;
#ifdef ENABLE_L3_CACHE
    case STATE_DELIVER_TO_L3: {
      logging::println(logging::DEBUG, "[vcpu {}] mb{}: deliver {} to L3 ({} cycles left) (iter {})", mb.threadid, i,
                       (void*)(Waddr)mb.addr, mb.cycles, machine.iterations);
      mb.cycles--;
      if unlikely (!mb.cycles) {
        hierarchy->L3.validate(mb.addr);
        mb.cycles = L3_LATENCY;
        mb.state = STATE_DELIVER_TO_L2;
        stats.dcache.missbuf.deliver.mem_to_L3++;
      }
      break;
    }
#endif
    case STATE_DELIVER_TO_L2: {
      logging::println(logging::DEBUG, "[vcpu {}] mb{}: deliver {} to L2 ({} cycles left) (iter {})", mb.threadid, i,
                       (void*)(Waddr)mb.addr, mb.cycles, machine.iterations);
      mb.cycles--;
      if unlikely (!mb.cycles) {
        logging::println(logging::DEBUG, "[vcpu {}] mb{}: delivered to L2", mb.threadid, i);
        hierarchy->L2.validate(mb.addr);
        mb.cycles = L2_LATENCY;
        mb.state = STATE_DELIVER_TO_L1;
        stats.dcache.missbuf.deliver.L3_to_L2++;
      }
      break;
    }
    case STATE_DELIVER_TO_L1: {
      logging::println(logging::DEBUG, "[vcpu {}] mb{}: deliver {} to L1 ({} cycles left) (iter {})", mb.threadid, i,
                       (void*)(Waddr)mb.addr, mb.cycles, machine.iterations);
      mb.cycles--;
      if unlikely (!mb.cycles) {
        logging::println(logging::DEBUG, "[vcpu {}] mb{}: delivered to L1 switch", mb.threadid, i);

        if likely (mb.dcache) {
          logging::println(logging::DEBUG, "[vcpu {}] mb{}: delivered {} to L1 dcache", mb.threadid, i,
                           (void*)(Waddr)mb.addr);
          // If the L2 line size is bigger than the L1 line size, this will validate multiple lines in the L1 when an L2 line arrives:
          // foreach (i, L2_LINE_SIZE / L1_LINE_SIZE) L1.validate(mb.addr + i*L1_LINE_SIZE, std::bitset<L1_LINE_SIZE>().set());
          hierarchy->L1.validate(mb.addr, std::bitset<L1_LINE_SIZE>().set());
          stats.dcache.missbuf.deliver.L2_to_L1D++;
          hierarchy->lfrq.wakeup(mb.addr, mb.lfrqmap);
        }
        if unlikely (mb.icache) {
          // Sometimes we can initiate an icache miss on an existing dcache line in the missbuf
          logging::println(logging::DEBUG, "[vcpu {}] mb{}: delivered {} to L1 icache", mb.threadid, i,
                           (void*)(Waddr)mb.addr);
          // If the L2 line size is bigger than the L1 line size, this will validate multiple lines in the L1 when an L2 line arrives:
          // foreach (i, L2_LINE_SIZE / L1I_LINE_SIZE) L1I.validate(mb.addr + i*L1I_LINE_SIZE, std::bitset<L1I_LINE_SIZE>().set());
          hierarchy->L1I.validate(mb.addr, std::bitset<L1I_LINE_SIZE>().set());
          stats.dcache.missbuf.deliver.L2_to_L1I++;
          LoadStoreInfo lsi = 0;
          lsi.rob = mb.rob;
          lsi.threadid = mb.threadid; /* FIXME: can mb.threadid be 0xfe at this point, and does that matter ? */
          if likely (hierarchy->callback)
            hierarchy->callback->icache_wakeup(lsi, mb.addr);
        }

        assert(!freemap[i]);
        freemap[i] = 1;
        mb.reset();
        count--;
        assert(count >= 0);
      }
      break;
    }
    }
  }
}

template<int SIZE>
void MissBuffer<SIZE>::annul_lfrq(int slot) {
  foreach (i, SIZE) {
    Entry& mb = missbufs[i];
    mb.lfrqmap[slot] = 0; // which LFRQ entries should this load wake up?
  }
}

} // namespace x86sim

template<int linesize>
auto std::formatter<std::pair<x86sim::CacheSubsystem::CacheLine<linesize>, x86sim::W64>>::format(
    const std::pair<x86sim::CacheSubsystem::CacheLine<linesize>, x86sim::W64>& p, std::format_context& ctx) const {
  using namespace x86sim;
  auto out = ctx.out();
  const x86sim::byte* data = (const x86sim::byte*)&p.first;
  foreach (i, linesize / 8) {
    out = std::format_to(out, "    {} \n", bytemaskstring(data + i * 8, (W64)-1LL, 8, 8));
  }
  return out;
}


template<int linesize>
auto std::formatter<std::pair<x86sim::CacheSubsystem::CacheLineWithValidMask<linesize>, x86sim::W64>>::format(
    const std::pair<x86sim::CacheSubsystem::CacheLineWithValidMask<linesize>, x86sim::W64>& p,
    std::format_context& ctx) const {
  using namespace x86sim;
  auto out = ctx.out();
  const x86sim::byte* data = (const x86sim::byte*)&p.first;
  foreach (i, linesize / 8) {
    out = std::format_to(out, "    {} \n", bytemaskstring(data + i * 8, p.first.valid(i * 8, 8).integer(), 8, 8));
  }
  return out;
}

namespace x86sim {

int CacheHierarchy::issueload_slowpath(Waddr physaddr, SFR& sfra, LoadStoreInfo lsi, bool& L2hit) {
  starttimer(load_slowpath_timer);

  L1CacheLine* L1line = L1.probe(physaddr);

  //
  // Loads and stores that also miss the L2 Stores that
  // miss both the L1 and L2 do not require this since
  // there could not possibly be a previous load or
  // store within the current trace that accessed that
  // line (otherwise it would already have been allocated
  // and locked in the L2). In this case, allocate a
  // fresh L2 line and wait for the data to arrive.
  //

  if (L1line)
    logging::println(logging::TRACE, "issue_load_slowpath: L1line for {} = {} validmask ({})", (void*)(Waddr)physaddr,
                     (void*)L1line, L1line->valid.to_string());
  else
    logging::println(logging::TRACE, "issue_load_slowpath: L1line for {} = {} validmask (not found)",
                     (void*)(Waddr)physaddr, (void*)L1line);

  if likely (!L1line) {
    //L1line = L1.select(physaddr);
    stats.dcache.load.transfer.L2_to_L1_full++;
  } else {
    stats.dcache.load.transfer.L2_to_L1_partial++;
  }

  L2hit = 0;

  L2CacheLine* L2line = L2.probe(physaddr);

  if likely (L2line) {
    //
    // We had at least a partial L2 hit, but is the requested data actually mapped into the line?
    //
    std::bitset<L2_LINE_SIZE> sframask, reqmask;
    prep_L2_sframask_and_reqmask((lsi.sfrused) ? &sfra : null, physaddr, lsi.sizeshift, sframask, reqmask);
    L2hit =
        (lsi.sfrused) ? ((reqmask & (sframask | L2line->valid)) == reqmask) : ((reqmask & L2line->valid) == reqmask);
#ifdef ISSUE_LOAD_STORE_DEBUG
    logging::println(logging::DEBUG, "L2hit = {}", L2hit);
    logging::println(logging::DEBUG, "  cachemask {}", L2line->valid);
    logging::println(logging::DEBUG, "  sframask  {}", sframask);
    logging::println(logging::DEBUG, "  reqmask   {}", reqmask);
#endif
  }

#ifdef CACHE_ALWAYS_HITS
  L1line = L1.select(physaddr);
  L1line->tag = L1.tagof(physaddr);
  L1line->valid.set();
  L2line->tag = L2.tagof(physaddr);
  L2line->valid.set();
  L2hit = 1;
#endif

#ifdef L2_ALWAYS_HITS
  L2line = L2.select(physaddr);
  L2line->tag = L2.tagof(physaddr);
  L2line->valid.set();
  L2line->lru = machine.sim_cycle;
  L2hit = 1;
#endif

  //
  // Regardless of whether or not we had a hit somewhere,
  // L1 and L2 lines have been allocated by this point.
  // Slap a lock on the L2 line it so it can't get evicted.
  // Once it's locked up, we can move it into the L1 later.
  //
  // If we did have a hit, but either the L1 or L2 lines
  // were still missing bytes, initiate prefetches to fill
  // them in.
  //

  LoadFillReq req(physaddr, lsi.sfrused ? sfra.data : 0, lsi.sfrused ? sfra.bytemask : 0, lsi, machine.sim_cycle);

  int lfrqslot = missbuf.initiate_miss(req, L2hit, lsi.rob);

  if unlikely (lfrqslot < 0) {
    logging::println(logging::DEBUG,
                     "iteration {}: LFRQ or MB has no free entries for L2->L1: forcing LFRQFull exception",
                     machine.iterations);
    stoptimer(load_slowpath_timer);
    return -1;
  }

  stoptimer(load_slowpath_timer);

  return lfrqslot;
}

int CacheHierarchy::get_lfrq_mb(int lfrqslot) const {
  assert(inrange(lfrqslot, 0, LFRQ_SIZE - 1));

  const LoadFillReq& req = lfrq.reqs[lfrqslot];
  return req.mbidx;
}

int CacheHierarchy::get_lfrq_mb_state(int lfrqslot) const {
  assert(inrange(lfrqslot, 0, LFRQ_SIZE - 1));

  const LoadFillReq& req = lfrq.reqs[lfrqslot];
  if unlikely (req.mbidx < 0)
    return -1;

  assert(!missbuf.freemap[req.mbidx]);
  return missbuf.missbufs[req.mbidx].state;
}

bool CacheHierarchy::covered_by_sfr(W64 addr, SFR* sfr, int sizeshift) {
  std::bitset<L1_LINE_SIZE> sframask, reqmask;
  prep_sframask_and_reqmask(sfr, addr, sizeshift, sframask, reqmask);
  return ((sframask & reqmask) == reqmask);
}

bool CacheHierarchy::probe_cache_and_sfr(W64 addr, const SFR* sfr, int sizeshift) {
  std::bitset<L1_LINE_SIZE> sframask, reqmask;
  prep_sframask_and_reqmask(sfr, addr, sizeshift, sframask, reqmask);

  //
  // Short circuit if the SFR covers the entire load: no need for cache probe
  //
  if unlikely ((sframask & reqmask) == reqmask)
    return true;

  L1CacheLine* L1line = L1.probe(addr);

  if unlikely (!L1line)
    return false;

  //
  // We have a hit on the L1 line itself, but still need to make
  // sure all the data can be filled by some combination of
  // bytes from sfra or the cache data.
  //
  // If not, put this request on the LFRQ and mark it as waiting.
  //

  bool hit = ((reqmask & (sframask | L1line->valid)) == reqmask);

  return hit;
}

void CacheHierarchy::annul_lfrq_slot(int lfrqslot) {
  lfrq.annul(lfrqslot);
}

//
// NOTE: lsi should specify destination of REG_null for prefetches!
//
static const int PREFETCH_STOPS_AT_L2 = 0;

void CacheHierarchy::initiate_prefetch(W64 addr, int cachelevel, int threadid) {
  addr = floor(addr, L1_LINE_SIZE);

  L1CacheLine* L1line = L1.probe(addr);

  if unlikely (L1line) {
    stats.dcache.prefetch.in_L1++;
    return;
  }

  L2CacheLine* L2line = L2.probe(addr);

  if unlikely (L2line) {
    stats.dcache.prefetch.in_L2++;
    if (PREFETCH_STOPS_AT_L2)
      return; // only move up to L2 level, and it's already there
  }

  logging::println(logging::DEBUG, "Prefetch requested for {} to cache level {}", (void*)(Waddr)addr, cachelevel);

  missbuf.initiate_miss(
      addr, L2line,
      false /*, 0xffff, threadid*/); // NB: no threadid -> default bogus threadid -> not flushed on pipeline flush
  stats.dcache.prefetch.required++;
}

//
// Instruction cache
//

bool CacheHierarchy::probe_icache(Waddr virtaddr, Waddr physaddr) {
  L1ICacheLine* L1line = L1I.probe(physaddr);
  bool hit = (L1line != null);

  return hit;
}

int CacheHierarchy::initiate_icache_miss(W64 addr, int rob, int threadid) {
  addr = floor(addr, L1I_LINE_SIZE);
  bool line_in_L2 = (L2.probe(addr) != null);
  int mb = missbuf.initiate_miss(addr, L2.probe(addr), true, rob, threadid);

  logging::println(logging::DEBUG, "[vcpu {}] Initiate icache miss on {} to missbuf {} ({})", threadid,
                   (void*)(Waddr)addr, mb, (line_in_L2 ? "in L2" : "not in L2"));

  return mb;
}

//
// Commit one store from an SFR to the L2 cache without locking
// any cache lines. The store must have already been checked
// to have no exceptions.
//
W64 CacheHierarchy::commitstore(const SFR& sfr, int threadid, bool perform_actual_write) {
  if unlikely (sfr.invalid | (sfr.bytemask == 0))
    return 0;

  static const bool DEBUG = 0;

  starttimer(store_flush_timer);

  W64 addr = sfr.physaddr << 3;

  L2CacheLine* L2line = L2.select(addr);

  if likely (perform_actual_write)
    storemask(addr, sfr.data, sfr.bytemask);

  L1CacheLine* L1line = L1.select(addr);

  L1line->valid |= ((W64)sfr.bytemask << lowbits(addr, 6));
  L2line->valid |= ((W64)sfr.bytemask << lowbits(addr, 6));

  if unlikely (!L1line->valid.all()) {
    per_context_dcache_stats_update(threadid, store.prefetches++);
    missbuf.initiate_miss(addr, L2line->valid.all(), false, 0xffff, threadid);
  }

  stoptimer(store_flush_timer);

  return 0;
}

//
// Submit a speculative store that marks the relevant bytes as valid
// so they can be immediately forwarded to loads, but do not actually
// write to the cache itself.
//
W64 CacheHierarchy::speculative_store(const SFR& sfr, int threadid) {
  return commitstore(sfr, threadid, false);
}

void CacheHierarchy::clock() {
  if unlikely ((machine.sim_cycle & 0x7fffffff) == 0x7fffffff) {
    // Clear any 32-bit cycle-related counters in the cache to prevent wraparound:
    L1.clearstats();
    L1I.clearstats();
    L2.clearstats();
#ifdef ENABLE_L3_CACHE
    L3.clearstats();
#endif
    logging::println(logging::INFO, "Clearing cache statistics to prevent wraparound...");
    logging::flush();
  }

  lfrq.clock();
  missbuf.clock();
}

void CacheHierarchy::complete() {
  lfrq.restart();
  missbuf.restart();
}

void CacheHierarchy::complete(int threadid) {
  lfrq.reset(threadid);
  missbuf.reset(threadid);
}

void CacheHierarchy::reset() {
  lfrq.reset();
  missbuf.reset();
#ifdef ENABLE_L3_CACHE
  L3.reset();
#endif
  L2.reset();
  L1.reset();
  L1I.reset();
  itlb.reset();
  dtlb.reset();
}

} // namespace x86sim

auto std::formatter<x86sim::CacheSubsystem::CacheHierarchy>::format(const x86sim::CacheSubsystem::CacheHierarchy& ch,
                                                                    std::format_context& ctx) const {
  using namespace x86sim;
  auto out = ctx.out();
  out = std::format_to(out, "Data Cache Subsystem:\n");
  out = std::format_to(out, "{}", ch.lfrq);
  out = std::format_to(out, "{}", ch.missbuf);
  out = std::format_to(out, "L1: {} sets x {} ways, line {} bytes\n", CacheSubsystem::L1_SET_COUNT,
                       CacheSubsystem::L1_WAY_COUNT, CacheSubsystem::L1_LINE_SIZE);
  out = std::format_to(out, "L2: {} sets x {} ways, line {} bytes\n", CacheSubsystem::L2_SET_COUNT,
                       CacheSubsystem::L2_WAY_COUNT, CacheSubsystem::L2_LINE_SIZE);
  return out;
}

namespace x86sim {

//
// Make sure the templates and vtables get instantiated:
//
void PerCoreCacheCallbacks::dcache_wakeup(LoadStoreInfo lsi, W64 physaddr) {}
void PerCoreCacheCallbacks::icache_wakeup(LoadStoreInfo lsi, W64 physaddr) {}

template struct CacheSubsystem::LoadFillReqQueue<LFRQ_SIZE>;
template struct CacheSubsystem::MissBuffer<MISSBUF_COUNT>;

/*
// Generator for expand_8bit_to_64bit_lut:

foreach (i, 256) {
byte* m = (byte*)(&expand_8bit_to_64bit_lut[i]);
m[0] = (bit(i, 0) ? 0xff : 0x00);
m[1] = (bit(i, 1) ? 0xff : 0x00);
m[2] = (bit(i, 2) ? 0xff : 0x00);
m[3] = (bit(i, 3) ? 0xff : 0x00);
m[4] = (bit(i, 4) ? 0xff : 0x00);
m[5] = (bit(i, 5) ? 0xff : 0x00);
m[6] = (bit(i, 6) ? 0xff : 0x00);
m[7] = (bit(i, 7) ? 0xff : 0x00);
logging::print(logging::INFO, "  0x{:016x}, ", expand_8bit_to_64bit_lut[i]);
if ((i & 3) == 3) logging::println(logging::INFO, "")
}
*/

} // namespace x86sim
