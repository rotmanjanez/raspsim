// -*- c++ -*-
//
// PTLsim: Cycle Accurate x86-64 Simulator
// Data Cache
//
// Copyright 2003-2008 Matt T. Yourst <yourst@yourst.com>
//

#ifndef _DCACHE_H_
#define _DCACHE_H_

#include "ptlsim.h"
#include "x86sim/logging.h"

namespace x86sim {

struct LoadStoreInfo {
  W16 rob;
  W8 threadid;
  W8 sizeshift : 2, aligntype : 2, sfrused : 1, internal : 1, signext : 1, pad1 : 1;
  W32 pad32;
  RawDataAccessors(LoadStoreInfo, W64);
};

#define per_context_dcache_stats_ref(vcpuid) (*(((PerContextDataCacheStats*)&stats.dcache.vcpu0) + (vcpuid)))
#define per_context_dcache_stats_update(vcpuid, expr) stats.dcache.total.expr, per_context_dcache_stats_ref(vcpuid).expr

namespace CacheSubsystem {
// How many load wakeups can be driven into the core each cycle:
const int MAX_WAKEUPS_PER_CYCLE = 2;

#ifndef STATS_ONLY

// non-debugging only:
//#define __RELEASE__
#ifdef __RELEASE__
#undef assert
#define assert(x) (x)
#endif

//#define CACHE_ALWAYS_HITS
//#define L2_ALWAYS_HITS

// 16 KB L1 at 2 cycles       // increase to 32 KB to match Core 2
const int L1_LINE_SIZE = 64;
const int L1_SET_COUNT = 64;
const int L1_WAY_COUNT = 4;
// #define ENFORCE_L1_DCACHE_BANK_CONFLICTS
const int L1_DCACHE_BANKS = 8; // 8 banks x 8 bytes/bank = 64 bytes/line

// 32 KB L1I
const int L1I_LINE_SIZE = 64;
const int L1I_SET_COUNT = 128;
const int L1I_WAY_COUNT = 4;

// 256 KB L2 at 6 cycles
const int L2_LINE_SIZE = 64;
const int L2_SET_COUNT = 256; // 256 KB
const int L2_WAY_COUNT = 16;
const int L2_LATENCY = 5; // don't include the extra wakeup cycle (waiting->ready state transition) in the LFRQ

#define ENABLE_L3_CACHE
#ifdef ENABLE_L3_CACHE
// 4 MB L3 cache (2048 sets, 32 ways) with 64-byte lines, latency 16 cycles
const int L3_SET_COUNT = 2048;
const int L3_WAY_COUNT = 32;
const int L3_LINE_SIZE = 64;
const int L3_LATENCY = 8; // Core 2 Duo 2.0 GHz has 14 cycle total L2 latency
#endif
// Load Fill Request Queue (maximum number of missed loads)
// const int LFRQ_SIZE = 63;
const int LFRQ_SIZE = 64;

// Allow up to 32 outstanding lines in the L2 awaiting service:
const int MISSBUF_COUNT = 64;
// const int MISSBUF_COUNT = 4;

// Main memory latency
const int MAIN_MEM_LATENCY = 140; // Core 2 Duo 2.4 GHz has 160 cycle total L2 latency

// TLBs
const int ITLB_SIZE = 32;
const int DTLB_SIZE = 32;

//#define ISSUE_LOAD_STORE_DEBUG
//#define CHECK_LOADS_AND_STORES

// Line Usage Statistics

//#define TRACK_LINE_USAGE

#ifdef TRACK_LINE_USAGE
#define DCACHE_L1_LINE_LIFETIME_INTERVAL 1
#define DCACHE_L1_LINE_DEADTIME_INTERVAL 1
#define DCACHE_L1_LINE_HITCOUNT_INTERVAL 1
#define DCACHE_L1_LINE_LIFETIME_SLOTS 8192
#define DCACHE_L1_LINE_DEADTIME_SLOTS 8192
#define DCACHE_L1_LINE_HITCOUNT_SLOTS 64

#define DCACHE_L1I_LINE_LIFETIME_INTERVAL 16
#define DCACHE_L1I_LINE_DEADTIME_INTERVAL 16
#define DCACHE_L1I_LINE_HITCOUNT_INTERVAL 1
#define DCACHE_L1I_LINE_LIFETIME_SLOTS 8192
#define DCACHE_L1I_LINE_DEADTIME_SLOTS 8192
#define DCACHE_L1I_LINE_HITCOUNT_SLOTS 1024

#define DCACHE_L2_LINE_LIFETIME_INTERVAL 4
#define DCACHE_L2_LINE_DEADTIME_INTERVAL 4
#define DCACHE_L2_LINE_HITCOUNT_INTERVAL 1
#define DCACHE_L2_LINE_LIFETIME_SLOTS 65536
#define DCACHE_L2_LINE_DEADTIME_SLOTS 65536
#define DCACHE_L2_LINE_HITCOUNT_SLOTS 256

#define DCACHE_L3_LINE_LIFETIME_INTERVAL 64
#define DCACHE_L3_LINE_DEADTIME_INTERVAL 64
#define DCACHE_L3_LINE_HITCOUNT_INTERVAL 1
#define DCACHE_L3_LINE_LIFETIME_SLOTS 16384
#define DCACHE_L3_LINE_DEADTIME_SLOTS 16384
#define DCACHE_L3_LINE_HITCOUNT_SLOTS 256
#endif

//
// Cache Line Types
//
template<int linesize>
struct CacheLine {
#ifdef TRACK_LINE_USAGE
  W32 filltime;
  W32 lasttime;
  W32 hitcount;
#else
  byte dummy;
#endif
  void reset() { clearstats(); }
  void invalidate() { reset(); }
  void fill(W64 tag, const std::bitset<linesize>& valid) {}

  void clearstats() {
#ifdef TRACK_LINE_USAGE
    filltime = sim_cycle;
    lasttime = sim_cycle;
    hitcount = 0;
#endif
  }
};

template<int linesize>
struct CacheLineWithValidMask {
  std::bitset<linesize> valid;
#ifdef TRACK_LINE_USAGE
  W32 filltime;
  W32 lasttime;
  W32 hitcount;
#endif

  void clearstats() {
#ifdef TRACK_LINE_USAGE
    filltime = sim_cycle;
    lasttime = sim_cycle;
    hitcount = 0;
#endif
  }

  void reset() {
    valid = 0;
    clearstats();
  }
  void invalidate() { reset(); }
  void fill(W64 tag, const std::bitset<linesize>& valid) { this->valid |= valid; }
};

typedef CacheLineWithValidMask<L1_LINE_SIZE> L1CacheLine;
typedef CacheLine<L1I_LINE_SIZE> L1ICacheLine;
typedef CacheLineWithValidMask<L2_LINE_SIZE> L2CacheLine;
#ifdef ENABLE_L3_CACHE
typedef CacheLine<L3_LINE_SIZE> L3CacheLine;
#endif

//
// L1 data cache
//
#ifdef TRACK_LINE_USAGE
static const char* cache_names[4] = {"L1", "I1", "L2", "L3"};

template<int uniq, typename V, int LIFETIME_INTERVAL, int LIFETIME_SLOTS, int DEADTIME_INTERVAL, int DEADTIME_SLOTS,
         int HITCOUNT_INTERVAL, int HITCOUNT_SLOTS>
struct HistogramAssociativeArrayStatisticsCollector {
  static W64 line_lifetime_histogram[LIFETIME_SLOTS];
  static W64 line_deadtime_histogram[DEADTIME_SLOTS];
  static W64 line_hitcount_histogram[HITCOUNT_SLOTS];

  HistogramAssociativeArrayStatisticsCollector() { reset(); }

  static void reset() {
    setzero(line_lifetime_histogram);
    setzero(line_deadtime_histogram);
    setzero(line_hitcount_histogram);
  }

  static void evicted(const V& line, W64 tag) {
    // Line has been evicted: update statistics
    W64s lifetime = line.lasttime - line.filltime;
    assert(lifetime >= 0);
    int lifetimeslot = std::clamp(lifetime / LIFETIME_INTERVAL, 0, LIFETIME_SLOTS - 1);
    line_lifetime_histogram[lifetimeslot]++;

    W64s deadtime = sim_cycle - line.lasttime;
    int deadtimeslot = std::clamp(deadtime / DEADTIME_INTERVAL, 0, DEADTIME_SLOTS - 1);
    line_deadtime_histogram[deadtimeslot]++;

    W64 hitcount = line.hitcount;
    int hitcountslot = std::clamp(hitcount / HITCOUNT_INTERVAL, 0, HITCOUNT_SLOTS - 1);
    line_hitcount_histogram[hitcountslot]++;

    logging::println("[{}] {}: evicted({}): lifetime {}, deadtime {}, hitcount {} (line addr {})", cache_names[uniq],
                     sim_cycle, (void*)tag, lifetime, deadtime, hitcount, (void*)&line);
  }

  static void filled(V& line, W64 tag) {
    line.filltime = sim_cycle;
    line.lasttime = sim_cycle;
    line.hitcount = 1;

    logging::println("[{}] {}: filled({}) (line addr {})", cache_names[uniq], sim_cycle, (void*)tag, (void*)&line);
  }

  static void inserted(V& line, W64 newtag, int way) { filled(line, newtag); }

  static void replaced(V& line, W64 oldtag, W64 newtag, int way) {
    evicted(line, oldtag);
    filled(line, newtag);
  }

  static void probed(V& line, W64 tag, int way, bool hit) {
    logging::println("[{}] {}: probe({}): {} way {}: hitcount {}, filltime {}, lasttime {} (line addr {})",
                     cache_names[uniq], sim_cycle, (void*)tag, (hit ? "HIT" : "miss"), way, line.hitcount,
                     line.filltime, line.lasttime, (void*)&line);
    if (hit) {
      line.hitcount++;
      line.lasttime = sim_cycle;
    }
  }

  static void overflow(W64 tag) {}

  static void locked(V& slot, W64 tag, int way) {}
  static void unlocked(V& slot, W64 tag, int way) {}

  static void invalidated(V& line, W64 oldtag, int way) { evicted(line, oldtag); }
};

typedef HistogramAssociativeArrayStatisticsCollector<
    0, L1CacheLine, DCACHE_L1_LINE_LIFETIME_INTERVAL, DCACHE_L1_LINE_LIFETIME_SLOTS, DCACHE_L1_LINE_DEADTIME_INTERVAL,
    DCACHE_L1_LINE_DEADTIME_SLOTS, DCACHE_L1_LINE_HITCOUNT_INTERVAL, DCACHE_L1_LINE_HITCOUNT_SLOTS>
    L1StatsCollectorBase;

typedef HistogramAssociativeArrayStatisticsCollector<1, L1ICacheLine, DCACHE_L1I_LINE_LIFETIME_INTERVAL,
                                                     DCACHE_L1I_LINE_LIFETIME_SLOTS, DCACHE_L1I_LINE_DEADTIME_INTERVAL,
                                                     DCACHE_L1I_LINE_DEADTIME_SLOTS, DCACHE_L1I_LINE_HITCOUNT_INTERVAL,
                                                     DCACHE_L1I_LINE_HITCOUNT_SLOTS>
    L1IStatsCollectorBase;

typedef HistogramAssociativeArrayStatisticsCollector<
    2, L2CacheLine, DCACHE_L2_LINE_LIFETIME_INTERVAL, DCACHE_L2_LINE_LIFETIME_SLOTS, DCACHE_L2_LINE_DEADTIME_INTERVAL,
    DCACHE_L2_LINE_DEADTIME_SLOTS, DCACHE_L2_LINE_HITCOUNT_INTERVAL, DCACHE_L2_LINE_HITCOUNT_SLOTS>
    L2StatsCollectorBase;

#ifdef ENABLE_L3_CACHE
typedef HistogramAssociativeArrayStatisticsCollector<
    3, L3CacheLine, DCACHE_L3_LINE_LIFETIME_INTERVAL, DCACHE_L3_LINE_LIFETIME_SLOTS, DCACHE_L3_LINE_DEADTIME_INTERVAL,
    DCACHE_L3_LINE_DEADTIME_SLOTS, DCACHE_L3_LINE_HITCOUNT_INTERVAL, DCACHE_L3_LINE_HITCOUNT_SLOTS>
    L3StatsCollectorBase;
#endif

struct L1StatsCollector : public L1StatsCollectorBase {};
struct L1IStatsCollector : public L1IStatsCollectorBase {};
struct L2StatsCollector : public L2StatsCollectorBase {};
#ifdef ENABLE_L3_CACHE
struct L3StatsCollector : public L3StatsCollectorBase {};
#endif

#else
typedef NullAssociativeArrayStatisticsCollector<W64, L1CacheLine> L1StatsCollector;
typedef NullAssociativeArrayStatisticsCollector<W64, L1ICacheLine> L1IStatsCollector;
typedef NullAssociativeArrayStatisticsCollector<W64, L2CacheLine> L2StatsCollector;
#ifdef ENABLE_L3_CACHE
typedef NullAssociativeArrayStatisticsCollector<W64, L3CacheLine> L3StatsCollector;
#endif
#endif

template<typename V, int setcount, int waycount, int linesize,
         typename stats = NullAssociativeArrayStatisticsCollector<W64, V>>
struct DataCache : public AssociativeArray<W64, V, setcount, waycount, linesize, stats> {
  typedef AssociativeArray<W64, V, setcount, waycount, linesize, stats> base_t;
  void clearstats() {
#ifdef TRACK_LINE_USAGE
    foreach (set, L1_SET_COUNT) {
      foreach (way, waycount) {
        base_t::sets[set][way].clearstats();
      }
    }
#endif
  }
};

struct L1Cache : public DataCache<L1CacheLine, L1_SET_COUNT, L1_WAY_COUNT, L1_LINE_SIZE, L1StatsCollector> {
  L1CacheLine* validate(W64 addr, const std::bitset<L1_LINE_SIZE>& valid) {
    addr = tagof(addr);
    L1CacheLine* line = select(addr);
    line->fill(addr, valid);
    return line;
  }
};

//
// L1 instruction cache
//

struct L1ICache : public DataCache<L1ICacheLine, L1I_SET_COUNT, L1I_WAY_COUNT, L1I_LINE_SIZE, L1IStatsCollector> {
  L1ICacheLine* validate(W64 addr, const std::bitset<L1I_LINE_SIZE>& valid) {
    addr = tagof(addr);
    L1ICacheLine* line = select(addr);
    line->fill(addr, valid);
    return line;
  }
};


//
// L2 cache
//

typedef DataCache<L2CacheLine, L2_SET_COUNT, L2_WAY_COUNT, L2_LINE_SIZE, L2StatsCollector> L2CacheBase;

struct L2Cache : public L2CacheBase {
  void validate(W64 addr) {
    L2CacheLine* line = select(addr);
    if (!line)
      return;
    line->valid.set();
  }

  void deliver(W64 address);
};

//
// L3 cache
//
#ifdef ENABLE_L3_CACHE
struct L3Cache : public DataCache<L3CacheLine, L3_SET_COUNT, L3_WAY_COUNT, L3_LINE_SIZE, L3StatsCollector> {
  L3CacheLine* validate(W64 addr) {
    W64 oldaddr;
    L3CacheLine* line = select(addr, oldaddr);
    return line;
  }
};
#endif

static inline void prep_sframask_and_reqmask(const SFR* sfr, W64 addr, int sizeshift,
                                             std::bitset<L1_LINE_SIZE>& sframask, std::bitset<L1_LINE_SIZE>& reqmask) {
  sframask =
      (sfr) ? (std::bitset<L1_LINE_SIZE>(sfr->bytemask) << 8 * lowbits(sfr->physaddr, log2(L1_LINE_SIZE) - 3)) : 0;
  reqmask = std::bitset<L1_LINE_SIZE>(bitmask(1 << sizeshift)) << lowbits(addr, log2(L1_LINE_SIZE));
}

static inline void prep_L2_sframask_and_reqmask(const SFR* sfr, W64 addr, int sizeshift,
                                                std::bitset<L2_LINE_SIZE>& sframask,
                                                std::bitset<L2_LINE_SIZE>& reqmask) {
  sframask =
      (sfr) ? (std::bitset<L2_LINE_SIZE>(sfr->bytemask) << 8 * lowbits(sfr->physaddr, log2(L2_LINE_SIZE) - 3)) : 0;
  reqmask = std::bitset<L2_LINE_SIZE>(bitmask(1 << sizeshift)) << lowbits(addr, log2(L2_LINE_SIZE));
}

//
// TLB class with one-hot semantics. 36 bit tags are required since
// virtual addresses are 48 bits, so 48 - 12 (2^12 bytes per page)
// is 36 bits.
//
template<int tlbid, int size>
struct TranslationLookasideBuffer : public x86sim::FullyAssociativeTagsNbitOneHot<size, 40> {
  typedef x86sim::FullyAssociativeTagsNbitOneHot<size, 40> base_t;
  TranslationLookasideBuffer() : base_t() {}

  void reset() { base_t::reset(); }

  // Get the 40-bit TLB tag (36 bit virtual page ID plus 4 bit threadid)
  static W64 tagof(W64 addr, W64 threadid) { return bits(addr, 12, 36) | (threadid << 36); }

  bool probe(W64 addr, int threadid = 0) {
    W64 tag = tagof(addr, threadid);
    return (base_t::probe(tag) >= 0);
  }

  bool insert(W64 addr, int threadid = 0) {
    addr = floor(addr, PAGE_SIZE);
    W64 tag = tagof(addr, threadid);
    W64 oldtag;
    int way = base_t::select(tag, oldtag);
    W64 oldaddr = lowbits(oldtag, 36) << 12;
    logging::println(logging::TRACE, "TLB insertion of virt page {} (virt addr {}) into way {}: {}", (void*)(Waddr)addr,
                     (void*)(Waddr)(addr), way, ((oldtag != tag) ? "evicted old entry" : "already present"));
    return (oldtag != tag);
  }

  int flush_all() {
    reset();
    return size;
  }

  int flush_thread(W64 threadid) {
    W64 tag = threadid << 36;
    W64 tagmask = 0xfULL << 36;
    std::bitset<size> slotmask = base_t::masked_match(tag, tagmask);
    int n = slotmask.count();
    base_t::masked_invalidate(slotmask);
    return n;
  }

  int flush_virt(Waddr virtaddr, W64 threadid) { return this->invalidate(tagof(virtaddr, threadid)); }
};

typedef TranslationLookasideBuffer<0, DTLB_SIZE> DTLB;
typedef TranslationLookasideBuffer<1, ITLB_SIZE> ITLB;

struct CacheHierarchy;

//
// Load fill request queue (LFRQ) contains any requests for outstanding
// loads from both the L2 or L1.
//
struct LoadFillReq {
  W64 addr; // physical address
  W64 data; // data already known so far (e.g. from SFR)
  LoadStoreInfo lsi;
  W32 initcycle;
  byte mask;
  byte fillL1 : 1, fillL2 : 1;
  W8s mbidx;

  inline LoadFillReq() {}

  LoadFillReq(W64 addr, W64 data, byte mask, LoadStoreInfo lsi, W64 initcycle);
};


template<int size>
struct LoadFillReqQueue {
  MachineImpl& machine;
  CacheHierarchy* hierarchy = nullptr;
  std::bitset<size> freemap; // Slot is free
  std::bitset<size> waiting; // Waiting for the line to arrive in the L1
  std::bitset<size> ready;   // Wait to extract/signext and write into register
  LoadFillReq reqs[size];
  int count;

  static const int SIZE = size;

  LoadFillReqQueue(MachineImpl& machine_, CacheHierarchy& hierarchy_) : machine(machine_), hierarchy(&hierarchy_) { reset(); }

  // Clear entries belonging to one thread
  void reset(int threadid);

  // Reset all threads
  void reset() {
    freemap.set();
    ready = 0;
    waiting = 0;
    count = 0;
  }

  void changestate(int idx, std::bitset<size>& oldstate, std::bitset<size>& newstate) {
    oldstate[idx] = 0;
    newstate[idx] = 1;
  }

  void free(int lfrqslot) {
    assert(waiting[lfrqslot]);
    changestate(lfrqslot, waiting, freemap);
    assert(count > 0);
    count--;
  }

  bool full() const { return freemap.none(); }

  int remaining() const { return (size - count); }

  void annul(int lfrqslot);

  void restart();

  int add(const LoadFillReq& req);

  void wakeup(W64 address, const std::bitset<LFRQ_SIZE>& lfrqmask);

  void clock();

  LoadFillReq& operator[](int idx) { return reqs[idx]; }
  const LoadFillReq& operator[](int idx) const { return reqs[idx]; }
};

enum { STATE_IDLE, STATE_DELIVER_TO_L3, STATE_DELIVER_TO_L2, STATE_DELIVER_TO_L1 };
static const char* missbuf_state_names[] = {"idle", "mem->L3", "L3->L2", "L2->L1"};

template<int SIZE>
struct MissBuffer {
  struct Entry {
    W64 addr; // physical line address we are waiting for
    W16 state;
    W16 dcache : 1, icache : 1; // L1I vs L1D
    W32 cycles;
    W16 rob;
    W8 threadid;

    std::bitset<LFRQ_SIZE> lfrqmap; // which LFRQ entries should this load wake up?
    void reset() {
      lfrqmap = 0;
      addr = 0xffffffffffffffffULL;
      state = STATE_IDLE;
      cycles = 0;
      icache = 0;
      dcache = 0;
      rob = 0xffff;
      threadid = 0xff;
    }
  };

  MachineImpl& machine;
  CacheHierarchy* hierarchy = nullptr;

  MissBuffer(MachineImpl& machine_, CacheHierarchy& hierarchy_) : machine(machine_), hierarchy(&hierarchy_) { reset(); }
  Entry missbufs[SIZE];
  std::bitset<SIZE> freemap;
  int count;

  void reset();
  void reset(int threadid);
  void restart();
  bool full() const { return freemap.none(); }
  int remaining() const { return (SIZE - count); }
  int find(W64 addr);
  int initiate_miss(W64 addr, bool hit_in_L2, bool icache = 0, int rob = 0xffff, int threadid = 0xfe);
  int initiate_miss(LoadFillReq& req, bool hit_in_L2, int rob = 0xffff);
  void annul_lfrq(int slot);
  void annul_lfrq(int slot, int threadid);
  void clock();
};

struct PerCoreCacheCallbacks {
  virtual void dcache_wakeup(LoadStoreInfo lsi, W64 physaddr);
  virtual void icache_wakeup(LoadStoreInfo lsi, W64 physaddr);
};

struct CacheHierarchy {
  MachineImpl& machine;
  LoadFillReqQueue<LFRQ_SIZE> lfrq;
  MissBuffer<MISSBUF_COUNT> missbuf;
  L1Cache L1;
  L1ICache L1I;
  L2Cache L2;
#ifdef ENABLE_L3_CACHE
  L3Cache L3;
#endif
  DTLB dtlb;
  ITLB itlb;

  PerCoreCacheCallbacks* callback;

  explicit CacheHierarchy(MachineImpl& machine_) : machine(machine_), lfrq(machine_, *this), missbuf(machine_, *this) { callback = null; }

  bool probe_cache_and_sfr(W64 addr, const SFR* sfra, int sizeshift);
  bool covered_by_sfr(W64 addr, SFR* sfr, int sizeshift);
  void annul_lfrq_slot(int lfrqslot);
  int issueload_slowpath(Waddr physaddr, SFR& sfra, LoadStoreInfo lsi, bool& L2hit);

  int issueload_slowpath(Waddr physaddr, SFR& sfra, LoadStoreInfo lsi) {
    bool L2hit = 0;
    return issueload_slowpath(physaddr, sfra, lsi, L2hit);
  }

  int get_lfrq_mb(int lfrqslot) const;
  int get_lfrq_mb_state(int lfrqslot) const;
  bool lfrq_or_missbuf_full() const { return lfrq.full() | missbuf.full(); }

  W64 commitstore(const SFR& sfr, int threadid = 0xff, bool perform_actual_write = true);
  W64 speculative_store(const SFR& sfr, int threadid = 0xff);

  void initiate_prefetch(W64 addr, int cachelevel, int threadid = 0xfe);

  bool probe_icache(Waddr virtaddr, Waddr physaddr);
  int initiate_icache_miss(W64 addr, int rob = 0xffff, int threadid = 0xff);

  void reset();
  void clock();
  void complete();
  void complete(int threadid);
};
#endif // STATS_ONLY
}; // namespace CacheSubsystem

struct PerContextDataCacheStats {
  struct load {
    struct hit {
      W64 L1;
      W64 L2;
      W64 L3;
      W64 mem;
    } hit;

    struct dtlb {
      W64 hits;
      W64 misses;
    } dtlb;

    struct tlbwalk {
      W64 L1_dcache_hit;
      W64 L1_dcache_miss;
      W64 no_lfrq_mb;
    } tlbwalk;
  } load;

  struct fetch {
    struct hit {
      W64 L1;
      W64 L2;
      W64 L3;
      W64 mem;
    } hit;

    struct itlb {
      W64 hits;
      W64 misses;
    } itlb;

    struct tlbwalk {
      W64 L1_dcache_hit;
      W64 L1_dcache_miss;
      W64 no_lfrq_mb;
    } tlbwalk;
  } fetch;

  struct store {
    W64 prefetches;
  } store;
};

struct DataCacheStats {
  struct load {
    struct transfer {
      W64 L2_to_L1_full;
      W64 L2_to_L1_partial;
      W64 L2_L1I_full;
    } transfer;
  } load;

  struct missbuf {
    W64 inserts;
    struct deliver {
      W64 mem_to_L3;
      W64 L3_to_L2;
      W64 L2_to_L1D;
      W64 L2_to_L1I;
    } deliver;
  } missbuf;

  struct prefetch {
    W64 in_L1;
    W64 in_L2;
    W64 required;
  } prefetch;

  struct lfrq {
    W64 inserts;
    W64 wakeups;
    W64 annuls;
    W64 resets;
    W64 total_latency;
    double average_latency;
    W64 width[CacheSubsystem::MAX_WAKEUPS_PER_CYCLE + 1];
  } lfrq;

  PerContextDataCacheStats total;
  // IMPORTANT: This list MUST be equal in length to the number of active VCPUs (at most MAX_CONTEXTS):
  PerContextDataCacheStats vcpu0;
  PerContextDataCacheStats vcpu1;
  PerContextDataCacheStats vcpu2;
  PerContextDataCacheStats vcpu3;
  PerContextDataCacheStats vcpu4;
  PerContextDataCacheStats vcpu5;
  PerContextDataCacheStats vcpu6;
  PerContextDataCacheStats vcpu7;
  PerContextDataCacheStats vcpu8;
  PerContextDataCacheStats vcpu9;
  PerContextDataCacheStats vcpu10;
  PerContextDataCacheStats vcpu11;
  PerContextDataCacheStats vcpu12;
  PerContextDataCacheStats vcpu13;
  PerContextDataCacheStats vcpu14;
  PerContextDataCacheStats vcpu15;
  PerContextDataCacheStats vcpu16;
  PerContextDataCacheStats vcpu17;
  PerContextDataCacheStats vcpu18;
  PerContextDataCacheStats vcpu19;
  PerContextDataCacheStats vcpu20;
  PerContextDataCacheStats vcpu21;
  PerContextDataCacheStats vcpu22;
  PerContextDataCacheStats vcpu23;
  PerContextDataCacheStats vcpu24;
  PerContextDataCacheStats vcpu25;
  PerContextDataCacheStats vcpu26;
  PerContextDataCacheStats vcpu27;
  PerContextDataCacheStats vcpu28;
  PerContextDataCacheStats vcpu29;
  PerContextDataCacheStats vcpu30;
  PerContextDataCacheStats vcpu31;
};

#ifndef STATS_ONLY
//
// std::formatter specializations for C++23 std::print/std::format support
//
} // namespace x86sim

namespace std {

template<>
struct formatter<x86sim::CacheSubsystem::LoadFillReq> {
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }
  auto format(const x86sim::CacheSubsystem::LoadFillReq& req, format_context& ctx) const;
};

template<int Size>
struct formatter<x86sim::CacheSubsystem::LoadFillReqQueue<Size>> {
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }
  auto format(const x86sim::CacheSubsystem::LoadFillReqQueue<Size>& lfrq, format_context& ctx) const {
    auto out = ctx.out();
    out = std::format_to(out, "LoadFillReqQueue<{}>: {} of {} entries ({} free)\n", Size, lfrq.count, Size,
                         Size - lfrq.count);
    out = std::format_to(out, "  Free:   {}\n", lfrq.freemap.to_string());
    out = std::format_to(out, "  Wait:   {}\n", lfrq.waiting.to_string());
    out = std::format_to(out, "  Ready:  {}\n", lfrq.ready.to_string());
    foreach (i, Size) {
      if (!lfrq.freemap[i]) {
        out = std::format_to(out, "  slot {:>2}: {}\n", i, lfrq.reqs[i]);
      }
    }
    return out;
  }
};

template<int Size>
struct formatter<x86sim::CacheSubsystem::MissBuffer<Size>> {
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }
  auto format(const x86sim::CacheSubsystem::MissBuffer<Size>& mb, format_context& ctx) const {
    auto out = ctx.out();
    out = std::format_to(out, "MissBuffer<{}>:\n", Size);
    foreach (i, Size) {
      if likely (mb.freemap[i])
        continue;
      const auto& entry = mb.missbufs[i];
      out =
          std::format_to(out, "slot {:>2}: vcpu {}, addr {} state {:<8} {} {} on {} cycles -> lfrq {}\n", i,
                         entry.threadid, (void*)(x86sim::Waddr)entry.addr,
                         x86sim::CacheSubsystem::missbuf_state_names[entry.state],
                         (entry.dcache ? "dcache" : " "), (entry.icache ? "icache" : " "), entry.cycles,
                         entry.lfrqmap.to_string());
    }
    return out;
  }
};

template<>
struct formatter<x86sim::CacheSubsystem::CacheHierarchy> {
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }
  auto format(const x86sim::CacheSubsystem::CacheHierarchy& ch, format_context& ctx) const;
};

template<int Type, int Size>
struct formatter<x86sim::CacheSubsystem::TranslationLookasideBuffer<Type, Size>> {
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }
  auto format(const x86sim::CacheSubsystem::TranslationLookasideBuffer<Type, Size>& tlb, format_context& ctx) const {
    using base_t = x86sim::FullyAssociativeTagsNbitOneHot<Size, 40>;
    return formatter<base_t>{}.format(static_cast<const base_t&>(tlb), ctx);
  }
};

template<int linesize>
struct formatter<std::pair<x86sim::CacheSubsystem::CacheLine<linesize>, x86sim::W64>> {
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }
  auto format(const std::pair<x86sim::CacheSubsystem::CacheLine<linesize>, x86sim::W64>& p, format_context& ctx) const;
};

template<int linesize>
struct formatter<std::pair<x86sim::CacheSubsystem::CacheLineWithValidMask<linesize>, x86sim::W64>> {
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }
  auto format(const std::pair<x86sim::CacheSubsystem::CacheLineWithValidMask<linesize>, x86sim::W64>& p, format_context& ctx) const;
};

} // namespace std
#else
} // namespace x86sim
#endif // STATS_ONLY

#endif // _DCACHE_H_
