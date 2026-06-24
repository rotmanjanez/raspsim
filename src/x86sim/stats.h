// -*- c++ -*-
//
// PTLsim: Cycle Accurate x86-64 Simulator
// Statistics data store tree
//
// Copyright 2005-2008 Matt T. Yourst <yourst@yourst.com>
//

#ifndef _STATS_H_
#define _STATS_H_

#include <ctime>
#include <string_view>

#include "globals.h"
#include "superstl.h"
#include "ptlsim.h"

#define STATS_ONLY
#include "decode.h"
#include "ooocore.h"
#include "dcache.h"
#include "branchpred.h"
namespace x86sim {

#undef STATS_ONLY


#define increment_clipped_histogram(h, slot, incr) h[std::clamp(W64(slot), W64(0), W64(lengthof(h) - 1))] += incr;

struct PTLsimStats {
  struct summary {
    W64 cycles;
    W64 insns;
    W64 uops;
    W64 basicblocks;
  } summary;

  // Run-time information, captured when the global stats object is
  // constructed; compile-time build provenance lives in build_info
  struct SimulatorInfo {
    W64 run_timestamp = static_cast<W64>(std::time(nullptr));
    std::string_view platform = host_platform_name();
  } simulator;

  //
  // Decoder and basic block cache
  //
  struct decoder {
    struct throughput {
      W64 basic_blocks;
      W64 x86_insns;
      W64 uops;
      W64 bytes;
    } throughput;

    W64 x86_decode_type[DECODE_TYPE_COUNT];

    struct bb_decode_type {
      W64 all_insns_fast;
      W64 some_complex_insns;
    } bb_decode_type;

    // Alignment of instructions within pages
    struct page_crossings {
      W64 within_page;
      W64 crosses_page;
    } page_crossings;

    // Basic block cache
    struct bbcache {
      W64 count;
      W64 inserts;
      W64 invalidates[INVALIDATE_REASON_COUNT];
    } bbcache;

    // Page cache
    struct pagecache {
      W64 count;
      W64 inserts;
      W64 invalidates[INVALIDATE_REASON_COUNT];
    } pagecache;

    W64 reclaim_rounds;
  } decoder;

  OutOfOrderCoreStats ooocore;
  DataCacheStats dcache;


  struct external {
    W64 assists[ASSIST_COUNT];
    W64 traps[256];
  } external;
};

extern struct PTLsimStats stats;


} // namespace x86sim

#endif // _STATS_H_
