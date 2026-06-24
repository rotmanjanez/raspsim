//
// PTLsim: Cycle Accurate x86-64 Simulator
// Branch Prediction
//
// Copyright 2003-2008 Matt T. Yourst <yourst@yourst.com>
//
// This program is free software; it is licensed under the
// GNU General Public License, Version 2.
//

#include "branchpred.h"
#include "stats.h"
#include "x86sim/logging.h"

namespace x86sim {

template<int SIZE>
struct BimodalPredictor {
  std::array<byte, SIZE> table;

  void reset() {
    foreach (i, SIZE)
      table[i] = 1;
  }

  inline int hash(W64 branchaddr) { return lowbits((branchaddr >> 16) ^ branchaddr, log2(SIZE)); }

  byte* predict(W64 branchaddr) { return &table[hash(branchaddr)]; }
};

template<int L1SIZE, int L2SIZE, int SHIFTWIDTH, bool HISTORYXOR>
struct TwoLevelPredictor {
  std::array<int, L1SIZE> shiftregs; // L1 history shift register(s)
  std::array<byte, L2SIZE> L2table;  // L2 prediction state table

  void reset() {
    // initialize counters uniformly
    foreach (i, L2SIZE)
      L2table[i] = 1;
  }

  byte* predict(W64 branchaddr) {
    int L1index = lowbits(branchaddr, log2(L1SIZE));
    int L2index = shiftregs[L1index];

    if (HISTORYXOR) {
      L2index ^= branchaddr;
    } else {
      L2index |= branchaddr << SHIFTWIDTH;
    }

    L2index = lowbits(L2index, log2(L2SIZE));

    return &L2table[L2index];
  }
};

struct BTBEntry {
  W64 target; // last destination of branch when taken

  void reset() { target = 0; }
};

template<int SETCOUNT, int WAYCOUNT>
struct BranchTargetBuffer : public AssociativeArray<W64, BTBEntry, SETCOUNT, WAYCOUNT, 1> {};

template<int SIZE>
struct ReturnAddressStack;


// Enable to debug the return address stack (RAS) predictor mechanism
// #define DEBUG_RAS

template<int SIZE>
struct ReturnAddressStack : public Queue<ReturnAddressStackEntry, SIZE> {
  typedef Queue<ReturnAddressStackEntry, SIZE> base_t;

  void push(W64 uuid, W64 rip, ReturnAddressStackEntry& old) {
    logging::println(logging::DEBUG, "ReturnAddressStack::push(uuid {}, rip {}):", uuid, (void*)(Waddr)rip);
    if (base_t::full()) {
      logging::println(logging::DEBUG, "  Return address stack overflow: removing oldest entry to make space");
      stats.ooocore.branchpred.ras.overflows++;
      base_t::pophead();
    }

    ReturnAddressStackEntry& e = *base_t::push();
    assert(&e);

    old = e;
    logging::println(logging::DEBUG, "  Old entry: {}", old);

    e.uuid = uuid;
    e.rip = rip;

    stats.ooocore.branchpred.ras.pushes++;
    logging::println(logging::DEBUG, "{}", *this);
  }

  ReturnAddressStackEntry& pop(ReturnAddressStackEntry& old) {
    logging::println(logging::DEBUG, "ReturnAddressStack::pop():");
    if (base_t::empty()) {
      stats.ooocore.branchpred.ras.underflows++;
      logging::println(logging::DEBUG, "  Return address stack underflow: returning entry with zero fields");
      old.idx = -1;
      old.uuid = 0;
      old.rip = 0;
      return old;
    }

    ReturnAddressStackEntry& e = *base_t::pop();
    assert(&e);
    old = e;

    logging::println(logging::DEBUG, "  Old entry: {}", old);
    logging::println(logging::DEBUG, "{}", *this);

    stats.ooocore.branchpred.ras.pops++;

    return e;
  }

  W64 peek() {
    logging::println(logging::DEBUG, "ReturnAddressStack::peek():");
    if (base_t::empty()) {
      logging::println(logging::DEBUG, "  Return address stack is empty: returning bogus rip 0");
      return 0;
    }

    logging::print(logging::INFO, "  Peeking entry {}", (*base_t::peektail()));

    return base_t::peektail()->rip;
  }

  //
  // Pop a speculative push from the stack
  //
  void annulpush(const ReturnAddressStackEntry& old) {
    logging::println(logging::DEBUG, "ReturnAddressStack::annulpush(old index {}, uuid {}, rip {}):", old.idx, old.uuid,
                     (void*)(Waddr)old.rip);

    if (base_t::empty()) {
      logging::println(logging::DEBUG, "  Cannot annul: return address stack is empty");
      return;
    }

    ReturnAddressStackEntry& e = *base_t::peektail();
    e.uuid = old.uuid;
    e.rip = old.rip;

    ReturnAddressStackEntry dummy;
    pop(dummy);
    logging::println(logging::DEBUG, "  Popped speculative push; e.index = {} vs tail {}", e.index(), base_t::tail);
    assert(e.index() == base_t::tail);

    stats.ooocore.branchpred.ras.annuls++;
  }

  //
  // Push the old data back on the stack
  //
  void annulpop(const ReturnAddressStackEntry& old) {
    logging::println(logging::DEBUG, "ReturnAddressStack::annulpop(old index {}, uuid {}, rip {}):", old.idx, old.uuid,
                     (void*)(Waddr)old.rip);

    if (base_t::full()) {
      logging::println(logging::DEBUG, "  Cannot annul: stack is full");
      return;
    }
    ReturnAddressStackEntry dummy;

    logging::println(logging::DEBUG, "  Pushed speculative pop; old.index = {} vs tail {}", old.index(), base_t::tail);
    assert(old.index() == base_t::tail);
    push(old.uuid, old.rip, dummy);
    stats.ooocore.branchpred.ras.annuls++;
  }
};

// std::formatter for ReturnAddressStack<SIZE> - inherits from Queue formatter
} // namespace x86sim

namespace std {
template<int SIZE>
struct formatter<x86sim::ReturnAddressStack<SIZE>> : formatter<x86sim::Queue<x86sim::ReturnAddressStackEntry, SIZE>> {};
} // namespace std

namespace x86sim {

template<int METASIZE, int BIMODSIZE, int L1SIZE, int L2SIZE, int SHIFTWIDTH, bool HISTORYXOR, int BTBSETS, int BTBWAYS,
         int RASSIZE>
struct CombinedPredictor {
  Options& config;

  TwoLevelPredictor<L1SIZE, L2SIZE, SHIFTWIDTH, HISTORYXOR> twolevel;
  BimodalPredictor<BIMODSIZE> bimodal;
  BimodalPredictor<METASIZE> meta;

  BranchTargetBuffer<BTBSETS, BTBWAYS> btb;
  ReturnAddressStack<RASSIZE> ras;

  explicit CombinedPredictor(Options& config_) : config(config_) {}

  void reset() {
    twolevel.reset();
    bimodal.reset();
    meta.reset();
    btb.reset();
    ras.reset();
  }

  void updateras(PredictorUpdate& predinfo, W64 rip) {
    if unlikely (predinfo.flags & BRANCH_HINT_RET) {
      predinfo.ras_push = 0;
      ras.pop(predinfo.ras_old);
    } else if likely (predinfo.flags & BRANCH_HINT_CALL) {
      predinfo.ras_push = 1;
      ras.push(predinfo.uuid, rip, predinfo.ras_old);
    }
  }

  //
  // NOTE: branchaddr should point to first byte *after* branching insn,
  // since x86 has variable length instructions.
  //
  W64 predict(PredictorUpdate& update, int type, W64 branchaddr, W64 target) {
    update.cp1 = null;
    update.cp2 = null;
    update.cpmeta = null;
    update.flags = type;

    if unlikely ((type & (BRANCH_HINT_COND | BRANCH_HINT_INDIRECT)) == 0) {
      // Unconditional: always return target
      return target;
    }

    if likely (type & BRANCH_HINT_COND) {
      byte& bimodalctr = *bimodal.predict(branchaddr);
      byte& twolevelctr = *twolevel.predict(branchaddr);
      byte& metactr = *meta.predict(branchaddr);
      update.cpmeta = &metactr;
      update.meta = (metactr >= 2);
      update.bimodal = (bimodalctr >= 2);
      update.twolevel = (twolevelctr >= 2);
      if (metactr >= 2) {
        update.cp1 = &twolevelctr;
        update.cp2 = &bimodalctr;
      } else {
        update.cp1 = &bimodalctr;
        update.cp2 = &twolevelctr;
      }
    }

    //
    // If this is a return, find next entry that would be popped
    // Caller is responsible for using updateras() to update the
    // RAS once annulable resources have been allocated for this
    // return insn.
    //
    if unlikely (type & BRANCH_HINT_RET) {
      logging::println(logging::DEBUG, "Peeking RAS for uuid {}:", update.uuid);
      return ras.peek();
    }

    BTBEntry* pbtb = btb.probe(branchaddr);

    // if this is a jump, ignore predicted direction; we know it's taken.
    if unlikely (!(type & BRANCH_HINT_COND)) {
      return (pbtb ? pbtb->target : target);
    }

    // Predict conditional branch. If this is the first time, predict backward
    // jumps as taken and forward jumps as not-taken.
    if unlikely (!pbtb && config.debug.static_branchpred) {
      return target < branchaddr ? target : branchaddr;
    }
    return (*(update.cp1) >= 2) ? target : branchaddr;
  }

  void update(PredictorUpdate& update, W64 branchaddr, W64 target) {
    int type = update.flags;

    bool taken = (target != branchaddr);

    //
    // keep stats about JMPs; also, but don't change any pred state for JMPs
    // which are returns.
    //
    if unlikely (type & BRANCH_HINT_INDIRECT) {
      if unlikely (type & BRANCH_HINT_RET)
        return;
    }

    //
    // L1 table is updated unconditionally for combining predictor too:
    //
    if likely (type & BRANCH_HINT_COND) {
      int l1index = lowbits(branchaddr, log2(L1SIZE));
      twolevel.shiftregs[l1index] = lowbits((twolevel.shiftregs[l1index] << 1) | taken, SHIFTWIDTH);
    }

    //
    // Find BTB entry. Allocate always to detect first use of conditional branch
    //
    BTBEntry* pbtb = btb.select(branchaddr);

    //
    // Now p is a possibly null pointer into the direction prediction table,
    // and pbtb is a possibly null pointer into the BTB (either to a
    // matched-on entry or a victim which was LRU in its set)
    //

    //
    // update state (but not for jumps)
    //
    if likely (update.cp1) {
      byte& counter = *update.cp1;
      counter = std::clamp(counter + (taken ? +1 : -1), 0, 3);
    }

    //
    // combining predictor also updates second predictor and meta predictor
    // second direction predictor
    //
    if likely (update.cp2) {
      byte& counter = *update.cp2;
      counter = std::clamp(counter + (taken ? +1 : -1), 0, 3);
    }

    //
    // Update meta predictor
    //
    if likely (update.cpmeta) {
      if (update.bimodal != update.twolevel) {
        //
        // We only update meta predictor if directions were different.
        // We increment the counter if the twolevel predictor was correct;
        // if the bimodal predictor was correct, we decrement it.
        //
        byte& counter = *update.cpmeta;
        bool twolevel_or_bimodal = (update.twolevel == taken);
        counter = std::clamp(counter + (twolevel_or_bimodal ? +1 : -1), 0, 3);
      }
    }

    //
    // update BTB (but only for taken branches)
    //
    if likely (pbtb) {
      // Update either the entry selected above, or if not found, use the LRU entry:
      pbtb->target = target;
    }
  }

  //
  // Speculative execution can corrupt the RAS, since entries will be pushed
  // as call insns are fetched. If those call insns were along an incorrect
  // branch path, they must be annulled.
  //
  void annulras(const PredictorUpdate& predinfo) {
    logging::println(logging::DEBUG, "Update RAS for uuid {}:", predinfo.uuid);
    if (predinfo.ras_push)
      ras.annulpush(predinfo.ras_old);
    else
      ras.annulpop(predinfo.ras_old);
  }
};

// template <int METASIZE, int BIMODSIZE, int L1SIZE, int L2SIZE, int SHIFTWIDTH, bool HISTORYXOR, int BTBSETS, int BTBWAYS, int RASSIZE>
// G-share constraints: METASIZE, BIMODSIZE, 1, L2SIZE, log2(L2SIZE), (HISTORYXOR = true), BTBSETS, BTBWAYS, RASSIZE
struct BranchPredictorImplementation : public CombinedPredictor<65536, 65536, 1, 65536, 16, 1, 1024, 4, 1024> {
  explicit BranchPredictorImplementation(Options& config)
      : CombinedPredictor<65536, 65536, 1, 65536, 16, 1, 1024, 4, 1024>(config) {}
};

void BranchPredictorInterface::destroy() {
  if (impl)
    delete impl;
  impl = null;
}

void BranchPredictorInterface::reset() {
  impl->reset();
}

void BranchPredictorInterface::init() {
  destroy();
  impl = new BranchPredictorImplementation(config);
  reset();
}

W64 BranchPredictorInterface::predict(PredictorUpdate& update, int type, W64 branchaddr, W64 target) {
  return impl->predict(update, type, branchaddr, target);
}

void BranchPredictorInterface::update(PredictorUpdate& update, W64 branchaddr, W64 target) {
  impl->update(update, branchaddr, target);
}

void BranchPredictorInterface::updateras(PredictorUpdate& predinfo, W64 branchaddr) {
  impl->updateras(predinfo, branchaddr);
};

void BranchPredictorInterface::annulras(const PredictorUpdate& predinfo) {
  impl->annulras(predinfo);
};

void BranchPredictorInterface::flush() {}

} // namespace x86sim
