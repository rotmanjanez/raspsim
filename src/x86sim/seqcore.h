// -*- c++ -*-
//
// PTLsim: Cycle Accurate x86-64 Simulator
// Sequential Core Simulator Configuration
//
// Copyright 2004-2008 Matt T. Yourst <yourst@yourst.com>
//

#ifndef _SEQCORE_H_
#define _SEQCORE_H_

#include "ptlsim.h"

#include <memory>
#include <vector>

namespace x86sim {

struct SequentialCore;

struct SequentialMachine : public MachineImpl {
  std::unique_ptr<SequentialCore> core;

  SequentialMachine(Machine& machine, Options config);
  ~SequentialMachine() override;
  std::string_view name() const override;
  Context& cpu_context() override;
  const Context& cpu_context() const override;

  //
  // Run the processor model, until a stopping point
  // is hit (as configured elsewhere in config).
  //
  int run() override;
  void flush_tlb(Context& ctx) override;
  void flush_tlb_virt(Context& ctx, Waddr virtaddr) override;

  virtual void dump_state();

  //
  // Update any statistics in stats in preparation
  // for writing it somewhere. The model may also
  // directly update the global stats structure
  // while it runs; this is only for cleanup tasks
  // or computing derived values.
  //
  virtual void update_stats(PTLsimStats& stats) override;
};

//
// Free-standing sequential execution of one basic block
//
int execute_sequential(Context& ctx);

enum {
  SEQEXEC_OK = 0,
  SEQEXEC_EARLY_EXIT,
  SEQEXEC_SMC,
  SEQEXEC_CHECK,
  SEQEXEC_UNALIGNED,
  SEQEXEC_EXCEPTION,
  SEQEXEC_INVALIDRIP,
  SEQEXEC_SKIPBLOCK,
  SEQEXEC_BARRIER,
  SEQEXEC_INTERRUPT,
  SEQEXEC_RESULT_COUNT,
};

extern const char* seqexec_result_names[SEQEXEC_RESULT_COUNT];

//
// Execute N basic blocks and capture state (but do not commit)
//
static const int MAX_STORES_IN_COMMIT_RECORD = 4096;

struct CommitRecord : public Context {
  int exit_reason;

  int store_count;
  SFR stores[MAX_STORES_IN_COMMIT_RECORD];

  int pte_update_count;
  // PTEUpdate pte_update_list[MAX_STORES_IN_COMMIT_RECORD];
  // Waddr pte_update_virt[MAX_STORES_IN_COMMIT_RECORD];

  void reset() {
    store_count = 0;
    pte_update_count = 0;
    exit_reason = SEQEXEC_OK;
  }
};

int execute_sequential(Context& ctx, CommitRecord* cmtrec = null, W64 bbcount = 1,
                       W64 insncount = std::numeric_limits<W64>::max());

extern W64 suppress_total_user_insn_count_updates_in_seqcore;

} // namespace x86sim

//
// std::formatter specialization
template<>
struct std::formatter<x86sim::CommitRecord> {
  constexpr auto parse(std::format_parse_context& ctx) { return ctx.begin(); }
  auto format(const x86sim::CommitRecord& cr, std::format_context& ctx) const;
};

#endif // _SEQCORE_H_
