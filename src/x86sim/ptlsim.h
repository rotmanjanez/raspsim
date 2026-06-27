// -*- c++ -*-
//
// PTLsim: Cycle Accurate x86-64 Simulator
// Simulator Structures
//
// Copyright 2003-2008 Matt T. Yourst <yourst@yourst.com>
//

#ifndef _PTLSIM_H_
#define _PTLSIM_H_

#include <cstddef>
#include <format>
#include <memory>
#include <string_view>

#include "globals.h"
#include "ptlhwdef.h"
#include "x86sim/x86sim.hpp"

namespace x86sim {

void user_process_terminated(int rc);

// Compile-time host platform name; replaces runtime uname(), which is not
// available on all targets (e.g. wasm)
constexpr const char* host_platform_name() {
#if defined(__EMSCRIPTEN__)
  return "emscripten";
#elif defined(__wasi__)
  return "wasi";
#elif defined(__APPLE__)
  return "darwin";
#elif defined(__linux__)
  return "linux";
#else
  return "unknown";
#endif
}

#ifndef BUILDHOST
#define BUILDHOST unknown
#endif
#ifndef GIT_REVISION
#define GIT_REVISION unknown
#endif

// Compile-time build provenance; single source of truth for the startup
// banner and PTLsimStats
namespace build_info {
inline constexpr std::string_view timestamp = __DATE__ " " __TIME__;
inline constexpr std::string_view hostname = stringify(BUILDHOST);
inline constexpr std::string_view git_revision = stringify(GIT_REVISION);
inline constexpr std::string_view compiler =
#if defined(__clang__)
    "clang-" stringify(__clang_major__) "." stringify(__clang_minor__) "." stringify(__clang_patchlevel__);
#elif defined(__GNUC__)
    "gcc-" stringify(__GNUC__) "." stringify(__GNUC_MINOR__) "." stringify(__GNUC_PATCHLEVEL__);
#else
    "unknown";
#endif
} // namespace build_info

static const int MAX_TRANSOP_BUFFER_SIZE = 4;

struct PTLsimStats;
struct BasicBlockCache;

struct MachineImpl {
  Machine& owner;
  // The host's syscall/cpuid callbacks (owner.callbacks_), reachable directly
  // from the core's dispatch path without a public Machine method.
  HostCallbacks& callbacks;
  Options config;
  std::unique_ptr<BasicBlockCache> bbcache;

  W64 sim_cycle = 0;
  W64 unhalted_cycle_count = 0;
  W64 iterations = 0;
  W64 total_uops_executed = 0;
  W64 total_uops_committed = 0;
  W64 total_user_insns_committed = 0;
  W64 total_basic_blocks_committed = 0;

  // The architectural state and address space the current run executes against.
  // Machine::run points these at the caller's CpuState/AddressSpace before
  // delegating to run(); each core loads its single internal context from them
  // at the start of run() and stores the committed state back at the end.
  CpuState* state = nullptr;
  AddressSpace* address_space = nullptr;

  // Set by the core's syscall/cpuid dispatch when the host asks execution to
  // stop; Machine::run consumes it and prefers it over a plain instruction
  // limit. Reset at the start of each run.
  std::optional<RunResult> pending_stop;

  explicit MachineImpl(Machine& owner_, Options config_);
  virtual ~MachineImpl();
  virtual std::string_view name() const = 0;

  // The core's single internal execution context: scratch holding the live
  // register file, the REG_ctx self-pointer and the uop helper methods.
  [[nodiscard]] virtual Context& cpu_context() = 0;
  [[nodiscard]] virtual const Context& cpu_context() const = 0;

  // Run the model until a stopping point, executing *state against
  // *address_space and writing the committed architectural state back to *state.
  virtual int run() = 0;
  virtual void update_stats(PTLsimStats& stats) = 0;
  virtual void flush_tlb(Context& ctx) = 0;
  virtual void flush_tlb_virt(Context& ctx, Waddr virtaddr) = 0;
};

struct TransOpBuffer {
  TransOp uops[MAX_TRANSOP_BUFFER_SIZE];
  UopImpl synthops[MAX_TRANSOP_BUFFER_SIZE];
  int index;
  int count;

  bool get(TransOp& uop, UopImpl& synthop) {
    if (!count)
      return false;
    uop = uops[index];
    synthop = synthops[index];
    index++;
    if (index >= count) {
      count = 0;
      index = 0;
    }
    return true;
  }

  void reset() {
    index = 0;
    count = 0;
  }

  int put() { return count++; }

  bool empty() const { return (count == 0); }

  TransOpBuffer() { reset(); }
};

void split_unaligned(const TransOp& transop, TransOpBuffer& buf);


void backup_and_reopen_logfile(const Options& options);
void shutdown_subsystems();

void update_progress();

struct PTLsimBanner {
  int argc = 0;
  char** argv = null;
};

extern "C" void switch_to_sim();

//
// uop implementations
//

struct AddrPair {
  byte* start;
  byte* end;
};

void init_uops();
void shutdown_uops();
UopImpl get_synthcode_for_uop(int op, int size, bool setflags, int cond, int extshift, bool except, bool internal);
UopImpl get_synthcode_for_cond_branch(int opcode, int cond, int size, bool except);
void synth_uops_for_bb(BasicBlock& bb);
struct PTLsimStats;

#define INVALIDRIP 0xffffffffffffffffULL

extern bool logenable;
void force_logging_enabled(Options& options);
bool handle_config_change(Options& options, int argc, char** argv);

} // namespace x86sim

namespace std {

template<>
struct formatter<x86sim::PTLsimBanner> {
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }
  auto format(const x86sim::PTLsimBanner& banner, format_context& ctx) const {
    auto out = std::format_to(ctx.out(), "//  \n");
    out = std::format_to(out, "//  PTLsim: Cycle Accurate x86-64 Simulator\n");
    out = std::format_to(out, "//  Copyright 1999-2007 Matt T. Yourst <yourst@yourst.com>\n");
    out = std::format_to(out, "// \n");
    out = std::format_to(out, "//  Built {} on {} using {}\n", x86sim::build_info::timestamp,
                         x86sim::build_info::hostname, x86sim::build_info::compiler);
    out = std::format_to(out, "//  Revision {}\n", x86sim::build_info::git_revision);
    out = std::format_to(out, "//  Running on {}\n", x86sim::host_platform_name());
    out = std::format_to(out, "//  \n");

    out = std::format_to(out, "//  Arguments:\n");
    for (int i = 0; i < banner.argc; i++) {
      const char* arg = (banner.argv && banner.argv[i]) ? banner.argv[i] : "";
      out = std::format_to(out, "{} \n", arg);
    }

    return std::format_to(out, "\n");
  }
};
} // namespace std

#endif // _PTLSIM_H_
