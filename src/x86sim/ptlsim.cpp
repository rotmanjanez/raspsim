//
// PTLsim: Cycle Accurate x86-64 Simulator
// Shared Functions and Structures
//
// Copyright 2000-2008 Matt T. Yourst <yourst@yourst.com>
//

#include <string>
#include <unordered_map>
#include <fcntl.h>
#include "globals.h"
#include "ptlsim.h"
#include "decode.h"
#include "x86sim/logging.h"
#include "stats.h"

namespace x86sim {

#ifndef CONFIG_ONLY
//
// Global variables
//
PTLsimStats stats;

bool logenable = 0;
#endif

MachineImpl::MachineImpl(Machine& owner_, Options config_)
    : owner(owner_), config(std::move(config_)), bbcache(std::make_unique<BasicBlockCache>()) {}

MachineImpl::~MachineImpl() = default;

int MachineImpl::run() { return 0; }

Context::Context(const Options& config, MachineImpl& core, int vcpuid_) : Context() {
  machine = &core.owner;
  machine_impl = &core;
  no_sse = !config.sse;
  no_x87 = !config.x87;
  vcpuid = vcpuid_;
}

#ifndef CONFIG_ONLY


std::string current_log_filename;
std::string current_bbcache_dump_filename;

void backup_and_reopen_logfile(const Options& options) {
  if (!options.log_filename.empty()) {
    std::string log_filename = options.log_filename.string();

    // Close existing log
    logging::flush();

    // Backup old log file
    std::string oldname = std::format("{}.backup", log_filename);
    std::remove(oldname.c_str());
    std::rename(log_filename.c_str(), oldname.c_str());

    // Open new log file
    logging::set_file_sink(log_filename.c_str());
  }
}

void force_logging_enabled(Options& options) {
  logenable = 1;
  options.start_log_at_iteration = 0;
  options.loglevel = static_cast<int>(logging::Level::VERBOSE); // Maximum verbosity
  logging::set_level(logging::VERBOSE);
  options.flush_event_log_every_cycle = 1;
}

bool handle_config_change(Options& options, int argc, char** argv) {
  static bool first_time = true;

  std::string log_filename = options.log_filename.string();
  if (!options.log_filename.empty() && (log_filename != current_log_filename)) {
    // Can also use "-logfile /dev/fd/1" to send to stdout (or /dev/fd/2 for stderr):
    backup_and_reopen_logfile(options);
    current_log_filename = log_filename;
  }

  // Set log level in the new logging system
  logging::set_level(static_cast<int>(options.loglevel));

  if ((options.loglevel > 0) & (options.start_log_at_rip == INVALIDRIP) &
      (options.start_log_at_iteration == std::numeric_limits<W64s>::max())) {
    options.start_log_at_iteration = 0;
  }

  // Force printing every cycle if loglevel <= TRACE (5):
  // (Note: new system has inverted semantics - lower numbers = more verbose)
  if (options.loglevel <= static_cast<int>(logging::Level::TRACE)) {
    options.event_log_enabled = 1;
    options.flush_event_log_every_cycle = 1;
  }

  //
  // Fix up parameter defaults:
  //
  if (options.start_log_at_rip != INVALIDRIP) {
    options.start_log_at_iteration = std::numeric_limits<W64s>::max();
    logenable = 0;
  } else if (options.start_log_at_iteration != std::numeric_limits<W64s>::max()) {
    options.start_log_at_rip = INVALIDRIP;
    logenable = 0;
  }

  logenable = 1;

  if (!options.bbcache_dump_filename.empty() && (options.bbcache_dump_filename != current_bbcache_dump_filename)) {
    // Can also use "-logfile /dev/fd/1" to send to stdout (or /dev/fd/2 for stderr):
    if (bbcache_dump_file)
      std::fclose(bbcache_dump_file);
    bbcache_dump_file = std::fopen(options.bbcache_dump_filename.c_str(), "wb");
    if (!bbcache_dump_file) {
      logging::println(logging::WARNING, "Cannot open bb dump file '{}'", options.bbcache_dump_filename);
    }
    current_bbcache_dump_filename = options.bbcache_dump_filename;
  }

  options.start_log_at_rip = signext64(options.start_log_at_rip, 48);
  options.log_backwards_from_trigger_rip = signext64(options.log_backwards_from_trigger_rip, 48);
  options.start_at_rip = signext64(options.start_at_rip, 48);
  options.stop_at_rip = signext64(options.stop_at_rip, 48);

  if (first_time) {
    if (!options.quiet) {
      logging::eprint("{}", PTLsimBanner{argc, argv});
    }
    logging::println(logging::INFO, "{}", options);
    first_time = false;
  }

  return true;
}

extern void shutdown_uops();

void shutdown_subsystems() {
  //
  // Let the subsystems close any special files or buffers
  // they may have open:
  //
  shutdown_uops();
  shutdown_decode();
}

#endif // CONFIG_ONLY

} // namespace x86sim
