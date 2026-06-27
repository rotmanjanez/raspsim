//
// PTLsim: Cycle Accurate x86-64 Simulator
// Shared Functions and Structures
//
// Copyright 2000-2008 Matt T. Yourst <yourst@yourst.com>
//

#include <cstdio>
#include <format>
#include <limits>
#include <memory>
#include <mutex>
#include <string>
#include <fcntl.h>
#include <utility>
#include "globals.h"
#include "ptlsim.h"
#include "decode.h"
#include "ptlhwdef.h"
#include "typedefs.h"
#include "x86sim/logging.hpp"
#include "stats.h"
#include "x86sim/x86sim.hpp"

namespace x86sim {

#ifndef CONFIG_ONLY
//
// Global variables
//
PTLsimStats stats;
#endif

MachineImpl::MachineImpl(Machine& owner_, Options config_)
    : owner(owner_), callbacks(owner_.callbacks_), config(std::move(config_)),
      bbcache(std::make_unique<BasicBlockCache>()) {
  bbcache->machine = this;
}

MachineImpl::~MachineImpl() = default;

Context::Context(const Options& config, MachineImpl& core, int vcpuid_) : Context() {
  machine = &core.owner;
  machine_impl = &core;
  no_sse = !config.sse;
  no_x87 = !config.x87;
  vcpuid = vcpuid_;
}

#ifndef CONFIG_ONLY


std::string current_bbcache_dump_filename;

void backup_and_reopen_logfile(logging::Logger& logger, const Options& options) {
  if (!options.log.log_filename.empty()) {
    std::string log_filename = options.log.log_filename.string();

    // Close existing log
    logger.flush();

    // Backup old log file
    std::string oldname = std::format("{}.backup", log_filename);
    std::remove(oldname.c_str());
    std::rename(log_filename.c_str(), oldname.c_str());

    // Open new log file
    logger.set_file_sink(log_filename.c_str());
  }
}

void force_logging_enabled(logging::Logger& logger, Options& options) {
  logger.set_enabled(true);
  options.log.start_log_at_iteration = 0;
  options.log.loglevel = static_cast<int>(logging::Level::VERBOSE); // Maximum verbosity
  logger.set_level(logging::VERBOSE);
  options.log.flush_event_log_every_cycle = 1;
}

bool handle_config_change(logging::Logger& logger, Options& options, int argc, char** argv) {
  static bool first_time = true;

  // Each Machine owns its logger, so configure this logger's sink directly from
  // its own options. There is deliberately no process-global "file already
  // opened" cache here: that leftover from the single-global-logger design made
  // every machine after the first skip sink setup and fall back to its default
  // (silent) state, while a stale machine's file stayed bound. An empty filename
  // means "no log output at all", which leaves the logger sink-less and silent.
  // Can also use "-logfile /dev/fd/1" to send to stdout (or /dev/fd/2 for stderr).
  if (!options.log.log_filename.empty()) {
    backup_and_reopen_logfile(logger, options);
  }

  // Set log level in the new logging system
  logger.set_level(static_cast<int>(options.log.loglevel));

  if ((options.log.loglevel > 0) & (options.log.start_log_at_rip == INVALIDRIP) &
      (options.log.start_log_at_iteration == std::numeric_limits<W64s>::max())) {
    options.log.start_log_at_iteration = 0;
  }

  // Force printing every cycle if loglevel <= TRACE (5):
  // (Note: new system has inverted semantics - lower numbers = more verbose)
  if (options.log.loglevel <= static_cast<int>(logging::Level::TRACE)) {
    options.log.event_log_enabled = 1;
    options.log.flush_event_log_every_cycle = 1;
  }

  //
  // Fix up parameter defaults:
  //
  if (options.log.start_log_at_rip != INVALIDRIP) {
    options.log.start_log_at_iteration = std::numeric_limits<W64s>::max();
    logger.set_enabled(false);
  } else if (options.log.start_log_at_iteration != std::numeric_limits<W64s>::max()) {
    options.log.start_log_at_rip = INVALIDRIP;
    logger.set_enabled(false);
  }

  logger.set_enabled(true);

  // Unlike the logger, the bb dump file is a single process-global resource
  // shared by every Machine, and current_bbcache_dump_filename caches which file
  // is open so a second machine with the same -bbdump path reuses it instead of
  // truncating it. The std::string cache cannot be atomic, so the compound
  // compare/open/swap is serialized under a mutex; the FILE* itself is
  // std::atomic, keeping the decoder's lock-free hot-path read race-free against
  // this construction-time reconfigure.
  {
    static std::mutex bbcache_mutex;
    const std::scoped_lock lock(bbcache_mutex);
    if (!options.debug.bbcache_dump_filename.empty() &&
        (options.debug.bbcache_dump_filename != current_bbcache_dump_filename)) {
      FILE* dump = std::fopen(options.debug.bbcache_dump_filename.c_str(), "wb");
      if (!dump) {
        logger.println(logging::WARNING, "Cannot open bb dump file '{}'", options.debug.bbcache_dump_filename);
      }
      if (FILE* old = bbcache_dump_file.exchange(dump, std::memory_order_acq_rel)) {
        std::fclose(old);
      }
      current_bbcache_dump_filename = options.debug.bbcache_dump_filename;
    }
  }

  options.log.start_log_at_rip = signext64(options.log.start_log_at_rip, 48);
  options.log.log_backwards_from_trigger_rip = signext64(options.log.log_backwards_from_trigger_rip, 48);
  options.debug.start_at_rip = signext64(options.debug.start_at_rip, 48);
  options.stop_at_rip = signext64(options.stop_at_rip, 48);

  if (first_time) {
    if (!options.quiet) {
      logging::eprint("{}", PTLsimBanner{argc, argv});
    }
    logger.println(logging::INFO, "{}", options);
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
