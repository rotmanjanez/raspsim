//
// PTLsim: Cycle Accurate x86-64 Simulator
// Configuration Management
//
// Copyright 2000-2008 Matt T. Yourst <yourst@yourst.com>
//

#include "ptlsim.h"
#include "config.h"
#include "x86sim/logging.h"

namespace logging = x86sim::logging;

ConfigurationParser<x86sim::Options> configparser;

std::string ConfigurationParserBase::format_to_string_usage(const void* baseptr) const {
  std::string result = "Options are:\n";
  ConfigurationOption* option = options;
  int maxlength = 0;
  while (option) {
    if likely (option->type != OPTION_TYPE_SECTION)
      maxlength = std::max(maxlength, int(strlen(option->name)));
    option = option->next;
  }

  option = options;
  while (option) {
    void* variable = (baseptr) ? ((void*)((x86sim::Waddr)baseptr + option->offset)) : nullptr;
    if (option->type == OPTION_TYPE_SECTION) {
      result += std::format("{}:\n", option->description);
      option = option->next;
      continue;
    }

    result += std::format("  -{:<{}} {} [", option->name, maxlength, option->description);

    switch (option->type) {
    case OPTION_TYPE_NONE:
      break;
    case OPTION_TYPE_W64: {
      x86sim::W64 v = *((x86sim::W64*)(variable));
      if (v == infinity)
        result += "inf";
      else
        result += std::format("{}", v);
      break;
    }
    case OPTION_TYPE_FLOAT:
      result += std::format("{}", *((double*)(variable)));
      break;
    case OPTION_TYPE_STRING:
      result += *(std::string*)variable;
      break;
    case OPTION_TYPE_PATH:
      result += ((std::filesystem::path*)variable)->string();
      break;
    case OPTION_TYPE_BOOL:
      result += ((*(bool*)variable) ? "enabled" : "disabled");
      break;
    case OPTION_TYPE_CORE_MODEL:
      result += std::format("{}", *(x86sim::CoreModel*)variable);
      break;
    default:
      assert(false);
    }
    result += "]\n";
    option = option->next;
  }
  result += "\n";

  return result;
}

int ConfigurationParserBase::parse(void* baseptr, int argc, char* argv[]) {
  int i = 0;

  while (i < argc) {
    if ((argv[i][0] == '-') && strlen(argv[i]) > 1) {
      char* name = &argv[i][1];
      i++;
      bool found = false;

      if (0 == strcmp(name, "-"))
        return i;

      ConfigurationOption* option = options;
      while (option) {
        if (option->type == OPTION_TYPE_SECTION) {
          option = option->next;
          continue;
        }
        if (x86sim::strequal(name, option->name)) {
          found = true;
          void* variable = (void*)((x86sim::Waddr)baseptr + option->offset);
          if ((option->type != OPTION_TYPE_NONE) && (option->type != OPTION_TYPE_BOOL) && (i == (argc + 1))) {
            logging::eprintln("Warning: missing value for option '{}'", argv[i - 1]);
            break;
          }
          switch (option->type) {
          case OPTION_TYPE_NONE:
            break;
          case OPTION_TYPE_W64: {
            char* p = (i < argc) ? argv[i] : nullptr;
            int len = (p) ? strlen(p) : 0;
            if (!len) {
              logging::eprintln("Warning: option {} had no argument; ignoring", argv[i - 1]);
              break;
            }

            char buf[256];
            strncpy(buf, p, sizeof(buf));
            p = buf;

            x86sim::W64 multiplier = 1;
            char* endp = p;
            bool isinf = (strncmp(p, "inf", 3) == 0);
            if (len > 1) {
              char& c = p[len - 1];
              switch (c) {
              case 'k':
              case 'K':
                multiplier = 1000LL;
                c = 0;
                break;
              case 'm':
              case 'M':
                multiplier = 1000000LL;
                c = 0;
                break;
              case 'g':
              case 'G':
                multiplier = 1000000000LL;
                c = 0;
                break;
              case 't':
              case 'T':
                multiplier = 1000000000000LL;
                c = 0;
                break;
              }
            }
            x86sim::W64 v = (isinf) ? infinity : strtoull(p, &endp, 0);
            if ((!isinf) && (endp[0] != 0)) {
              logging::eprintln("Warning: invalid value '{}' for option {}; ignoring", p, argv[i - 1]);
            }
            v *= multiplier;
            *((x86sim::W64*)variable) = v;
            i++;

            break;
          }
          case OPTION_TYPE_FLOAT:
            if (i >= argc) {
              logging::eprintln("Warning: option {} had no argument; ignoring", argv[i - 1]);
              break;
            }
            *((double*)variable) = atof(argv[i++]);
            break;
          case OPTION_TYPE_BOOL:
            *((bool*)variable) = (!(*((bool*)variable)));
            break;
          case OPTION_TYPE_STRING: {
            if (i >= argc) {
              logging::eprintln("Warning: option {} had no argument; ignoring", argv[i - 1]);
              break;
            }
            std::string& s = *((std::string*)variable);
            s = argv[i++];
            break;
          }
          case OPTION_TYPE_PATH: {
            if (i >= argc) {
              logging::eprintln("Warning: option {} had no argument; ignoring", argv[i - 1]);
              break;
            }
            std::filesystem::path& path = *((std::filesystem::path*)variable);
            path = argv[i++];
            break;
          }
          case OPTION_TYPE_CORE_MODEL: {
            if (i >= argc) {
              logging::eprintln("Warning: option {} had no argument; ignoring", argv[i - 1]);
              break;
            }
            std::string_view model = argv[i++];
            if (model == "seq" || model == "sequential") {
              *((x86sim::CoreModel*)variable) = x86sim::CoreModel::sequential;
            } else if (model == "ooo" || model == "out_of_order") {
              *((x86sim::CoreModel*)variable) = x86sim::CoreModel::out_of_order;
            } else {
              logging::eprintln("Warning: invalid core model '{}'; ignoring", model);
            }
            break;
          }
          default:
            assert(false);
          }
          break;
        }

        option = option->next;
      }
      if (!found) {
        logging::eprintln("Warning: invalid option '{}'",
                          (x86sim::inrange(i - 1, 0, argc - 1) ? argv[i - 1] : "<missing>"));
        i++;
      }
    } else {
      return i; // trailing arguments, if any
    }
  }

  // no trailing arguments
  return -1;
}

std::string ConfigurationParserBase::format_to_string_config(const void* baseptr) const {
  std::string result = "Active parameters:\n";

  ConfigurationOption* option = options;
  while (option) {
    void* variable = (baseptr) ? ((void*)((x86sim::Waddr)baseptr + option->offset)) : nullptr;

    if (option->type == OPTION_TYPE_SECTION) {
      option = option->next;
      continue;
    }
    result += std::format("  -{:<12} ", option->name);
    switch (option->type) {
    case OPTION_TYPE_NONE:
    case OPTION_TYPE_SECTION:
      break;
    case OPTION_TYPE_W64: {
      x86sim::W64 v = *((x86sim::W64*)(variable));
      if (v == 0) {
        result += "0";
      } else if (v == infinity) {
        result += "infinity";
      } else if ((v % 1000000000LL) == 0) {
        result += std::format("{} G", (v / 1000000000LL));
      } else if ((v % 1000000LL) == 0) {
        result += std::format("{} M", (v / 1000000LL));
      } else {
        result += std::format("{}", v);
      }
      break;
    }
    case OPTION_TYPE_FLOAT:
      result += std::format("{}", *((double*)(variable)));
      break;
    case OPTION_TYPE_BOOL:
      result += (*((bool*)(variable)) ? "enabled" : "disabled");
      break;
    case OPTION_TYPE_STRING:
      result += *((std::string*)(variable));
      break;
    case OPTION_TYPE_PATH:
      result += ((std::filesystem::path*)variable)->string();
      break;
    case OPTION_TYPE_CORE_MODEL:
      result += std::format("{}", *((x86sim::CoreModel*)(variable)));
      break;
    default:
      break;
    }
    result += "\n";

    option = option->next;
  }

  return result;
}



template<>
void ConfigurationParser<x86sim::Options>::setup() {
  section("Simulation Control");

  add(core, "core", "Run using specified core (-core <corename>)");
  add(core_count, "cores", "Number of simulated cores/contexts");

  section("General Logging Control");
  add(quiet, "quiet", "Do not print PTLsim system information banner");
  add(log.log_filename, "logfile", "Log filename (use /dev/fd/1 for stdout, /dev/fd/2 for stderr)");
  add(log.loglevel, "loglevel", "Log level (0 to 99)");
  add(log.start_log_at_iteration, "startlog", "Start logging after iteration <startlog>");
  add(log.start_log_at_rip, "startlogrip", "Start logging after first translation of basic block starting at rip");
  add(log.log_on_console, "consolelog", "Replicate log file messages to console");
  add(log.log_ptlsim_boot, "bootlog", "Log PTLsim early boot and injection process (for debugging)");
  add(log.log_buffer_size, "logbufsize", "Size of PTLsim logfile buffer (not related to -ringbuf)");
  add(debug.dump_state_now, "dump-state-now", "Dump the event log ring buffer and internal state of the active core");
  add(debug.abort_at_end, "abort-at-end", "Abort current simulation after next command (don't wait for next x86 boundary)");
  add(log.mm_logfile, "mm-logfile", "Log PTLsim memory manager requests (alloc, free) to this file (use with ptlmmlog)");
  add(log.mm_log_buffer_size, "mm-logbuf-size", "Size of PTLsim memory manager log buffer (in events, not bytes)");
  add(log.enable_inline_mm_logging, "mm-log-inline", "Print every memory manager request in the main log file");
  add(log.enable_mm_validate, "mm-validate", "Validate every memory manager request against internal structures (slow)");

  section("Event Ring Buffer Logging Control");
  add(log.event_log_enabled, "ringbuf", "Log all core events to the ring buffer for backwards-in-time debugging");
  add(log.event_log_ring_buffer_size, "ringbuf-size", "Core event log ring buffer size: only save last <ringbuf> entries");
  add(log.flush_event_log_every_cycle, "flush-events", "Flush event log ring buffer to logfile after every cycle");
  add(log.log_backwards_from_trigger_rip, "ringbuf-trigger-rip",
      "Print event ring buffer when first uop in this rip is committed");

  // Userspace only
  section("Start Point");
  add(debug.start_at_rip, "startrip", "Start at rip <startrip>");
  add(debug.include_dyn_linker, "excludeld", "Exclude dynamic linker execution");
  add(debug.trigger_mode, "trigger", "Trigger mode: wait for user process to do simcall before entering PTL mode");
  add(debug.pause_at_startup, "pause-at-startup", "Pause for N seconds after starting up (to allow debugger to attach)");

  section("Trace Stop Point");
  add(stop_at_user_insns, "stopinsns", "Stop after executing <stopinsns> user instructions");
  add(stop_at_iteration, "stopiter", "Stop after <stop> iterations (does not apply to cycle-accurate cores)");
  add(stop_at_rip, "stoprip", "Stop before rip <stoprip> is translated for the first time");

  section("Out of Order Core (ooocore)");
  add(debug.perfect_cache, "perfect-cache", "Perfect cache performance: all loads and stores hit in L1");
  add(debug.static_branchpred, "static-branchpred", "Use static branch prediction when there is no BTB entry");

  section("Miscellaneous");
  add(debug.bbcache_dump_filename, "bbdump", "Basic block cache dump filename");
  add(debug.sequential_mode_insns, "seq", "Run in sequential mode for <seq> instructions before switching to out of order");
  add(debug.exit_after_fullsim, "exitend", "Kill the thread after full simulation completes rather than going native");
};

void print_usage(int argc, char** argv) {
  logging::eprintln("Syntax: ptlsim <executable> <arguments...>");
  logging::eprintln("All other options come from file /home/<username>/.ptlsim/path/to/executable");
  logging::eprintln("");
  x86sim::Options defaults;
  logging::eprint("{}", configparser.options.format_to_string_usage(&defaults));
}
