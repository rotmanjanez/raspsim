//
// PTLsim: Cycle Accurate x86-64 Simulator
// RASPsim application
//
// Copyright 2020-2020 Alexis Engelke <engelke@in.tum.de>
//

#include <globals.h>
#include <superstl.h>
#include <mm.h>

#include <elf.h>

#include <ptlsim.h>
#include <ptlsim-api.h>
#include <ptlhwdef.h>
#include <config.h>
#include <stats.h>
#include <raspsim-hwsetup.h>

struct PTLsimConfig;
extern PTLsimConfig config;

extern ConfigurationParser<PTLsimConfig> configparser;



extern "C" void assert_fail(const char *__assertion, const char *__file, unsigned int __line, const char *__function) {
  stringbuf sb;
  sb << "Assert ", __assertion, " failed in ", __file, ":", __line, " (", __function, ") at ", sim_cycle, " cycles, ", iterations, " iterations, ", total_user_insns_committed, " user commits", endl;

  cerr << sb, flush;

  if (logfile) {
    logfile << sb, flush;
    PTLsimMachine* machine = PTLsimMachine::getcurrent();
    if (machine) machine->dump_state(logfile);
    logfile.close();
  }

  sys_exit(1); // Well, we don't want core dumps.

  // Crash and make a core dump:
  asm("ud2a");
  abort();
}

// Saved and restored by asm code:
FXSAVEStruct x87state;
W16 saved_cs;
W16 saved_ss;
W16 saved_ds;
W16 saved_es;
W16 saved_fs;
W16 saved_gs;

void Raspsim::propagate_x86_exception(byte exception, W32 errorcode, Waddr virtaddr) {
  Context& ctx{Raspsim::getContext()};

  const char* msg = Raspsim::formatException(exception, errorcode, virtaddr);
  logfile << msg, endl, flush;
  cerr << msg, endl, flush;
  free((void*)msg);

  cerr << "End state:", endl;
  cerr << ctx, endl;
  sys_exit(1);
}

#ifdef __x86_64__
//
// SYSCALL instruction from x86-64 mode
//
void Raspsim::handle_syscall_64bit() {
  bool DEBUG = 1; //analyze_in_detail();
  //
  // Handle an x86-64 syscall:
  // (This is called from the assist_syscall ucode assist)
  //
  Context& ctx{Raspsim::getContext()};

  size_t syscallid = ctx.commitarf[REG_rax];
  W64 arg1 = ctx.commitarf[REG_rdi];
  W64 arg2 = ctx.commitarf[REG_rsi];
  W64 arg3 = ctx.commitarf[REG_rdx];
  W64 arg4 = ctx.commitarf[REG_r10];
  W64 arg5 = ctx.commitarf[REG_r8];
  W64 arg6 = ctx.commitarf[REG_r9];

  if (DEBUG)
    logfile << "handle_syscall -> (#", syscallid, " ", ((syscallid < lengthofSyscallNames()) ? syscall_names_64bit[syscallid] : "???"),
      ") from ", (void*)ctx.commitarf[REG_rcx], " args ", " (", (void*)arg1, ", ", (void*)arg2, ", ", (void*)arg3, ", ", (void*)arg4, ", ",
      (void*)arg5, ", ", (void*)arg6, ") at iteration ", iterations, endl, flush;

  ctx.commitarf[REG_rax] = -ENOSYS;
  ctx.commitarf[REG_rip] = ctx.commitarf[REG_rcx];

  if (DEBUG) logfile << "handle_syscall: result ", ctx.commitarf[REG_rax], " (", (void*)ctx.commitarf[REG_rax], "); returning to ", (void*)ctx.commitarf[REG_rip], endl, flush;
}

#endif // __x86_64__

void Raspsim::handle_syscall_32bit(int semantics) {
  Context& ctx{Raspsim::getContext()};

  bool DEBUG = 1; //analyze_in_detail();
  //
  // Handle a 32-bit syscall:
  // (This is called from the assist_syscall ucode assist)
  //
  // All 32-bit syscalls are unsupported in simulator-only mode
  ctx.commitarf[REG_rax] = -ENOSYS;

  ctx.commitarf[REG_rip] = ctx.commitarf[REG_nextrip];
}


bool handle_config_arg(Raspsim& sim, char* line, dynarray<Waddr>* dump_pages) {
  if (*line == '\0') return false;
  dynarray<char*> toks;
  toks.tokenize(line, " ");
  if (toks.empty())
    return false;

  if (toks[0][0] == '#') {
    return false;
  }

  if (toks[0][0] == 'M') { // allocate page M<addr> <prot>
    if (toks.size() != 2) {
      cerr << "Error: option ", line, " has wrong number of arguments", endl;
      return true;
    }
    char* endp;
    W64 addr = strtoull(toks[0] + 1, &endp, 16);
    if (*endp != '\0' || lowbits(addr, 12)) {
      cerr << "Error: invalid value ", toks[0], " ", endp, endl;
      return true;
    }
    int prot = 0;
    if (!strcmp(toks[1], "ro")) prot = PROT_READ;
    else if (!strcmp(toks[1], "rw")) prot = PROT_READ | PROT_WRITE;
    else if (!strcmp(toks[1], "rx")) prot = PROT_READ | PROT_EXEC;
    else if (!strcmp(toks[1], "rwx")) prot = PROT_READ | PROT_WRITE | PROT_EXEC;
    else {
      cerr << "Error: invalid mem prot ", toks[1], endl;
      return true;
    }
    sim.map(addr, Raspsim::getPageSize(), prot);
  } else if (toks[0][0] == 'W') { // write to mem W<addr> <hexbytes>, may not cross page boundaries
    if (toks.size() != 2) {
      cerr << "Error: option ", line, " has wrong number of arguments", endl;
      return true;
    }
    char* endp;
    W64 addr = strtoull(toks[0] + 1, &endp, 16);
    if (*endp != '\0') {
      cerr << "Error: invalid value ", toks[0], endl;
      return true;
    }
    byte* mapped = sim.getMappedPage(addr);
    if (!mapped) {
      cerr << "Error: page not mapped ", (void*) addr, endl;
      return true;
    }
    Waddr arglen = strlen(toks[1]);
    if ((arglen & 1) || arglen/2 > 4096-lowbits(addr, 12)) {
      cerr << "Error: arg has odd size or crosses page boundary", (void*) addr, endl;
      return true;
    }
    unsigned n = min((Waddr)(4096 - lowbits(addr, 12)), arglen/2);
    foreach (i, n) {
      char hex_byte[3] = {toks[1][i*2],toks[1][i*2+1], 0};
      mapped[i] = strtoul(hex_byte, NULL, 16);
    }
  } else if (toks[0][0] == 'D') { // dump page D<page>
    if (toks.size() != 1) {
      cerr << "Error: option ", line, " has wrong number of arguments", endl;
      return true;
    }
    char* endp;
    W64 addr = strtoull(toks[0] + 1, &endp, 16);
    if (*endp != '\0') {
      cerr << "Error: invalid value ", toks[0], endl;
      return true;
    }
    dump_pages->push(floor(addr, PAGE_SIZE));
  } else if (!strcmp(toks[0], "Fnox87")) {
    sim.disableX87();
  } else if (!strcmp(toks[0], "Fnosse")) {
    sim.disableSSE();
  } else if (!strcmp(toks[0], "Fnocache")) {
    sim.enablePerfectCache();
  } else if (!strcmp(toks[0], "Fstbrpred")) {
    sim.enableStaticBranchPrediction();
  } else {
    if (toks.size() != 2) {
      cerr << "Error: option ", line, " has wrong number of arguments", endl;
      return true;
    }
    int reg = sim.getRegisterIndex(toks[0]);
    if (reg < 0) {
      cerr << "Error: invalid register ", toks[0], endl;
      return true;
    }
    char* endp;
    W64 v = strtoull(toks[1], &endp, 0);
    if (*endp != '\0') {
      cerr << "Error: invalid value ", toks[1], endl;
      return true;
    }
    sim.setRegisterValue(reg, v);
  }

  return false;
}

//
// PTLsim main: called after ptlsim_preinit() brings up boot subsystems
//
int main(int argc, char** argv) {    
 configparser.setup();
 config.reset();

 int ptlsim_arg_count = 1 + configparser.parse(config, argc-1, argv+1);
 if (ptlsim_arg_count == 0) ptlsim_arg_count = argc;
 handle_config_change(config, ptlsim_arg_count - 1, argv+1);

 Raspsim sim{};
 dynarray<Waddr> dump_pages;
 // TODO(AE): set seccomp filter before parsing arguments
 bool parse_err = false;
 for (unsigned i = ptlsim_arg_count; i < argc; i++) {
  if (argv[i][0] == '@') {
    stringbuf line;
    istream is(argv[i] + 1);
    if (!is) {
      cerr << "Warning: cannot open command list file '", argv[i], "'", endl;
      continue;
    }  
    for (;;) {
      line.reset();
      if (!is) break;
      is >> line;  
      char* p = strchr(line, '#');
      if (p) *p = 0;
      parse_err |= handle_config_arg(sim, line, &dump_pages);
    }
  } else {
    parse_err |= handle_config_arg(sim, argv[i], &dump_pages);
  }
}  
 if (parse_err) {
   cerr << "Error: could not parse all arguments", endl, flush;
   sys_exit(1);
 }

 logfile << endl,  "=== Switching to simulation mode at rip ", (void*)(Waddr) sim.getRegisterValue(REG_rip), " ===", endl, endl, flush;
 logfile << "Baseline state:", endl;
 logfile << sim.getContext();

 sim.run();

 cerr << "End state:", endl;
 cerr << sim.getContext(), endl;

 foreach (i, dump_pages.length) {
   Waddr addr = dump_pages[i];
   byte* mapped = sim.getMappedPage(addr);
   if (!mapped) {
     cerr << "Error dumping memory: page not mapped ", (void*) addr, endl;
   } else {
     cerr << "Dump of memory at ", (void*) addr, ": ", endl;
     cerr << bytestring(mapped, PAGE_SIZE), endl;
   }
 }
 
 cerr << "Decoder stats:";
 foreach(i, DECODE_TYPE_COUNT) {
   cerr << " ", decode_type_names[i], "=", stats.decoder.x86_decode_type[i];
 }
 cerr << endl;
 cerr << flush;

 cerr << endl, "=== Exiting after full simulation on tid ", sys_gettid(), " at rip ", (void*)(Waddr) sim.getRegisterValue(REG_rip), " (",
   sim_cycle, " cycles, ", total_user_insns_committed, " user commits, ", iterations, " iterations) ===", endl, endl;

 Raspsim::stutdown();

 sys_exit(0);
}
