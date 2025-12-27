#include <raspsim-hwsetup.h>
#include <addrspace.h>
#include <ptlsim.h>

struct PTLsimConfig;
extern PTLsimConfig config;

static Context ctx alignto(4096) insection(".ctx"){};
static AddressSpace asp{};

// Raspsim class implementation

Raspsim::Raspsim() {
  config.reset();

  init_uops();
// Set up initial context:
  ctx.reset();
  asp.reset();
  ctx.use32 = 1;
  ctx.use64 = 1;
  ctx.commitarf[REG_rsp] = 0;
  ctx.commitarf[REG_rip] = 0x100000;
  ctx.commitarf[REG_flags] = 0;
  ctx.internal_eflags = 0;

  ctx.seg[SEGID_CS].selector = 0x33;
  ctx.seg[SEGID_SS].selector = 0x2b;
  ctx.seg[SEGID_DS].selector = 0x00;
  ctx.seg[SEGID_ES].selector = 0x00;
  ctx.seg[SEGID_FS].selector = 0x00;
  ctx.seg[SEGID_GS].selector = 0x00;
  ctx.update_shadow_segment_descriptors();

  // ctx.fxrstor(x87state);

  ctx.vcpuid = 0;
  ctx.running = 1;
  ctx.commitarf[REG_ctx] = (Waddr)&ctx;
  ctx.commitarf[REG_fpstack] = (Waddr)&ctx.fpstack;

  //
  // Swap the FP control registers to the user process version, so FP uopimpls
  // can use the real rounding control bits.
  //
  x86_set_mxcsr(ctx.mxcsr | MXCSR_EXCEPTION_DISABLE_MASK);
}

Raspsim::~Raspsim() {
    sim_cycle = 0;
    unhalted_cycle_count = 0;
    iterations = 0;
    total_uops_executed = 0;
    total_uops_committed = 0;
    total_user_insns_committed = 0;
    total_basic_blocks_committed = 0;

    getMachine()->reset();
}

PTLsimMachine* Raspsim::getMachine() {
  return PTLsimMachine::getmachine(config.core_name);
}

const char* Raspsim::getCoreName() { return config.core_name; }

void Raspsim::setLogfile(const char* filename) {
  config.log_filename = filename;
  backup_and_reopen_logfile();
}

Waddr Raspsim::getPageSize() { return PAGE_SIZE; }

AddressSpace& Raspsim::getAddrspace() { return asp; }

Context& Raspsim::getContext() { return ctx; }

W64 Raspsim::getRegisterValue(int reg) { return ctx.commitarf[reg]; }

void Raspsim::map(Waddr start, W64 size, int prot) {
  asp.map(start, size, prot);
}

void Raspsim::unmap(Waddr start, W64 size) {
  asp.unmap(start, size);
}

void* Raspsim::page_virt_to_mapped(Waddr start) {
  return (byte*)asp.page_virt_to_mapped(start);
}

void Raspsim::disableSSE() { ctx.no_sse = 1; }
void Raspsim::disableX87() { ctx.no_x87 = 1; }
void Raspsim::enablePerfectCache() { config.perfect_cache = 1; }
void Raspsim::enableStaticBranchPrediction() { config.static_branchpred = 1; }
void Raspsim::setTimeout(W64 ninstr) { config.stop_at_user_insns = ninstr; }

int Raspsim::getRegisterIndex(const char* regname) {
  int reg = -1;
  foreach (j, sizeof(arch_reg_names) / sizeof(arch_reg_names[0])) {
    if (!strcmp(regname, arch_reg_names[j])) {
      reg = j; break;
    }
  }
  return reg;
}

byte* Raspsim::getMappedPage(Waddr addr) {
  assert(addr % PAGE_SIZE == 0);
  return (byte*)Raspsim::getAddrspace().page_virt_to_mapped(addr);
}

W64 Raspsim::cycles() { return sim_cycle; }
W64 Raspsim::instructions() { return total_user_insns_committed; }

void Raspsim::setRegisterValue(int reg, W64 value) { ctx.commitarf[reg] = value; }

void Raspsim::stutdown() {
  shutdown_subsystems();
  logfile.flush();
  cerr.flush();
}

void Raspsim::run() {
  simulate(config.core_name);
}

char* Raspsim::formatException(byte exception, W32 errorcode, Waddr virtaddr) {
  stringbuf buf{};

  buf << "Exception ", exception, " (", x86_exception_names[exception], ") code=", errorcode, " addr=", (void*)virtaddr, " @ rip ", (void*)(Waddr)ctx.commitarf[REG_rip], " (", total_user_insns_committed, " commits, ", sim_cycle, " cycles)", endl;

  // PF
  if (exception == 14) {
    // PF Flags
    W8 p    = errorcode & 0x00000001;
    W8 wr   = errorcode & 0x00000002;
    W8 us   = errorcode & 0x00000004;
    W8 rsvd = errorcode & 0x00000008;
    W8 id   = errorcode & 0x00000010;
    W8 pk   = errorcode & 0x00000020;

    buf << "PageFault error code: 0x", hexstring(errorcode, 32), ", Flags: ", (pk ? "PK " : ""), (id ? "I " : "D "), (rsvd ? "RSVD " : ""), (us ? "U " : "S "), (wr ? "W " : "R "), (p ? "P" : ""), endl;
  }
  return strdup((char*)buf);
}

char* Raspsim::formatContext(const Context& ctx) {
  stringbuf buf{};

  buf << ctx;

  return strdup((char*)buf);
}


int Raspsim::getPageProtection(void* addr) {
  return asp.getattr(addr);
}

// Begin Virtual Hardware Setup for RASPsim

// Userspace PTLsim only supports one VCPU:
int current_vcpuid() { return 0; }

bool asp_check_exec(void* addr) { return Raspsim::getAddrspace().fastcheck(addr, Raspsim::getAddrspace().execmap); }

bool smc_isdirty(Waddr mfn) { return Raspsim::getAddrspace().isdirty(mfn); }
void smc_setdirty(Waddr mfn) { Raspsim::getAddrspace().setdirty(mfn); }
void smc_cleardirty(Waddr mfn) { Raspsim::getAddrspace().cleardirty(mfn); }

bool check_for_async_sim_break() { return iterations >= config.stop_at_iteration || total_user_insns_committed >= config.stop_at_user_insns; }

int inject_events() { return 0; }
void print_sysinfo(ostream& os) {}

// Only one VCPU in userspace PTLsim:
Context& contextof(int vcpu) { return Raspsim::getContext(); }

W64 loadphys(Waddr addr) {
  W64& data = *(W64*)addr;
  return data;
}

W64 storemask(Waddr addr, W64 data, byte bytemask) {
  W64& mem = *(W64*)addr;
  mem = mux64(expand_8bit_to_64bit_lut[bytemask], mem, data);
  return data;
}

int Context::copy_from_user(void* target, Waddr addr, int bytes, PageFaultErrorCode& pfec, Waddr& faultaddr, bool forexec, Level1PTE& ptelo, Level1PTE& ptehi) {
  logfile << "VMEM: Read from user ", (void*)addr, " (", bytes, ")", endl, flush;

  bool readable;
  bool executable;

  int n = 0;
  pfec = 0;

  ptelo = 0;
  ptehi = 0;

  readable = Raspsim::getAddrspace().fastcheck((byte*)addr, asp.readmap);
  if likely (forexec) executable = asp.fastcheck((byte*)addr, asp.execmap);
  if unlikely ((!readable) | (forexec & !executable)) {
    faultaddr = addr;
    pfec.p = readable;
    pfec.nx = (forexec & (!executable));
    pfec.us = 1;
    return n;
  }

  n = min((Waddr)(4096 - lowbits(addr, 12)), (Waddr)bytes);

  void* mapped_addr = asp.page_virt_to_mapped(addr);
  assert(mapped_addr);
  logfile << "VMEM: Read ", mapped_addr, " = ", *(W8*)mapped_addr, "[", hexstring(*(W8*)mapped_addr, 64), "]", endl, flush;
  memcpy(target, mapped_addr, n);

  // All the bytes were on the first page
  if likely (n == bytes) return n;

  // Go on to second page, if present
  readable = asp.fastcheck((byte*)(addr + n), asp.readmap);
  if likely (forexec) executable = asp.fastcheck((byte*)(addr + n), asp.execmap);
  if unlikely ((!readable) | (forexec & !executable)) {
    faultaddr = addr + n;
    pfec.p = readable;
    pfec.nx = (forexec & (!executable));
    pfec.us = 1;
    return n;
  }

  memcpy((byte*)target + n, asp.page_virt_to_mapped(addr + n), bytes - n);
  return bytes;
}

int Context::copy_to_user(Waddr target, void* source, int bytes, PageFaultErrorCode& pfec, Waddr& faultaddr) {
  logfile << "VMEM: Write to user ", (void*)target, " (", bytes, ")", endl, flush;

  pfec = 0;
  bool writable = asp.fastcheck((byte*)target, asp.writemap);
  if unlikely (!writable) {
    faultaddr = target;
    pfec.p = asp.fastcheck((byte*)target, asp.readmap);
    pfec.rw = 1;
    return 0;
  }

  byte* targetlo = (byte*)asp.page_virt_to_mapped(target);
  int nlo = min((Waddr)(4096 - lowbits(target, 12)), (Waddr)bytes);

  smc_setdirty(target >> 12);

  // All the bytes were on the first page
  if likely (nlo == bytes) {
    memcpy(targetlo, source, nlo);
    return bytes;
  }

  // Go on to second page, if present
  writable = asp.fastcheck((byte*)(target + nlo), asp.writemap);
  if unlikely (!writable) {
    faultaddr = target + nlo;
    pfec.p = asp.fastcheck((byte*)(target + nlo), asp.readmap);
    pfec.rw = 1;
    pfec.us = 1;
    return nlo;
  }

  memcpy(asp.page_virt_to_mapped(target + nlo), (byte*)source + nlo, bytes - nlo);
  memcpy(targetlo, source, nlo);

  smc_setdirty((target + nlo) >> 12);

  return bytes;
}

Waddr Context::check_and_translate(Waddr virtaddr, int sizeshift, bool store, bool internal, int& exception, PageFaultErrorCode& pfec, PTEUpdate& pteupdate, Level1PTE& pteused) {
  exception = 0;
  pteupdate = 0;
  pteused = 0;
  pfec = 0;

  if unlikely (lowbits(virtaddr, sizeshift)) {
    exception = EXCEPTION_UnalignedAccess;
    return INVALID_PHYSADDR;
  }

  if unlikely (internal) {
    // Directly mapped to PTL space:
    return virtaddr;
  }

  AddressSpace::spat_t top = (store) ? asp.writemap : asp.readmap;

  if unlikely (!asp.fastcheck(virtaddr, top)) {
    exception = (store) ? EXCEPTION_PageFaultOnWrite : EXCEPTION_PageFaultOnRead;
    pfec.p = asp.fastcheck(virtaddr, asp.readmap);
    pfec.rw = store;
    pfec.us = 1;
    return 0;
  }

  return (Waddr) asp.page_virt_to_mapped(floor(signext64(virtaddr, 48), 8));
}

int Context::write_segreg(unsigned int segid, W16 selector) {
  // Well, we don't want to play with the fire...
  return EXCEPTION_x86_gp_fault;
}

void Context::update_shadow_segment_descriptors() {
  W64 limit = (use64) ? 0xffffffffffffffffULL : 0xffffffffULL;

  SegmentDescriptorCache& cs = seg[SEGID_CS];
  cs.present = 1;
  cs.base = 0;
  cs.limit = limit;

  virt_addr_mask = limit;

  SegmentDescriptorCache& ss = seg[SEGID_SS];
  ss.present = 1;
  ss.base = 0;
  ss.limit = limit;

  SegmentDescriptorCache& ds = seg[SEGID_DS];
  ds.present = 1;
  ds.base = 0;
  ds.limit = limit;

  SegmentDescriptorCache& es = seg[SEGID_ES];
  es.present = 1;
  es.base = 0;
  es.limit = limit;

  SegmentDescriptorCache& fs = seg[SEGID_FS];
  fs.present = 1;
  fs.base = 0;
  fs.limit = limit;

  SegmentDescriptorCache& gs = seg[SEGID_GS];
  gs.present = 1;
  gs.base = 0;
  gs.limit = limit;
}

// In userspace PTLsim, virtual == physical:
// FIXME(AE): software virtual memory
RIPVirtPhys& RIPVirtPhys::update(Context& ctx, int bytes) {
  use64 = ctx.use64;
  kernel = 0;
  df = ((ctx.internal_eflags & FLAG_DF) != 0);
  padlo = 0;
  padhi = 0;
  mfnlo = rip >> 12;
  mfnhi = (rip + (bytes-1)) >> 12;
  return *this;
}


void Context::propagate_x86_exception(byte exception, W32 errorcode, Waddr virtaddr) {
    Raspsim::propagate_x86_exception(exception, errorcode, virtaddr);
}

#ifdef __x86_64__

void handle_syscall_64bit() {
  Raspsim::handle_syscall_64bit();
}

#endif // __x86_64__
  
void handle_syscall_32bit(int semantics) {
  Raspsim::handle_syscall_32bit(semantics);
}


const char* const syscall_names_64bit[] = {
  "read", "write", "open", "close", "stat", "fstat", "lstat", "poll", "lseek", "mapPage", "mprotect", "munmap", "brk", "rt_sigaction", "rt_sigprocmask", "rt_sigreturn", "ioctl", "pread64", "pwrite64", "readv", "writev", "access", "pipe", "select", "sched_yield", "mremap", "msync", "mincore", "madvise", "shmget", "shmat", "shmctl", "dup", "dup2", "pause", "nanosleep", "getitimer", "alarm", "setitimer", "getpid", "sendfile", "socket", "connect", "accept", "sendto", "recvfrom", "sendmsg", "recvmsg", "shutdown", "bind", "listen", "getsockname", "getpeername", "socketpair", "setsockopt", "getsockopt", "clone", "fork", "vfork", "execve", "exit", "wait4", "kill", "uname", "semget", "semop", "semctl", "shmdt", "msgget", "msgsnd", "msgrcv", "msgctl", "fcntl", "flock", "fsync", "fdatasync", "truncate", "ftruncate", "getdents", "getcwd", "chdir", "fchdir", "rename", "mkdir", "rmdir", "creat", "link", "unlink", "symlink", "readlink", "chmod", "fchmod", "chown", "fchown", "lchown", "umask", "gettimeofday", "getrlimit", "getrusage", "sysinfo", "times", "ptrace", "getuid", "syslog", "getgid", "setuid", "setgid", "geteuid", "getegid", "setpgid", "getppid", "getpgrp", "setsid", "setreuid", "setregid", "getgroups", "setgroups", "setresuid", "getresuid", "setresgid", "getresgid", "getpgid", "setfsuid", "setfsgid", "getsid", "capget", "capset", "rt_sigpending", "rt_sigtimedwait", "rt_sigqueueinfo", "rt_sigsuspend", "sigaltstack", "utime", "mknod", "uselib", "personality", "ustat", "statfs", "fstatfs", "sysfs", "getpriority", "setpriority", "sched_setparam", "sched_getparam", "sched_setscheduler", "sched_getscheduler", "sched_get_priority_max", "sched_get_priority_min", "sched_rr_get_interval", "mlock", "munlock", "mlockall", "munlockall", "vhangup", "modify_ldt", "pivot_root", "_sysctl", "prctl", "arch_prctl", "adjtimex", "setrlimit", "chroot", "sync", "acct", "settimeofday", "mount", "umount2", "swapon", "swapoff", "reboot", "sethostname", "setdomainname", "iopl", "ioperm", "create_module", "init_module", "delete_module", "get_kernel_syms", "query_module", "quotactl", "nfsservctl", "getpmsg", "putpmsg", "afs_syscall", "tuxcall", "security", "gettid", "readahead", "setxattr", "lsetxattr", "fsetxattr", "getxattr", "lgetxattr", "fgetxattr", "listxattr", "llistxattr", "flistxattr", "removexattr", "lremovexattr", "fremovexattr", "tkill", "time", "futex", "sched_setaffinity", "sched_getaffinity", "set_thread_area", "io_setup", "io_destroy", "io_getevents", "io_submit", "io_cancel", "get_thread_area", "lookup_dcookie", "epoll_create", "epoll_ctl_old", "epoll_wait_old", "remap_file_pages", "getdents64", "set_tid_address", "restart_syscall", "semtimedop", "fadvise64", "timer_create", "timer_settime", "timer_gettime", "timer_getoverrun", "timer_delete", "clock_settime", "clock_gettime", "clock_getres", "clock_nanosleep", "exit_group", "epoll_wait", "epoll_ctl", "tgkill", "utimes", "vserver", "vserver", "mbind", "set_mempolicy", "get_mempolicy", "mq_open", "mq_unlink", "mq_timedsend", "mq_timedreceive", "mq_notify", "mq_getsetattr", "kexec_load", "waitid"};

W64 lengthofSyscallNames() { return lengthof(syscall_names_64bit); }