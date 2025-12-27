#ifndef _RASPSIM_HWSETUP_H_
#define _RASPSIM_HWSETUP_H_

#include <typedefs.h>

class AddressSpace;
class Context;
class PTLsimMachine;

extern const char* const syscall_names_64bit[];

W64 lengthofSyscallNames();

// Class for setting up user space simulation mode and common settup code for the simulator
// Wrapper to static or global state for better binding support
// TODO: (JR) add lock so only a single instance of raspsim can be used at a time
class Raspsim {
public:

  Raspsim();
  ~Raspsim();

  void run();

protected:
  PTLsimMachine* getMachine();
  const char* getCoreName();
public:

  int getRegisterIndex(const char* regname);
  void setRegisterValue(int reg, W64 value);
  W64 getRegisterValue(int reg);
  
  byte* getMappedPage(Waddr addr);
  int getPageProtection(void* addr);

  void disableSSE();
  void disableX87();
  void enablePerfectCache();
  void enableStaticBranchPrediction();
  void setLogfile(const char* filename);
  void setTimeout(W64 ninstr);

  void map(Waddr start, W64 size, int prot);
  void unmap(Waddr start, W64 size);
  void* page_virt_to_mapped(Waddr start);

  static W64 cycles();
  static W64 instructions();

  static Waddr getPageSize();

  static void stutdown();
  static AddressSpace& getAddrspace();
  static Context& getContext();
  
  static char* formatException(byte exception, W32 errorcode, Waddr virtaddr);
  static char* formatContext(const Context& ctx);
  

  // Implement as desired depending on the bindings
  static void propagate_x86_exception(byte exception, W32 errorcode, Waddr virtaddr);

#ifdef __x86_64__
  static void handle_syscall_64bit();
#endif // __x86_64__
  static void handle_syscall_32bit(int semantics);
};

#endif
