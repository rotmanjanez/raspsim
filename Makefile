# -*- makefile -*-
#
# PTLsim: Cycle Accurate x86-64 Simulator
# Makefile
#
# Copyright 2000-2008 Matt T. Yourst <yourst@yourst.com>
#

#
# If you are running on a 64-bit distro but want to build
# a 32-bit PTLsim binary, and your distro doesn't provide
# the "linux32" or "32bit" uname-changing commands, you
# will need to manually override the checks below:
#
ifndef MACHTYPE
	MACHTYPE = "$(shell uname -m)"
endif

ifneq (,$(findstring x86_64,"$(MACHTYPE)"))
	__x86_64__=1
endif

# For GCC versions > 4.2 install version 4.2 and uncomment the following line:
# CC = g++-4.2
CC = g++

GCCVER_SPECIFIC =

#SVNREV=$(shell svn info | grep "Last Changed Rev" | cut -d " " -f4)
#SVNDATE=$(shell svn info | grep "Last Changed Date" | cut -d " " -f4)

#ifeq (,$(SVNREV))
# Subversion is either not installed or the current directory isn't a PTLsim repository:
SVNREV=0
SVNDATE=unknown
#endif

INCFLAGS = -Isrc -DBUILDHOST="`hostname -f`" -DSVNREV="$(SVNREV)" -DSVNDATE="$(SVNDATE)"

# Cross-platform C++23 flags - no architecture-specific optimizations
CFLAGS = -std=c++23 -O2 -fomit-frame-pointer -pipe -falign-functions=16 -funroll-loops
CFLAGS += -fno-trapping-math -fno-stack-protector -fno-strict-aliasing -Wreturn-type $(GCCVER_SPECIFIC)
# Note: Removed -fno-exceptions -fno-rtti for C++23 standard library compatibility



BASEOBJS = superstl.o config.o syscalls.o
COMMONOBJS = ptlsim.o mm.o ptlhwdef.o decode-core.o decode-fast.o decode-complex.o decode-x87.o decode-sse.o uopimpl.o seqcore.o

OOOOBJS = branchpred.o dcache.o ooocore.o ooopipe.o oooexec.o
RASPSIMOBJS = raspsim-hwsetup.o addrspace.o

COMMONINCLUDES = logic.h ptlhwdef.h decode.h seqexec.h dcache.h dcache-amd-k8.h config.h ptlsim.h superstl.h globals.h ptlsim-api.h mm.h syscalls.h stats.h typedefs.h registers.def
OOOINCLUDES = branchpred.h ooocore.h ooocore-amd-k8.h

COMMONCPPFILES = ptlsim.cpp raspsim.cpp mm.cpp superstl.cpp ptlhwdef.cpp decode-core.cpp decode-fast.cpp decode-complex.cpp decode-x87.cpp decode-sse.cpp uopimpl.cpp dcache.cpp config.cpp syscalls.cpp

OOOCPPFILES = ooocore.cpp ooopipe.cpp oooexec.cpp seqcore.cpp branchpred.cpp

OBJS = $(addprefix src/, $(BASEOBJS) $(COMMONOBJS) $(OOOOBJS) $(RASPSIMOBJS))
INCLUDEFILES = $(addprefix src/, $(COMMONINCLUDES) $(OOOINCLUDES))
CPPFILES = $(addprefix src/, $(COMMONCPPFILES) $(OOOCPPFILES))

CFLAGS += -D__PTLSIM_OOO_ONLY__

TOPLEVEL = raspsim

all: $(TOPLEVEL)
	@echo "Compiled successfully..."

# Cross-platform build targets (no architecture restrictions)
raspsim: src/raspsim.o $(OBJS) Makefile
	$(CXX) $< $(OBJS) -o $@

PYRASPSIM = pyraspsim/raspsim/core$(shell python3-config --extension-suffix)

$(PYRASPSIM): CFLAGS += -fPIC
$(PYRASPSIM): pyraspsim/raspsim/pyraspsim.cpp $(OBJS) Makefile
	@python3 -c "import pybind11" || (echo "pybind11 is not installed. Please install it using 'pip3 install pybind11'"; exit 1)
	$(CXX) $(INCFLAGS) -O3 -Wall -shared -std=c++23 -fPIC $(shell python3 -m pybind11 --includes) $< -o $@ $(OBJS)

pyraspsim/raspsim/core.pyi: $(PYRASPSIM) Makefile
	@which pybind11-stubgen > /dev/null 2>&1 || (echo "pybind11-stubgen is not installed. Please install it using 'pip install pybind11-stubgen'"; exit 1)
	cd pyraspsim/raspsim && env PYTHONPATH=. pybind11-stubgen -o . core

pyraspsim: $(PYRASPSIM) pyraspsim/raspsim/core.pyi

src/%.o: src/%.cpp
	$(CC) $(CFLAGS) $(INCFLAGS) -o $@ -c $<

src/%.o: src/%.S
	$(CC) $(CFLAGS) $(INCFLAGS) -o $@ -c $<

src/%.o: src/%.c
	$(CC) $(CFLAGS) $(INCFLAGS) -o $@ -c $<

clean:
	rm -fvrd raspsim src/*.o core core.[0-9]* .depend *.gch *.so *.log.* pyraspsim/raspsim/*.so pyraspsim/raspsim/*.pyi pyraspsim/*.egg-info **/__pycache__  **/**/__pycache__

INCLUDEFILES += $(PT2XINCLUDES)
CPPFILES += $(PT2XCPPFILES)

#
# Miscellaneous:
#

DISTFILES = $(CPPFILES) $(INCLUDEFILES) Makefile COPYING README

dist: $(DISTFILES)
	tar zcvf ptlsim-`date "+%Y%m%d%H%M%S"`.tar.gz $(DISTFILES)

backup: dist

distfiles: $(DISTFILES)
	@echo $(DISTFILES)

.depend:
	$(CC) $(CFLAGS) $(INCFLAGS) -MM $(CPPFILES) $(ASMFILES) > .depend

-include .depend

.PHONY: all clean dist distfiles backup pyraspsim(venv) 