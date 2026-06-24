# Thin wrapper around CMake presets.
#
# This file is named GNUmakefile on purpose: GNU make resolves makefiles in the
# order GNUmakefile > makefile > Makefile, so `make` always picks this wrapper
# even when CMake generates its own `Makefile` in the source tree (in-source
# builds). CMake never writes a GNUmakefile, so the two never collide.

CMAKE ?= cmake

# Default generator is Ninja. Override per-invocation with `make GEN=...` or via
# the environment (e.g. `CMAKE_GENERATOR=...`). CMake honors CMAKE_GENERATOR
# because the presets deliberately do not pin a generator.
GEN ?= Ninja
export CMAKE_GENERATOR := $(GEN)

PRESETS := debug release relwithdeb tsan
CMAKE_TARGETS := raspsim x86sim x86sim_defaults x86sim-linux

.PHONY: all $(PRESETS) $(CMAKE_TARGETS) test clean distclean format format-check hooks

all: debug

$(PRESETS):
	$(CMAKE) --preset $@
	$(CMAKE) --build --preset $@
	@ln -sfn build/$@/compile_commands.json compile_commands.json

$(CMAKE_TARGETS):
	$(CMAKE) --preset debug
	$(CMAKE) --build --preset debug --target $@
	@ln -sfn build/debug/compile_commands.json compile_commands.json

test: debug
	ctest --test-dir build/debug --output-on-failure

tsan-test: tsan
	ctest --preset tsan

# Run the configured build's clean target for every preset that exists.
clean:
	@for p in $(PRESETS); do \
	  if [ -d build/$$p ]; then $(CMAKE) --build --preset $$p --target clean; fi; \
	done

# Wipe all build artifacts.
distclean:
	rm -rf build compile_commands.json

# Reformat all project C/C++ sources in place with clang-format.
format:
	tools/clang-format.sh fix

# Fail if any source is not clang-format clean (same check CI runs).
format-check:
	tools/clang-format.sh check

# Point this clone's git hooks at the tracked .githooks directory so the
# clang-format pre-commit hook runs.
hooks:
	git config core.hooksPath .githooks
	@echo "core.hooksPath set to .githooks"
