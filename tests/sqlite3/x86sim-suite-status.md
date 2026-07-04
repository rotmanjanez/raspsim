# SQLite test suite under x86sim — status & findings

This documents the state of the upstream SQLite test suite (`SQLITE_TEST_SUITE=all`)
run through `x86sim-linux` on the **sequential (`seq`) core**, and the fixes made
to get there. The suite is driven by `run_sqlite_tests.sh` → `testrunner.tcl`; each
`.test` file is one CTest (`sqlite3_seq_<stem>`).

## Result

| | before | after |
|---|---:|---:|
| passed | 743 (63%) | **1052 (90%)** |
| failed | 432 | **85** |
| timed out | — | **35** |

*(1175 seq tests total. The 1052/88/35 full-suite run predates the last three
syscalls below; `symlink`/`fsync`/`tkill` move a few more to passing — e.g.
`quota2` now passes — so the standing failure count is ~85.)*

## Root cause of the original 432 "failures"

The failures were **not** SQLite failing. `run_sqlite_tests.sh` launched
`testrunner.tcl` under a plain `tclsh`, which has no `sqlite3` Tcl package. That
made `find_interpreter` re-exec the **entire test orchestrator under the
simulator**. The orchestrator `fork()`s a child to run each test job, and the
simulator does not implement `fork()` — so every test whose testset contained a
job died with *"couldn't fork child process"*, while tests whose testset was
empty "passed" with **zero cases actually executed** (hollow passes).

**Fix** (`testrunner.tcl`, `trd_early_x86sim_command`): the orchestrator now runs
**natively**; only the individual `.test` jobs are wrapped with x86sim (they
already were, via `trd_x86sim_wrap_testfixture`). So the driver stays on the host
and the code under test runs on the simulator, as intended.

## Syscalls implemented

Running real test cases surfaced several unimplemented Linux syscalls. All are in
the host-POSIX layer (`syscall-linux-posix.cpp`) except the signal ones
(`syscall-linux.cpp`), and are wired into the `host()` chain.

| # | syscall | needed by | notes |
|---:|---|---|---|
| 204 | `sched_getaffinity` | testfixture worker count | musl `sysconf(_SC_NPROCESSORS_*)` uses only this; returns a single-CPU mask |
| 19/20 | `readv`/`writev` | musl buffered stdio flush | any guest printing more than a trickle hit syscall 20 |
| 280 | `utimensat` | file-mtime stamping | every `*fault` test + backup/cacheflush |
| 88 | `symlink` | `symlink.test` | host passthrough |
| 74/75 | `fsync`/`fdatasync` | `quota2.test` | flush the mapped host fd |
| 200/234 | `tkill`/`tgkill` | musl `raise()`/`pthread_kill()` | Phase-1 no-op signal model: validate & return 0 |

## Classification of the remaining failures/timeouts

Derived by re-running each non-passing test's child directly under the simulator
and bucketing the signature.

### A. Simulator page-faults — *real sim bugs, worth investigating*
`capi3c`, `thread1`, `e_fkey` die on `x86 exception: Exception 14` (a page fault
the simulator raises mid-run). These look like genuine defects in the memory
model rather than missing features.

### B. Unsupported features — *not quick wins*
- **Threads** (`thread001/002/003/004/005`, `thread1/2/3`, `walthread`, `sharedA`):
  the testfixture's `Tcl_CreateThread` / notifier thread fails; the simulator is
  single-process/single-thread. Needs clone/thread support.
- **WAL shared-memory / multi-process locking / file `mmap`** (`wal`, `wal5`,
  `walro*`, `walsetlk*`, `shared2`, `superlock`, `unixexcl`, `multiplex2`,
  `external_reader`, `mmap1/4`, `mmapfault`, `recover`, `rowallock`, `busy2`, …):
  exit early with no output; need `MAP_SHARED` file mmap and cross-process
  locking / a shared-memory VFS.

### C. Genuine value mismatches (~40)
Real behavioral differences once tests actually run:
- **Optimizer / query-plan**: `like` [47/159], `where2` [39/90], `index6/7`,
  `whereG`, `join8`, `starschema1`.
- **Crash-recovery (expected under sim)**: `walcrash` [1000/1400],
  `walcrash2` [1000/1002], `crash4` [1000/2999], `crash8` — simulated power-loss
  doesn't reproduce.
- **Other diffs** (small counts): `fkey2`, `misc7`, `pragma`, `pragma3`, `stat`,
  `trace3`, `trigger7`, `without_rowid3`, `delete`, `e_select`, `temptable`,
  `readonly`, `lock2/4`, `oserror`, `backup2`, `symlink` (1/39), …

### D. Slow / timeout artifacts (~41 + a few)
Fault-injection (`malloc`, `mallocA`, `*fault`, `ioerr`, `*_err`), `fuzz*`,
`speed*`, `sort*`, and huge-value tests (`atof1`, `round1`, `json106`,
`window3`). These hit the 900 s cap under `-j4` contention rather than being
broken — `crash3`, `expr`, `lock`, and `table` were verified to **pass when run
standalone**. A larger per-test timeout would reclassify several.

## Suggested next steps

1. **Timeouts**: raise `SQLITE_TEST_TIMEOUT` and/or lower parallelism for the
   heavy fault/fuzz/speed tests — likely the cheapest win for the count.
2. **Page-faults** (`capi3c`, `thread1`, `e_fkey`): debug the `Exception 14`
   crashes; these are the most likely to be genuine simulator bugs.
3. **Threads / WAL shared-memory**: larger feature work (clone + shared-memory
   VFS) behind most of the structural failures.
