# Generate one CTest per SQLite .test file, per x86sim core.
#
# Each test runs a single upstream SQLite .test script through the x86sim-linux
# guest executor via run_sqlite_tests.sh. Tests are labelled by core ("seq" /
# "ooo") so the two groups can be selected with `ctest -L seq` / `ctest -L ooo`.
#
# The set of .test files is configurable via SQLITE_TEST_SUITE so a local run
# can do only the fastest subset while CI/main runs everything.

set(SQLITE_TEST_SUITE "veryquick" CACHE STRING
    "SQLite test scope run through x86sim-linux: extraquick | veryquick | all")
set_property(CACHE SQLITE_TEST_SUITE PROPERTY STRINGS extraquick veryquick all)

set(SQLITE_TEST_CORES "seq;ooo" CACHE STRING
    "x86sim cores to generate test groups for (semicolon-separated, e.g. seq;ooo)")

set(SQLITE_TEST_TIMEOUT 900 CACHE STRING
    "Per-test timeout in seconds for SQLite-through-x86sim CTests")

find_program(TCLSH_EXECUTABLE NAMES tclsh REQUIRED)

set(_runner  "${CMAKE_CURRENT_SOURCE_DIR}/run_sqlite_tests.sh")
set(_fixture "${CMAKE_CURRENT_SOURCE_DIR}/sqlite3-test")

file(GLOB _tests CONFIGURE_DEPENDS
     RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/*.test")

# Permutation / meta drivers — not standalone leaf tests. Excluded for every
# suite (mirrors the $alltests exclusions in permutations.test).
set(_meta_excludes
  all quick veryquick extraquick full permutations memleak soak
  fts3 mallocAll rtree session rbu releasetest)
foreach(_m IN LISTS _meta_excludes)
  list(FILTER _tests EXCLUDE REGEX "^${_m}\\.test$")
endforeach()

if(SQLITE_TEST_SUITE STREQUAL "all")
  set(_perm full)
else()
  set(_perm veryquick)

  # veryquick glob excludes (substring match), from permutations.test:
  foreach(_pat malloc ioerr fault bigfile _err fts5corrupt fts5big fts5aj rbucrash)
    list(FILTER _tests EXCLUDE REGEX "${_pat}")
  endforeach()

  # $allquicktests literal exclusions: crash / thread / fuzz / speed / slow.
  set(_slow_excludes
    backup_ioerr corrupt corruptC crash crash2 crash3 crash4 crash5 crash6 crash7
    delete3 e_fts3 fts3rnd fkey_malloc fuzz fuzz3 fuzz_malloc in2 loadext misc7
    mutex2 onefile pagerfault2 savepoint4 savepoint6 select9
    speed1 speed1p speed2 speed3 speed4 speed4p sqllimits1 tkt2686
    thread001 thread002 thread003 thread004 thread005 trans2 vacuum3
    incrvacuum_ioerr autovacuum_crash btree8 shared_err vtab_err walslow
    walcrash walcrash3 walthread rtree3 indexfault securedel2 sort3 sort4
    fts4growth fts4growth2 bigsort walprotocol mmap4 fuzzer2 walcrash2 e_fkey backup
    fts4merge fts4merge2 fts4merge4 fts4check fts4merge5 fts3cov fts3snippet
    fts3corrupt2 fts3an fts3defer fts4langid fts3sort fts5unicode recovercorrupt
    rtree4 sessionbig writecrash view3 fts5dlidx fts5ac fts4merge3 fts5prefix sessionB)
  if(SQLITE_TEST_SUITE STREQUAL "extraquick")
    # extraquick drops a few more slow files on top of veryquick.
    list(APPEND _slow_excludes
      wal3 sort2 mmap1 percentile where8m walcksum savepoint3
      fuzzer1 fuzzer3 fts3expr3)
  endif()
  foreach(_s IN LISTS _slow_excludes)
    list(FILTER _tests EXCLUDE REGEX "^${_s}\\.test$")
  endforeach()
endif()

list(SORT _tests)
list(LENGTH _tests _n)
message(STATUS
  "SQLite/x86sim CTests: ${_n} .test files x [${SQLITE_TEST_CORES}] "
  "(suite=${SQLITE_TEST_SUITE}, perm=${_perm})")

foreach(_core IN LISTS SQLITE_TEST_CORES)
  foreach(_tf IN LISTS _tests)
    cmake_path(GET _tf STEM _stem)
    set(_name "sqlite3_${_core}_${_stem}")
    # Per-test workdir so tests can run under `ctest -j` without colliding on
    # the symlinked testfixture inside the shared default workdir.
    set(_wd "${CMAKE_CURRENT_BINARY_DIR}/run/${_core}/${_stem}")
    add_test(NAME "${_name}"
      COMMAND "${_runner}"
        --core=${_core}
        --x86sim      "$<TARGET_FILE:x86sim-linux>"
        --testfixture "${_fixture}"
        --testdir     "${CMAKE_CURRENT_SOURCE_DIR}"
        --workdir     "${_wd}"
        -- --jobs 1 --stop-on-error ${_perm} ${_tf})
    set_tests_properties("${_name}" PROPERTIES
      LABELS "sqlite3;${_core}"
      TIMEOUT ${SQLITE_TEST_TIMEOUT}
      ENVIRONMENT "TCLSH=${TCLSH_EXECUTABLE}")
  endforeach()
endforeach()
