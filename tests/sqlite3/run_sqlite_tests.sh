#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: tests/sqlite3/run_sqlite_tests.sh [RUNNER_OPTIONS] [--] [TESTRUNNER_ARGS...]

Runner options:
  --core=seq|ooo          x86sim core to use for testfixture (default: seq)
  --seq                   shorthand for --core=seq
  --ooo, --ooocore        shorthand for --core=ooo
  --x86sim PATH           x86sim-linux binary
  --x86sim-arg ARG        extra argument passed to x86sim-linux; may repeat
  --testfixture PATH      sqlite3-test/testfixture binary
  --testdir PATH          SQLite test script directory
  --workdir PATH          scratch directory for testrunner state
  --help                  show this help

Examples:
  tests/sqlite3/run_sqlite_tests.sh --core=seq veryquick
  tests/sqlite3/run_sqlite_tests.sh --core=ooo all select1.test
  tests/sqlite3/run_sqlite_tests.sh --seq -- --dryrun --jobs 1 all select1.test
USAGE
}

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "$script_dir/../.." && pwd)"

x86sim="${X86SIM:-$repo_root/build/release/x86sim-linux}"
core="${X86SIM_CORE:-seq}"
testfixture="${SQLITE_TESTFIXTURE:-$script_dir/sqlite3-test}"
testdir="${SQLITE_TESTDIR:-$script_dir}"
workdir="${SQLITE_TEST_WORKDIR:-$repo_root/build/sqlite-test-run}"
tclsh="${TCLSH:-tclsh}"
x86sim_args=()
testrunner_args=()

while (($#)); do
  case "$1" in
    --help|-h)
      usage
      exit 0
      ;;
    --core)
      if (($# < 2)); then
        echo "missing value for --core" >&2
        exit 2
      fi
      core="$2"
      shift 2
      ;;
    --core=*)
      core="${1#*=}"
      shift
      ;;
    --x86sim-core)
      if (($# < 2)); then
        echo "missing value for --x86sim-core" >&2
        exit 2
      fi
      core="$2"
      shift 2
      ;;
    --x86sim-core=*)
      core="${1#*=}"
      shift
      ;;
    --seq)
      core="seq"
      shift
      ;;
    --ooo|--ooocore)
      core="ooo"
      shift
      ;;
    --x86sim)
      if (($# < 2)); then
        echo "missing value for --x86sim" >&2
        exit 2
      fi
      x86sim="$2"
      shift 2
      ;;
    --x86sim=*)
      x86sim="${1#*=}"
      shift
      ;;
    --x86sim-arg)
      if (($# < 2)); then
        echo "missing value for --x86sim-arg" >&2
        exit 2
      fi
      x86sim_args+=("$2")
      shift 2
      ;;
    --x86sim-arg=*)
      x86sim_args+=("${1#*=}")
      shift
      ;;
    --testfixture)
      if (($# < 2)); then
        echo "missing value for --testfixture" >&2
        exit 2
      fi
      testfixture="$2"
      shift 2
      ;;
    --testfixture=*)
      testfixture="${1#*=}"
      shift
      ;;
    --testdir)
      if (($# < 2)); then
        echo "missing value for --testdir" >&2
        exit 2
      fi
      testdir="$2"
      shift 2
      ;;
    --testdir=*)
      testdir="${1#*=}"
      shift
      ;;
    --workdir)
      if (($# < 2)); then
        echo "missing value for --workdir" >&2
        exit 2
      fi
      workdir="$2"
      shift 2
      ;;
    --workdir=*)
      workdir="${1#*=}"
      shift
      ;;
    --)
      shift
      testrunner_args+=("$@")
      break
      ;;
    *)
      testrunner_args+=("$1")
      shift
      ;;
  esac
done

case "$core" in
  seq|ooo) ;;
  *)
    echo "unsupported x86sim core '$core' (expected seq or ooo)" >&2
    exit 2
    ;;
esac

if [[ ! -x "$x86sim" ]]; then
  echo "x86sim binary is not executable: $x86sim" >&2
  exit 1
fi

if [[ ! -x "$testfixture" ]]; then
  echo "sqlite testfixture is not executable: $testfixture" >&2
  exit 1
fi

if [[ ! -d "$testdir" ]]; then
  echo "SQLite test directory does not exist: $testdir" >&2
  exit 1
fi

mkdir -p "$workdir"
rm -f "$workdir/testfixture"
ln -s "$testfixture" "$workdir/testfixture"

wrapped_args=()
for arg in "${x86sim_args[@]}"; do
  wrapped_args+=(--x86sim-arg "$arg")
done

cd "$workdir"
export SQLITE_TESTDIR="$testdir"

exec "$tclsh" "$script_dir/testrunner.tcl" \
  --x86sim "$x86sim" \
  --x86sim-core "$core" \
  "${wrapped_args[@]}" \
  "${testrunner_args[@]}"
