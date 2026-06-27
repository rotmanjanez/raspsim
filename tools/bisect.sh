#!/usr/bin/env bash
set -euo pipefail
shopt -s extglob

sim="./build/debug/x86sim-linux"
backend="auto"
core="seq"
max_instr=""
regs=(rip)

usage() {
  cat <<'EOF'
usage: tools/diffcpu.sh [OPTIONS] -- PROGRAM [ARG...]

Compare x86sim-linux against a real x86-64 CPU under a debugger and bisect the
first divergent instruction count.

Options:
  --sim PATH           x86sim-linux path (default: ./build/debug/x86sim-linux)
  --backend NAME      auto, gdb, or lldb (default: auto)
  --core NAME         seq or ooo (default: seq)
  --max-instr N       upper instruction count to probe (required)
  --regs "LIST"       space-separated registers to compare (default: rip)
  -h, --help          show this help

The native side requires Linux on real x86-64 hardware. V1 is intended for
static/non-PIE x86-64 Linux ELFs and deterministic single-thread repros.
EOF
}

die() {
  echo "diffcpu: $*" >&2
  exit 2
}

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"
}

canon_hex() {
  local value="$1"
  value="${value%%,*}"
  value="$(printf '%s' "$value" | tr '[:upper:]' '[:lower:]')"
  [[ "$value" == 0x* ]] || die "expected hex register value, got: $1"
  value="${value#0x}"
  value="${value##+(0)}"
  [[ -n "$value" ]] || value="0"
  printf '0x%s' "$value"
}

parse_reg_line() {
  local prefix="$1"
  local line="$2"
  local -A found=()
  local word key value

  for word in $line; do
    [[ "$word" == "$prefix" ]] && continue
    [[ "$word" == *=* ]] || continue
    key="${word%%=*}"
    value="${word#*=}"
    found["$key"]="$(canon_hex "$value")"
  done

  for key in "${regs[@]}"; do
    [[ -n "${found[$key]:-}" ]] || die "missing $key in $prefix output: $line"
    printf '%s=%s\n' "$key" "${found[$key]}"
  done
}

parse_debugger_regs() {
  local output="$1"
  local -A found=()
  local line name key value

  while IFS= read -r line; do
    name=""
    if [[ "$line" =~ ^[[:space:]]*([A-Za-z0-9_]+)[[:space:]]*= ]]; then
      name="${BASH_REMATCH[1]}"
    elif [[ "$line" =~ ^[[:space:]]*([A-Za-z0-9_]+)[[:space:]]+ ]]; then
      name="${BASH_REMATCH[1]}"
    fi
    [[ -n "$name" ]] || continue
    if [[ "$line" =~ (0x[0-9A-Fa-f]+) ]]; then
      value="${BASH_REMATCH[1]}"
    else
      continue
    fi

    case "$name" in
      eflags|rflags) key="flags" ;;
      *) key="$name" ;;
    esac
    found["$key"]="$(canon_hex "$value")"
  done <<<"$output"

  for key in "${regs[@]}"; do
    [[ -n "${found[$key]:-}" ]] || die "missing $key in debugger output"
    printf '%s=%s\n' "$key" "${found[$key]}"
  done
}

show_state() {
  local state="$1"
  while IFS= read -r line; do
    echo "  $line"
  done <<<"$state"
}

quote_words() {
  local word
  for word in "$@"; do
    printf '%q ' "$word"
  done
}

show_commands() {
  echo "sim command template: $(quote_words "$sim" "--core=$core" --stopinsns N --dump-regs -- "${program_args[@]}")"
  case "$backend" in
    gdb)
      echo "native command template: gdb -q -nx -batch -ex 'set disable-randomization on' -ex 'break *$entry' -ex run -ex 'si N' -ex 'info registers $(native_reg_names)' --args $(quote_words "${program_args[@]}")"
      ;;
    lldb)
      echo "native command template: lldb --batch -o 'settings set target.disable-aslr true' -o 'breakpoint set --address $entry' -o run -o 'thread step-inst --count N' -o 'register read $(native_reg_names)' -- $(quote_words "${program_args[@]}")"
      ;;
  esac
}

diff_state() {
  local sim_state="$1"
  local native_state="$2"
  local key sim_value native_value

  for key in "${regs[@]}"; do
    sim_value="$(awk -F= -v key="$key" '$1 == key { print $2 }' <<<"$sim_state")"
    native_value="$(awk -F= -v key="$key" '$1 == key { print $2 }' <<<"$native_state")"
    if [[ "$sim_value" != "$native_value" ]]; then
      printf '  %s: sim=%s native=%s\n' "$key" "$sim_value" "$native_value"
    fi
  done
}

states_equal() {
  [[ "$1" == "$2" ]]
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --sim)
        [[ $# -ge 2 ]] || die "missing value for --sim"
        sim="$2"
        shift 2
        ;;
      --backend)
        [[ $# -ge 2 ]] || die "missing value for --backend"
        backend="$2"
        shift 2
        ;;
      --core)
        [[ $# -ge 2 ]] || die "missing value for --core"
        core="$2"
        shift 2
        ;;
      --max-instr)
        [[ $# -ge 2 ]] || die "missing value for --max-instr"
        max_instr="$2"
        shift 2
        ;;
      --regs)
        [[ $# -ge 2 ]] || die "missing value for --regs"
        read -r -a regs <<<"$2"
        shift 2
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      --)
        shift
        program_args=("$@")
        return
        ;;
      *)
        die "unknown option or missing -- before program: $1"
        ;;
    esac
  done
  die "missing PROGRAM"
}

choose_backend() {
  case "$backend" in
    gdb|lldb) ;;
    auto)
      if command -v gdb >/dev/null 2>&1; then
        backend="gdb"
      elif command -v lldb >/dev/null 2>&1; then
        backend="lldb"
      else
        die "neither gdb nor lldb found"
      fi
      ;;
    *) die "unknown backend: $backend" ;;
  esac
}

native_reg_names() {
  local names=()
  local reg
  for reg in "${regs[@]}"; do
    if [[ "$reg" == "flags" ]]; then
      if [[ "$backend" == "gdb" ]]; then
        names+=(eflags)
      else
        names+=(rflags)
      fi
    else
      names+=("$reg")
    fi
  done
  printf '%s ' "${names[@]}"
}

run_sim() {
  local count="$1"
  local output status line
  set +e
  output=$("$sim" "--core=$core" --stopinsns "$count" --dump-regs -- "${program_args[@]}" 2>&1 >/dev/null)
  status=$?
  set -e
  [[ $status -eq 0 ]] || die "simulator failed at count $count: $output"
  line="$(grep '^x86sim-regs ' <<<"$output" | tail -n 1 || true)"
  [[ -n "$line" ]] || die "simulator did not print x86sim-regs at count $count"
  parse_reg_line "x86sim-regs" "$line"
}

run_native_gdb() {
  local count="$1"
  local reg_names output status
  reg_names="$(native_reg_names)"
  local cmd=(gdb -q -nx -batch
    -ex 'set pagination off'
    -ex 'set confirm off'
    -ex 'set disable-randomization on'
    -ex "break *$entry"
    -ex run)
  if (( count > 0 )); then
    cmd+=(-ex "si $count")
  fi
  cmd+=(-ex "info registers $reg_names" --args "${program_args[@]}")

  set +e
  output=$("${cmd[@]}" 2>&1)
  status=$?
  set -e
  [[ $status -eq 0 ]] || die "gdb failed at count $count: $output"
  parse_debugger_regs "$output"
}

run_native_lldb() {
  local count="$1"
  local reg_names output status
  reg_names="$(native_reg_names)"
  local cmd=(lldb --batch
    -o 'settings set target.disable-aslr true'
    -o "breakpoint set --address $entry"
    -o run)
  if (( count > 0 )); then
    cmd+=(-o "thread step-inst --count $count")
  fi
  cmd+=(-o "register read $reg_names" -- "${program_args[@]}")

  set +e
  output=$("${cmd[@]}" 2>&1)
  status=$?
  set -e
  [[ $status -eq 0 ]] || die "lldb failed at count $count: $output"
  parse_debugger_regs "$output"
}

run_native() {
  case "$backend" in
    gdb) run_native_gdb "$1" ;;
    lldb) run_native_lldb "$1" ;;
  esac
}

probe() {
  local count="$1"
  sim_state="$(run_sim "$count")"
  native_state="$(run_native "$count")"
}

program_args=()
parse_args "$@"

[[ ${#program_args[@]} -gt 0 ]] || die "missing PROGRAM"
[[ -n "$max_instr" ]] || die "--max-instr is required"
[[ "$max_instr" =~ ^[0-9]+$ ]] || die "--max-instr must be a decimal integer"
[[ ${#regs[@]} -gt 0 ]] || die "--regs must name at least one register"
[[ -x "$sim" ]] || die "simulator is not executable: $sim"

case "$(uname -s):$(uname -m)" in
  Linux:x86_64|Linux:amd64) ;;
  *) die "native comparison requires Linux on real x86-64 hardware" ;;
esac

need_cmd readelf
choose_backend
need_cmd "$backend"

program="${program_args[0]}"
[[ -f "$program" ]] || die "program not found: $program"

elf_type="$(readelf -h "$program" | awk -F: '/Type:/ { gsub(/^[ \t]+/, "", $2); print $2; exit }')"
[[ "$elf_type" == EXEC* ]] || die "expected static/non-PIE ET_EXEC ELF, got: $elf_type"
entry="$(readelf -h "$program" | awk '/Entry point address:/ { print $4; exit }')"
[[ "$entry" =~ ^0x[0-9A-Fa-f]+$ ]] || die "could not determine ELF entry point"

echo "diffcpu: backend=$backend core=$core max_instr=$max_instr entry=$entry regs=${regs[*]}" >&2

probe 0
if ! states_equal "$sim_state" "$native_state"; then
  echo "first divergent instruction count: 0"
  echo "differing registers:"
  diff_state "$sim_state" "$native_state"
  echo "sim state:"
  show_state "$sim_state"
  echo "native state:"
  show_state "$native_state"
  show_commands
  exit 1
fi

probe "$max_instr"
if states_equal "$sim_state" "$native_state"; then
  echo "no divergence found through instruction count $max_instr"
  echo "state:"
  show_state "$sim_state"
  exit 0
fi

lo=0
hi="$max_instr"
last_equal_state=""
first_bad_sim="$sim_state"
first_bad_native="$native_state"

while (( hi - lo > 1 )); do
  mid=$(((lo + hi) / 2))
  probe "$mid"
  if states_equal "$sim_state" "$native_state"; then
    lo="$mid"
    last_equal_state="$sim_state"
  else
    hi="$mid"
    first_bad_sim="$sim_state"
    first_bad_native="$native_state"
  fi
done

echo "previous equal instruction count: $lo"
echo "first divergent instruction count: $hi"
echo "differing registers:"
diff_state "$first_bad_sim" "$first_bad_native"
if [[ -n "$last_equal_state" ]]; then
  echo "previous equal state:"
  show_state "$last_equal_state"
fi
echo "sim state at $hi:"
show_state "$first_bad_sim"
echo "native state at $hi:"
show_state "$first_bad_native"
show_commands
exit 1
