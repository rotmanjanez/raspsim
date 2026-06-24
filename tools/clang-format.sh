#!/usr/bin/env sh
# Format (or check) the project's C/C++ sources with clang-format.
#
# This is the single source of truth for *which* files are subject to the
# project clang-format style; CI and the pre-commit hook both use the same set.
# The vendored SQLite amalgamation (under build/ and tests/) and the Python
# bindings are intentionally excluded.
#
# Usage:
#   tools/clang-format.sh check   # fail if any file is not formatted (default)
#   tools/clang-format.sh fix     # reformat files in place
#
# Override the binary with CLANG_FORMAT=/path/to/clang-format.
set -eu

mode="${1:-check}"
CLANG_FORMAT="${CLANG_FORMAT:-clang-format}"

cd "$(git rev-parse --show-toplevel)"

if ! command -v "$CLANG_FORMAT" >/dev/null 2>&1; then
  echo "error: '$CLANG_FORMAT' not found (set CLANG_FORMAT to override)" >&2
  exit 1
fi

# Tracked C/C++ sources owned by this project.
files=$(git ls-files \
  'src/*.cpp' 'src/*.cc' 'src/*.c' 'src/*.h' 'src/*.hpp' \
  'include/*.h' 'include/*.hpp' 'include/*.cpp')

if [ -z "$files" ]; then
  echo "no source files found" >&2
  exit 0
fi

case "$mode" in
  fix)
    echo "$files" | xargs "$CLANG_FORMAT" -i
    echo "Formatted $(echo "$files" | wc -l | tr -d ' ') files."
    ;;
  check)
    echo "$files" | xargs "$CLANG_FORMAT" --dry-run -Werror
    ;;
  *)
    echo "usage: $0 [check|fix]" >&2
    exit 2
    ;;
esac
