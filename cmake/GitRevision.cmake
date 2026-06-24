# Determines the current git revision and exposes it as GIT_REVISION.
# Falls back to "unknown" when git is unavailable or the source tree is not a repo.

find_package(Git QUIET)
if(Git_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} -C ${CMAKE_CURRENT_SOURCE_DIR} describe --tags --always --dirty
    OUTPUT_VARIABLE GIT_REVISION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
endif()

if(NOT GIT_REVISION)
  set(GIT_REVISION unknown)
endif()
