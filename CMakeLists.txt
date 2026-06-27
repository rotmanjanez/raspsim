cmake_minimum_required(VERSION 3.25)

project(x86sim
  DESCRIPTION "Cycle-accurate x86-64 virtual CPU core and command-line simulators"
  LANGUAGES CXX
)

option(X86SIM_BUILD_DEFAULTS "Build the Linux syscall/cpuid support library" ON)
option(X86SIM_BUILD_RASPSIM  "Build the raspsim command-line simulator"      ON)
option(X86SIM_BUILD_LINUX    "Build the x86sim-linux command-line simulator" ON)
option(X86SIM_BUILD_TESTS    "Build the CTest suite (SQLite through x86sim-linux)" ON)
option(X86SIM_BUILD_PYTHON   "Build the pybind11 Python extension module"   OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(GNUInstallDirs)
include(GitRevision)

# Build-time provenance, stamped into the core library and the simulators.
cmake_host_system_information(RESULT buildhost QUERY FQDN)
set(build_info BUILDHOST=${buildhost} GIT_REVISION=${GIT_REVISION})

set(install_targets)

# --- x86sim: portable virtual CPU core (always built) ---------------------
file(GLOB vcore_sources CONFIGURE_DEPENDS src/x86sim/*.cpp)
add_library(x86sim)
add_library(x86sim::x86sim ALIAS x86sim)
target_sources(x86sim
  PRIVATE
    ${vcore_sources}
  PUBLIC
    FILE_SET HEADERS
    BASE_DIRS include
    FILES
      include/x86sim/logging.hpp
      include/x86sim/registerfile.hpp
      include/x86sim/x86sim.hpp
)
target_include_directories(x86sim PRIVATE src/x86sim)
target_compile_definitions(x86sim PRIVATE ${build_info})
target_compile_features(x86sim PUBLIC cxx_std_26)
list(APPEND install_targets x86sim)

# --- x86sim_defaults: Linux syscall + cpuid support library ---------------
if(X86SIM_BUILD_DEFAULTS)
  add_library(x86sim_defaults)
  add_library(x86sim::defaults ALIAS x86sim_defaults)
  set_target_properties(x86sim_defaults PROPERTIES EXPORT_NAME defaults)
  target_sources(x86sim_defaults
    PRIVATE
      src/x86sim-support/cpuid.cpp
      src/x86sim-support/syscall-linux.cpp
    PUBLIC
      FILE_SET HEADERS
      BASE_DIRS include
      FILES
        include/x86sim-support/cpuid.hpp
        include/x86sim-support/syscall-linux.hpp
  )
  target_link_libraries(x86sim_defaults PUBLIC x86sim::x86sim)
  target_include_directories(x86sim_defaults PRIVATE src/x86sim)
  list(APPEND install_targets x86sim_defaults)
endif()

# --- raspsim: command-line simulator --------------------------------------
if(X86SIM_BUILD_RASPSIM)
  add_executable(raspsim
    src/raspsim/config.cpp
    src/raspsim/raspsim.cpp
  )
  target_link_libraries(raspsim PRIVATE x86sim::defaults)
  target_include_directories(raspsim PRIVATE src/x86sim)
  list(APPEND install_targets raspsim)
endif()

# --- x86sim-linux: command-line simulator ---------------------------------
if(X86SIM_BUILD_LINUX)
  add_executable(x86sim-linux
    src/x86sim-linux/x86sim-linux.cpp
  )
  target_link_libraries(x86sim-linux PRIVATE x86sim::defaults)
  list(APPEND install_targets x86sim-linux)
endif()

# --- install / export -----------------------------------------------------
install(TARGETS ${install_targets}
  EXPORT x86simTargets
  FILE_SET HEADERS
)
install(EXPORT x86simTargets
  NAMESPACE x86sim::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/x86sim"
)

# --- python bindings ------------------------------------------------------
# The extension only links the portable core (plus cpuid), never the Linux
# support library, so it stays buildable off-Linux. Driven by pyproject.toml,
# which points scikit-build-core at this file and disables the Linux tools.
if(X86SIM_BUILD_PYTHON)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(core MODULE
    pyraspsim/src/binding.cpp
    src/x86sim-support/cpuid.cpp
  )
  target_link_libraries(core PRIVATE x86sim::x86sim)
  target_include_directories(core PRIVATE include)
  install(TARGETS core LIBRARY DESTINATION x86sim)
endif()

# --- tests ----------------------------------------------------------------
if(X86SIM_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests/sqlite3)
endif()
