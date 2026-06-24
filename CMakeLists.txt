cmake_minimum_required(VERSION 3.25)

project(x86sim
  DESCRIPTION "Cycle-accurate x86-64 virtual CPU core and command-line simulators"
  LANGUAGES CXX
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

cmake_host_system_information(RESULT buildhost QUERY FQDN)
include(GNUInstallDirs)
include(GitRevision)

file(GLOB vcore_sources CONFIGURE_DEPENDS src/x86sim/*.cpp)

set(build_info_definitions BUILDHOST=${buildhost} GIT_REVISION=${GIT_REVISION})

add_library(x86sim)
add_library(x86sim::x86sim ALIAS x86sim)
target_sources(x86sim
  PRIVATE
    ${vcore_sources}
  PUBLIC
    FILE_SET HEADERS
    BASE_DIRS include
    FILES
      include/x86sim/logging.h
      include/x86sim/registerfile.hpp
      include/x86sim/x86sim.hpp
)
target_include_directories(x86sim
  PRIVATE
    src/x86sim
)
target_compile_definitions(x86sim PRIVATE ${build_info_definitions})

add_library(x86sim_defaults)
add_library(x86sim::defaults ALIAS x86sim_defaults)
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
set_target_properties(x86sim_defaults PROPERTIES EXPORT_NAME defaults)

add_executable(raspsim
  src/raspsim/config.cpp
  src/raspsim/raspsim.cpp
)
target_link_libraries(raspsim PRIVATE x86sim::defaults)
target_include_directories(raspsim PRIVATE src/x86sim)
target_compile_definitions(raspsim PRIVATE ${build_info_definitions})

add_executable(x86sim-linux
  src/x86sim-linux/x86sim-linux.cpp
)
target_link_libraries(x86sim-linux PRIVATE x86sim::defaults)
target_compile_definitions(x86sim-linux PRIVATE ${build_info_definitions})

install(TARGETS x86sim x86sim_defaults raspsim x86sim-linux
  EXPORT x86simTargets
  FILE_SET HEADERS
)
install(EXPORT x86simTargets
  NAMESPACE x86sim::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/x86sim"
)

option(X86SIM_BUILD_TESTS "Build the CTest suite (SQLite through x86sim-linux)" ON)
if(X86SIM_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests/sqlite3)
endif()
