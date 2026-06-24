cmake_minimum_required(VERSION 3.25)

project(raspsim
  DESCRIPTION "Cycle-accurate x86-64 virtual CPU core and RASPsim command-line simulator"
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

add_library(vcore_defaults)
add_library(x86sim::defaults ALIAS vcore_defaults)
target_sources(vcore_defaults
  PRIVATE
    src/x86sim-support/raspsim-defaults.cpp
    src/x86sim-support/syscall-linux.cpp
  PUBLIC
    FILE_SET HEADERS
    BASE_DIRS include
    FILES
      include/x86sim-support/defaults.hpp
      include/x86sim-support/syscall-linux.hpp
)
target_link_libraries(vcore_defaults PUBLIC x86sim::x86sim)
target_include_directories(vcore_defaults PRIVATE src/x86sim)
set_target_properties(vcore_defaults PROPERTIES EXPORT_NAME defaults)

add_executable(raspsim
  src/raspsim/config.cpp
  src/raspsim/raspsim.cpp
)
target_link_libraries(raspsim PRIVATE x86sim::defaults)
target_include_directories(raspsim PRIVATE src/x86sim)
target_compile_definitions(raspsim PRIVATE ${build_info_definitions})

install(TARGETS x86sim vcore_defaults raspsim
  EXPORT vcoreTargets
  FILE_SET HEADERS
)
install(EXPORT vcoreTargets
  NAMESPACE x86sim::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/x86sim"
)
