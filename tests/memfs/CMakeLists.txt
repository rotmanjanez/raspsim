# Unit tests for the portable in-memory filesystem. Plain assert-style
# executable (the repo carries no test framework dependency).
add_executable(memfs-unit memfs_unit.cpp)
target_link_libraries(memfs-unit PRIVATE x86sim::defaults)

add_test(NAME memfs-unit COMMAND memfs-unit)
