// Unit tests for the portable in-memory filesystem (x86sim-support/memfs.hpp)
// and the fd bookkeeping of the memfs syscall state
// (x86sim-support/syscall-linux-memfs.hpp). Errno expectations are the raw
// Linux values (ENOENT=2, EEXIST=17, ENOTDIR=20, EISDIR=21, EINVAL=22,
// EFBIG=27, ENOTEMPTY=39, ELOOP=40, EPERM=1).
#include "x86sim-support/memfs.hpp"
#include "x86sim-support/syscall-linux-memfs.hpp"

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace {

int failures = 0;

#define CHECK(condition)                                                                                               \
  do {                                                                                                                 \
    if (!(condition)) {                                                                                                \
      std::fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #condition);                                        \
      ++failures;                                                                                                      \
    }                                                                                                                  \
  } while (0)

using x86sim::memfs::Filesystem;
using x86sim::memfs::Node;

std::span<const std::byte> bytes_of(const char* text) {
  return std::as_bytes(std::span(text, std::strlen(text)));
}

void test_tree_and_resolution() {
  Filesystem fs;
  const auto& root = fs.root();
  CHECK(root->is_dir());
  CHECK(root->nlink == 2);

  auto data = fs.make_directory(*root, "data", 0755);
  CHECK(data.has_value());
  CHECK(root->nlink == 3);
  auto file = fs.create_file(**data, "in.txt", 0644);
  CHECK(file.has_value());

  // Plain resolution, "." and "..", clamping at the base.
  CHECK(fs.resolve("/data/in.txt", root, "/").has_value());
  CHECK(fs.resolve("data/in.txt", root, "/").has_value());
  CHECK(fs.resolve("in.txt", root, "/data").has_value());
  CHECK(fs.resolve("./in.txt", root, "/data").has_value());
  CHECK(fs.resolve("../data/in.txt", root, "/data").has_value());
  CHECK(fs.resolve("/../../data/in.txt", root, "/").has_value()); // ".." clamps at the root
  CHECK(fs.resolve("/", root, "/").value() == root);
  CHECK(fs.resolve(".", root, "/data").value() == *data);

  // Failure modes.
  CHECK(fs.resolve("/missing", root, "/").error() == 2);           // ENOENT
  CHECK(fs.resolve("/data/in.txt/x", root, "/").error() == 20);    // ENOTDIR: file used as directory
  CHECK(fs.resolve("/data/in.txt/", root, "/").error() == 20);     // ENOTDIR: trailing slash on a file
  CHECK(fs.resolve("/data/in.txt/.", root, "/").error() == 20);    // ENOTDIR: trailing "."
  CHECK(fs.resolve("", root, "/").error() == 2);                   // ENOENT: empty path
  CHECK(fs.resolve(std::string(5000, 'a'), root, "/").error() == 36); // ENAMETOOLONG

  // Deterministic listing order (std::map).
  (void)fs.create_file(**data, "b", 0644);
  (void)fs.create_file(**data, "a", 0644);
  std::vector<std::string> names;
  for (const auto& [name, node] : (*data)->dir().entries)
    names.push_back(name);
  CHECK((names == std::vector<std::string>{"a", "b", "in.txt"}));
}

void test_symlinks() {
  Filesystem fs;
  const auto& root = fs.root();
  auto dir = fs.make_directory(*root, "dir", 0755);
  (void)fs.create_file(**dir, "target", 0644);

  CHECK(fs.make_symlink(*root, "abs", "/dir/target").has_value());
  CHECK(fs.make_symlink(*root, "rel", "dir/target").has_value());
  CHECK(fs.make_symlink(*root, "todir", "dir/").has_value());
  CHECK(fs.make_symlink(*root, "dangling", "nowhere").has_value());

  CHECK(fs.resolve("/abs", root, "/").value()->is_file());
  CHECK(fs.resolve("/rel", root, "/").value()->is_file());
  CHECK(fs.resolve("/todir/target", root, "/").value()->is_file());
  CHECK(fs.resolve("/dangling", root, "/").error() == 2); // ENOENT through the dangling link

  // No-follow returns the link node itself.
  auto link_node = fs.resolve("/abs", root, "/", /*follow_final_symlink=*/false);
  CHECK(link_node.has_value() && (*link_node)->is_symlink());
  CHECK((*link_node)->symlink().target == "/dir/target");

  // Loops fail with ELOOP.
  CHECK(fs.make_symlink(*root, "loop_a", "/loop_b").has_value());
  CHECK(fs.make_symlink(*root, "loop_b", "/loop_a").has_value());
  CHECK(fs.resolve("/loop_a", root, "/").error() == 40); // ELOOP

  // Symlinked directories canonicalize to their real path.
  auto canonical = fs.canonicalize_directory("/todir", root, "/");
  CHECK(canonical.has_value() && *canonical == "/dir");
}

void test_mutations() {
  Filesystem fs;
  const auto& root = fs.root();
  auto a = fs.make_directory(*root, "a", 0755);
  auto b = fs.make_directory(*root, "b", 0755);
  auto file = fs.create_file(**a, "f", 0644);

  CHECK(fs.create_file(**a, "f", 0644).error() == 17);     // EEXIST
  CHECK(fs.make_directory(*root, "a", 0755).error() == 17); // EEXIST
  CHECK(fs.unlink(*root, "a").error() == 21);               // EISDIR: unlink a directory
  CHECK(fs.remove_directory(*root, "missing").error() == 2);
  CHECK(fs.remove_directory(**a, "f").error() == 20); // ENOTDIR: rmdir a file
  CHECK(fs.remove_directory(*root, "a").error() == 39); // ENOTEMPTY

  // rename: file over nothing, file over file, directory rules.
  CHECK(fs.rename(**a, "f", **b, "g").has_value());
  CHECK(fs.resolve("/b/g", root, "/").has_value());
  CHECK(fs.resolve("/a/f", root, "/").error() == 2);

  (void)fs.create_file(**b, "h", 0644);
  CHECK(fs.rename(**b, "g", **b, "h").has_value()); // replaces existing file
  CHECK((*b)->dir().entries.size() == 1);

  // Moving a directory into its own subtree is EINVAL.
  auto outer = fs.make_directory(*root, "outer", 0755);
  auto inner = fs.make_directory(**outer, "inner", 0755);
  CHECK(fs.rename(*root, "outer", **inner, "nested").error() == 22);

  // Directory over non-empty directory is ENOTEMPTY.
  auto c = fs.make_directory(*root, "c", 0755);
  (void)fs.create_file(**c, "keep", 0644);
  CHECK(fs.rename(*root, "outer", *root, "c").error() == 39);

  // rename onto the same node is a silent success.
  CHECK(fs.rename(*root, "c", *root, "c").has_value());

  // Directory rename updates nlink accounting.
  const auto root_links_before = root->nlink;
  CHECK(fs.rename(*root, "outer", **c, "moved").has_value());
  CHECK(root->nlink == root_links_before - 1);
  CHECK((*c)->nlink == 3);

  // Hard links: nlink counts, no links to directories.
  auto target = fs.create_file(*root, "t", 0644);
  CHECK(fs.link(*root, "t2", *target).has_value());
  CHECK((*target)->nlink == 2);
  CHECK(fs.link(*root, "cc", *c).error() == 1); // EPERM
  CHECK(fs.unlink(*root, "t2").has_value());
  CHECK((*target)->nlink == 1);
}

void test_file_contents() {
  Filesystem fs;
  const auto& root = fs.root();
  auto file = fs.create_file(*root, "f", 0644);
  Node& node = **file;

  CHECK(fs.write(node, 0, bytes_of("hello")).value() == 5);
  std::vector<std::byte> buffer(16);
  CHECK(fs.read(node, 0, buffer) == 5);
  CHECK(std::memcmp(buffer.data(), "hello", 5) == 0);
  CHECK(fs.read(node, 5, buffer) == 0);  // at EOF
  CHECK(fs.read(node, 99, buffer) == 0); // past EOF
  CHECK(fs.read(node, 3, buffer) == 2);  // short read

  // Sparse write zero-fills the hole.
  CHECK(fs.write(node, 8, bytes_of("x")).value() == 1);
  CHECK(node.size() == 9);
  CHECK(fs.read(node, 5, buffer) == 4);
  CHECK(buffer[0] == std::byte{0} && buffer[1] == std::byte{0} && buffer[2] == std::byte{0});

  // Truncate down and back up zero-fills.
  CHECK(fs.truncate(node, 2).has_value());
  CHECK(node.size() == 2);
  CHECK(fs.truncate(node, 4).has_value());
  CHECK(fs.read(node, 0, buffer) == 4);
  CHECK(std::memcmp(buffer.data(), "he\0\0", 4) == 0);

  // Size cap.
  fs.set_max_file_size(1024);
  CHECK(fs.truncate(node, 2048).error() == 27);       // EFBIG
  CHECK(fs.write(node, 1020, bytes_of("abcde")).error() == 27);
  CHECK(fs.write(node, 1019, bytes_of("abcde")).value() == 5);

  // Directories reject content operations.
  auto dir = fs.make_directory(*root, "d", 0755);
  CHECK(fs.write(**dir, 0, bytes_of("x")).error() == 21); // EISDIR
  CHECK(fs.truncate(**dir, 0).error() == 21);

  // Timestamps advance deterministically with mutations.
  const auto before = node.mtime;
  CHECK(fs.write(node, 0, bytes_of("y")).has_value());
  CHECK(node.mtime > before);
}

void test_fd_state() {
  using x86sim::linux_syscalls::MemFsState;
  using x86sim::linux_syscalls::make_memfs_state;

  auto state = make_memfs_state();
  state->reserved_fds = {0, 1, 2};
  auto& proc = state->proc(1);

  // Allocation skips reserved fds and fds already in the table.
  CHECK(state->allocate_fd(proc) == 3);
  state->install_at(1, 3, MemFsState::OpenFile{});
  CHECK(state->allocate_fd(proc) == 4);
  CHECK(state->allocate_fd(proc, 10) == 10);

  // A tombstoned reserved fd is free for reuse.
  proc.closed_reserved.insert(1);
  CHECK(state->allocate_fd(proc) == 1);

  // install_at may shadow any fd, including reserved ones.
  state->install_at(1, 0, MemFsState::OpenFile{});
  CHECK(proc.fds.contains(0));

  // clone_proc shares open-file descriptions (shared offsets) and namespaces.
  auto file = state->fs->create_file(*state->fs->root(), "f", 0644);
  state->install_at(1, 7, MemFsState::OpenFile{.node = *file});
  state->clone_proc(1, 2);
  auto& child = state->proc(2);
  CHECK(child.fds.at(7).file == state->proc(1).fds.at(7).file);
  CHECK(child.ns.mount == state->proc(1).ns.mount);

  // Distinct states never share filesystems.
  auto other = make_memfs_state();
  CHECK(other->fs != state->fs);
}

} // namespace

int main() {
  test_tree_and_resolution();
  test_symlinks();
  test_mutations();
  test_file_contents();
  test_fd_state();

  if (failures != 0) {
    std::fprintf(stderr, "%d check(s) failed\n", failures);
    return 1;
  }
  std::printf("all memfs checks passed\n");
  return 0;
}
