// Portable in-memory filesystem. Like syscall-linux.cpp, this translation unit
// must not include any host header: everything is deterministic, in-process
// state so the memfs works identically on every platform the wheel builds on.
// Errno values are reported as plain positive ints matching the Linux ABI
// (shared constants live in syscall-linux-detail.hpp).
#include "x86sim-support/memfs.hpp"

#include "syscall-linux-detail.hpp"

#include <algorithm>
#include <deque>
#include <utility>

namespace x86sim::memfs {

namespace errno_values = x86sim::linux_syscalls::detail;

namespace {

// Splits a path into components, dropping empty segments ("//" collapses).
// "." components are kept: the walk uses them to require directory-ness of
// the node they apply to ("a/." fails with ENOTDIR when a is a file).
[[nodiscard]] std::deque<std::string> split_components(std::string_view path) {
  std::deque<std::string> components;
  std::size_t begin = 0;
  while (begin < path.size()) {
    const std::size_t end = std::min(path.find('/', begin), path.size());
    if (end != begin)
      components.emplace_back(path.substr(begin, end - begin));
    begin = end + 1;
  }
  return components;
}

// Walk state: parallel stacks of nodes and their names relative to the base.
// stack.front() is always the base and never popped (".." clamps there, like
// ".." in the real root directory).
struct Walker {
  std::vector<std::shared_ptr<Node>> nodes;
  std::vector<std::string> names;

  explicit Walker(std::shared_ptr<Node> base) { reset(std::move(base)); }

  void reset(std::shared_ptr<Node> base) {
    nodes.clear();
    names.clear();
    nodes.push_back(std::move(base));
  }

  [[nodiscard]] const std::shared_ptr<Node>& current() const { return nodes.back(); }

  void push(std::shared_ptr<Node> node, std::string name) {
    nodes.push_back(std::move(node));
    names.push_back(std::move(name));
  }

  void pop() {
    if (nodes.size() > 1) {
      nodes.pop_back();
      names.pop_back();
    }
  }
};

[[nodiscard]] bool node_contains(const Node* directory, const Node* candidate) {
  if (directory == candidate)
    return true;
  if (!directory->is_dir())
    return false;
  for (const auto& [name, child] : directory->dir().entries)
    if (node_contains(child.get(), candidate))
      return true;
  return false;
}

} // namespace

std::uint64_t Node::size() const noexcept {
  if (is_file())
    return file().bytes.size();
  if (is_symlink())
    return symlink().target.size();
  return 4096; // conventional directory size, matching tmpfs
}

Filesystem::Filesystem() : root_(make_node(mode_directory | default_directory_permissions)) {
  root_->nlink = 2;
}

std::shared_ptr<Node> Filesystem::make_node(std::uint32_t mode) {
  auto node = std::make_shared<Node>();
  node->ino = next_ino_++;
  node->mode = mode;
  node->nlink = 1;
  node->atime = node->mtime = node->ctime = next_tick();
  switch (mode & mode_type_mask) {
  case mode_directory:
    node->data = DirData{};
    node->nlink = 2;
    break;
  case mode_symlink:
    node->data = SymlinkData{};
    break;
  default:
    node->data = FileData{};
    break;
  }
  return node;
}

// Shared resolution loop. Consumes `work` against the walker; expands symlinks
// (intermediate always, final per follow_final_symlink). Returns a Linux errno
// on failure.
namespace {

[[nodiscard]] std::expected<void, int> walk(Walker& walker, std::deque<std::string>& work,
                                            const std::shared_ptr<Node>& base, bool follow_final_symlink) {
  std::size_t symlink_hops = 0;
  for (;;) {
    // Expand a symlink on top of the stack whenever the walk must continue
    // through it (more components pending) or the caller wants the target.
    if (walker.current()->is_symlink() && (!work.empty() || follow_final_symlink)) {
      if (++symlink_hops > max_symlink_hops)
        return std::unexpected(errno_values::linux_eloop);
      const std::string target = walker.current()->symlink().target;
      if (target.empty())
        return std::unexpected(errno_values::linux_enoent);
      walker.pop();
      auto expanded = split_components(target);
      if (!target.empty() && target.back() == '/')
        expanded.emplace_back(".");
      work.insert(work.begin(), expanded.begin(), expanded.end());
      if (target.front() == '/')
        walker.reset(base);
      continue;
    }

    if (work.empty())
      return {};

    const std::string component = std::move(work.front());
    work.pop_front();

    const std::shared_ptr<Node>& current = walker.current();
    if (!current->is_dir())
      return std::unexpected(errno_values::linux_enotdir);

    if (component == ".")
      continue;
    if (component == "..") {
      walker.pop();
      continue;
    }

    auto it = current->dir().entries.find(component);
    if (it == current->dir().entries.end())
      return std::unexpected(errno_values::linux_enoent);
    walker.push(it->second, component);
  }
}

// Builds the walker start state and work list for `path` under base/cwd.
[[nodiscard]] std::expected<std::pair<Walker, std::deque<std::string>>, int>
prepare_walk(std::string_view path, const std::shared_ptr<Node>& base, std::string_view cwd) {
  if (path.empty())
    return std::unexpected(errno_values::linux_enoent);
  if (path.size() >= max_path_length)
    return std::unexpected(errno_values::linux_enametoolong);

  Walker walker(base);
  std::deque<std::string> work;
  if (path.front() != '/') {
    // Relative paths start at the cwd, itself a canonical path under base.
    auto cwd_components = split_components(cwd);
    work.insert(work.end(), cwd_components.begin(), cwd_components.end());
  }
  auto components = split_components(path);
  work.insert(work.end(), components.begin(), components.end());
  // A trailing slash requires the result to be a directory; an explicit "."
  // component enforces that through the walk (a lone "/" becomes just ".").
  if (path.back() == '/')
    work.emplace_back(".");
  return std::make_pair(std::move(walker), std::move(work));
}

} // namespace

std::expected<std::shared_ptr<Node>, int> Filesystem::resolve(std::string_view path,
                                                              const std::shared_ptr<Node>& base, std::string_view cwd,
                                                              bool follow_final_symlink) const {
  auto prepared = prepare_walk(path, base, cwd);
  if (!prepared)
    return std::unexpected(prepared.error());
  auto& [walker, work] = *prepared;
  if (auto walked = walk(walker, work, base, follow_final_symlink); !walked)
    return std::unexpected(walked.error());
  return walker.current();
}

std::expected<Filesystem::ParentRef, int> Filesystem::resolve_parent(std::string_view path,
                                                                     const std::shared_ptr<Node>& base,
                                                                     std::string_view cwd) const {
  auto prepared = prepare_walk(path, base, cwd);
  if (!prepared)
    return std::unexpected(prepared.error());
  auto& [walker, work] = *prepared;

  // Strip the leaf so the walk stops at the parent. A leaf of "." or ".." (or
  // a path naming the base itself) is degenerate: resolve the whole path and
  // report an empty leaf, letting each caller pick its errno.
  std::string leaf;
  while (!work.empty() && work.back() == ".")
    work.pop_back();
  if (!work.empty() && work.back() != "..") {
    leaf = std::move(work.back());
    work.pop_back();
  }

  if (auto walked = walk(walker, work, base, /*follow_final_symlink=*/true); !walked)
    return std::unexpected(walked.error());
  if (!walker.current()->is_dir())
    return std::unexpected(errno_values::linux_enotdir);
  return ParentRef{.dir = walker.current(), .leaf = std::move(leaf)};
}

std::expected<std::string, int> Filesystem::canonicalize_directory(std::string_view path,
                                                                   const std::shared_ptr<Node>& base,
                                                                   std::string_view cwd) const {
  auto prepared = prepare_walk(path, base, cwd);
  if (!prepared)
    return std::unexpected(prepared.error());
  auto& [walker, work] = *prepared;
  if (auto walked = walk(walker, work, base, /*follow_final_symlink=*/true); !walked)
    return std::unexpected(walked.error());
  if (!walker.current()->is_dir())
    return std::unexpected(errno_values::linux_enotdir);

  std::string canonical = "/";
  for (std::size_t i = 0; i < walker.names.size(); ++i) {
    if (i > 0)
      canonical += '/';
    canonical += walker.names[i];
  }
  return canonical;
}

std::expected<std::shared_ptr<Node>, int> Filesystem::create_file(Node& parent, const std::string& leaf,
                                                                  std::uint32_t permissions) {
  if (leaf.empty())
    return std::unexpected(errno_values::linux_eexist);
  if (!parent.is_dir())
    return std::unexpected(errno_values::linux_enotdir);
  if (parent.dir().entries.contains(leaf))
    return std::unexpected(errno_values::linux_eexist);
  auto node = make_node(mode_regular | (permissions & mode_permission_mask));
  parent.dir().entries.emplace(leaf, node);
  parent.mtime = parent.ctime = next_tick();
  return node;
}

std::expected<std::shared_ptr<Node>, int> Filesystem::make_directory(Node& parent, const std::string& leaf,
                                                                     std::uint32_t permissions) {
  if (leaf.empty())
    return std::unexpected(errno_values::linux_eexist);
  if (!parent.is_dir())
    return std::unexpected(errno_values::linux_enotdir);
  if (parent.dir().entries.contains(leaf))
    return std::unexpected(errno_values::linux_eexist);
  auto node = make_node(mode_directory | (permissions & mode_permission_mask));
  parent.dir().entries.emplace(leaf, node);
  ++parent.nlink; // the child's ".." entry
  parent.mtime = parent.ctime = next_tick();
  return node;
}

std::expected<std::shared_ptr<Node>, int> Filesystem::make_symlink(Node& parent, const std::string& leaf,
                                                                   std::string target) {
  if (leaf.empty())
    return std::unexpected(errno_values::linux_eexist);
  if (!parent.is_dir())
    return std::unexpected(errno_values::linux_enotdir);
  if (parent.dir().entries.contains(leaf))
    return std::unexpected(errno_values::linux_eexist);
  if (target.empty())
    return std::unexpected(errno_values::linux_enoent);
  auto node = make_node(mode_symlink | 0777);
  node->symlink().target = std::move(target);
  parent.dir().entries.emplace(leaf, node);
  parent.mtime = parent.ctime = next_tick();
  return node;
}

std::expected<void, int> Filesystem::link(Node& parent, const std::string& leaf, std::shared_ptr<Node> target) {
  if (leaf.empty())
    return std::unexpected(errno_values::linux_eexist);
  if (!parent.is_dir())
    return std::unexpected(errno_values::linux_enotdir);
  if (!target || target->is_dir())
    return std::unexpected(errno_values::linux_eperm); // hard links to directories are forbidden
  if (parent.dir().entries.contains(leaf))
    return std::unexpected(errno_values::linux_eexist);
  ++target->nlink;
  target->ctime = next_tick();
  parent.dir().entries.emplace(leaf, std::move(target));
  parent.mtime = parent.ctime = tick_;
  return {};
}

std::expected<void, int> Filesystem::unlink(Node& parent, const std::string& leaf) {
  if (leaf.empty())
    return std::unexpected(errno_values::linux_eisdir);
  if (!parent.is_dir())
    return std::unexpected(errno_values::linux_enotdir);
  auto it = parent.dir().entries.find(leaf);
  if (it == parent.dir().entries.end())
    return std::unexpected(errno_values::linux_enoent);
  if (it->second->is_dir())
    return std::unexpected(errno_values::linux_eisdir);
  if (it->second->nlink > 0)
    --it->second->nlink;
  it->second->ctime = next_tick();
  parent.dir().entries.erase(it);
  parent.mtime = parent.ctime = tick_;
  return {};
}

std::expected<void, int> Filesystem::remove_directory(Node& parent, const std::string& leaf) {
  if (leaf.empty())
    return std::unexpected(errno_values::linux_einval); // rmdir "." / the namespace root
  if (!parent.is_dir())
    return std::unexpected(errno_values::linux_enotdir);
  auto it = parent.dir().entries.find(leaf);
  if (it == parent.dir().entries.end())
    return std::unexpected(errno_values::linux_enoent);
  if (!it->second->is_dir())
    return std::unexpected(errno_values::linux_enotdir);
  if (!it->second->dir().entries.empty())
    return std::unexpected(errno_values::linux_enotempty);
  parent.dir().entries.erase(it);
  --parent.nlink;
  parent.mtime = parent.ctime = next_tick();
  return {};
}

std::expected<void, int> Filesystem::rename(Node& source_parent, const std::string& source_leaf, Node& target_parent,
                                            const std::string& target_leaf) {
  if (source_leaf.empty() || target_leaf.empty())
    return std::unexpected(errno_values::linux_einval);
  if (!source_parent.is_dir() || !target_parent.is_dir())
    return std::unexpected(errno_values::linux_enotdir);

  auto source_it = source_parent.dir().entries.find(source_leaf);
  if (source_it == source_parent.dir().entries.end())
    return std::unexpected(errno_values::linux_enoent);
  std::shared_ptr<Node> node = source_it->second;

  auto target_it = target_parent.dir().entries.find(target_leaf);
  if (target_it != target_parent.dir().entries.end() && target_it->second == node)
    return {}; // same file: rename(2) succeeds without doing anything

  if (node->is_dir() && node_contains(node.get(), &target_parent))
    return std::unexpected(errno_values::linux_einval); // moving a directory into its own subtree

  if (target_it != target_parent.dir().entries.end()) {
    const std::shared_ptr<Node>& existing = target_it->second;
    if (node->is_dir()) {
      if (!existing->is_dir())
        return std::unexpected(errno_values::linux_enotdir);
      if (!existing->dir().entries.empty())
        return std::unexpected(errno_values::linux_enotempty);
      --target_parent.nlink;
    } else if (existing->is_dir()) {
      return std::unexpected(errno_values::linux_eisdir);
    } else if (existing->nlink > 0) {
      --existing->nlink;
    }
    target_parent.dir().entries.erase(target_it);
  }

  source_parent.dir().entries.erase(source_it);
  target_parent.dir().entries.emplace(target_leaf, node);
  if (node->is_dir()) {
    --source_parent.nlink;
    ++target_parent.nlink;
  }
  node->ctime = next_tick();
  source_parent.mtime = source_parent.ctime = tick_;
  target_parent.mtime = target_parent.ctime = tick_;
  return {};
}

std::size_t Filesystem::read(const Node& node, std::uint64_t offset, std::span<std::byte> out) const {
  if (!node.is_file())
    return 0;
  const auto& bytes = node.file().bytes;
  if (offset >= bytes.size())
    return 0;
  const std::size_t available = bytes.size() - static_cast<std::size_t>(offset);
  const std::size_t count = std::min(out.size(), available);
  std::copy_n(bytes.begin() + static_cast<std::ptrdiff_t>(offset), count, out.begin());
  return count;
}

std::expected<std::size_t, int> Filesystem::write(Node& node, std::uint64_t offset, std::span<const std::byte> in) {
  if (node.is_dir())
    return std::unexpected(errno_values::linux_eisdir);
  if (!node.is_file())
    return std::unexpected(errno_values::linux_einval);
  if (offset > max_file_size_ || in.size() > max_file_size_ - offset)
    return std::unexpected(errno_values::linux_efbig);

  auto& bytes = node.file().bytes;
  const std::uint64_t end = offset + in.size();
  if (end > bytes.size())
    bytes.resize(static_cast<std::size_t>(end)); // zero-fills any hole before offset
  std::copy(in.begin(), in.end(), bytes.begin() + static_cast<std::ptrdiff_t>(offset));
  node.mtime = node.ctime = next_tick();
  return in.size();
}

std::expected<void, int> Filesystem::truncate(Node& node, std::uint64_t length) {
  if (node.is_dir())
    return std::unexpected(errno_values::linux_eisdir);
  if (!node.is_file())
    return std::unexpected(errno_values::linux_einval);
  if (length > max_file_size_)
    return std::unexpected(errno_values::linux_efbig);
  node.file().bytes.resize(static_cast<std::size_t>(length)); // growth zero-fills
  node.mtime = node.ctime = next_tick();
  return {};
}

StatInfo Filesystem::stat(const Node& node) const noexcept {
  return StatInfo{
      .ino = node.ino,
      .mode = node.mode,
      .nlink = node.nlink,
      .size = node.size(),
      .atime = node.atime,
      .mtime = node.mtime,
      .ctime = node.ctime,
  };
}

void Filesystem::set_permissions(Node& node, std::uint32_t permissions) noexcept {
  node.mode = (node.mode & mode_type_mask) | (permissions & mode_permission_mask);
  node.ctime = next_tick();
}

void Filesystem::set_times(Node& node, std::uint64_t atime, std::uint64_t mtime) noexcept {
  node.atime = atime;
  node.mtime = mtime;
  node.ctime = next_tick();
}

void Filesystem::touch(Node& node) noexcept {
  node.atime = node.mtime = node.ctime = next_tick();
}

bool Filesystem::is_ancestor_of(const std::shared_ptr<Node>& node, const std::shared_ptr<Node>& candidate) const {
  if (!node || !candidate)
    return false;
  return node_contains(node.get(), candidate.get());
}

} // namespace x86sim::memfs
