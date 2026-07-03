// Portable, deterministic in-memory filesystem with Linux semantics.
//
// This is the backing store for the opt-in memfs syscall personality
// (x86sim-support/syscall-linux-memfs.hpp) and the host-side query API exposed
// through the Python bindings. It is deliberately independent of the syscall
// layer: operations report failures as plain Linux errno values (positive
// ints) inside std::expected, and nothing here touches a host OS header, the
// host clock, or host entropy — a fixed tick counter that advances on every
// mutation stands in for time, so simulations stay reproducible on any
// platform.
#ifndef X86SIM_SUPPORT_MEMFS_HPP
#define X86SIM_SUPPORT_MEMFS_HPP

#include <cstddef>
#include <cstdint>
#include <expected>
#include <map>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace x86sim::memfs {

// Linux st_mode file-type bits (values match the Linux ABI).
inline constexpr std::uint32_t mode_type_mask = 0170000;
inline constexpr std::uint32_t mode_regular = 0100000;
inline constexpr std::uint32_t mode_directory = 0040000;
inline constexpr std::uint32_t mode_symlink = 0120000;
inline constexpr std::uint32_t mode_permission_mask = 07777;

inline constexpr std::uint32_t default_file_permissions = 0644;
inline constexpr std::uint32_t default_directory_permissions = 0755;

inline constexpr std::size_t max_symlink_hops = 40;
inline constexpr std::size_t max_path_length = 4096;

struct Node;

struct FileData {
  std::vector<std::byte> bytes;
};

// std::map keeps directory listings in a deterministic order and gives
// getdents64 a cursor (upper_bound on the last returned name) that stays valid
// across concurrent directory mutation.
struct DirData {
  std::map<std::string, std::shared_ptr<Node>> entries;
};

struct SymlinkData {
  std::string target;
};

// One inode. Children are owned by their parent directory only (no parent
// back-pointers), so shared_ptr ownership is acyclic. An open-file handle that
// outlives an unlink keeps the node alive, matching POSIX unlink semantics.
struct Node {
  std::uint64_t ino = 0;
  std::uint32_t mode = 0;
  std::uint32_t nlink = 0;
  std::uint64_t atime = 0;
  std::uint64_t mtime = 0;
  std::uint64_t ctime = 0;
  std::variant<FileData, DirData, SymlinkData> data;

  [[nodiscard]] bool is_file() const noexcept { return std::holds_alternative<FileData>(data); }
  [[nodiscard]] bool is_dir() const noexcept { return std::holds_alternative<DirData>(data); }
  [[nodiscard]] bool is_symlink() const noexcept { return std::holds_alternative<SymlinkData>(data); }

  [[nodiscard]] FileData& file() { return std::get<FileData>(data); }
  [[nodiscard]] const FileData& file() const { return std::get<FileData>(data); }
  [[nodiscard]] DirData& dir() { return std::get<DirData>(data); }
  [[nodiscard]] const DirData& dir() const { return std::get<DirData>(data); }
  [[nodiscard]] SymlinkData& symlink() { return std::get<SymlinkData>(data); }
  [[nodiscard]] const SymlinkData& symlink() const { return std::get<SymlinkData>(data); }

  [[nodiscard]] std::uint64_t size() const noexcept;
};

struct StatInfo {
  std::uint64_t ino = 0;
  std::uint32_t mode = 0;
  std::uint32_t nlink = 0;
  std::uint64_t size = 0;
  std::uint64_t atime = 0;
  std::uint64_t mtime = 0;
  std::uint64_t ctime = 0;
};

class Filesystem {
public:
  Filesystem();

  Filesystem(const Filesystem&) = delete;
  Filesystem& operator=(const Filesystem&) = delete;

  [[nodiscard]] const std::shared_ptr<Node>& root() const noexcept { return root_; }

  // --- path resolution ------------------------------------------------------
  // Paths resolve relative to `base` (the namespace root; "/" and ".." clamp
  // there) and `cwd`, a canonical absolute path interpreted under `base`.
  // Intermediate symlinks always resolve; the final component follows
  // `follow_final_symlink`. Errors: ENOENT, ENOTDIR, ELOOP, ENAMETOOLONG.
  [[nodiscard]] std::expected<std::shared_ptr<Node>, int> resolve(std::string_view path,
                                                                  const std::shared_ptr<Node>& base,
                                                                  std::string_view cwd,
                                                                  bool follow_final_symlink = true) const;

  // Resolves everything but the last component. `leaf` is empty when the path
  // names the root itself.
  struct ParentRef {
    std::shared_ptr<Node> dir;
    std::string leaf;
  };
  [[nodiscard]] std::expected<ParentRef, int> resolve_parent(std::string_view path, const std::shared_ptr<Node>& base,
                                                             std::string_view cwd) const;

  // Canonical absolute form of `path` ("/a/b", no ".", "..", or duplicate
  // slashes), for chdir/getcwd bookkeeping. Fails if the path does not name an
  // existing directory.
  [[nodiscard]] std::expected<std::string, int> canonicalize_directory(std::string_view path,
                                                                       const std::shared_ptr<Node>& base,
                                                                       std::string_view cwd) const;

  // --- namespace mutations ---------------------------------------------------
  [[nodiscard]] std::expected<std::shared_ptr<Node>, int> create_file(Node& parent, const std::string& leaf,
                                                                      std::uint32_t permissions);
  [[nodiscard]] std::expected<std::shared_ptr<Node>, int> make_directory(Node& parent, const std::string& leaf,
                                                                         std::uint32_t permissions);
  [[nodiscard]] std::expected<std::shared_ptr<Node>, int> make_symlink(Node& parent, const std::string& leaf,
                                                                       std::string target);
  [[nodiscard]] std::expected<void, int> link(Node& parent, const std::string& leaf, std::shared_ptr<Node> target);
  [[nodiscard]] std::expected<void, int> unlink(Node& parent, const std::string& leaf);
  [[nodiscard]] std::expected<void, int> remove_directory(Node& parent, const std::string& leaf);
  [[nodiscard]] std::expected<void, int> rename(Node& source_parent, const std::string& source_leaf,
                                                Node& target_parent, const std::string& target_leaf);

  // --- file content -----------------------------------------------------------
  // Reads never fail: they return 0 at or past EOF (short reads otherwise).
  [[nodiscard]] std::size_t read(const Node& node, std::uint64_t offset, std::span<std::byte> out) const;
  // Writing past EOF zero-fills the hole. EFBIG past the size cap, EISDIR on
  // directories.
  [[nodiscard]] std::expected<std::size_t, int> write(Node& node, std::uint64_t offset,
                                                      std::span<const std::byte> in);
  [[nodiscard]] std::expected<void, int> truncate(Node& node, std::uint64_t length);

  // --- metadata ---------------------------------------------------------------
  [[nodiscard]] StatInfo stat(const Node& node) const noexcept;
  void set_permissions(Node& node, std::uint32_t permissions) noexcept;
  void set_times(Node& node, std::uint64_t atime, std::uint64_t mtime) noexcept;
  void touch(Node& node) noexcept;

  // True when `node` is `candidate` or one of its ancestors (used by rename to
  // reject moving a directory into its own subtree).
  [[nodiscard]] bool is_ancestor_of(const std::shared_ptr<Node>& node, const std::shared_ptr<Node>& candidate) const;

  [[nodiscard]] std::uint64_t max_file_size() const noexcept { return max_file_size_; }
  void set_max_file_size(std::uint64_t limit) noexcept { max_file_size_ = limit; }

  [[nodiscard]] std::uint64_t current_tick() const noexcept { return tick_; }

private:
  [[nodiscard]] std::uint64_t next_tick() noexcept { return ++tick_; }
  [[nodiscard]] std::shared_ptr<Node> make_node(std::uint32_t mode);

  std::shared_ptr<Node> root_;
  std::uint64_t next_ino_ = 1;
  std::uint64_t tick_ = 0;
  std::uint64_t max_file_size_ = std::uint64_t{1} << 30; // 1 GiB
};

} // namespace x86sim::memfs

#endif
