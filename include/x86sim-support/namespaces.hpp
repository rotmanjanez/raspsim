// Per-process namespace attachments for the Linux syscall support library.
//
// Everything a syscall handler may derive from process context alone hangs off
// ProcessNamespaces. Today that is only the mount namespace (which memfs tree a
// process sees, and which node acts as its root); the struct is deliberately a
// seam for future per-process policy — capability sets, allowed-syscall
// filters, pid/net namespaces — which attach here without new plumbing because
// every handler already receives a ProcessId and can look its namespaces up.
//
// Sharing model: clone/fork shares the parent's ProcessNamespaces members by
// shared_ptr copy (the behavior of fork() without CLONE_NEW*); an
// unshare()-style deep copy is the future hook for namespace isolation.
#ifndef X86SIM_SUPPORT_NAMESPACES_HPP
#define X86SIM_SUPPORT_NAMESPACES_HPP

#include "x86sim-support/memfs.hpp"

#include <memory>

namespace x86sim::linux_syscalls {

// Which filesystem view a process has. `root` is the node "/" resolves to and
// ".." clamps at — always fs->root() until chroot/pivot_root support lands.
struct MountNamespace {
  std::shared_ptr<memfs::Filesystem> fs;
  std::shared_ptr<memfs::Node> root;
};

[[nodiscard]] inline std::shared_ptr<MountNamespace>
make_mount_namespace(std::shared_ptr<memfs::Filesystem> fs) {
  auto ns = std::make_shared<MountNamespace>();
  ns->root = fs->root();
  ns->fs = std::move(fs);
  return ns;
}

struct ProcessNamespaces {
  std::shared_ptr<MountNamespace> mount;
  // future: std::shared_ptr<Capabilities> capabilities;
  // future: std::shared_ptr<SyscallFilter> allowed_syscalls;
};

} // namespace x86sim::linux_syscalls

#endif
