from __future__ import annotations

from typing import IO, Any, Callable, Literal
from pathlib import Path

from . import bindings
from .bindings import Machine as _Machine, RaspsimException
from .elf import ELF
from .fs import MemPath


class Machine(_Machine):
    def __init__(
        self,
        *,
        glibc: bool = False,
        stdin: IO[bytes] | Any | None = None,
        stdout: IO[bytes] | Any | None = None,
        stderr: IO[bytes] | Any | None = None,
        readlink: Callable[[str], str | None] | None = None,
        core: Literal["ooo", "out_of_order", "seq", "sequential"] = "ooo",
        memfs: bool | bindings.Fs | MemPath = False,
    ):
        """
        Create a new Machine instance.

        Args:
            glibc: Enable the portable Linux glibc-startup syscalls: the
                malloc/free heap (brk + anonymous mmap/munmap/mremap) plus
                arch_prctl, set_tid_address, set_robust_list, rseq, prlimit64,
                uname, futex and the synthetic signal syscalls. Off by default.
            stdin: A file-like object exposing ``.read(n)`` that the guest's
                ``read(0, ...)`` is routed to (e.g. ``io.BytesIO``). Host fds are
                never used; ``None`` leaves fd 0 unconfigured.
            stdout: A file-like object exposing ``.write(bytes)`` that the
                guest's ``write(1, ...)`` is routed to. ``None`` leaves fd 1
                unconfigured.
            stderr: As ``stdout`` but for the guest's ``write(2, ...)``.
            readlink: A callable ``readlink(path: str) -> str | None`` that
                resolves guest ``readlink``/``readlinkat`` calls (e.g. for
                ``/proc/self/exe``). Returning ``None`` yields ``-ENOENT``. When
                unset, the guest's readlink returns ``-ENOSYS``.
            core: CPU model, ``"ooo"`` (out-of-order, default) or ``"seq"``
                (sequential). Both execute unaligned memory accesses correctly
                and can run glibc; the sequential core is simpler and slower,
                the out-of-order core is the cycle-accurate default.
            memfs: Opt-in in-memory Linux-like filesystem for the guest.
                ``True`` creates a fresh sandbox; an existing ``bindings.Fs``
                or ``MemPath`` shares its filesystem (prepopulation, or
                sharing between Machines). The Unix file syscalls (open, read,
                write, stat, getdents64, mkdir, ...) then work against the
                sandbox, which the host can inspect through :attr:`fs`. No fd
                number is special: stdio stays with the ``stdin``/``stdout``/
                ``stderr`` streams unless the guest ``dup2()``s over them or
                :meth:`map_fd` binds a memfs file there. Mutually exclusive
                with ``readlink``.
        """
        if isinstance(memfs, MemPath):
            memfs = memfs.filesystem
        super().__init__(
            glibc=glibc,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            readlink=readlink,
            core=core,
            memfs=memfs,
        )
        self._current_elf: ELF | None = None

    @property
    def fs(self) -> MemPath:
        """Root of the guest's in-memory filesystem as a pathlib-like object.

        Only available when the Machine was constructed with ``memfs``.
        """
        return MemPath(super().fs)

    def map_fd(self, fd: int, path: MemPath | str, mode: str = "r") -> int:
        """Bind a memfs file to an arbitrary guest fd number.

        This is how stdio becomes a memfs file: setting up fds 0/1/2 (or any
        other number) is entirely the caller's choice. ``mode`` follows
        ``open()``: ``"r"``, ``"w"``, ``"a"``, optionally with ``"+"``.
        """
        return super().map_fd(fd, str(path), mode)

    def load_elf(self, elf: ELF, abi: Literal["sysv"] = "sysv") -> None:
        """
        Populate the simulator with the ELF binary data.

        Args:
            sim: Machine instance.
            elf: ELF binary data.
            abi: ABI to use. Currently only supports sysv.
        """
        del abi  # unused

        self._current_elf = elf
        self.registers.rip = elf.entry

        if elf.stack is None:
            raise RaspsimException("stack must be set")
        self.registers.rsp = elf.stack

        for segment in elf.segments:
            self.memmap(
                segment.vaddr,
                segment.prot,
                length=segment.size,
                data=segment.data or b"",
            )

    def run(self, ninstr: int | None = None) -> None:
        """
        Run the simulator.
        """
        try:
            if ninstr is not None:
                super().run(ninstr)
            else:
                super().run()
        except Exception as e:
            path: str | None = None
            lineno: int | None = None
            line: str | None = None
            if self._current_elf is not None:
                path, lineno = self._current_elf.loc(self.registers.rip)
                try:
                    lines = Path(path).read_text().splitlines()
                    line = lines[lineno - 1]
                except Exception as ex:
                    line = "Internal error: " + str(ex)

            setattr(e, "file", path)
            setattr(e, "lineno", lineno)
            setattr(e, "line", line)
            raise e
