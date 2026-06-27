from __future__ import annotations

from typing import IO, Any, Literal
from pathlib import Path

from .bindings import Machine as _Machine, RaspsimException
from .elf import ELF


class Machine(_Machine):
    def __init__(
        self,
        *,
        heap: bool = False,
        stdin: IO[bytes] | Any | None = None,
        stdout: IO[bytes] | Any | None = None,
        stderr: IO[bytes] | Any | None = None,
    ):
        """
        Create a new Machine instance.

        Args:
            heap: Enable the Linux malloc/free heap (brk + anonymous
                mmap/munmap/mremap) for the guest. Off by default.
            stdin: A file-like object exposing ``.read(n)`` that the guest's
                ``read(0, ...)`` is routed to (e.g. ``io.BytesIO``). Host fds are
                never used; ``None`` leaves fd 0 unconfigured.
            stdout: A file-like object exposing ``.write(bytes)`` that the
                guest's ``write(1, ...)`` is routed to. ``None`` leaves fd 1
                unconfigured.
            stderr: As ``stdout`` but for the guest's ``write(2, ...)``.
        """
        super().__init__(heap=heap, stdin=stdin, stdout=stdout, stderr=stderr)
        self._current_elf: ELF | None = None

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
