from __future__ import annotations

from typing import Literal
from pathlib import Path

from .core import Core, RaspsimException
from .elf import ELF


class Raspsim(Core):
    def __init__(self):
        super().__init__()
        self._current_elf: ELF | None = None

    def load_elf(self, elf: ELF, abi: Literal["sysv"] = "sysv") -> None:
        """
        Populate the simulator with the ELF binary data.

        Args:
            sim: Raspsim instance.
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
