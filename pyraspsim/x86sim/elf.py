from __future__ import annotations
import logging
from dataclasses import dataclass, field
from typing import IO

from elftools.elf.elffile import ELFFile

from .core import Prot, getProtFromELFSegment, RaspsimException


@dataclass
class Segment:
    """
    A class representing a memory segment.
    Attributes:
        vaddr (int): The virtual address of the segment.
        prot (Prot): The protection attributes of the segment.
        size (int): The size of the segment.
        data (bytes | None): Optional data contained in the segment.
    Methods:
        __repr__(): Returns a string representation of the Segment instance.
    Note:
        Length of `data` shall not exceed the size of the segment.
        If size is greater than the length of the data, the rest of the segment denotes zeroed memory.
    """

    vaddr: int
    prot: Prot
    size: int
    data: bytes | None = None
    stack: bool = False

    def __repr__(self):
        data_repr = f" data={self.data!r}" if self.data else ""
        return f"Segment(vaddr={self.vaddr:#x}, prot={self.prot}, size={self.size:#x}{data_repr})"


class InvalidELFException(Exception):
    """
    Exception raised when an invalid ELF file is encountered.
    """


@dataclass
class ELF:
    """
    ELF binary data.

    Attributes:
        entry (int): The entry point of the ELF binary.
        stack (int | None): The stack pointer, if available.
        segments (list[Segment]): A list of segments loaded into memory.
    Methods:
        __repr__(): Returns a string representation of the ELF object.
    Note:
        Only contains the segments that are loaded into memory. Requires x86_64 executable ELF.
    """

    entry: int
    stack: int | None = None
    segments: list[Segment] = field(default_factory=list)
    line_map: list[tuple[int, str, int]] = field(default_factory=list)

    def loc(self, addr: int) -> tuple[str, int]:
        """
        Get the filename and line number of the address.

        Args:
            addr (int): The address to look up.

        Returns:
            tuple[str, int]: A tuple containing the filename and line number.
        """
        for a, filename, line in self.line_map:
            if a < addr:
                return filename, line

        return "<unknown>", 0

    def __repr__(self):
        line_map_repr = ""
        if self.line_map:
            line_map_repr = (
                ", line_map=["
                + ", ".join(
                    f"({hex(addr)}, {filename}:{line})"
                    for addr, filename, line in self.line_map
                )
                + "]"
            )

        return f"ELF(entry={self.entry:#x}, stack={self.stack:#x}, segments={self.segments}{line_map_repr})"

    @staticmethod
    def _extract_linemap(elffile):
        if not elffile.has_dwarf_info():
            return []

        linemap: list[tuple[int, str, int]] = []
        dwarfinfo = elffile.get_dwarf_info()
        for cu in dwarfinfo.iter_CUs():
            line_program = dwarfinfo.line_program_for_CU(cu)
            if not line_program:
                continue
            file_entries = line_program.header["file_entry"]
            for lpe in line_program.get_entries():
                if not lpe.state:
                    continue
                if lpe.state.file:
                    filename = file_entries[lpe.state.file - 1].name.decode()
                else:
                    filename = "<unknown>"
                if linemap and linemap[-1][0] == lpe.state.address:
                    continue
                linemap.append((lpe.state.address, filename, lpe.state.line))

        linemap = sorted(linemap, key=lambda x: x[0])
        linemap.reverse()
        return linemap

    @classmethod
    def from_file(
        cls,
        f: IO[bytes],
        gnu_stack_size: int = 0x1000,
        gnu_stack_addr: int = 0x7FE000000000,
    ) -> ELF:
        """
        Load an ELF file and parse its segments.

        Args:
            elf (IO[bytes]): A byte stream representing the ELF file.
            gnu_stack_size (int, optional): The size of the stack when a GNU Stack Segment is found. Defaults to 0x1000.
            gnu_stack_addr (int, optional): The address of the segment when a GNU Stack Segment is found. Defaults to 0x7FE000000000.

        Returns:
            ELF: An ELF object containing the parsed segments and entry point.

        Raises:
            InvalidELFException: If the ELF file is not of type ET_EXEC or machine type EM_X86_64.
            InvalidELFException: If a segment's file size is greater than its memory size.
            ValueError: If the segment data length does not match the segment's file size.
            InvalidELFException: If an unhandled segment type is encountered.
        """
        # new stream from bytes
        elf = ELFFile(f)

        if elf["e_type"] != "ET_EXEC":
            raise InvalidELFException(f"Invalid ELF type {elf['e_type']}")

        if elf["e_machine"] != "EM_X86_64":
            raise InvalidELFException(f"Invalid ELF machine {elf['e_machine']}")

        loaded_elf = cls(entry=elf["e_entry"])

        for segment in elf.iter_segments():
            logging.debug(
                "[.] ELF segment: %s V=%#x FSZ=%#x MSZ=%#x flags=%d",
                segment["p_type"],
                segment["p_vaddr"],
                segment["p_filesz"],
                segment["p_memsz"],
                segment["p_flags"],
            )
            if segment["p_type"] == "PT_LOAD":
                if segment["p_filesz"] > segment["p_memsz"]:
                    raise InvalidELFException(
                        "Invalid ELF: corrupt binary. p_filesz > p_memsz"
                    )

                # Try to read segment data
                segment_data = segment.data()
                if len(segment_data) != segment["p_filesz"]:
                    raise ValueError("Invalid ELF: corrupt binary")

                logging.debug("[.] Read %d bytes of code/data", len(segment_data))
                loaded_elf.segments.append(
                    Segment(
                        vaddr=segment["p_vaddr"],
                        size=segment["p_memsz"],
                        data=segment_data,
                        prot=getProtFromELFSegment(segment["p_flags"]),
                    )
                )
            elif segment["p_type"] == "PT_GNU_STACK":
                loaded_elf.stack = gnu_stack_addr
                loaded_elf.segments.append(
                    Segment(
                        vaddr=gnu_stack_addr - gnu_stack_size,
                        size=gnu_stack_size,
                        prot=Prot.RW,
                        stack=True,
                    )
                )
            elif segment["p_type"] == "PT_NULL":
                pass
            else:
                raise InvalidELFException(
                    f"Invalid ELF: unhandled segment {segment['p_type']}"
                )

        loaded_elf.line_map = ELF._extract_linemap(elf)

        return loaded_elf

    def add_trampoline(self):
        """
        Adds a trampoline to the ELF binary.

        The trampoline is a call instruction to the entry point followed by an terminate instruction.
        If reaching the terminator indicates that the program has gracefully returned from the call.

        Parameters:
            elf (ELF): The ELF object to which the trampoline will be added.
        Raises:
            RaspsimException: If the entry point is not page-aligned or if the entry
                              point is less than 0x3000.

        Modifies:
            - Adjusts the entry point of the ELF object.
            - Appends a new segment containing the trampoline code to the ELF object's segments.
        """
        code_segments = [s for s in self.segments if s.prot in [Prot.RX, Prot.RWX]]
        if not all(s.vaddr & 0xFFF == 0 for s in code_segments):
            raise RaspsimException("code segments must be page-aligned")

        if not all(s.vaddr >= 0x3000 for s in code_segments):
            raise RaspsimException("code segments must be above address 0x3000")

        trampoline_bias = 0x2000
        addr_start = max(0x1000, self.entry - trampoline_bias)
        reloff = self.entry - addr_start - 5
        trampoline = b"\xe8" + reloff.to_bytes(4, "little") + b"\xcd\x80"
        self.segments.append(
            Segment(
                vaddr=addr_start, prot=Prot.RX, size=len(trampoline), data=trampoline
            )
        )
        self.entry = addr_start
