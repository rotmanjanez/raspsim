"""
This module provides utilities for compiling code to ELF binaries, handling ELF files, and populating a simulator with ELF binary data.

Classes:
    ELF: Represents ELF binary data.
    Segment: Represents a memory segment.
    InvalidELFException: Exception raised for invalid ELF files.

Functions:
    rscompile: Context manager to compile code to an ELF binary.
    elf_add_trampoline: Adds a trampoline to the ELF binary.
    asm_preamble: Generates an assembly preamble.
    asm_stop_sim: Generates an assembly instruction to stop the simulator.
    populate_from_elf: Populates the simulator with ELF binary data.

Most of the code is taken from work done by Alexis Engelke for the rasptest project.
"""

import subprocess
import tempfile
from typing import Literal, IO, Generator, cast, Any
from contextlib import contextmanager
from pathlib import Path
from inspect import cleandoc


@contextmanager
def rscompile(
    *,
    code: str | None = None,
    file: str | Path | None = None,
    entry_label: str = "_start",
    lang: Literal["c", "assembler"] = "assembler",
    compiler: str = "cc",
    std: str = "c17",
    text_segment_addr: int = 0x4000,
    additional_flags: list[str] | None = None,
    debug_symbols: bool = False,
) -> Generator[IO[bytes], None, None]:
    """
    Compile Code to an ELF binary.

    Args:
        code: Code to compile.
        file: File containing the code to compile.
        entry_label: Entry point label.
        lang: Language of the code.
        compiler: Compiler to use.
        std: Standard to use.
        text_segment_addr: Address where the text segment should be placed.
        additional_flags: Additional flags to pass to the compiler.
        debug_symbols: Whether to include debug symbols in the binary.

    Yields:
        A file-like object containing the compiled ELF binary.
    """
    if code is None and file is None:
        raise ValueError("Either code or file must be provided")

    if code is not None and file is not None:
        raise ValueError("Only one of code or file must be provided")

    cmd = [
        compiler,
        f"-x{lang}",
        "-pipe",
        f"-std={std}",
        "-static",
        "-march=k8",
        "-mno-mmx",
        "-mno-80387",
        "-mtune=generic",  # prevent excessive loop unrolling
        "-ffreestanding",
        "-nostdlib",
        "-D_MM_MALLOC_H_INCLUDED",
        "-fcf-protection=none",
        # Assembler arguments
        "-Wa,-march=k8+cmov+nommx",
        "-Wa,-mx86-used-note=no",
        "-Wa,--fatal-warnings",
        # Linker arguments
        "-Wl,-n",
        "-Wl,--fatal-warnings",  # fail if entry point doesn't exist
        "-Wl,--no-dynamic-linker",
        "-Wl,--build-id=none",
        "-Wl,-z,defs",
        "-Wl,-z,noexecstack",
        "-Wl,-z,norelro",
        "-Wl,-z,noseparate-code",
        f"-Wl,-Ttext={text_segment_addr:#x}",
        f"-Wl,-e,{entry_label}",
    ]
    if debug_symbols:
        cmd += ["-g", "-gdwarf-4"]
    cmd += additional_flags or []

    with tempfile.NamedTemporaryFile("w+b", delete=False) as f:
        cmd += ["-o", f.name]
        kwargs: dict[str, Any] = {"text": True, "capture_output": True, "check": True}

        if code is not None:
            cmd += ["-"]
            kwargs["input"] = code
        else:
            if not Path(cast(str | Path, file)).exists():
                raise FileNotFoundError(f"File not found: {file}")
            cmd += [str(file)]

        subprocess.run(args=cmd, **kwargs)

        f.seek(0)
        yield f


def asm_preamble(label: str = "_start") -> str:
    return cleandoc(
        f"""
        .intel_syntax noprefix
        .global {label}
        .section .text
        {label}:
        """
    )


def asm_stop_sim() -> str:
    return "int 0x80\n"
