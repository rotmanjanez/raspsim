# x86sim

`x86sim` is a cycle-accurate x86-64 simulator, originally developed by Alexis Engelke as a fork of PTLsim by Matt Yourst. This package provides Python bindings to the `x86sim` C++ library, allowing you to configure the virtual address space, set initial register values, and run simulations of ELF binaries directly from Python.

## Installation

The extension is built from the in-tree C++ core via CMake (driven by
[scikit-build-core](https://scikit-build-core.readthedocs.io/)). From the
repository root:

```bash
pip install .
```

A C++26-capable compiler and CMake >= 3.25 are required.

## Features

- **Cycle-accurate simulation** of x86-64 instructions.
- **ELF Binary Loading**: Parse and load ELF binaries, including segments and debugging information.
- **Custom Memory Management**: Map memory segments with specific permissions.
- **Register Access**: Read and write CPU registers during simulation.
- **Exception Handling**: Detailed exceptions for various fault conditions.
- **Python Binding**: Interface directly with the simulator using Python.

## Key Components

### `Machine` Class

The main class used to interact with the simulator (`x86sim.Machine`, a thin
wrapper over the compiled `x86sim.bindings.Machine`). It provides methods to:

- Load ELF binaries into the simulator.
- Run the simulation.
- Access and modify CPU registers.
- Handle exceptions and retrieve simulation statistics.

### `ELF` Class

Represents an ELF binary and provides methods to:

- Parse ELF files from bytes and extract segments.
- Handle debugging information like line mappings.
- Add trampolines to binaries for controlled execution flow.

### `Segment` Class

Represents a memory segment with attributes:

- `vaddr`: Virtual address.
- `prot`: Memory protection (read, write, execute).
- `size`: Size of the segment.
- `data`: Optional data contained in the segment.
- `stack`: Indicates if the segment is a stack.

### `iN` Classes (`i8`, `i16`, `i32`, `i64`)

Represents N-bit integers with proper wrapping and arithmetic operations. Useful for working with register values and ensuring correct bit-width behavior.

### Utility Functions

- `simcompile`: Context manager to compile code into an ELF binary using appropriate flags.
- `asm_preamble`: Generates assembly code preamble with an entry label.
- `asm_stop_sim`: Generates assembly instruction to stop the simulator (`int 0x80`).

## Usage

### Compiling and Loading an ELF Binary

```python
from x86sim.simcompile import simcompile
from x86sim.elf import ELF
from x86sim import Machine

# Write your assembly code
code = """
.global _start
.intel_syntax noprefix
.text
_start:
    mov rax, -1
    int 0x80
"""

# Compile the code into an ELF binary
with simcompile(code=code) as f:
    elf = ELF.from_file(f)

# Create a simulator instance and load the ELF binary
sim = Machine()
sim.load_elf(elf)

# Run the simulation
sim.run()

# Access register values
print(f"RAX: {sim.registers.rax:#x}")
```

### Accessing Memory
```python
# Map a memory segment
address = sim.memmap(start=0x1000, prot=Prot.RW, length=0x1000)

# Write to memory
address[0x0] = b'\x01\x02\x03\x04'

# Read from memory
data = address.read(size=4)
print(f"Data: {data.hex()}")
```

### glibc syscalls and Stream I/O

`Machine` is a register/memory simulator by default: an `int 0x80` stops it and
every other syscall raises. Several opt-in features extend it, all off by
default:

- `glibc=True` enables the portable Linux glibc-startup syscalls: the malloc/free
  heap (`brk` and anonymous `mmap`/`munmap`/`mremap`) plus `arch_prctl`,
  `set_tid_address`, `set_robust_list`, `rseq`, `prlimit64`, `uname`,
  `getpid`/`getuid`/…, `futex` and the synthetic signal syscalls. With
  `glibc=False` those syscalls raise.
- `stdin`/`stdout`/`stderr` route the guest's `read(0, ...)` / `write(1, ...)` /
  `write(2, ...)` to Python file-like objects (anything exposing `.read(n)` /
  `.write(bytes)`, such as `io.BytesIO`). No host file descriptors are ever
  touched — Python fully controls the mapping. An fd left as `None` is
  unconfigured and a guest I/O on it raises.
- `readlink=callable` resolves the guest's `readlink`/`readlinkat` calls via a
  Python `readlink(path: str) -> str | None` (e.g. for `/proc/self/exe`, which
  glibc reads at startup). Returning `None` yields `-ENOENT`; without a callback
  the syscall returns `-ENOSYS`. No real filesystem is touched.
- `core="seq"` selects the sequential CPU model instead of the default
  out-of-order one (`core="ooo"`). Both cores execute unaligned memory accesses
  correctly and can run glibc; the sequential core is simpler and slower, the
  out-of-order core is the cycle-accurate default.

These features require the guest to use the 64-bit `syscall` instruction (the
Linux syscall ABI: number in `rax`, arguments in `rdi`, `rsi`, `rdx`, ...).
`int 0x80` remains the guest-exit sentinel.

#### Capturing guest output

```python
import io
from x86sim import Machine

out = io.BytesIO()
sim = Machine(stdout=out)
# ... load a guest that does write(1, msg, len) then exit_group(0) via syscall ...
sim.run()
print(out.getvalue())  # the bytes the guest wrote; the real terminal stays silent
```

#### Feeding guest input

```python
import io
from x86sim import Machine

sim = Machine(stdin=io.BytesIO(b"hello"))
# ... load a guest that does read(0, buf, n) via syscall ...
sim.run()
# the guest buffer now contains b"hello"
```

#### Enabling the glibc-startup syscalls

```python
from x86sim import Machine

sim = Machine(glibc=True)
# ... load a guest that calls brk / anonymous mmap (e.g. via malloc) ...
sim.run()  # glibc syscalls succeed; with Machine() (glibc=False) they would raise
```

#### Running a real glibc binary

A static, non-PIE glibc executable can run end to end once you build its initial
stack (argc/argv/envp and the auxiliary vector — see
`tests/python/elfload.py` for a reference loader). It runs on either core:

```python
import io
from x86sim import Machine

out = io.BytesIO()
sim = Machine(glibc=True, stdout=out, readlink=lambda p: "/hello")
# ... map the ELF's PT_LOAD segments and lay out the argc/argv/envp/auxv stack ...
sim.run()  # __libc_start_main -> main -> write(1, ...) -> exit_group
print(out.getvalue())
```

### Working with Registers
```python
# Set register values
sim.registers.rax = 0xdeadbeef

# Get register values
rip = sim.registers["rip"]
print(f"RIP: {rip:#x}")
```

### Handling Exceptions

The simulator provides detailed exceptions for various fault conditions, such as PageFaultException, InvalidOpcodeException, etc. These can be caught and handled appropriately. In case that the address space of the simulator was populated with an ELF binary that contains debugging information, the exception will contain the `file` and `lineno`and `line` attributes that can be used to provide more context about the fault. Either all three attributes are present or none of them.
```python
try:
    sim.run()
except RaspsimException as e:
    print(f"Simulation error: {e}")
    if hasattr(e, 'lineno'):
        print(f"At {e.file}:{e.lineno}: {e.line}")
```

### Command-Line Interface

You can use x86sim as a command-line tool to simulate code directly from the terminal.

```bash
python -m x86sim << EOF
.global _start
.intel_syntax noprefix
.text
_start:
    mov rax, 1
    int 0x80
EOF
```

## Contributing
Contributions are welcome! Please submit pull requests or open issues for any bugs or feature requests at the [GitHub repository](https://github.com/rotmanjanez/x86sim).