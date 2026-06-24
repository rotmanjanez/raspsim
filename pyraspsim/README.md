# RASPsim

`RASPsim` is a cycle-accurate x86-64 simulator, originally developed by Alexis Engelke as a fork of PTLsim by Matt Yourst. This package provides Python bindings to `RASPsim`, allowing you to configure the virtual address space, set initial register values, and run simulations of ELF binaries directly from Python.

## Features

- **Cycle-accurate simulation** of x86-64 instructions.
- **ELF Binary Loading**: Parse and load ELF binaries, including segments and debugging information.
- **Custom Memory Management**: Map memory segments with specific permissions.
- **Register Access**: Read and write CPU registers during simulation.
- **Exception Handling**: Detailed exceptions for various fault conditions.
- **Python Binding**: Interface directly with the simulator using Python.

## Key Components

### `Raspsim` Class

The main class used to interact with the simulator. It provides methods to:

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

- `rscompile`: Context manager to compile code into an ELF binary using appropriate flags.
- `asm_preamble`: Generates assembly code preamble with an entry label.
- `asm_stop_sim`: Generates assembly instruction to stop the simulator (`int 0x80`).

## Usage

### Compiling and Loading an ELF Binary

```python
from raspsim.utils import rscompile, asm_preamble, asm_stop_sim
from raspsim.elf import ELF
from raspsim.raspsim import Raspsim

# Write your assembly code
code = asm_preamble() + """
    mov rax, -1
""" + asm_stop_sim()

# Compile the code into an ELF binary
with rscompile(code=code) as f:
    elf = ELF.from_file(f)

# Create a simulator instance and load the ELF binary
sim = Raspsim()
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

You can use raspsim as a command-line tool to simulate code directly from the terminal.

```bash
python -m raspsim << EOF
.global _start
.intel_syntax noprefix
.text
_start:
    mov rax, 1
    int 0x80
EOF
```

## Contributing
Contributions are welcome! Please submit pull requests or open issues for any bugs or feature requests at the [GitHub repository](https://github.com/Joshy-R/raspsim).