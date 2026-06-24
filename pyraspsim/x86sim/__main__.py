import sys
import argparse

from .utils import rscompile
from .raspsim import Raspsim
from .core import RaspsimException
from .elf import ELF


def main() -> None:
    parser = argparse.ArgumentParser(description="Raspsim simulator")
    parser.add_argument("file", type=argparse.FileType("r"), help="Source file", nargs="?")
    parser.add_argument("--debug", action="store_true", help="Enable debug mode")
    parser.add_argument(
        "-x",
        "--lang",
        choices=["c", "assembler"],
        default="assembler",
        help="Language of the source file",
    )
    parser.add_argument(
        "--entry-label",
        default="_start",
        help="Entry point label",
    )
    parser.add_argument(
        "--compiler",
        default="cc",
        help="Compiler to use",
    )
    parser.add_argument(
        "--std",
        default="c17",
        help="Standard to use",
    )

    args = parser.parse_args()

    if args.file:
        program = args.file.read()
    else:
        program = "\n".join(sys.stdin.readlines())

    with rscompile(code=program, entry_label=args.entry_label, lang=args.lang, compiler=args.compiler, std=args.std, debug_symbols=True) as f:
        elf = ELF.from_file(f)

    # Create a new Raspsim instance
    sim = Raspsim()
    sim.load_elf(elf)

    # Run the simulation
    try:
        sim.run()

        # Print the simulation state
        print(sim)
    except RaspsimException as e:
        if hasattr(e, "file") and hasattr(e, "lineno") and hasattr(e, "line"):
            print(f"Error in {e.file}:{e.lineno}: {e}", file=sys.stderr)
            raise e



if __name__ == "__main__":
    main()
