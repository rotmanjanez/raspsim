try:
    from .raspsim import *
    from .core import *
    from .elf import *
except ImportError as e:
    raise ImportError("""Failed to import core module. This is most likely due to an unsupoorted python version, plattform or architecture.
    Supported versions are 3.7 to 3.15 on linux x86_64""") from e
