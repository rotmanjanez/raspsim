"""Pathlib-flavored access to a Machine's in-memory guest filesystem (memfs).

``MemPath`` wraps the compiled ``bindings.Fs`` handle so host code can inspect
and prepopulate the sandbox the guest sees, with the ergonomics of
``pathlib.Path``::

    m = x86sim.Machine(memfs=True)
    (m.fs / "data").mkdir()
    (m.fs / "data" / "in.txt").write_text("hello")
    m.run()
    print((m.fs / "data" / "out.txt").read_text())
    print([p.name for p in (m.fs / "data").iterdir()])

Note that Python's built-in ``open()`` always dispatches to the host OS and
cannot reach a virtual filesystem; use ``read_bytes``/``write_text`` and
friends instead.
"""

from __future__ import annotations

from pathlib import PurePosixPath
from typing import Iterator, NamedTuple, Union

from . import bindings

__all__ = ["MemPath", "MemStat"]


class MemStat(NamedTuple):
    """Stat result of a memfs node.

    Timestamps are deterministic mutation ticks (the memfs never reads the
    host clock, keeping simulations reproducible), not wall-clock seconds.
    """

    st_mode: int
    st_ino: int
    st_nlink: int
    st_size: int
    st_atime: int
    st_mtime: int
    st_ctime: int


class MemPath:
    """A path inside an in-memory guest filesystem.

    Deliberately *not* ``os.PathLike``: passing it to ``open()`` or ``os.*``
    would silently hit the host filesystem instead of the sandbox.
    """

    __slots__ = ("_fs", "_p")

    def __init__(self, fs: "bindings.Fs", path: Union[str, PurePosixPath] = "/") -> None:
        self._fs = fs
        p = PurePosixPath(path)
        if not p.is_absolute():
            p = PurePosixPath("/") / p
        self._p = p

    # -- path arithmetic ---------------------------------------------------

    def __truediv__(self, part: Union[str, PurePosixPath, "MemPath"]) -> "MemPath":
        if isinstance(part, MemPath):
            part = part._p
        return MemPath(self._fs, self._p / part)

    def joinpath(self, *parts: Union[str, PurePosixPath]) -> "MemPath":
        return MemPath(self._fs, self._p.joinpath(*parts))

    @property
    def name(self) -> str:
        return self._p.name

    @property
    def parent(self) -> "MemPath":
        return MemPath(self._fs, self._p.parent)

    @property
    def parts(self) -> tuple[str, ...]:
        return self._p.parts

    @property
    def filesystem(self) -> "bindings.Fs":
        """The underlying compiled Fs handle."""
        return self._fs

    def __str__(self) -> str:
        return str(self._p)

    def __repr__(self) -> str:
        return f"MemPath({str(self._p)!r})"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, MemPath) and self._p == other._p and self._fs == other._fs

    def __hash__(self) -> int:
        return hash(self._p)

    # -- queries -------------------------------------------------------------

    def exists(self) -> bool:
        return self._fs.exists(str(self._p))

    def is_file(self) -> bool:
        return self._fs.is_file(str(self._p))

    def is_dir(self) -> bool:
        return self._fs.is_dir(str(self._p))

    def is_symlink(self) -> bool:
        return self._fs.is_symlink(str(self._p))

    def stat(self, *, follow_symlinks: bool = True) -> MemStat:
        raw = self._fs.stat(str(self._p), follow_symlinks)
        return MemStat(
            st_mode=raw["st_mode"],
            st_ino=raw["st_ino"],
            st_nlink=raw["st_nlink"],
            st_size=raw["st_size"],
            st_atime=raw["st_atime"],
            st_mtime=raw["st_mtime"],
            st_ctime=raw["st_ctime"],
        )

    def iterdir(self) -> Iterator["MemPath"]:
        for name in self._fs.listdir(str(self._p)):
            yield self / name

    # -- content -------------------------------------------------------------

    def read_bytes(self) -> bytes:
        return self._fs.read_bytes(str(self._p))

    def read_text(self, encoding: str = "utf-8", errors: str = "strict") -> str:
        return self.read_bytes().decode(encoding, errors)

    def write_bytes(self, data: bytes) -> int:
        self._fs.write_bytes(str(self._p), bytes(data))
        return len(data)

    def write_text(self, data: str, encoding: str = "utf-8", errors: str = "strict") -> int:
        return self.write_bytes(data.encode(encoding, errors))

    def truncate(self, size: int = 0) -> None:
        self._fs.truncate(str(self._p), size)

    # -- mutations -----------------------------------------------------------

    def mkdir(self, *, parents: bool = False, exist_ok: bool = False) -> None:
        try:
            self._fs.mkdir(str(self._p))
        except FileNotFoundError:
            if not parents or self._p.parent == self._p:
                raise
            self.parent.mkdir(parents=True, exist_ok=True)
            self._fs.mkdir(str(self._p))
        except FileExistsError:
            if not exist_ok or not self.is_dir():
                raise

    def unlink(self, *, missing_ok: bool = False) -> None:
        try:
            self._fs.unlink(str(self._p))
        except FileNotFoundError:
            if not missing_ok:
                raise

    def rmdir(self) -> None:
        self._fs.rmdir(str(self._p))

    def rename(self, target: Union[str, PurePosixPath, "MemPath"]) -> "MemPath":
        if isinstance(target, MemPath):
            target_path = str(target._p)
        else:
            target_path = str(target)
        self._fs.rename(str(self._p), target_path)
        return MemPath(self._fs, target_path)

    def symlink_to(self, target: Union[str, PurePosixPath, "MemPath"]) -> None:
        if isinstance(target, MemPath):
            target = str(target._p)
        self._fs.symlink(str(target), str(self._p))

    def readlink(self) -> PurePosixPath:
        return PurePosixPath(self._fs.readlink(str(self._p)))
