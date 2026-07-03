"""Host-side tests for the in-memory guest filesystem (memfs).

Everything here runs without executing a single guest instruction: it
exercises the MemPath/Fs surface the host uses to prepopulate and inspect the
sandbox.
"""

import pytest

import x86sim


@pytest.fixture
def machine() -> x86sim.Machine:
    return x86sim.Machine(memfs=True)


def test_memfs_is_opt_in() -> None:
    plain = x86sim.Machine()
    with pytest.raises(ValueError):
        plain.fs  # noqa: B018


def test_memfs_accepts_bool_or_fs_only() -> None:
    with pytest.raises(ValueError):
        x86sim.Machine(memfs="yes")  # type: ignore[arg-type]


def test_memfs_and_readlink_are_mutually_exclusive() -> None:
    with pytest.raises(ValueError):
        x86sim.Machine(memfs=True, readlink=lambda path: None)


def test_write_read_roundtrip(machine: x86sim.Machine) -> None:
    path = machine.fs / "hello.txt"
    assert path.write_text("hello world") == 11
    assert path.read_text() == "hello world"
    assert path.read_bytes() == b"hello world"
    assert path.exists() and path.is_file() and not path.is_dir()


def test_path_arithmetic(machine: x86sim.Machine) -> None:
    p = machine.fs / "a" / "b" / "c.txt"
    assert str(p) == "/a/b/c.txt"
    assert p.name == "c.txt"
    assert str(p.parent) == "/a/b"
    assert p.parts == ("/", "a", "b", "c.txt")
    assert p == machine.fs / "a" / "b" / "c.txt"
    assert p != machine.fs / "a"
    assert len({p, machine.fs / "a" / "b" / "c.txt"}) == 1


def test_iterdir_is_sorted(machine: x86sim.Machine) -> None:
    d = machine.fs / "dir"
    d.mkdir()
    for name in ("zeta", "alpha", "mid"):
        (d / name).write_bytes(b"")
    assert [p.name for p in d.iterdir()] == ["alpha", "mid", "zeta"]


def test_mkdir_parents_and_exist_ok(machine: x86sim.Machine) -> None:
    nested = machine.fs / "x" / "y" / "z"
    with pytest.raises(FileNotFoundError):
        nested.mkdir()
    nested.mkdir(parents=True)
    assert nested.is_dir()
    with pytest.raises(FileExistsError):
        nested.mkdir()
    nested.mkdir(exist_ok=True)


def test_oserror_subclasses(machine: x86sim.Machine) -> None:
    root = machine.fs
    (root / "d").mkdir()
    (root / "d" / "f").write_bytes(b"x")

    with pytest.raises(FileNotFoundError) as excinfo:
        (root / "missing").read_bytes()
    assert excinfo.value.errno == 2
    assert excinfo.value.filename == "/missing"

    with pytest.raises(IsADirectoryError):
        (root / "d").read_bytes()
    with pytest.raises(NotADirectoryError):
        list((root / "d" / "f").iterdir())
    with pytest.raises(OSError) as excinfo:
        (root / "d").rmdir()
    assert excinfo.value.errno == 39  # ENOTEMPTY


def test_unlink_and_rmdir(machine: x86sim.Machine) -> None:
    root = machine.fs
    (root / "f").write_bytes(b"x")
    (root / "f").unlink()
    assert not (root / "f").exists()
    with pytest.raises(FileNotFoundError):
        (root / "f").unlink()
    (root / "f").unlink(missing_ok=True)

    (root / "d").mkdir()
    (root / "d").rmdir()
    assert not (root / "d").exists()


def test_rename(machine: x86sim.Machine) -> None:
    root = machine.fs
    (root / "src").write_bytes(b"payload")
    target = (root / "src").rename(root / "dst")
    assert isinstance(target, x86sim.MemPath)
    assert target.read_bytes() == b"payload"
    assert not (root / "src").exists()


def test_symlinks(machine: x86sim.Machine) -> None:
    root = machine.fs
    (root / "target.txt").write_text("via link")
    (root / "link").symlink_to("/target.txt")
    assert (root / "link").is_symlink()
    assert (root / "link").read_text() == "via link"
    assert str((root / "link").readlink()) == "/target.txt"
    # stat follows by default; lstat via follow_symlinks=False
    assert (root / "link").stat().st_size == 8
    assert (root / "link").stat(follow_symlinks=False).st_size == len("/target.txt")


def test_stat_fields(machine: x86sim.Machine) -> None:
    f = machine.fs / "f"
    f.write_bytes(b"12345")
    st = f.stat()
    assert st.st_size == 5
    assert st.st_nlink == 1
    assert st.st_mode & 0o170000 == 0o100000  # S_IFREG
    assert st.st_ino > 0

    d = machine.fs / "d"
    d.mkdir()
    assert d.stat().st_mode & 0o170000 == 0o040000  # S_IFDIR
    assert d.stat().st_nlink == 2


def test_deterministic_timestamps(machine: x86sim.Machine) -> None:
    f = machine.fs / "f"
    f.write_bytes(b"a")
    first = f.stat().st_mtime
    f.write_bytes(b"b")
    assert f.stat().st_mtime > first


def test_truncate(machine: x86sim.Machine) -> None:
    f = machine.fs / "f"
    f.write_bytes(b"abcdef")
    f.truncate(3)
    assert f.read_bytes() == b"abc"
    f.truncate(6)
    assert f.read_bytes() == b"abc\0\0\0"


def test_machines_do_not_share_state() -> None:
    a = x86sim.Machine(memfs=True)
    b = x86sim.Machine(memfs=True)
    (a.fs / "only-a").write_bytes(b"x")
    assert not (b.fs / "only-a").exists()


def test_shared_fs_prepopulation() -> None:
    fs = x86sim.Fs()
    fs.mkdir("/etc")
    fs.write_bytes("/etc/hosts", b"127.0.0.1 localhost\n")

    m = x86sim.Machine(memfs=fs)
    assert (m.fs / "etc" / "hosts").read_bytes() == b"127.0.0.1 localhost\n"
    # Same filesystem, both directions.
    (m.fs / "added-by-machine").write_bytes(b"y")
    assert fs.read_bytes("/added-by-machine") == b"y"
    assert m.fs.filesystem == fs

    # A MemPath is accepted too and shares the same tree.
    m2 = x86sim.Machine(memfs=m.fs)
    assert (m2.fs / "etc" / "hosts").exists()


def test_write_requires_existing_parent(machine: x86sim.Machine) -> None:
    with pytest.raises(FileNotFoundError):
        (machine.fs / "no" / "such" / "dir.txt").write_bytes(b"x")


def test_map_fd_requires_memfs() -> None:
    plain = x86sim.Machine()
    with pytest.raises(ValueError):
        plain.map_fd(1, "/log", "w")


def test_map_fd_modes(machine: x86sim.Machine) -> None:
    (machine.fs / "existing").write_bytes(b"data")
    machine.map_fd(3, machine.fs / "existing", "r")
    machine.map_fd(4, machine.fs / "created", "w")
    assert (machine.fs / "created").exists()
    with pytest.raises(FileNotFoundError):
        machine.map_fd(5, machine.fs / "missing", "r")
    with pytest.raises(ValueError):
        machine.map_fd(6, machine.fs / "existing", "q")
