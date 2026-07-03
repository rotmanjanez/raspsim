"""Guest-side tests for the in-memory filesystem (memfs).

Each test compiles a small freestanding assembly guest with simcompile() and
runs it against a Machine with memfs enabled. Guests communicate results
either through the memfs itself (the host inspects files afterwards) or by
writing to fd 1, which is routed to an io.BytesIO unless a test rebinds it.
"""

import io
import struct

import x86sim
from x86sim.simcompile import simcompile


def make_machine(code: str, **kwargs) -> tuple[x86sim.Machine, io.BytesIO]:
    out = io.BytesIO()
    kwargs.setdefault("memfs", True)
    kwargs.setdefault("stdout", out)
    machine = x86sim.Machine(**kwargs)
    with simcompile(code=code) as elf_file:
        machine.load_elf(x86sim.ELF.from_file(elf_file))
    return machine, out


EXIT = """
    xor %rdi, %rdi
    mov $60, %rax
    syscall
"""


def test_open_write_close_roundtrip() -> None:
    machine, out = make_machine(
        r"""
.globl _start
_start:
    sub $64, %rsp

    # fd = open("/out.txt", O_WRONLY|O_CREAT, 0644)
    lea out_path(%rip), %rdi
    mov $65, %rsi
    mov $420, %rdx
    mov $2, %rax
    syscall
    mov %rax, %r12

    # write(fd, "hello", 5); write(1, "stdout\n", 7); close(fd)
    mov %r12, %rdi
    lea msg(%rip), %rsi
    mov $5, %rdx
    mov $1, %rax
    syscall
    mov $1, %rdi
    lea stream_msg(%rip), %rsi
    mov $7, %rdx
    mov $1, %rax
    syscall
    mov %r12, %rdi
    mov $3, %rax
    syscall

    # cat /in.txt to fd 1
    lea in_path(%rip), %rdi
    xor %rsi, %rsi
    mov $2, %rax
    syscall
    mov %rax, %rdi
    mov %rsp, %rsi
    mov $32, %rdx
    xor %rax, %rax
    syscall
    mov %rax, %rdx
    mov $1, %rdi
    mov %rsp, %rsi
    mov $1, %rax
    syscall
"""
        + EXIT
        + r"""
out_path:   .asciz "/out.txt"
in_path:    .asciz "/in.txt"
msg:        .ascii "hello"
stream_msg: .ascii "stdout\n"
"""
    )
    (machine.fs / "in.txt").write_text("from-host")
    machine.run()

    # Guest-created file visible to the host; fd 1 stayed a Python stream.
    assert (machine.fs / "out.txt").read_bytes() == b"hello"
    assert out.getvalue() == b"stdout\nfrom-host"
    # First guest open() picked the lowest non-reserved fd (0/1/2 are streams).
    assert (machine.fs / "out.txt").stat().st_nlink == 1


def test_dup2_shadows_python_stream() -> None:
    machine, out = make_machine(
        r"""
.globl _start
_start:
    # fd = open("/captured.txt", O_WRONLY|O_CREAT, 0644); dup2(fd, 1)
    lea path(%rip), %rdi
    mov $65, %rsi
    mov $420, %rdx
    mov $2, %rax
    syscall
    mov %rax, %rdi
    mov $1, %rsi
    mov $33, %rax
    syscall

    # write(1, ...) now lands in the memfs file, not the Python stream
    mov $1, %rdi
    lea msg(%rip), %rsi
    mov $8, %rdx
    mov $1, %rax
    syscall
"""
        + EXIT
        + r"""
path: .asciz "/captured.txt"
msg:  .ascii "captured"
"""
    )
    machine.run()
    assert (machine.fs / "captured.txt").read_bytes() == b"captured"
    assert out.getvalue() == b""


def test_map_fd_binds_stdout_to_memfs() -> None:
    code = (
        r"""
.globl _start
_start:
    mov $1, %rdi
    lea msg(%rip), %rsi
    mov $6, %rdx
    mov $1, %rax
    syscall
"""
        + EXIT
        + r"""
msg: .ascii "logged"
"""
    )
    # No Python stdout stream at all: fd 1 is a memfs file, set up by the host.
    machine = x86sim.Machine(memfs=True)
    with simcompile(code=code) as elf_file:
        machine.load_elf(x86sim.ELF.from_file(elf_file))
    machine.map_fd(1, machine.fs / "stdout.log", "w")
    machine.run()
    assert (machine.fs / "stdout.log").read_bytes() == b"logged"


def test_chdir_getcwd() -> None:
    machine, out = make_machine(
        r"""
.globl _start
_start:
    sub $128, %rsp

    # chdir("/data/sub")
    lea dir_path(%rip), %rdi
    mov $80, %rax
    syscall

    # getcwd(buf, 128) returns strlen+1; write(1, buf, rax-1)
    mov %rsp, %rdi
    mov $128, %rsi
    mov $79, %rax
    syscall
    lea -1(%rax), %rdx
    mov $1, %rdi
    mov %rsp, %rsi
    mov $1, %rax
    syscall

    # open a file relative to the cwd and write to it
    lea rel_path(%rip), %rdi
    mov $65, %rsi
    mov $420, %rdx
    mov $2, %rax
    syscall
    mov %rax, %rdi
    lea msg(%rip), %rsi
    mov $3, %rdx
    mov $1, %rax
    syscall
"""
        + EXIT
        + r"""
dir_path: .asciz "/data/sub"
rel_path: .asciz "rel.txt"
msg:      .ascii "rel"
"""
    )
    (machine.fs / "data" / "sub").mkdir(parents=True)
    machine.run()
    assert out.getvalue() == b"/data/sub"
    assert (machine.fs / "data" / "sub" / "rel.txt").read_bytes() == b"rel"


def test_getdents64() -> None:
    machine, out = make_machine(
        r"""
.globl _start
_start:
    sub $512, %rsp

    # fd = open("/d", O_RDONLY|O_DIRECTORY)
    lea path(%rip), %rdi
    mov $65536, %rsi
    xor %rdx, %rdx
    mov $2, %rax
    syscall

    # n = getdents64(fd, buf, 512); write(1, buf, n)
    mov %rax, %rdi
    mov %rsp, %rsi
    mov $512, %rdx
    mov $217, %rax
    syscall
    mov %rax, %rdx
    mov $1, %rdi
    mov %rsp, %rsi
    mov $1, %rax
    syscall
"""
        + EXIT
        + r"""
path: .asciz "/d"
"""
    )
    d = machine.fs / "d"
    d.mkdir()
    (d / "zeta").write_bytes(b"")
    (d / "alpha").mkdir()
    (d / "link").symlink_to("/d/zeta")
    machine.run()

    # Parse the raw linux_dirent64 records the guest streamed to fd 1.
    raw = out.getvalue()
    entries: dict[str, int] = {}
    offset = 0
    while offset < len(raw):
        ino, _, reclen, dtype = struct.unpack_from("<QQHB", raw, offset)
        name = raw[offset + 19 : offset + reclen].split(b"\0", 1)[0].decode()
        entries[name] = dtype
        assert ino > 0
        offset += reclen
    # "." and ".." are synthesized first; entries come in map (sorted) order.
    assert list(entries) == [".", "..", "alpha", "link", "zeta"]
    assert entries["alpha"] == 4  # DT_DIR
    assert entries["zeta"] == 8   # DT_REG
    assert entries["link"] == 10  # DT_LNK


def test_lseek_and_pread() -> None:
    machine, out = make_machine(
        r"""
.globl _start
_start:
    sub $64, %rsp

    # fd = open("/in.txt", O_RDONLY)
    lea path(%rip), %rdi
    xor %rsi, %rsi
    mov $2, %rax
    syscall
    mov %rax, %r12

    # lseek(fd, 3, SEEK_SET); read 3 -> "def"
    mov %r12, %rdi
    mov $3, %rsi
    xor %rdx, %rdx
    mov $8, %rax
    syscall
    mov %r12, %rdi
    mov %rsp, %rsi
    mov $3, %rdx
    xor %rax, %rax
    syscall
    mov $1, %rdi
    mov %rsp, %rsi
    mov $3, %rdx
    mov $1, %rax
    syscall

    # pread64(fd, buf, 2, 1) -> "bc" (offset unaffected by lseek)
    mov %r12, %rdi
    mov %rsp, %rsi
    mov $2, %rdx
    mov $1, %r10
    mov $17, %rax
    syscall
    mov $1, %rdi
    mov %rsp, %rsi
    mov $2, %rdx
    mov $1, %rax
    syscall
"""
        + EXIT
        + r"""
path: .asciz "/in.txt"
"""
    )
    (machine.fs / "in.txt").write_text("abcdef")
    machine.run()
    assert out.getvalue() == b"defbc"


def test_stat_reports_size() -> None:
    machine, out = make_machine(
        r"""
.globl _start
_start:
    sub $256, %rsp

    # stat("/in.txt", buf); write(1, buf+48, 8) -> little-endian st_size
    lea path(%rip), %rdi
    mov %rsp, %rsi
    mov $4, %rax
    syscall
    mov $1, %rdi
    lea 48(%rsp), %rsi
    mov $8, %rdx
    mov $1, %rax
    syscall
"""
        + EXIT
        + r"""
path: .asciz "/in.txt"
"""
    )
    (machine.fs / "in.txt").write_bytes(b"x" * 1234)
    machine.run()
    assert struct.unpack("<Q", out.getvalue())[0] == 1234


def test_guest_rename_and_unlink_visible_to_host() -> None:
    machine, _ = make_machine(
        r"""
.globl _start
_start:
    # rename("/a.txt", "/b.txt"); unlink("/c.txt"); mkdir("/made", 0755)
    lea a_path(%rip), %rdi
    lea b_path(%rip), %rsi
    mov $82, %rax
    syscall
    lea c_path(%rip), %rdi
    mov $87, %rax
    syscall
    lea d_path(%rip), %rdi
    mov $493, %rsi
    mov $83, %rax
    syscall
"""
        + EXIT
        + r"""
a_path: .asciz "/a.txt"
b_path: .asciz "/b.txt"
c_path: .asciz "/c.txt"
d_path: .asciz "/made"
"""
    )
    (machine.fs / "a.txt").write_bytes(b"payload")
    (machine.fs / "c.txt").write_bytes(b"doomed")
    machine.run()
    assert not (machine.fs / "a.txt").exists()
    assert (machine.fs / "b.txt").read_bytes() == b"payload"
    assert not (machine.fs / "c.txt").exists()
    assert (machine.fs / "made").is_dir()


def test_open_missing_file_returns_enoent() -> None:
    machine, out = make_machine(
        r"""
.globl _start
_start:
    # open("/missing", O_RDONLY) must return -ENOENT (-2)
    lea path(%rip), %rdi
    xor %rsi, %rsi
    mov $2, %rax
    syscall
    cmp $-2, %rax
    jne fail
    mov $1, %rdi
    lea ok_msg(%rip), %rsi
    mov $2, %rdx
    mov $1, %rax
    syscall
fail:
"""
        + EXIT
        + r"""
path:   .asciz "/missing"
ok_msg: .ascii "ok"
"""
    )
    machine.run()
    assert out.getvalue() == b"ok"


def test_file_backed_mmap() -> None:
    machine, out = make_machine(
        r"""
.globl _start
_start:
    # fd = open("/blob", O_RDONLY)
    lea path(%rip), %rdi
    xor %rsi, %rsi
    mov $2, %rax
    syscall

    # addr = mmap(NULL, 4096, PROT_READ, MAP_PRIVATE, fd, 0)
    mov %rax, %r8
    xor %rdi, %rdi
    mov $4096, %rsi
    mov $1, %rdx
    mov $2, %r10
    xor %r9, %r9
    mov $9, %rax
    syscall

    # write(1, addr, 10)
    mov %rax, %rsi
    mov $1, %rdi
    mov $10, %rdx
    mov $1, %rax
    syscall
"""
        + EXIT
        + r"""
path: .asciz "/blob"
""",
        glibc=True,  # mmap itself comes from the portable glibc chain
    )
    (machine.fs / "blob").write_bytes(b"0123456789abcdef")
    machine.run()
    assert out.getvalue() == b"0123456789"


# The syscalls glibc's stdio startup and file paths lean on, exercised
# directly. (A full static-glibc guest is currently blocked by the Python ELF
# loader, which predates memfs and does not handle PT_NOTE/PT_TLS segments.)
def test_glibc_hot_path_syscalls() -> None:
    machine, out = make_machine(
        r"""
.globl _start
_start:
    sub $512, %rsp

    # fstat(1, buf) on the Python stream: synthetic char device (glibc stdio
    # startup does this before its first write)
    mov $1, %rdi
    mov %rsp, %rsi
    mov $5, %rax
    syscall
    test %rax, %rax
    jnz fail
    mov 24(%rsp), %eax        # st_mode
    and $0170000, %eax
    cmp $0020000, %eax        # S_IFCHR
    jne fail

    # ioctl(1, TCGETS, buf) -> -ENOTTY (-25)
    mov $1, %rdi
    mov $0x5401, %rsi
    mov %rsp, %rdx
    mov $16, %rax
    syscall
    cmp $-25, %rax
    jne fail

    # statx(AT_FDCWD, "/f", 0, STATX_BASIC_STATS, buf)
    mov $-100, %rdi
    lea path(%rip), %rsi
    xor %rdx, %rdx
    mov $0x7ff, %r10
    mov %rsp, %r8
    mov $332, %rax
    syscall
    test %rax, %rax
    jnz fail
    mov 40(%rsp), %rax        # stx_size
    cmp $6, %rax
    jne fail

    # access("/f", R_OK)
    lea path(%rip), %rdi
    mov $4, %rsi
    mov $21, %rax
    syscall
    test %rax, %rax
    jnz fail

    # fd = open("/f", O_RDONLY); dup(fd) shares the offset
    lea path(%rip), %rdi
    xor %rsi, %rsi
    mov $2, %rax
    syscall
    mov %rax, %r12
    mov %r12, %rdi
    mov $32, %rax
    syscall
    mov %rax, %r13

    # read 3 through the original, then 3 through the dup: "abc" + "def"
    mov %r12, %rdi
    mov %rsp, %rsi
    mov $3, %rdx
    xor %rax, %rax
    syscall
    mov %r13, %rdi
    lea 3(%rsp), %rsi
    mov $3, %rdx
    xor %rax, %rax
    syscall
    mov $1, %rdi
    mov %rsp, %rsi
    mov $6, %rdx
    mov $1, %rax
    syscall

    # fcntl(fd, F_GETFL) -> O_RDONLY (0)
    mov %r12, %rdi
    mov $3, %rsi
    mov $72, %rax
    syscall
    test %rax, %rax
    jnz fail

    mov $1, %rdi
    lea ok_msg(%rip), %rsi
    mov $3, %rdx
    mov $1, %rax
    syscall
fail:
"""
        + EXIT
        + r"""
path:   .asciz "/f"
ok_msg: .ascii " ok"
"""
    )
    (machine.fs / "f").write_text("abcdef")
    machine.run()
    assert out.getvalue() == b"abcdef ok"
