.section .rodata
dev_null:
  .asciz "/dev/null"

.section .text
.global _start
.type _start, @function
_start:
  mov $2, %rax
  lea dev_null(%rip), %rdi
  xor %rsi, %rsi
  xor %rdx, %rdx
  syscall

  mov %rax, %rdi
  mov $3, %rax
  syscall

  mov $9, %rax
  xor %rdi, %rdi
  mov $4096, %rsi
  mov $3, %rdx
  mov $0x22, %r10
  mov $-1, %r8
  xor %r9, %r9
  syscall

  mov %rax, %r12
  movabs $0x57202c6f6c6c6548, %rax
  mov %rax, (%r12)
  movabs $0x0a21646c726f, %rax
  mov %rax, 8(%r12)

  mov $1, %rax
  mov $1, %rdi
  mov %r12, %rsi
  mov $14, %rdx
  syscall

  mov $11, %rax
  mov %r12, %rdi
  mov $4096, %rsi
  syscall

  mov $231, %rax
  xor %rdi, %rdi
  syscall

.section .note.GNU-stack,"",@progbits
