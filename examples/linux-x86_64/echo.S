.section .rodata
space:
  .ascii " "
newline:
  .ascii "\n"

.section .text
.global _start
.type _start, @function
_start:
  mov (%rsp), %r12
  lea 8(%rsp), %r13
  mov $1, %r14

arg_loop:
  cmp %r12, %r14
  jae done_args

  mov (%r13,%r14,8), %rsi
  mov %rsi, %r15

strlen_loop:
  cmpb $0, (%r15)
  je write_arg
  inc %r15
  jmp strlen_loop

write_arg:
  mov %r15, %rdx
  sub %rsi, %rdx
  mov $1, %rax
  mov $1, %rdi
  syscall

  inc %r14
  cmp %r12, %r14
  jae done_args

  mov $1, %rax
  mov $1, %rdi
  lea space(%rip), %rsi
  mov $1, %rdx
  syscall
  jmp arg_loop

done_args:
  mov $1, %rax
  mov $1, %rdi
  lea newline(%rip), %rsi
  mov $1, %rdx
  syscall

  mov $231, %rax
  xor %rdi, %rdi
  syscall

.section .note.GNU-stack,"",@progbits
