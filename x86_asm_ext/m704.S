/* m704.S -- (C) Geoffrey Reynolds, December 2008.

   These routines are based on Pierrick Gaudrys REDC code in GMP-ECM.


   void mulmod704(uint64_t *R, const uint64_t *A, const uint64_t *B,
                  const uint64_t *N, uint64_t inv);

   Assign R <-- A*B (mod N) in Montgomery form, where A,B,N < 2^704.
   Assumes that A,B are in Montgomery form, N is odd, and N[0]*inv = -1.


   void mulmod704_proth0(uint64_t *R, const uint64_t *A, const uint64_t *B,
                         const uint64_t *N);

   As mulmod704, but assumes that N[8-1] = 0, N[0] = 1.


   void mulmod704_proth1(uint64_t *R, const uint64_t *A, const uint64_t *B,
                         const uint64_t *N);

   As mulmod704, but assumes that N[9-1] = 0, N[0] = 1.


   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

	.text
	.globl	_mulmod704
	.globl	mulmod704
	.align	16

_mulmod704:
mulmod704:
	push	%rbx
	push	%rbp
#ifdef _WIN64
	push	%rdi
	push	%rsi
#endif
	push	%r12
	push	%r13
	push	%r14
	sub	$96, %rsp

#ifdef _WIN64
	mov	%rdx, %r13	/* A */
	mov	%r9, %r10	/* N */
	mov	%r8, %r9	/* B */
	mov	%rcx, %rdi	/* R */
	mov	192(%rsp), %r8	/* inv */
#else
	mov	%rsi, %r13	/* A */
	mov	%rdx, %r9	/* B */
	mov	%rcx, %r10	/* N */
#endif

	/* i = 0, j = 0 */

	movq	(%r13), %r14
	movq	(%r9), %rax

	xor	%ecx, %ecx
	lea	(%rsp), %rbp
	xor	%r12, %r12

	mulq	%r14
	addq	$1, %r12

	movq 	%rax, %rsi
	movq	%rdx, %rbx

	imulq	%r8, %rax
	movq	%rax, %r11

	mulq	(%r10)
	addq	%rax, %rsi
	movq	8(%r9), %rax
	adcq	%rdx, %rbx
	setc	%cl

	/* j = 1 */

	movq	%rcx, %rsi

	mulq	%r14
	addq	%rax, %rbx
	movq	8(%r10), %rax
	adcq	%rdx, %rsi

	mulq	%r11
	addq	%rbx, %rax
	movq	%rax, 0(%rbp)
	movq	16(%r9), %rax
	adcq	%rdx, %rsi
	setc	%cl

	/* j = 2 */

	movq	%rcx, %rbx

	mulq	%r14
	addq	%rax, %rsi
	movq	16(%r10), %rax
	adcq	%rdx, %rbx

	mulq	%r11
	addq	%rsi, %rax
	movq	%rax, 8(%rbp)
	movq	24(%r9), %rax
	adcq	%rdx, %rbx
	setc	%cl

	/* j = 3 */

	movq	%rcx, %rsi

	mulq	%r14
	addq	%rax, %rbx
	movq	24(%r10), %rax
	adcq	%rdx, %rsi

	mulq	%r11
	addq	%rbx, %rax
	movq	%rax, 16(%rbp)
	movq	32(%r9), %rax
	adcq	%rdx, %rsi
	setc	%cl

	/* j = 4 */

	movq	%rcx, %rbx

	mulq	%r14
	addq	%rax, %rsi
	movq	32(%r10), %rax
	adcq	%rdx, %rbx

	mulq	%r11
	addq	%rsi, %rax
	movq	%rax, 24(%rbp)
	movq	40(%r9), %rax
	adcq	%rdx, %rbx
	setc	%cl

	/* j = 5 */

	movq	%rcx, %rsi

	mulq	%r14
	addq	%rax, %rbx
	movq	40(%r10), %rax
	adcq	%rdx, %rsi

	mulq	%r11
	addq	%rbx, %rax
	movq	%rax, 32(%rbp)
	movq	48(%r9), %rax
	adcq	%rdx, %rsi
	setc	%cl

	/* j = 6 */

	movq	%rcx, %rbx

	mulq	%r14
	addq	%rax, %rsi
	movq	48(%r10), %rax
	adcq	%rdx, %rbx

	mulq	%r11
	addq	%rsi, %rax
	movq	%rax, 40(%rbp)
	movq	56(%r9), %rax
	adcq	%rdx, %rbx
	setc	%cl

	/* j = 7 */

	movq	%rcx, %rsi

	mulq	%r14
	addq	%rax, %rbx
	movq	56(%r10), %rax
	adcq	%rdx, %rsi

	mulq	%r11
	addq	%rbx, %rax
	movq	%rax, 48(%rbp)
	movq	64(%r9), %rax
	adcq	%rdx, %rsi
	setc	%cl

	/* j = 8 */

	movq	%rcx, %rbx

	mulq	%r14
	addq	%rax, %rsi
	movq	64(%r10), %rax
	adcq	%rdx, %rbx

	mulq	%r11
	addq	%rsi, %rax
	movq	%rax, 56(%rbp)
	movq	72(%r9), %rax
	adcq	%rdx, %rbx
	setc	%cl

	/* j = 9 */

	movq	%rcx, %rsi

	mulq	%r14
	addq	%rax, %rbx
	movq	72(%r10), %rax
	adcq	%rdx, %rsi

	mulq	%r11
	addq	%rbx, %rax
	movq	%rax, 64(%rbp)
	movq	80(%r9), %rax
	adcq	%rdx, %rsi
	setc	%cl

	/* j = 10 */

	movq	%rcx, %rbx

	mulq	%r14
	addq	%rax, %rsi
	movq	80(%r10), %rax
	adcq	%rdx, %rbx
	mulq    %r11
	addq	%rax, %rsi
	movq	%rsi, 72(%rbp)
	adcq	%rdx, %rbx
	movq	%rbx, 80(%rbp)
	setc	%cl
	movq	%rcx, 88(%rbp)

	.align 16
1:
	/* i > 0, j = 0 */

	movq	(%r13,%r12,8), %r14
	movq	(%r9), %rax
        movq	(%rbp), %rsi
	movq	8(%rbp), %rbx

	mulq	%r14
	addq	$1, %r12

	addq	%rsi, %rax
	adcq	%rdx, %rbx
	setc	%cl

	movq 	%rax, %rsi
	imulq	%r8, %rax
	movq	%rax, %r11

	mulq	(%r10)
	addq	%rax, %rsi
	adcq	%rdx, %rbx

	movq	8(%r9), %rax

	/* j = 1 */

	movq	%rcx, %rsi
	adcq	16(%rbp), %rsi
	setc	%cl

	mulq	%r14
	addq	%rax, %rbx
	movq	8(%r10), %rax
	adcq	%rdx, %rsi
	adcb	$0, %cl

	mulq	%r11
	addq	%rbx, %rax
	movq	%rax, 0(%rbp)
	adcq	%rdx, %rsi
	movq	16(%r9), %rax

	/* j = 2 */

	movq	%rcx, %rbx
	adcq	24(%rbp), %rbx
	setc	%cl

	mulq	%r14
	addq	%rax, %rsi
	movq	16(%r10), %rax
	adcq	%rdx, %rbx
	adcb	$0, %cl

	mulq	%r11
	addq	%rsi, %rax
	movq	%rax, 8(%rbp)
	adcq	%rdx, %rbx
	movq	24(%r9), %rax

	/* j = 3 */

	movq	%rcx, %rsi
	adcq	32(%rbp), %rsi
	setc	%cl

	mulq	%r14
	addq	%rax, %rbx
	movq	24(%r10), %rax
	adcq	%rdx, %rsi
	adcb	$0, %cl

	mulq	%r11
	addq	%rbx, %rax
	movq	%rax, 16(%rbp)
	adcq	%rdx, %rsi
	movq	32(%r9), %rax

	/* j = 4 */

	movq	%rcx, %rbx
	adcq	40(%rbp), %rbx
	setc	%cl

	mulq	%r14
	addq	%rax, %rsi
	movq	32(%r10), %rax
	adcq	%rdx, %rbx
	adcb	$0, %cl

	mulq	%r11
	addq	%rsi, %rax
	movq	%rax, 24(%rbp)
	adcq	%rdx, %rbx
	movq	40(%r9), %rax

	/* j = 5 */

	movq	%rcx, %rsi
	adcq	48(%rbp), %rsi
	setc	%cl

	mulq	%r14
	addq	%rax, %rbx
	movq	40(%r10), %rax
	adcq	%rdx, %rsi
	adcb	$0, %cl

	mulq	%r11
	addq	%rbx, %rax
	movq	%rax, 32(%rbp)
	adcq	%rdx, %rsi
	movq	48(%r9), %rax

	/* j = 6 */

	movq	%rcx, %rbx
	adcq	56(%rbp), %rbx
	setc	%cl

	mulq	%r14
	addq	%rax, %rsi
	movq	48(%r10), %rax
	adcq	%rdx, %rbx
	adcb	$0, %cl

	mulq	%r11
	addq	%rsi, %rax
	movq	%rax, 40(%rbp)
	adcq	%rdx, %rbx
	movq	56(%r9), %rax

	/* j = 7 */

	movq	%rcx, %rsi
	adcq	64(%rbp), %rsi
	setc	%cl

	mulq	%r14
	addq	%rax, %rbx
	movq	56(%r10), %rax
	adcq	%rdx, %rsi
	adcb	$0, %cl

	mulq	%r11
	addq	%rbx, %rax
	movq	%rax, 48(%rbp)
	adcq	%rdx, %rsi
	movq	64(%r9), %rax

	/* j = 8 */

	movq	%rcx, %rbx
	adcq	72(%rbp), %rbx
	setc	%cl

	mulq	%r14
	addq	%rax, %rsi
	movq	64(%r10), %rax
	adcq	%rdx, %rbx
	adcb	$0, %cl

	mulq	%r11
	addq	%rsi, %rax
	movq	%rax, 56(%rbp)
	adcq	%rdx, %rbx
	movq	72(%r9), %rax

	/* j = 9 */

	movq	%rcx, %rsi
	adcq	80(%rbp), %rsi
	setc	%cl

	mulq	%r14
	addq	%rax, %rbx
	movq	72(%r10), %rax
	adcq	%rdx, %rsi
	adcb	$0, %cl

	mulq	%r11
	addq	%rbx, %rax
	movq	%rax, 64(%rbp)
	adcq	%rdx, %rsi
	movq	80(%r9), %rax

	/* j = 10 */

	movq	%rcx, %rbx
	adcq	88(%rbp), %rbx
	setc	%cl

	mulq	%r14
	addq	%rax, %rsi
	movq	80(%r10), %rax
	adcq	%rdx, %rbx
	adcb	$0, %cl
	mulq    %r11
	addq	%rax, %rsi
	movq	%rsi, 72(%rbp)
	adcq	%rdx, %rbx
	movq	%rbx, 80(%rbp)
	adcb	$0, %cl
	movq	%rcx, 88(%rbp)

	cmpq	$11, %r12
	jb	1b


	mov	(%rbp), %rax
	mov	8(%rbp), %rbx
	mov	16(%rbp), %rdx
	mov	24(%rbp), %rsi
	mov	32(%rbp), %r8
	mov	40(%rbp), %r9
	mov	48(%rbp), %r11
	mov	56(%rbp), %r12
	mov	64(%rbp), %r13
	mov	72(%rbp), %r14
	mov	80(%rbp), %rbp
0:
	mov	%rax, (%rdi)
	mov	%rbx, 8(%rdi)
	mov	%rdx, 16(%rdi)
	mov	%rsi, 24(%rdi)
	mov	%r8, 32(%rdi)
	mov	%r9, 40(%rdi)
	mov	%r11, 48(%rdi)
	mov	%r12, 56(%rdi)
	mov	%r13, 64(%rdi)
	mov	%r14, 72(%rdi)
	mov	%rbp, 80(%rdi)

	sub	(%r10), %rax
	sbb	8(%r10), %rbx
	sbb	16(%r10), %rdx
	sbb	24(%r10), %rsi
	sbb	32(%r10), %r8
	sbb	40(%r10), %r9
	sbb	48(%r10), %r11
	sbb	56(%r10), %r12
	sbb	64(%r10), %r13
	sbb	72(%r10), %r14
	sbb	80(%r10), %rbp
	sbb	$0, %rcx
	jnc	0b

	add	$96, %rsp
	pop	%r14
	pop	%r13
	pop	%r12
#ifdef _WIN64
	pop	%rsi
	pop	%rdi
#endif
	pop	%rbp
	pop	%rbx
	ret




	.text
	.globl	_mulmod704_proth0
	.globl	mulmod704_proth0
	.align	16

_mulmod704_proth0:
mulmod704_proth0:
	push	%rbx
	push	%rbp
#ifdef _WIN64
	push	%rdi
	push	%rsi
#endif
	push	%r12
	push	%r13
	push	%r14
	sub	$96, %rsp

#ifdef _WIN64
	mov	%rdx, %r13	/* A */
	mov	%r9, %r10	/* N */
	mov	%r8, %r9	/* B */
	mov	%rcx, %rdi	/* R */
#else
	mov	%rsi, %r13	/* A */
	mov	%rdx, %r9	/* B */
	mov	%rcx, %r10	/* N */
#endif

	/* i = 0, j = 0 */

	movq	(%r13), %r8
	movq	(%r9), %rax

	xor	%ecx, %ecx
	lea	(%rsp), %rbp
	xor	%r12, %r12

	mulq	%r8
	addq	$1, %r12

	movq	%rdx, %rbx

	neg	%rax
	movq	%rax, %r11
	movq	8(%r9), %rax
	adcq	$0, %rbx

	/* j = 1 */

	mulq	%r8
	xor	%esi, %esi
	addq	%rax, %rbx
	adcq	%rdx, %rsi

	movq	%rbx, 0(%rbp)
	movq	16(%r9), %rax

	/* j = 2 */

	mulq	%r8
	xor	%ebx, %ebx
	addq	%rax, %rsi
	adcq	%rdx, %rbx

	movq	%rsi, 8(%rbp)
	movq	24(%r9), %rax

	/* j = 3 */

	mulq	%r8
	xor	%esi, %esi
	addq	%rax, %rbx
	adcq	%rdx, %rsi

	movq	%rbx, 16(%rbp)
	movq	32(%r9), %rax

	/* j = 4 */

	mulq	%r8
	xor	%ebx, %ebx
	addq	%rax, %rsi
	adcq	%rdx, %rbx

	movq	%rsi, 24(%rbp)
	movq	40(%r9), %rax

	/* j = 5 */

	mulq	%r8
	xor	%esi, %esi
	addq	%rax, %rbx
	adcq	%rdx, %rsi

	movq	%rbx, 32(%rbp)
	movq	48(%r9), %rax

	/* j = 6 */

	mulq	%r8
	xor	%ebx, %ebx
	addq	%rax, %rsi
	adcq	%rdx, %rbx

	movq	%rsi, 40(%rbp)
	movq	56(%r9), %rax

	/* j = 7 */

	mulq	%r8
	xor	%esi, %esi
	addq	%rax, %rbx
	adcq	%rdx, %rsi

	movq	%rbx, 48(%rbp)
	movq	64(%r9), %rax

	/* j = 8 */

	mulq	%r8
	xor	%ebx, %ebx
	addq	%rax, %rsi
	adcq	%rdx, %rbx

	movq	%rsi, 56(%rbp)
	movq	72(%r9), %rax

	/* j = 9 */

	mulq	%r8
	xor	%esi, %esi
	addq	%rax, %rbx
	movq	72(%r10), %rax
	adcq	%rdx, %rsi

	mulq	%r11
	addq	%rbx, %rax
	movq	%rax, 64(%rbp)
	movq	80(%r9), %rax
	adcq	%rdx, %rsi
	setc	%cl

	/* j = 10 */

	movq	%rcx, %rbx

	mulq	%r8
	addq	%rax, %rsi
	movq	80(%r10), %rax
	adcq	%rdx, %rbx
	mulq    %r11
	addq	%rax, %rsi
	movq	%rsi, 72(%rbp)
	adcq	%rdx, %rbx
	movq	%rbx, 80(%rbp)
	setc	%cl
	movq	%rcx, 88(%rbp)

	.align 16
1:
	/* i > 0, j = 0 */

	movq	(%r13,%r12,8), %r8
	movq	(%r9), %rax
        movq	(%rbp), %rsi
	movq	8(%rbp), %rbx

	mulq	%r8
	addq	$1, %r12

	addq	%rsi, %rax
	adcq	%rdx, %rbx
	setc	%cl

	neg	%rax
	movq	%rax, %r11
	adcq	$0, %rbx
	movq	8(%r9), %rax

	/* j = 1 */

	movq	%rcx, %rsi
	adcq	16(%rbp), %rsi
	setc	%cl

	mulq	%r8
	addq	%rax, %rbx
	adcq	%rdx, %rsi
	adcb	$0, %cl

	movq	%rbx, 0(%rbp)
	movq	16(%r9), %rax

	/* j = 2 */

	movq	%rcx, %rbx
	add	24(%rbp), %rbx
	setc	%cl

	mulq	%r8
	addq	%rax, %rsi
	adcq	%rdx, %rbx
	adcb	$0, %cl

	movq	%rsi, 8(%rbp)
	movq	24(%r9), %rax

	/* j = 3 */

	movq	%rcx, %rsi
	add	32(%rbp), %rsi
	setc	%cl

	mulq	%r8
	addq	%rax, %rbx
	adcq	%rdx, %rsi
	adcb	$0, %cl

	movq	%rbx, 16(%rbp)
	movq	32(%r9), %rax

	/* j = 4 */

	movq	%rcx, %rbx
	add	40(%rbp), %rbx
	setc	%cl

	mulq	%r8
	addq	%rax, %rsi
	adcq	%rdx, %rbx
	adcb	$0, %cl

	movq	%rsi, 24(%rbp)
	movq	40(%r9), %rax

	/* j = 5 */

	movq	%rcx, %rsi
	add	48(%rbp), %rsi
	setc	%cl

	mulq	%r8
	addq	%rax, %rbx
	adcq	%rdx, %rsi
	adcb	$0, %cl

	movq	%rbx, 32(%rbp)
	movq	48(%r9), %rax

	/* j = 6 */

	movq	%rcx, %rbx
	add	56(%rbp), %rbx
	setc	%cl

	mulq	%r8
	addq	%rax, %rsi
	adcq	%rdx, %rbx
	adcb	$0, %cl

	movq	%rsi, 40(%rbp)
	movq	56(%r9), %rax

	/* j = 7 */

	movq	%rcx, %rsi
	add	64(%rbp), %rsi
	setc	%cl

	mulq	%r8
	addq	%rax, %rbx
	adcq	%rdx, %rsi
	adcb	$0, %cl

	movq	%rbx, 48(%rbp)
	movq	64(%r9), %rax

	/* j = 8 */

	movq	%rcx, %rbx
	add	72(%rbp), %rbx
	setc	%cl

	mulq	%r8
	addq	%rax, %rsi
	adcq	%rdx, %rbx
	adcb	$0, %cl

	movq	%rsi, 56(%rbp)
	movq	72(%r9), %rax

	/* j = 9 */

	movq	%rcx, %rsi
	add	80(%rbp), %rsi
	setc	%cl

	mulq	%r8
	addq	%rax, %rbx
	movq	72(%r10), %rax
	adcq	%rdx, %rsi
	adcb	$0, %cl

	mulq	%r11
	addq	%rbx, %rax
	movq	%rax, 64(%rbp)
	adcq	%rdx, %rsi
	movq	80(%r9), %rax

	/* j = 10 */

	movq	%rcx, %rbx
	adcq	88(%rbp), %rbx
	setc	%cl

	mulq	%r8
	addq	%rax, %rsi
	movq	80(%r10), %rax
	adcq	%rdx, %rbx
	adcb	$0, %cl
	mulq    %r11
	addq	%rax, %rsi
	movq	%rsi, 72(%rbp)
	adcq	%rdx, %rbx
	movq	%rbx, 80(%rbp)
	adcb	$0, %cl
	movq	%rcx, 88(%rbp)

	cmpq	$11, %r12
	jb	1b


	mov	(%rbp), %rax
	mov	8(%rbp), %rbx
	mov	16(%rbp), %rdx
	mov	24(%rbp), %rsi
	mov	32(%rbp), %r8
	mov	40(%rbp), %r9
	mov	48(%rbp), %r11
	mov	56(%rbp), %r12
	mov	64(%rbp), %r13
	mov	72(%rbp), %r14
	mov	80(%rbp), %rbp
0:
	mov	%rax, (%rdi)
	mov	%rbx, 8(%rdi)
	mov	%rdx, 16(%rdi)
	mov	%rsi, 24(%rdi)
	mov	%r8, 32(%rdi)
	mov	%r9, 40(%rdi)
	mov	%r11, 48(%rdi)
	mov	%r12, 56(%rdi)
	mov	%r13, 64(%rdi)
	mov	%r14, 72(%rdi)
	mov	%rbp, 80(%rdi)

	sub	$1, %rax
	sbb	$0, %rbx
	sbb	$0, %rdx
	sbb	$0, %rsi
	sbb	$0, %r8
	sbb	$0, %r9
	sbb	$0, %r11
	sbb	$0, %r12
	sbb	$0, %r13
	sbb	72(%r10), %r14
	sbb	80(%r10), %rbp
	sbb	$0, %rcx
	jnc	0b

	add	$96, %rsp
	pop	%r14
	pop	%r13
	pop	%r12
#ifdef _WIN64
	pop	%rsi
	pop	%rdi
#endif
	pop	%rbp
	pop	%rbx
	ret




	.text
	.globl	_mulmod704_proth1
	.globl	mulmod704_proth1
	.align	16

_mulmod704_proth1:
mulmod704_proth1:
	push	%rbx
	push	%rbp
#ifdef _WIN64
	push	%rdi
	push	%rsi
#endif
	push	%r12
	push	%r13
	push	%r14
	sub	$96, %rsp

#ifdef _WIN64
	mov	%rdx, %r13	/* A */
	mov	80(%r9), %r10	/* N[10] */
	mov	%r8, %r9	/* B */
	mov	%rcx, %rdi	/* R */
#else
	mov	%rsi, %r13	/* A */
	mov	%rdx, %r9	/* B */
	mov	80(%rcx), %r10	/* N[10] */
#endif

	/* i = 0, j = 0*/

	movq	(%r13), %r8
	movq	(%r9), %rax

	xor	%ecx, %ecx
	lea	(%rsp), %rbp
	xor	%r12, %r12

	mulq	%r8
	addq	$1, %r12

	movq	%rdx, %rbx

	neg	%rax
	movq	%rax, %r11
	movq	8(%r9), %rax
	adcq	$0, %rbx

	/* j = 1 */

	mulq	%r8
	xor	%esi, %esi
	addq	%rax, %rbx
	adcq	%rdx, %rsi

	movq	%rbx, 0(%rbp)
	movq	16(%r9), %rax

	/* j = 2 */

	mulq	%r8
	xor	%ebx, %ebx
	addq	%rax, %rsi
	adcq	%rdx, %rbx

	movq	%rsi, 8(%rbp)
	movq	24(%r9), %rax

	/* j = 3 */

	mulq	%r8
	xor	%esi, %esi
	addq	%rax, %rbx
	adcq	%rdx, %rsi

	movq	%rbx, 16(%rbp)
	movq	32(%r9), %rax

	/* j = 4 */

	mulq	%r8
	xor	%ebx, %ebx
	addq	%rax, %rsi
	adcq	%rdx, %rbx

	movq	%rsi, 24(%rbp)
	movq	40(%r9), %rax

	/* j = 5 */

	mulq	%r8
	xor	%esi, %esi
	addq	%rax, %rbx
	adcq	%rdx, %rsi

	movq	%rbx, 32(%rbp)
	movq	48(%r9), %rax

	/* j = 6 */

	mulq	%r8
	xor	%ebx, %ebx
	addq	%rax, %rsi
	adcq	%rdx, %rbx

	movq	%rsi, 40(%rbp)
	movq	56(%r9), %rax

	/* j = 7 */

	mulq	%r8
	xor	%esi, %esi
	addq	%rax, %rbx
	adcq	%rdx, %rsi

	movq	%rbx, 48(%rbp)
	movq	64(%r9), %rax

	/* j = 8 */

	mulq	%r8
	xor	%ebx, %ebx
	addq	%rax, %rsi
	adcq	%rdx, %rbx

	movq	%rsi, 56(%rbp)
	movq	72(%r9), %rax

	/* j = 9 */

	mulq	%r8
	xor	%esi, %esi
	addq	%rax, %rbx
	adcq	%rdx, %rsi

	movq	%rbx, 64(%rbp)
	movq	80(%r9), %rax

	/* j = 10 */

	mulq	%r8
	xor	%ebx, %ebx
	addq	%rax, %rsi
	movq	%r10, %rax
	adcq	%rdx, %rbx
	mulq    %r11
	addq	%rax, %rsi
	movq	%rsi, 72(%rbp)
	adcq	%rdx, %rbx
	movq	%rbx, 80(%rbp)
	setc	%cl
	movq	%rcx, 88(%rbp)

	.align 16
1:
	/* i > 0, j = 0 */

	movq	(%r13,%r12,8), %r8
	movq	(%r9), %rax
        movq	(%rbp), %rsi
	movq	8(%rbp), %rbx

	mulq	%r8
	addq	$1, %r12

	addq	%rsi, %rax
	adcq	%rdx, %rbx
	setc	%cl

	neg	%rax
	movq	%rax, %r11
	adcq	$0, %rbx
	movq	8(%r9), %rax

	/* j = 1 */

	movq	%rcx, %rsi
	adcq	16(%rbp), %rsi
	setc	%cl

	mulq	%r8
	addq	%rax, %rbx
	adcq	%rdx, %rsi
	adcb	$0, %cl

	movq	%rbx, 0(%rbp)
	movq	16(%r9), %rax

	/* j = 2 */

	movq	%rcx, %rbx
	add	24(%rbp), %rbx
	setc	%cl

	mulq	%r8
	addq	%rax, %rsi
	adcq	%rdx, %rbx
	adcb	$0, %cl

	movq	%rsi, 8(%rbp)
	movq	24(%r9), %rax

	/* j = 3 */

	movq	%rcx, %rsi
	add	32(%rbp), %rsi
	setc	%cl

	mulq	%r8
	addq	%rax, %rbx
	adcq	%rdx, %rsi
	adcb	$0, %cl

	movq	%rbx, 16(%rbp)
	movq	32(%r9), %rax

	/* j = 4 */

	movq	%rcx, %rbx
	add	40(%rbp), %rbx
	setc	%cl

	mulq	%r8
	addq	%rax, %rsi
	adcq	%rdx, %rbx
	adcb	$0, %cl

	movq	%rsi, 24(%rbp)
	movq	40(%r9), %rax

	/* j = 5 */

	movq	%rcx, %rsi
	add	48(%rbp), %rsi
	setc	%cl

	mulq	%r8
	addq	%rax, %rbx
	adcq	%rdx, %rsi
	adcb	$0, %cl

	movq	%rbx, 32(%rbp)
	movq	48(%r9), %rax

	/* j = 6 */

	movq	%rcx, %rbx
	add	56(%rbp), %rbx
	setc	%cl

	mulq	%r8
	addq	%rax, %rsi
	adcq	%rdx, %rbx
	adcb	$0, %cl

	movq	%rsi, 40(%rbp)
	movq	56(%r9), %rax

	/* j = 7 */

	movq	%rcx, %rsi
	add	64(%rbp), %rsi
	setc	%cl

	mulq	%r8
	addq	%rax, %rbx
	adcq	%rdx, %rsi
	adcb	$0, %cl

	movq	%rbx, 48(%rbp)
	movq	64(%r9), %rax

	/* j = 8 */

	movq	%rcx, %rbx
	add	72(%rbp), %rbx
	setc	%cl

	mulq	%r8
	addq	%rax, %rsi
	adcq	%rdx, %rbx
	adcb	$0, %cl

	movq	%rsi, 56(%rbp)
	movq	72(%r9), %rax

	/* j = 9 */

	movq	%rcx, %rsi
	add	80(%rbp), %rsi
	setc	%cl

	mulq	%r8
	addq	%rax, %rbx
	adcq	%rdx, %rsi
	adcb	$0, %cl

	movq	%rbx, 64(%rbp)
	movq	80(%r9), %rax

	/* j = 10 */

	movq	%rcx, %rbx
	add	88(%rbp), %rbx
	setc	%cl

	mulq	%r8
	addq	%rax, %rsi
	movq	%r10, %rax
	adcq	%rdx, %rbx
	adcb	$0, %cl
	mulq    %r11
	addq	%rax, %rsi
	movq	%rsi, 72(%rbp)
	adcq	%rdx, %rbx
	movq	%rbx, 80(%rbp)
	adcb	$0, %cl
	movq	%rcx, 88(%rbp)

	cmpq	$11, %r12
	jb	1b


	mov	(%rbp), %rax
	mov	8(%rbp), %rbx
	mov	16(%rbp), %rdx
	mov	24(%rbp), %rsi
	mov	32(%rbp), %r8
	mov	40(%rbp), %r9
	mov	48(%rbp), %r11
	mov	56(%rbp), %r12
	mov	64(%rbp), %r13
	mov	72(%rbp), %r14
	mov	80(%rbp), %rbp
0:
	mov	%rax, (%rdi)
	mov	%rbx, 8(%rdi)
	mov	%rdx, 16(%rdi)
	mov	%rsi, 24(%rdi)
	mov	%r8, 32(%rdi)
	mov	%r9, 40(%rdi)
	mov	%r11, 48(%rdi)
	mov	%r12, 56(%rdi)
	mov	%r13, 64(%rdi)
	mov	%r14, 72(%rdi)
	mov	%rbp, 80(%rdi)

	sub	$1, %rax
	sbb	$0, %rbx
	sbb	$0, %rdx
	sbb	$0, %rsi
	sbb	$0, %r8
	sbb	$0, %r9
	sbb	$0, %r11
	sbb	$0, %r12
	sbb	$0, %r13
	sbb	$0, %r14
	sbb	%r10, %rbp
	sbb	$0, %rcx
	jnc	0b

	add	$96, %rsp
	pop	%r14
	pop	%r13
	pop	%r12
#ifdef _WIN64
	pop	%rsi
	pop	%rdi
#endif
	pop	%rbp
	pop	%rbx
	ret
