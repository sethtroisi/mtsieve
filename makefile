# Makefile for mtsieve applications. (C) Mark Rodenkirch, January 2017
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.

# Set DEBUG=yes to compile with debugging information and internal checks.

DEBUG=no
CC=g++

CPPFLAGS=-Isieve -m64 -Wall
LDFLAGS=-lstdc++
GPUCPPFLAGS=-DHAVE_GPU_WORKERS
GPULIBS=

ifeq ($(strip $(DEBUG)),yes)
   CPPFLAGS+=-g
else
   CPPFLAGS+=-O2
endif

ifeq ($(OS),Windows_NT)
   LDFLAGS+=-static
   CPPFLAGS+=-I"F:\stuff\gmp-6.1.2"
   GPUCPPFLAGS+=-I"C:\Program Files (x86)\AMD APP SDK\3.0\include"
   GMPLDFLAGS=-L"F:\stuff\gmp-6.1.2\.libs" -lgmp
   GPULDFLAGS="C:\Program Files (x86)\AMD APP SDK\3.0\lib\x86_64\OpenCl.lib"
else
   UNAME_S := $(shell uname -s)
   
   ifeq ($(UNAME_S),Darwin)
      CPPFLAGS+=-std=c++0x
      GMPLDFLAGS=-lgmp
      GPULDFLAGS+=-framework OpenCL
   else
      CPPFLAGS+=-std=c++11
      EXTRALDFLAGS+=-lpthread
      GMPLDFLAGS=-lgmp
      GPULDFLAGS+=-lOpenCL
   endif
endif
   
# No changes should be needed below here.

CPU_PROGS=afsieve cksieve dmdsieve gcwsieve gfndsieve fbncsieve fkbnsieve kbbsieve mfsieve pixsieve psieve srsieve2 twinsieve xyyxsieve

GPU_PROGS=afsievecl mfsievecl gfndsievecl pixsievecl xyyxsievecl

CPU_CORE_OBJS=core/App_cpu.o core/FactorApp.o core/AlgebraicFactorApp.o \
   core/Clock.o core/Parser.o core/Worker_cpu.o core/HashTable.o core/main.o core/SharedMemoryItem.o
   
GPU_CORE_OBJS=core/App_gpu.o core/FactorApp.o core/AlgebraicFactorApp.o \
   core/Clock.o core/Parser.o core/Worker_gpu.o core/HashTable.o core/main.o core/SharedMemoryItem.o \
   opencl/Device_gpu.o opencl/Kernel_gpu.o opencl/KernelArgument_gpu.o opencl/ErrorChecker_gpu.o

ASM_OBJS=x86_asm/fpu_mod_init_fini.o x86_asm/fpu_push_pop.o x86_asm/sse_mulmod.o \
   x86_asm/fpu_mulmod.o x86_asm/fpu_powmod.o x86_asm/fpu_powmod_4b_1n_4p.o \
   x86_asm/fpu_mulmod_iter.o x86_asm/fpu_mulmod_iter_4a.o x86_asm/fpu_mulmod_4a_4b_4p.o \
   x86_asm/sse_mod_init_fini.o  x86_asm/sse_powmod_4b_1n_4p.o x86_asm/sse_mulmod_4a_4b_4p.o \
   x86_asm/avx_set_a.o x86_asm/avx_set_b.o x86_asm/avx_get.o \
   x86_asm/avx_compute_reciprocal.o x86_asm/avx_compare.o \
   x86_asm/avx_mulmod.o x86_asm/avx_powmod.o x86_asm/sse_powmod_4b_1n_4p_mulmod_1k.o

ASM_EXT_OBJS=x86_asm_ext/m320.o x86_asm_ext/m384.o x86_asm_ext/m448.o x86_asm_ext/m512.o \
   x86_asm_ext/m576.o x86_asm_ext/m640.o x86_asm_ext/m704.o x86_asm_ext/m768.o \
   x86_asm_ext/mulmod128.o x86_asm_ext/mulmod192.o x86_asm_ext/mulmod256.o \
   x86_asm_ext/sqrmod128.o x86_asm_ext/sqrmod192.o x86_asm_ext/sqrmod256.o \
   x86_asm_ext/redc.o
   
PRIMESIEVE_OBJS=sieve/Erat.o sieve/EratBig.o sieve/EratMedium.o sieve/EratSmall.o sieve/PreSieve.o \
   sieve/CpuInfo.o sieve/MemoryPool.o sieve/PrimeGenerator.o sieve/PrimeSieve.o sieve/Wheel.o \
   sieve/IteratorHelper.o sieve/popcount.o sieve/nthPrime.o sieve/PrintPrimes.o \
   sieve/ParallelSieve.o sieve/iterator.o sieve/api.o sieve/SievingPrimes.o

AF_OBJS=alternating_factorial/AlternatingFactorialApp.o alternating_factorial/AlternatingFactorialWorker.o alternating_factorial/afsieve.o
CK_OBJS=carol_kynea/CarolKyneaApp.o carol_kynea/CarolKyneaWorker.o
DMD_OBJS=dm_divisor/DMDivisorApp.o dm_divisor/DMDivisorWorker.o
FBNC_OBJS=fixed_bnc/FixedBNCApp.o fixed_bnc/FixedBNCWorker.o
FKBN_OBJS=fixed_kbn/FixedKBNApp.o fixed_kbn/FixedKBNWorker.o
GCW_OBJS=cullen_woodall/CullenWoodallApp.o cullen_woodall/CullenWoodallWorker.o
GFND_OBJS=gfn_divisor/GFNDivisorApp.o gfn_divisor/GFNDivisorWorker.o
KBB_OBJS=kbb/KBBApp.o kbb/KBBWorker.o
MF_OBJS=multi_factorial/MultiFactorialApp.o multi_factorial/MultiFactorialWorker.o multi_factorial/mfsieve.o multi_factorial/multifactorial.o
PIX_OBJS=primes_in_x/PrimesInXApp.o primes_in_x/PrimesInXWorker.o primes_in_x/pixsieve.o
PRIM_OBJS=primorial/PrimorialApp.o primorial/PrimorialWorker.o primorial/primorial.o
TWIN_OBJS=twin/TwinApp.o twin/TwinWorker.o
SR2_OBJS=sierpinski_riesel/SierpinskiRieselApp.o sierpinski_riesel/AlgebraicFactorHelper.o \
   sierpinski_riesel/AbstractSubsequenceHelper.o sierpinski_riesel/AbstractWorker.o \
   sierpinski_riesel/GenericSubsequenceHelper.o sierpinski_riesel/GenericWorker.o
XYYX_OBJS=xyyx/XYYXApp.o xyyx/XYYXWorker.o

AF_GPU_OBJS=alternating_factorial/AlternatingFactorialGpuWorker.o
GCW_GPU_OBJS=cullen_woodall/CullenWoodallGpuWorker.o
GFND_GPU_OBJS=gfn_divisor/GFNDivisorGpuWorker.o
MF_GPU_OBJS=multi_factorial/MultiFactorialGpuWorker.o
PIX_GPU_OBJS=primes_in_x/PrimesInXGpuWorker.o
XYYX_GPU_OBJS=xyyx/XYYXGpuWorker.o

ALL_OBJS=$(PRIMESIEVE_OBJS) $(ASM_OBJS) $(ASM_EXT_OBJS) $(CPU_CORE_OBJS) $(GPU_CORE_OBJS) \
   $(AF_OBJS) $(MF_OBJS) $(FBNC_OBJS) $(FKBN_OBJS) $(GFND_OBJS) $(CK_OBJS) \
   $(PIX_OBJS) $(XYYX_OBJS) $(KBB_OBJS) $(GCW_OBJS) $(PRIM_OBJS) $(TWIN_OBJS) \
   $(DMD_OBJS) $(SR2_OBJS) \
   $(AF_GPU_OBJS) $(GCW_GPU_OBJS) $(GFND_GPU_OBJS) $(MF_GPU_OBJS) $(PIX_GPU_OBJS) $(XYYX_GPU_OBJS)

all: $(CPU_PROGS) $(GPU_PROGS)

cpu_all: $(CPU_PROGS)

gpu_all: $(GPU_PROGS)

%.o: %.S
	$(CC) $(CPPFLAGS) -c -o $@ $< 
   
%.o: %.cpp
	$(CC) $(CPPFLAGS) -c -o $@ $< 
   
%_cpu.o: %.cpp
	$(CC) $(CPPFLAGS) -c -o $@ $< 
    
%_gpu.o: %.cpp
	$(CC) $(CPPFLAGS) $(GPUCPPFLAGS) -c -o $@ $<

afsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(AF_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(EXTRALDFLAGS)

afsievecl: $(GPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(AF_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(GPULDFLAGS) $(EXTRALDFLAGS)
   
cksieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(CK_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(EXTRALDFLAGS)
   
dmdsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(ASM_EXT_OBJS) $(DMD_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(GMPLDFLAGS) $(EXTRALDFLAGS)
   
fbncsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(FBNC_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(EXTRALDFLAGS)

fkbnsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(FKBN_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(EXTRALDFLAGS)
   
gcwsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(GCW_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(GMPLDFLAGS) $(EXTRALDFLAGS)

gcwsievecl: $(GPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(GCW_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(GPULDFLAGS) $(GMPLDFLAGS) $(EXTRALDFLAGS)
   
gfndsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(ASM_EXT_OBJS) $(GFND_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^  $(GMPLDFLAGS) $(EXTRALDFLAGS)

gfndsievecl: $(GPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(ASM_EXT_OBJS) $(GFND_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^  $(GMPLDFLAGS) $(GPULDFLAGS) $(EXTRALDFLAGS)

kbbsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(KBB_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(EXTRALDFLAGS)

mfsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(MF_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(EXTRALDFLAGS)

mfsievecl: $(GPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(MF_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(GPULDFLAGS) $(EXTRALDFLAGS)
   
pixsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(PIX_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(EXTRALDFLAGS)

pixsievecl: $(GPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(PIX_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(GPULDFLAGS) $(EXTRALDFLAGS)

psieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(PRIM_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(EXTRALDFLAGS)
   
srsieve2: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(SR2_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(EXTRALDFLAGS)

twinsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(TWIN_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(EXTRALDFLAGS)

xyyxsieve: $(CPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(XYYX_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(EXTRALDFLAGS)

xyyxsievecl: $(GPU_CORE_OBJS) $(PRIMESIEVE_OBJS) $(ASM_OBJS) $(XYYX_OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(GPULDFLAGS) $(EXTRALDFLAGS)

clean_objs:
	rm -f $(ALL_OBJS) *.log *.pfgw
   
clean:
	rm -f $(CPU_PROGS) $(GPU_PROGS) $(ALL_OBJS) *.log *.pfgw
