/* mf_kernel.cl -- (C) Mark Rodenkirch, September 2016

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
   
   Thanks to Yves Gallot for this implementation based upon 
   Peter L. Montgomery, Modular multiplication without trial division, Math. Comp.44 (1985), 519–521.
 */

#pragma OPENCL EXTENSION cl_khr_int64_extended_atomics: enable

ulong  invmod(ulong a, ulong p);

void collect_factor(uint    n,
                    ulong   p,
  volatile __global uint   *factorCount,
           __global ulong2 *factors);

ulong mmmInvert(ulong _p);
ulong mmmOne(ulong _p);
ulong mmmR2(ulong _p, ulong _q, ulong _one);
ulong mmmAdd(ulong a, ulong b, ulong _p);
ulong mmmSub(ulong a, ulong b, ulong _p);
ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q);
ulong mmmNtoRes(ulong n, ulong _p, ulong _q, ulong _r2);
ulong mmmPowmod(ulong resB, ulong exp, ulong _p, ulong _q, ulong _one);
ulong mmmResToN(ulong res, ulong _p, ulong _q);

__kernel void sm_kernel7(__global const  ulong  *primes,
                         __global const  uint   *terms,
                volatile __global        uint   *factorCount,
                         __global        ulong2 *factors)
{
   int    gid = get_global_id(0);
   
   ulong  thePrime = primes[gid];

   if (thePrime == 0) return;

   ulong  two9sq = 99;
   ulong  three9sq = 999;
   ulong  four9sq = 9999;
   ulong  five9sq = 99999;
   ulong  six9sq = 999999;
   ulong  seven9sq = 9999999;

   ulong  m = ((ulong) terms[0] * 9999999) + 10000000;
   ulong  exp = 7*terms[0] - 6999993;
   ulong  invmod2, invmod3, invmod4, invmod5, invmod6;
   ulong  resTemp, resInvmod;
   
   two9sq *= two9sq;
   three9sq *= three9sq;
   four9sq *= four9sq;
   five9sq *= five9sq;
   six9sq *= six9sq;
   seven9sq *= seven9sq;
   
   invmod2 = invmod(two9sq, thePrime);
  
   if (thePrime < three9sq)
      invmod3 = invmod(three9sq%thePrime, thePrime);
   else
      invmod3 = invmod(three9sq, thePrime);
   
   if (thePrime < four9sq)
      invmod4 = invmod(four9sq%thePrime, thePrime);
   else
      invmod4 = invmod(four9sq, thePrime);
   
   if (thePrime < five9sq)
      invmod5 = invmod(five9sq%thePrime, thePrime);
   else
      invmod5 = invmod(five9sq, thePrime);
   
   if (thePrime < six9sq)
      invmod6 = invmod(six9sq%thePrime, thePrime);
   else
      invmod6 = invmod(six9sq, thePrime);
      
   ulong _q = mmmInvert(thePrime);
   ulong _one = mmmOne(thePrime);
   ulong _r2 = mmmR2(thePrime, _q, _one);

   ulong res10e1 = mmmNtoRes(10, thePrime, _q, _r2);
   ulong res10en = mmmNtoRes(10, thePrime, _q, _r2);
   
   // calculate a1x
   //    ax = 123456789%f;
   ulong resAx = mmmNtoRes(123456789, thePrime, _q, _r2);

   // calculate a2x
   //   ax = mulmod(ax, 99*99, f);
   //   ax = (ax + 991)%f;
   //    t = powmod(10,180,f);
   //   ax = mulmod(ax, t, f);
   //        submod(ax, 100);
   //        divmod(ax, 99*99);
   //        submod(ax, 1);
         
   resTemp = mmmNtoRes(two9sq, thePrime, _q, _r2);
   resAx = mmmMulmod(resAx, resTemp, thePrime, _q);
   resTemp = mmmNtoRes(991, thePrime, _q, _r2);
   resAx = mmmAdd(resAx, resTemp, thePrime);
   resTemp = mmmPowmod(res10e1, 180, thePrime, _q, _one);
   resAx = mmmMulmod(resAx, resTemp, thePrime, _q);
   res10en = mmmMulmod(res10e1, res10e1, thePrime, _q);
   resAx = mmmSub(resAx, res10en, thePrime);
   resInvmod = mmmNtoRes(invmod2, thePrime, _q, _r2);
   resAx = mmmMulmod(resAx, resInvmod, thePrime, _q);
   resAx = mmmSub(resAx, _one, thePrime);

   // calculate a3x
   //   ax = mulmod(ax, 999*999, f);
   //   ax = (ax + 99901)%f;
   //    t = powmod(10,2700,f);
   //   ax = mulmod(ax, t, f);
   //        submod(ax, 1000);
   //        divmod(ax, 999*999);
   //        submod(ax, 1);
               
   resTemp = mmmNtoRes(three9sq, thePrime, _q, _r2);
   resAx = mmmMulmod(resAx, resTemp, thePrime, _q);
   resTemp = mmmNtoRes(99901, thePrime, _q, _r2);
   resAx = mmmAdd(resAx, resTemp, thePrime);
   resTemp = mmmPowmod(res10e1, 2700, thePrime, _q, _one);
   resAx = mmmMulmod(resAx, resTemp, thePrime, _q);
   res10en = mmmMulmod(res10en, res10e1, thePrime, _q);
   resAx = mmmSub(resAx, res10en, thePrime);
   resInvmod = mmmNtoRes(invmod3, thePrime, _q, _r2);
   resAx = mmmMulmod(resAx, resInvmod, thePrime, _q);
   resAx = mmmSub(resAx, _one, thePrime);
  
    // calculate a4x
   //   ax = mulmod(ax, 9999*9999, f);
   //   ax = (ax + 9999001)%f;
   //    t = powmod(10,36000,f);
   //   ax = mulmod(ax, t, f);
   //        submod(ax, 10000);
   //        divmod(ax, 9999*9999);
   //        submod(ax, 1);

   resTemp = mmmNtoRes(four9sq, thePrime, _q, _r2);
   resAx = mmmMulmod(resAx, resTemp, thePrime, _q);
   resTemp = mmmNtoRes(9999001, thePrime, _q, _r2);
   resAx = mmmAdd(resAx, resTemp, thePrime);
   resTemp = mmmPowmod(res10e1, 36000, thePrime, _q, _one);
   resAx = mmmMulmod(resAx, resTemp, thePrime, _q);
   res10en = mmmMulmod(res10en, res10e1, thePrime, _q);
   resAx = mmmSub(resAx, res10en, thePrime);
   resInvmod = mmmNtoRes(invmod4, thePrime, _q, _r2);
   resAx = mmmMulmod(resAx, resInvmod, thePrime, _q);
   resAx = mmmSub(resAx, _one, thePrime);

   // calculate a5x
   //   ax = mulmod(ax, 99999ULL*99999ULL, f);
   //   ax = (ax + 999990001ULL)%f;
   //    t = powmod(10,450000,f);
   //   ax = mulmod(ax, t, f);
   //        submod(ax, 100000);
   //        divmod(ax, 99999ULL*99999ULL);
   //        submod(ax, 1);

   resTemp = mmmNtoRes(five9sq, thePrime, _q, _r2);
   resAx = mmmMulmod(resAx, resTemp, thePrime, _q);
   resTemp = mmmNtoRes(999990001, thePrime, _q, _r2);
   resAx = mmmAdd(resAx, resTemp, thePrime);
   resTemp = mmmPowmod(res10e1, 450000, thePrime, _q, _one);
   resAx = mmmMulmod(resAx, resTemp, thePrime, _q);
   res10en = mmmMulmod(res10en, res10e1, thePrime, _q);
   resAx = mmmSub(resAx, res10en, thePrime);
   resInvmod = mmmNtoRes(invmod5, thePrime, _q, _r2);
   resAx = mmmMulmod(resAx, resInvmod, thePrime, _q);
   resAx = mmmSub(resAx, _one, thePrime);

   // calculate a6x
   //   ax = mulmod(ax, 999999ULL*999999ULL, f);
   //   ax = (ax + 99999900001ULL)%f;
   //    t = powmod(10,5400000,f);
   //   ax = mulmod(ax, t, f);
   //        submod(ax, 1000000);
   //        divmod(ax, 999999ULL*999999ULL);
   //        submod(ax, 1);
      
   resTemp = mmmNtoRes(six9sq, thePrime, _q, _r2);
   resAx = mmmMulmod(resAx, resTemp, thePrime, _q);
   resTemp = mmmNtoRes((ulong) 99999900001, thePrime, _q, _r2);
   resAx = mmmAdd(resAx, resTemp, thePrime);
   resTemp = mmmPowmod(res10e1, 5400000, thePrime, _q, _one);
   resAx = mmmMulmod(resAx, resTemp, thePrime, _q);
   res10en = mmmMulmod(res10en, res10e1, thePrime, _q);
   resAx = mmmSub(resAx, res10en, thePrime);
   resInvmod = mmmNtoRes(invmod6, thePrime, _q, _r2);
   resAx = mmmMulmod(resAx, resInvmod, thePrime, _q);
   resAx = mmmSub(resAx, _one, thePrime);
   
   // calculate a7(n)
   //   ax = mulmod(ax, 9999999 ULL*9999999ULL, f);
   //   ax = (ax + 9999999000001ULL)%f;
   //    t = powmod(10, 7*v[0]-6999993, f);
   //   ax = mulmod(ax, t, f);
   
   resTemp = mmmNtoRes(seven9sq, thePrime, _q, _r2);
   resAx = mmmMulmod(resAx, resTemp, thePrime, _q);
   resTemp = mmmNtoRes((ulong) 9999999000001, thePrime, _q, _r2);
   resAx = mmmAdd(resAx, resTemp, thePrime);
   resTemp = mmmPowmod(res10e1, exp, thePrime, _q, _one);
   resAx = mmmMulmod(resAx, resTemp, thePrime, _q);

   //    M = (9999999ULL*v[0]+10000000)%f;
   ulong resM = mmmNtoRes(m, thePrime, _q, _r2);
   
   ulong res9s = mmmNtoRes(9999999, thePrime, _q, _r2);
   resTemp = mmmNtoRes((ulong) 10000000000, thePrime, _q, _r2);
   ulong resR = mmmNtoRes((ulong) 100000000000, thePrime, _q, _r2);
   ulong resT[7];

   //   r = mulmod(10000000000ULL, 100000000000ULL, f); // r = 10^21%f
   //   T[1] = mulmod(r, r, f);
   //   T[2] = mulmod(T[1], T[1], f);
   //   T[3] = mulmod(T[1], T[2], f);
   //   T[4] = mulmod(T[2], T[2], f);
   //   T[5] = mulmod(T[2], T[3], f);
   //   T[6] = mulmod(T[3], T[3], f);

   resR = mmmMulmod(resR, resTemp, thePrime, _q);
   resT[1] = mmmMulmod(resR, resR, thePrime, _q);
   resT[2] = mmmMulmod(resT[1], resT[1], thePrime, _q);
   resT[3] = mmmMulmod(resT[1], resT[2], thePrime, _q);
   resT[4] = mmmMulmod(resT[2], resT[2], thePrime, _q);
   resT[5] = mmmMulmod(resT[2], resT[3], thePrime, _q);
   resT[6] = mmmMulmod(resT[3], resT[3], thePrime, _q);

   // This is only for the first term    
   if (resAx == resM)
      collect_factor(terms[0], thePrime, factorCount, factors);

   // This is for the remaining terms
   uint idx = 1;
   while (terms[idx] > 0)
   {         
      uint dn = terms[idx] - terms[idx-1];
      
      //    if(dn<=36)
      //      C = mulmod(C, T[dn/6], f);
      //    else {
      //      t = powmod(r, dn/3, f);
      //      C = mulmod(C, t, f);
      //    }
      
      if (dn <= 36)
         resAx = mmmMulmod(resAx, resT[dn/6], thePrime, _q);
      else
      {
         resTemp = mmmPowmod(resT[1], dn/6, thePrime, _q, _one);
         resAx = mmmMulmod(resAx, resTemp, thePrime, _q);
      }
      
      // M = (M + 9999999ULL*dn) % f;
      
      resTemp = mmmNtoRes(dn, thePrime, _q, _r2);
      resTemp = mmmMulmod(resTemp, res9s, thePrime, _q);
      resM = mmmAdd(resM, resTemp, thePrime);
      
      if (resAx == resM)
         collect_factor(terms[idx], thePrime, factorCount, factors);
         
      idx++;
   }
}

ulong mmmInvert(ulong _p)
{
   ulong p_inv = 1;
   ulong prev = 0;
   
   while (p_inv != prev)
   {
      prev = p_inv;
      p_inv *= (2 - _p * p_inv);
   }
   
   return p_inv;
}

// Compute the residual of 1 (mod p)
ulong mmmOne(ulong _p)
{
   return ((-_p) % _p);
}

// Compute the residual of 2^64 (mod p)
ulong mmmR2(ulong _p, ulong _q, ulong _one)
{
	ulong t = mmmAdd(_one, _one, _p);
   
   t = mmmAdd(t, t, _p);   // 4
	for (size_t i=0; i<5; i++)
      t = mmmMulmod(t, t, _p, _q);   // 4^{2^5} = 2^64
      
	return t;
}

ulong mmmAdd(ulong a, ulong b, ulong _p)
{
   ulong c = (a >= _p - b) ? _p : 0;
   return a + b - c;
}

ulong mmmSub(ulong a, ulong b, ulong _p)
{
   ulong c = (a < b) ? _p : 0;
   return a - b + c;
}

ulong mmmMulmod(ulong a, ulong b, ulong _p, ulong _q)
{
   ulong lo, hi;

#ifdef __NV_CL_C_VERSION
   const uint a0 = (uint)(a), a1 = (uint)(a >> 32);
   const uint b0 = (uint)(b), b1 = (uint)(b >> 32);

   uint c0 = a0 * b0, c1 = mul_hi(a0, b0), c2, c3;

   asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c1) : "r" (a0), "r" (b1), "r" (c1));
   asm volatile ("madc.hi.u32 %0, %1, %2, 0;" : "=r" (c2) : "r" (a0), "r" (b1));

   asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a1), "r" (b1), "r" (c2));
   asm volatile ("madc.hi.u32 %0, %1, %2, 0;" : "=r" (c3) : "r" (a1), "r" (b1));

   asm volatile ("mad.lo.cc.u32 %0, %1, %2, %3;" : "=r" (c1) : "r" (a1), "r" (b0), "r" (c1));
   asm volatile ("madc.hi.cc.u32 %0, %1, %2, %3;" : "=r" (c2) : "r" (a1), "r" (b0), "r" (c2));
   asm volatile ("addc.u32 %0, %1, 0;" : "=r" (c3) : "r" (c3));

   lo = upsample(c1, c0); hi = upsample(c3, c2);
#else
   lo = a * b; hi = mul_hi(a, b);
#endif

   ulong m = lo * _q;
   ulong mp = mul_hi(m, _p);
   long r = (long)(hi - mp);
   return (r < 0) ? r + _p : r;
}

// Compute the residual of n (mod p)
ulong mmmNtoRes(ulong n, ulong _p, ulong _q, ulong _r2)
{
   return mmmMulmod(n, _r2, _p, _q);
}

// Convert a residual back to n
ulong mmmResToN(ulong res, ulong _p, ulong _q)
{
   return mmmMulmod(res, 1, _p, _q);
}

// Compute the residual of b ^ n (mod p)
ulong mmmPowmod(ulong resB, ulong exp, ulong _p, ulong _q, ulong _one)
{
   ulong x = resB;
   ulong y = _one;

   while (true)
   {
      if (exp & 1)
         y = mmmMulmod(x, y, _p, _q);

      exp >>= 1;

      if (!exp)
         return y;

      x = mmmMulmod(x, x, _p, _q);
   }

   // Should never get here
   return 0;
}

ulong  invmod(ulong a, ulong p)
{
   ulong ps1, ps2, q, r, t, dividend, divisor;
   uint parity;

   if (a < 3)
      return (a < 2) ? a : (p+1)/2;

   q = p / a;
   r = p % a;

   dividend = a;
   divisor = r;
   ps1 = q;
   ps2 = 1;
   parity = 0;

   while (divisor > 1)
   {
      r = dividend - divisor;
      t = r - divisor;
      if (r >= divisor) {
         q += ps1; r = t; t -= divisor;
         if (r >= divisor) {
            q += ps1; r = t; t -= divisor;
            if (r >= divisor) {
               q += ps1; r = t; t -= divisor;
               if (r >= divisor) {
                  q += ps1; r = t; t -= divisor;
                  if (r >= divisor) {
                     q += ps1; r = t; t -= divisor;
                     if (r >= divisor) {
                        q += ps1; r = t; t -= divisor;
                        if (r >= divisor) {
                           q += ps1; r = t; t -= divisor;
                           if (r >= divisor) {
                              q += ps1; r = t;
                              if (r >= divisor) {
                                 q = dividend / divisor;
                                 r = dividend % divisor;
                                 q *= ps1;
                              } } } } } } } } }
      q += ps2;
      parity = ~parity;
      dividend = divisor;
      divisor = r;
      ps2 = ps1;
      ps1 = q;
   }

   return (parity) ? ps1 : p - ps1;
}

void collect_factor(uint    n,
                    ulong   p,
  volatile __global uint   *factorsCount,
           __global ulong2 *factors)
{

   uint old = atomic_inc(factorsCount);

   // If we reach the end, stop adding to the buffer.  The CPU code will end
   // with an error as the buffer is not large enough to capture all factors.
   if (old >= D_MAX_FACTORS)
      return;
   
   factors[old].x = n;
   factors[old].y = p;
}