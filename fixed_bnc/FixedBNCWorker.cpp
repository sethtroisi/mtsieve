/* FixedBNCWorker.cpp -- (C) Mark Rodenkirch, November 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>
#include <stdint.h>

#include "FixedBNCWorker.h"
#include "../x86_asm/fpu-asm-x86.h"

FixedBNCWorker::FixedBNCWorker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   ip_FixedBNCApp = (FixedBNCApp *) theApp;

   il_MinK = ip_FixedBNCApp->GetMinK();
   il_MaxK = ip_FixedBNCApp->GetMaxK();
   ii_Base = ip_FixedBNCApp->GetBase();
   ii_N = ip_FixedBNCApp->GetN();
   ii_C = ip_FixedBNCApp->GetC();

   ii_BaseInverses = 0;

   DeterminePrimeTermRange();

   BuildBaseInverses();

   il_PrimeList = (uint64_t *) xmalloc((ip_FixedBNCApp->GetCpuWorkSize() + 10) * sizeof(uint64_t));
   ii_InverseList = (uint32_t *) xmalloc((ip_FixedBNCApp->GetCpuWorkSize() + 10) * sizeof(uint32_t));

   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  FixedBNCWorker::CleanUp(void)
{
}

void  FixedBNCWorker::TestMegaPrimeChunk(void)
{
   uint64_t p1 = 0, p2, p3, p4;
   uint64_t k1, k2, k3, k4, ks[4];
   uint64_t maxPrime = ip_App->GetMaxPrime();
   int32_t  svb = 0;
   int32_t  pmb, count, idx;

   // Evaluate primes in the vector to determine if can yield a factor.  Only
   // put primes that can yield a factor into an array for the second loop.
   count = 0;

   for (uint32_t pIdx=0; pIdx<ii_PrimesInList; pIdx++)
   {
      p1 = il_PrimeList[pIdx];

      pmb = (p1 % ii_Base);

      if (ii_BaseInverses[pmb] == 0)
         continue;

      if (p1 > maxPrime)
         break;

      svb = ii_BaseInverses[pmb];

      il_PrimeList[count] = p1;
      ii_InverseList[count] = svb;
      count++;
   }

   if (count == 0)
      return;

   // Duplicate the last few entries so that the
   // number of valid entries is divisible by 4.
   while (count % 4 != 0)
   {
      il_PrimeList[count] = p1;
      ii_InverseList[count] = svb;
      count++;
   }

   for (idx=0; idx<count; idx+=4)
   {
      p1 = il_PrimeList[idx+0];
      p2 = il_PrimeList[idx+1];
      p3 = il_PrimeList[idx+2];
      p4 = il_PrimeList[idx+3];

      ks[0] = k1 = (1+ii_InverseList[idx+0]*p1)/ii_Base;
      ks[1] = k2 = (1+ii_InverseList[idx+1]*p2)/ii_Base;
      ks[2] = k3 = (1+ii_InverseList[idx+2]*p3)/ii_Base;
      ks[3] = k4 = (1+ii_InverseList[idx+3]*p4)/ii_Base;

      // Starting with k*2^n = 1 (mod p)
      //           --> k = (1/2)^n (mod p)
      //           --> k = inverse^n (mod p)
      fpu_powmod_4b_1n_4p(ks, ii_N, &il_PrimeList[idx+0]);

      if (ii_C == +1)
      {
         k1 = p1 - ks[0];
         k2 = p2 - ks[1];
         k3 = p3 - ks[2];
         k4 = p4 - ks[3];
      }
      else
      {
         k1 = ks[0];
         k2 = ks[1];
         k3 = ks[2];
         k4 = ks[3];
      }

      if (p1 <= il_MaxK)
      {
         if (k1 <= il_MaxK) RemoveTermsSmallPrime(p1, k1);
         if (k2 <= il_MaxK) RemoveTermsSmallPrime(p2, k2);
         if (k3 <= il_MaxK) RemoveTermsSmallPrime(p3, k3);
         if (k4 <= il_MaxK) RemoveTermsSmallPrime(p4, k4);
      }
      else
      {
         if (k1 <= il_MaxK) RemoveTermsBigPrime(p1, k1);
         if (k2 <= il_MaxK) RemoveTermsBigPrime(p2, k2);
         if (k3 <= il_MaxK) RemoveTermsBigPrime(p3, k3);
         if (k4 <= il_MaxK) RemoveTermsBigPrime(p4, k4);
      }

      SetLargestPrimeTested(p4, 4);
   }

   // Adjust for the possibility that we tested the same prime
   // more than once at the end of the list.
   if (p4 == p3)
      SetLargestPrimeTested(p4, -1);
   if (p4 == p2)
      SetLargestPrimeTested(p4, -1);
   if (p4 == p1)
      SetLargestPrimeTested(p4, -1);
}

void  FixedBNCWorker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("FixedBNCWorker::TestMiniPrimeChunk not implemented");
}

// This must be used when prime <= il_MaxK.
void    FixedBNCWorker::RemoveTermsSmallPrime(uint64_t prime, uint64_t k)
{
   double    tooBig = 0.0;
   uint64_t  maxPrime;
   uint64_t  bPowNModP;

   if (il_BpowN > 0)
   {
      maxPrime = ip_FixedBNCApp->GetMaxPrime();

      if (maxPrime == 0)
         tooBig = (double) KMAX_MAX;
      else
         tooBig = (double) maxPrime;

      tooBig -= 1.0;
   }

   // Make sure that k >= il_MinK
   if (k < il_MinK)
   {
      if (prime >= il_MinK)
         k += prime;
      else
      {
         // Compute k such that il_MinK <= k < il_MinK + p
         k += prime * ((il_MinK - k + prime - 1)/prime);
      }
   }

   if (k > il_MaxK)
      return;

   fpu_push_1divp(prime);

   bPowNModP = fpu_powmod(ii_Base, ii_N, prime);

   do
	{
      if (il_BpowN > 0)
      {
         double kBpowN = (double) il_BpowN;
         kBpowN *= (double) k;

         if (kBpowN < tooBig)
         {
            // If k*b^n+c is prime, then we won't treat this as a factor
            if (k * il_BpowN != prime - ii_C)
               if (ip_FixedBNCApp->ReportFactor(prime, k))
                  VerifyFactor(prime, k, bPowNModP);
         }
         else
            if (ip_FixedBNCApp->ReportFactor(prime, k))
               VerifyFactor(prime, k, bPowNModP);

      }
      else
      {
         if (ip_FixedBNCApp->ReportFactor(prime, k))
            VerifyFactor(prime, k, bPowNModP);
      }

		k += prime;
	} while (k <= il_MaxK);

   fpu_pop();
}

// Using this bypasses a number of if checks that can be done when prime > il_MaxK.
void    FixedBNCWorker::RemoveTermsBigPrime(uint64_t prime, uint64_t k)
{
   uint64_t  bPowNModP;

   // Make sure that k >= il_MinK
   if (k < il_MinK)
      return;

   if (ip_FixedBNCApp->ReportFactor(prime, k))
   {
      fpu_push_1divp(prime);

      bPowNModP = fpu_powmod(ii_Base, ii_N, prime);

      VerifyFactor(prime, k, bPowNModP);

      fpu_pop();
   }
}

void  FixedBNCWorker::VerifyFactor(uint64_t prime, uint64_t k, uint64_t bPowNModP)
{
   uint64_t rem;

   rem = fpu_mulmod(bPowNModP, k, prime);

   if (ii_C == +1 && rem != prime - 1)
      FatalError("%" PRIu64"*%u^%u%+d mod %" PRIu64" = %" PRIu64"", k, ii_Base, ii_N, ii_C, prime, rem + 1);

   if (ii_C == -1 && rem != 1)
      FatalError("%" PRIu64"*%u^%u%+d mod %" PRIu64" = %" PRIu64"", k, ii_Base, ii_N, ii_C, prime, rem - 1);
}

void    FixedBNCWorker::BuildBaseInverses(void)
{
	ii_BaseInverses = (uint32_t *) xmalloc(ii_Base * sizeof(uint32_t));

	ii_BaseInverses[0] = 0;
	ii_BaseInverses[1] = ii_Base-1;

	for (uint32_t i = 2; i < ii_Base; i++)
		ii_BaseInverses[i] = EuclidExtendedGCD(i, ii_Base);
}

uint32_t FixedBNCWorker::EuclidExtendedGCD(uint32_t a, uint32_t base)
{
	int64_t u = 0, d = base, v1 = 1, v3 = a;
   int64_t q, t1, t3;

	while (v3)
	{
		q = d/v3;
		t3 = d - (q*v3);
		t1 = u - (q*v1);

		u = v1;
		d = v3;

		v1 = t1;
		v3 = t3;
	}

	return (d == 1) ? (uint32_t)((u>0)?base-u:-u) : 0;
}

// Determine the minimum prime between kmin*b^n+c and kmax*b^n+c so
// that we will avoid removing those terms in that range that are prime
void  FixedBNCWorker::DeterminePrimeTermRange(void)
{
   uint64_t  maxPrime = ip_FixedBNCApp->GetMaxPrime();
   double    tooBig = (double) KMAX_MAX;
   double    b = (double) ii_Base;
   double    bpown = 1.0;
   double    mink = (double) il_MinK;

   // KMAX_MAX is the same as PMAX_MAX, 2^62.
   if (maxPrime == 0)
      tooBig = (double) KMAX_MAX;
   else
      tooBig = (double) maxPrime;

   il_BpowN = 0;

   for (uint32_t i=0; i<ii_N; i++)
   {
      bpown *= b;

      // If mink*b^n-1.0 > 2^62, then we don't need to worry about removing terms
      // that are prime.
      if (mink*bpown-1.0 > tooBig)
         return;
   }

   il_BpowN = (uint64_t) bpown;
}

