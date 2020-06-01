/* K1B2Worker.cpp -- (C) Mark Rodenkirch, November 2017

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include <cinttypes>

#include "K1B2Worker.h"
#include "../x86_asm/fpu-asm-x86.h"

K1B2Worker::K1B2Worker(uint32_t myId, App *theApp) : Worker(myId, theApp)
{
   ip_K1B2App = (K1B2App *) theApp;
         
   ii_MinN = ip_K1B2App->GetMinN();
   ii_MaxN = ip_K1B2App->GetMaxN();
   
   il_MinC = ip_K1B2App->GetMinC();
   il_MaxC = ip_K1B2App->GetMaxC();
   
   // The thread can't start until initialization is done
   ib_Initialized = true;
}

void  K1B2Worker::CleanUp(void)
{
}

void  K1B2Worker::TestMegaPrimeChunk(void)
{
   uint64_t p1, p2, p3, p4, ps[4];
   uint64_t twoExpN[4];
   uint64_t maxPrime = ip_App->GetMaxPrime();
   uint32_t n;
   bool     useSmallPLogic = true;
   
   vector<uint64_t>::iterator it = iv_Primes.begin(); 
   
   while (it != iv_Primes.end())
   {
      p1 = ps[0] = *it;
      it++;
      
      p2 = ps[1] = *it;
      it++;
      
      p3 = ps[2] = *it;
      it++;
      
      p4 = ps[3] = *it;
      it++;

      if (useSmallPLogic)
      {
         uint64_t maxC = MAX(abs(il_MinC), abs(il_MaxC));
         if (p4 > maxC)
            useSmallPLogic = false;
      }

      twoExpN[0] = twoExpN[1] = twoExpN[2] = twoExpN[3] = 2;
      
      fpu_powmod_4b_1n_4p(twoExpN, ii_MinN, ps);

      n = ii_MinN;

      while (n <= ii_MaxN)
      {
         if (useSmallPLogic)
         {
            RemoveTermsSmallP(p1, n, twoExpN[0]);
            RemoveTermsSmallP(p2, n, twoExpN[1]);
            RemoveTermsSmallP(p3, n, twoExpN[2]);
            RemoveTermsSmallP(p4, n, twoExpN[3]);
         }
         else
         {
            RemoveTermsLargeP(p1, n, twoExpN[0]);
            RemoveTermsLargeP(p2, n, twoExpN[1]);
            RemoveTermsLargeP(p3, n, twoExpN[2]);
            RemoveTermsLargeP(p4, n, twoExpN[3]);
         }
         
         // Multiple each term by 2.
         twoExpN[0] <<= 1;
         twoExpN[1] <<= 1;
         twoExpN[2] <<= 1;
         twoExpN[3] <<= 1;
    
         // If we have exceeded p, then subtract p.
         if (twoExpN[0] >= ps[0]) twoExpN[0] -= ps[0];
         if (twoExpN[1] >= ps[1]) twoExpN[1] -= ps[1];
         if (twoExpN[2] >= ps[2]) twoExpN[2] -= ps[2];
         if (twoExpN[3] >= ps[3]) twoExpN[3] -= ps[3];
         
         n++;
      }

      SetLargestPrimeTested(p4, 4);
      
      if (p4 > maxPrime)
         break;
   }
}

void  K1B2Worker::TestMiniPrimeChunk(uint64_t *miniPrimeChunk)
{
   FatalError("K1B2Worker::TestMiniPrimeChunk not implemented");
}

void    K1B2Worker::RemoveTermsSmallP(uint64_t prime, uint32_t n, uint64_t twoExpN)
{
   int64_t  c, cc;
   int64_t  sPrime;
      
   c = prime - twoExpN;
   
   // Need this to be signed for the next calculation
   sPrime = (int64_t) prime;
   
   // Compute c such that il_MinC <= c < il_MinC + p
   cc = sPrime * ((il_MinC - c)/sPrime);

   c += cc;
   
   while (c > il_MaxC) c -= prime;
   while (c < il_MinC) c += prime;
   
   if (c > il_MaxC)
      return;

   fpu_push_1divp(prime);
      
   do 
   {
      if (ip_K1B2App->ReportFactor(prime, n, c))
         VerifyFactor(prime, n, c);

      c += prime;
   } while (c <= il_MaxC);
   
   fpu_pop();
}

void    K1B2Worker::RemoveTermsLargeP(uint64_t prime, uint32_t n, uint64_t twoExpN)
{
   int64_t  c, cc;
   int64_t  sPrime;
      
   c = prime - twoExpN;
   
   // Need this to be signed for the next calculation
   sPrime = (int64_t) prime;
   
   // Compute c such that il_MinC <= c < il_MinC + p
   cc = sPrime * ((il_MinC - c)/sPrime);

   while (c > il_MaxC) c -= prime;
   while (c < il_MinC) c += prime;
   
   c += cc;
   
   if (c > il_MaxC)
      return;

   if (ip_K1B2App->ReportFactor(prime, n, c))
   {
      fpu_push_1divp(prime);
      VerifyFactor(prime, n, c);
      fpu_pop();
   }
}

void  K1B2Worker::VerifyFactor(uint64_t prime, uint32_t n, int64_t c)
{
   int64_t  rem;
   int64_t  sPrime;

   rem = fpu_powmod(2, n, prime);
   
   rem += c;

   // Need this to be signed for the next calculation
   sPrime = (int64_t) prime;
      
   rem = (rem % sPrime);

   if (rem < 0)
      rem += sPrime;
   
   if (rem != 0)
      FatalError("2^%u%+" PRId64" mod %" PRIu64" = %" PRIu64"", n, c, prime, rem + 1);
}
