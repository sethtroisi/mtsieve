/* GenericWorker.h -- (C) Mark Rodenkirch, October 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _GenericWorker_H
#define _GenericWorker_H

#include "SierpinskiRieselApp.h"
#include "AbstractSequenceHelper.h"
#include "AbstractWorker.h"
#include "../core/HashTable.h"
#include "../core/MpArithVector.h"

using namespace std;

class GenericWorker : public AbstractWorker
{
public:
   GenericWorker(uint32_t myId, App *theApp, AbstractSequenceHelper *appHelper);

   ~GenericWorker(void) {};

   void              Prepare(uint64_t largestPrimeTested, uint32_t bestQ);

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:

private:
   void              InitializeWorker(void);

   uint64_t          ProcessSmallPrimes(void);

   void              SetupDiscreteLog(uint32_t *b, uint64_t *p, MpArithVec mp, MpResVec mb);

   void              DiscreteLogSmallPrimes(uint32_t *b, uint64_t *p);
   void              DiscreteLogLargePrimes(uint32_t *b, uint64_t *p);

   void              BabySteps(MpArithVec mp, MpResVec mb, uint32_t *orderOfB);

   bool              ib_CanUseCIsOneLogic;
   uint64_t          il_MaxK;

   uint32_t          ii_SieveLow;
   uint32_t          ii_SieveRange;

   uint64_t          il_SmallPrimeSieveLimit;

   uint32_t          ii_BabySteps;
   uint32_t          ii_GiantSteps;

   HashTable        *ip_HashTable[4];

   uint32_t         *ssHash;        // there is one per subsequence

   MpResVec         *mBD;           // there is one set of 4 per Q
   MpResVec         *mBDCK;         // there is one set of 4 per subsequence

   MpResVec          mBM;
};

#endif

