/* AFSieveApp.h -- (C) Mark Rodenkirch, November 2016

   This class inherits from App.h and has the implementation for this project
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _GFNDivisorApp_H
#define _GFNDivisorApp_H

#include <gmp.h>
#include "../core/FactorApp.h"

class GFNDivisorApp : public FactorApp
{
public:
   GFNDivisorApp(void);

   ~GFNDivisorApp(void) {};

   void              Help(void);
   void              AddCommandLineOptions(string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText);
   
   uint64_t          GetMinK(void) { return il_MinK; };
   uint64_t          GetMaxK(void) { return il_MaxK; };
   uint64_t          GetKCount(void) { return (il_MaxK - il_MinK + 1); };
   uint32_t          GetMinN(void) { return ii_MinN; };
   uint32_t          GetMaxN(void) { return ii_MaxN; };
   uint32_t          GetNCount(void) { return (ii_MaxN - ii_MinN + 1); };

   bool              ReportFactor(uint64_t p, uint64_t k, uint32_t n);

protected:
   void              PreSieveHook(void);
   bool              PostSieveHook(void);
   
   void              NotifyAppToRebuild(void) {};
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);
   
   void              ProcessInputTermsFile(bool haveBitMap);
   void              ProcessInputTermsFile(bool haveBitMap, FILE *fPtr, char *fileName, bool firstFile);
   void              WriteOutputTermsFile(uint64_t largestPrime);
   bool              IsWritingOutputTermsFile(void){ return !ib_TestTerms; };
   uint64_t          WriteABCDTermsFile(char *fileName, uint32_t minN, uint64_t maxPrime);

private:
   void              TestRemainingTerms(void);
   bool              IsFermatDivisor(uint64_t k, uint32_t n);
   void              CheckRedc(mp_limb_t *xp, uint32_t xn, uint32_t b, uint32_t m, uint64_t k, uint32_t n);
   
   vector<vector<bool>>  iv_Terms;
   string            is_OutputTermsFilePrefix;
   
   uint32_t          ii_NsPerFile;
   uint64_t          il_MinK;
   uint64_t          il_MaxK;
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   
   bool              ib_TestTerms;
   uint64_t          il_KPerChunk;
   uint32_t          ii_NPerChunk;
   
   uint64_t          il_MinKOriginal;
   uint64_t          il_MaxKOriginal;
   uint64_t          il_MinKInChunk;
   uint64_t          il_MaxKInChunk;
   
   uint32_t          ii_MinNOriginal;
   uint32_t          ii_MaxNOriginal;
   uint32_t          ii_MinNInChunk;
   uint32_t          ii_MaxNInChunk;
   
   uint64_t          il_StartSievingUS;
   uint64_t          il_TotalTerms;
   uint64_t          il_TotalTermsEvaluated;
   uint64_t          il_TotalTermsInChunk;
};

#endif

