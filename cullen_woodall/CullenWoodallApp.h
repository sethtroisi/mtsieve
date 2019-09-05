/* CullenWoodallApp.h -- (C) Mark Rodenkirch, May 2018

   This class inherits from App.h and has the implementation for this project
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _CullenWoodallApp_H
#define _CullenWoodallApp_H

#include <vector>

#include "../core/AlgebraicFactorApp.h"

class CullenWoodallApp : public AlgebraicFactorApp
{
public:
   CullenWoodallApp(void);

   ~CullenWoodallApp(void) {};

   void              Help(void);
   void              AddCommandLineOptions(string &shortOpts, struct option *longOpts);
   parse_t           ParseOption(int opt, char *arg, const char *source);
   void              ValidateOptions(void);
   bool              ApplyFactor(const char *term);
   void              GetExtraTextForSieveStartedMessage(char *extraText);
   
   bool              ReportFactor(uint64_t p, uint32_t n, int32_t c);
   
   bool              IsCullenSearch() { return ib_Cullen; };
   bool              IsWoodallSearch() { return ib_Woodall; };
   
   int32_t           GetBase(void) { return ii_Base; };
   int32_t           GetMinN(void) { return ii_MinN; };
   int32_t           GetMaxN(void) { return ii_MaxN; };
   int32_t           GetGpuSteps(void) { return ii_GpuSteps; };
   
   uint64_t          GetTermCount(void) { return il_TermCount; };
   
   uint32_t          GetTerms(uint32_t *terms, uint32_t minGroupSize, uint32_t maxGroupSize);
   
protected:
   void              PreSieveHook(void) {};
   bool              PostSieveHook(void) { return true; };
   
   void              NotifyAppToRebuild(void) {};
   
   Worker           *CreateWorker(uint32_t id, bool gpuWorker, uint64_t largestPrimeTested);
   
   void              ProcessInputTermsFile(bool haveBitMap);
   bool              IsWritingOutputTermsFile(void){ return true; };
   void              WriteOutputTermsFile(uint64_t largestPrime);

private:
   void              SetInitialTerms(void);
   
   void              EliminateGfnAndMersenneTerms(void);
   void              EliminateTermsWithAlgebraicFactors(void);
   bool              CheckAlgebraicFactor(uint32_t n, int32_t c, const char *fmt, ...);
   
   vector<bool>      iv_CullenTerms;
   vector<bool>      iv_WoodallTerms;
      
   uint32_t          ii_Base;
   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   bool              ib_Woodall;
   bool              ib_Cullen;
   
   uint32_t          ii_GpuSteps;
};

#endif

