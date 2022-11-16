/* MultiFactorialGpuWorker.h -- (C) Mark Rodenkirch, September 2016

   This is the header for the class sets up the call to the fsieve GPU function
   and parses the output from the GPU to determine if we have a factor.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _MultiFactorialGpuWorker_H
#define _MultiFactorialGpuWorker_H

#include "MultiFactorialApp.h"
#include "../core/Worker.h"

#include "../core/GpuKernel.h"

class MultiFactorialGpuWorker : public Worker
{
public:
   MultiFactorialGpuWorker(uint32_t myId, App *theApp);

   ~MultiFactorialGpuWorker() {};

   void              TestMegaPrimeChunk(void);
   void              TestMiniPrimeChunk(uint64_t *miniPrimeChunk);
   void              CleanUp(void);

protected:
   MultiFactorialApp *ip_MultiFactorialApp;

   uint32_t          ii_MinN;
   uint32_t          ii_MaxN;
   uint32_t          ii_MaxGpuSteps;
   uint32_t          ii_MultiFactorial;
   uint32_t          ii_MaxGpuFactors;

   // Variables used by the kernel
   uint64_t         *il_RemainderList;
   uint32_t         *ii_Parameters;
   uint32_t         *ii_FactorCount;
   int64_t          *il_FactorList;

   GpuKernel        *ip_Kernel;
};

#endif

