/* Worker.h -- (C) Mark Rodenkirch, February 2012

   This class provides the interface for functions called by each thread
   that is sieving.  There will be one of these for each thread.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _WORKERTHREAD_H
#define _WORKERTHREAD_H

#include "main.h"
#include <vector>

class Worker;

#ifndef WIN32
#include <pthread.h>
#endif

#include "App.h"
#include "SharedMemoryItem.h"

typedef enum { WS_INITIALIZING,
               WS_WAITING_FOR_WORK, // Indicates this thread is initialized and waiting for work
               WS_HAS_WORK_TO_DO,   // Indidates this thread has work and can start working on it
               WS_WORKING,          // Indicates this thread is working
               WS_STOPPING,         // Indicates this thread is requested to stop
               WS_STOPPED           // Indicates this thread is stopped
             } workerstatus_t;

class Worker
{  
public:
   Worker(void) {};
   Worker(uint32_t myId, App *theApp);

   virtual ~Worker(void) = 0;

   // Mini-chunking can be used by worker classes that test in groups of 4 or more.
   // Mega-chunking should be used by worker classes that test individual primes.
   // Use SetMiniChunkRange to switch between mini and mega chunking.
   // The Mega version is deprecated.
   
   virtual void      TestMiniPrimeChunk(uint64_t *miniPrimeChunk) = 0;
   
   virtual void      TestMegaPrimeChunk(void) = 0;
   
   virtual void      CleanUp(void) = 0;
   
   // These are only to be used if calling multiple getters in succession
   void              LockStats(void)    { ip_StatsLocker->Lock(); }
   void              ReleaseStats(void) { ip_StatsLocker->Release(); }

   uint64_t          GetWorkerCpuUS(void)  { return il_WorkerCpuUS; }
   uint64_t          GetPrimesTested(void)    { return il_PrimesTested; }
   uint64_t          GetLargestPrimeTested(void)  { return il_LargestPrimeTested; }
 
   bool              IsInitialized(void) { return ib_Initialized; };
   bool              HasWorkToDo(void) { return (((workerstatus_t) ip_WorkerStatus->GetValueNoLock()) == WS_HAS_WORK_TO_DO); };
   bool              IsWaitingForWork(void) { return (((workerstatus_t) ip_WorkerStatus->GetValueNoLock()) == WS_WAITING_FOR_WORK); };
   bool              IsWorking(void) { return (((workerstatus_t) ip_WorkerStatus->GetValueNoLock()) == WS_WORKING); };
   bool              IsStopped(void) { return (((workerstatus_t) ip_WorkerStatus->GetValueNoLock()) == WS_STOPPED); };
   bool              IsGpuWorker(void) { return ib_GpuWorker; };
   
   void              DetermineWorkSize(void);
   void              WaitForHandOff(void);
   
   void              StopASAP(void) { ip_WorkerStatus->SetValueNoLock(WS_STOPPING); };
   uint64_t          ProcessNextPrimeChunk(uint64_t startFrom, uint64_t maxPrimeForChunk);

protected:
   bool              IsQuadraticResidue(uint64_t n, uint64_t p);
   uint64_t          InvMod32(uint32_t a, uint64_t p);
   uint64_t          InvMod64(uint64_t a, uint64_t p);
   
   // This is called by GPU workers.  It will allocate memory for il_PrimeList
   // and return the number of primes per worker.
   void              AllocatePrimeList(uint32_t workGroupSize);

   // This function will return a boolean indicating if the CPU supports the instructions used
   // by the avx256_xxx.S assembler code which rely on the ymm registers.
   bool              CpuSupportsAvx(void) { return (__builtin_cpu_supports("avx") && __builtin_cpu_supports("fma")); };
   
   // This function will return a boolean indicating if the CPU supports the instructions used
   // by the avx512_xxx.S assembler code which rely on the zmm registers.
   bool              CpuSupportsAvx512(void) { return (__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512vl")); };
   
   void              SetMiniChunkRange(uint64_t minPrimeForMiniChunkMode, uint64_t maxPrimeForMiniChunkMode, uint32_t chunkSize);

   void              SetLargestPrimeTested(uint64_t largestPrimeTested, uint64_t primesTested)  { il_LargestPrimeTested = largestPrimeTested; il_PrimesTested += primesTested; }
   
   uint32_t          ii_MyId;
   bool              ib_Initialized;
   bool              ib_GpuWorker;
   
   vector<uint64_t>  iv_Primes;
   
   // This is only used by GPU workers
   uint64_t         *il_PrimeList;
   
   App              *ip_App;

   SharedMemoryItem *ip_StatsLocker;
   SharedMemoryItem *ip_WorkerStatus;

   uint32_t          ii_WorkSize;
#ifndef WIN32
   pthread_t         ih_Thread;
#endif

private:
   void              SetHasWorkToDo(void) { ip_WorkerStatus->SetValueNoLock(WS_HAS_WORK_TO_DO); };
   void              SetWorking(void) { ip_WorkerStatus->SetValueNoLock(WS_WORKING); };
   void              SetWaitingForWork(void) { ip_WorkerStatus->SetValueNoLock(WS_WAITING_FOR_WORK); };
   void              SetStopped(void) { ip_WorkerStatus->SetValueNoLock(WS_STOPPED); };
   
   void              TestWithMiniChunks(void);
   
   uint32_t          ii_MiniChunkSize;
   uint64_t          il_MinPrimeForMiniChunkMode;
   uint64_t          il_MaxPrimeForMiniChunkMode;
   
   vector<uint64_t>::iterator it_PrimeIterator;
   
   // Total number of milliseconds spent in the thread.
   uint64_t          il_WorkerCpuUS;
   uint64_t          il_PrimesTested;
   uint64_t          il_LargestPrimeTested;
};

#endif

