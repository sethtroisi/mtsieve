/* HashTable.h -- (C) Mark Rodenkirch, January 2018

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _HASHTABLE_H
#define _HASHTABLE_H

#define HASH_NOT_FOUND     UINT32_MAX
#define HASH_MASK1         (1<<15)
#define HASH_MASK2         (HASH_MASK1-1)
#define HASH_MAX_ELTS      HASH_MASK2-1
#define HASH_MINIMUM_SHIFT 10

class HashTable
{
public:
   HashTable(uint32_t elements);

   ~HashTable(void);

   void  Cleanup(void);

   inline void Clear(void)
   {
      for(uint32_t i=0; i < hsize; i++) {
          htable[i] = empty_slot;
      }
      // will never match a Lookup value
      BJ64[empty_slot] = UINT64_MAX;
   }

   inline uint64_t get(uint32_t x) {return BJ64[x]; };

   inline void Insert(uint64_t bj, uint32_t j)
   {
      uint32_t slot;

      BJ64[j] = bj;
      slot = bj & hsize_minus1;
      if (htable[slot] == empty_slot)
         htable[slot] = j;
      else
      {
         olist[j] = htable[slot];
         htable[slot] = (j | HASH_MASK1);
      }
   };

   inline uint32_t Lookup(uint64_t bj)
   {
      uint32_t slot;
      uint16_t elt, elt_low;

      slot = bj & hsize_minus1;
      elt = htable[slot];
      elt_low = elt & HASH_MASK2;

      //if (elt == empty_slot)
      //   return HASH_NOT_FOUND;

      if (BJ64[elt_low] == bj)
         return elt_low;

      while (elt != elt_low) {
         elt = olist[elt_low];
         elt_low = elt & HASH_MASK2;
         if (BJ64[elt_low] == bj)
            return elt_low;
      }

      return HASH_NOT_FOUND;
   };

private:
   /**
    * `hsize` is the size of the hashtable (always a power of two).
    * `hsize_minus1` is a convience constant for MOD hsize.
    * `empty_slot` is the empty slot marker (e.g. nothing here)
    *
    * `BJ64` stores the inserted values
    *
    * `htable` tracks status of each slot:
    *   low bits track where in BJ64 the associated value is held.
    *   upper bit tracks if any item was moved (AKA "keep probing").
    *
    * `olist` [o]pen addressing list, when a conflict happens in slot s
    *     store probing continue point, olist[j] = htable[s]
    *     updates htable[s] = j
    *
    *     At lookup time if probing needs to continue, look in olist[s] to find the old htable[s]
    *     Can keep doing elt = olist[elt] until you find an elemnt that doesn't need any more probing
    */
   uint16_t *htable;
   uint16_t *olist;
   uint64_t *BJ64;
   uint32_t  hsize;
   uint32_t  hsize_minus1;
   uint16_t  empty_slot;
};

#endif
