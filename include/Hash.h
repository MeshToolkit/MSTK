/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#ifndef _H_HASH
#define _H_HASH

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MSTK_types.h"
#include "MSTK_defines.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _H_Hash_Private
  typedef struct Hash {
    unsigned int nentdat;
    void **entry;
  } Hash, *Hash_ptr;

  void pvtHash_Get_Pars(Hash_ptr h, unsigned int *nent, unsigned int *p, int *t);
  unsigned int  pvtHash_Get_Nent(Hash_ptr h);
  unsigned int  pvtHash_Get_Size(Hash_ptr h);
  void pvtHash_Set_Pars(Hash_ptr h, unsigned int nent, unsigned int p, int t);
  unsigned int pvtHash_Function(unsigned int np, void* *p);
  unsigned int pvtHash_Enlarge(Hash_ptr h);
  int pvtHash_CheckKeys(unsigned int np1, void* *p1, unsigned int np2, void* *p2);
  void pvtHash_CheckSize(Hash_ptr h);
#else
  typedef void *Hash_ptr;
#endif

  Hash_ptr Hash_New(unsigned int inisize, int type);
  void Hash_Delete(Hash_ptr h);
  
  Hash_ptr Hash_Add(Hash_ptr h, void *entry, unsigned int np, void* *p);
  Hash_ptr Hash_ChknAdd(Hash_ptr h, void *entry, unsigned int np, void* *p);
  int      Hash_Rem(Hash_ptr h, void *entry, unsigned int np, void* *p);
  void    *Hash_Entry(Hash_ptr h, unsigned int np, void* *p);
  int      Hash_Num_Entries(Hash_ptr h);
  List_ptr Hash_Entries(Hash_ptr h);

  void     Hash_Print(Hash_ptr h);

  void Hash_Lock(int *plock);
  void Hash_UnLock(int *plock);
  int Hash_IsLocked(int lock);

  int  Hash_AutoRemove(Hash_ptr h);
  void Hash_Set_AutoRemove(Hash_ptr h, int t);
  unsigned int Hash_Remove_Unused(Hash_ptr h);

#ifdef __cplusplus
}
#endif

#endif

