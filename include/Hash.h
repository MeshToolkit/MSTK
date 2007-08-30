#ifndef _H_HASH
#define _H_HASH

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MSTK_util.h"
#include "MSTK_malloc.h"
#include "List.h"

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
  int  pvtHash_Get_Type(Hash_ptr h); /* 0 for edges, 1 for faces */
  void pvtHash_Set_Pars(Hash_ptr h, unsigned int nent, unsigned int p, int t);
  unsigned int pvtHash_Function(unsigned int np, void* *p);
  unsigned int pvtHash_Enlarge(Hash_ptr h);
  int pvtHash_CheckKeys(unsigned int np1, void* *p1, unsigned int np2, void* *p2);
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

#ifdef DEBUG
  void     Hash_Print(Hash_ptr h);
#endif

#ifdef __cplusplus
}
#endif

#endif

