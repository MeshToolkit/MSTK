#ifndef _H_LIST
#define _H_LIST

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MSTK_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _H_List_Private
  typedef struct List {
    unsigned int nentdat;
    unsigned int *remdat;
    void **entry;    
  } List, *List_ptr;

  void pvtList_Get_Pars(List_ptr l, int *nent, int *p, int *nrem, int *rem1);
  int  pvtList_Get_Nent(List_ptr l);
  void pvtList_Set_Pars(List_ptr l, int nent, int p, int nrem, int rem1);
#else
  typedef void *List_ptr;
#endif

  List_ptr List_New(int inisize);
  void List_Delete(List_ptr list);
  List_ptr List_Compress(List_ptr list);
  List_ptr List_Copy(List_ptr list);
  
  List_ptr List_Add(List_ptr l, void *entry);
  List_ptr List_ChknAdd(List_ptr l, void *entry);
  List_ptr List_Insert(List_ptr l, void *nuentry, void *b4entry);
  List_ptr List_Inserti(List_ptr l, void *nuentry, int i);
  int      List_Rem(List_ptr l, void *entry);
  int      List_Remi(List_ptr l, int i);
  int      List_Remi_Raw(List_ptr l, int i);
  int      List_RemSorted(List_ptr l, void *entry, int (*entry2int)(void *));
  int      List_Replace(List_ptr l, void *entry, void *nuentry);
  int      List_Replacei(List_ptr l, int i, void *nuentry);
  int      List_Contains(List_ptr l, void *entry);
  int      List_Locate(List_ptr l, void *entry);
  void    *List_Entry(List_ptr l, int i);
  void    *List_Entry_Raw(List_ptr l, int i);
  void    *List_Next_Entry(List_ptr l, int *i);
  int      List_Num_Entries(List_ptr l);
  int      List_Num_Entries_Raw(List_ptr l);
  List_ptr List_Cat(List_ptr dest, List_ptr src);
  void     List_Sort(List_ptr l, size_t num, size_t size,
		     int(*comp)(const void *,const void *));
  void    *List_Search(List_ptr l, const void *key, size_t num, size_t size,
                       int(*comp)(const void *,const void *));

#ifdef DEBUG
  void     List_Print(List_ptr l);
#endif

  /* Extra functionality for hash-tables */

  void*   *List_Entries(List_ptr l);


#ifdef __cplusplus
}
#endif

#endif
