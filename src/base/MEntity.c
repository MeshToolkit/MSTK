#define _H_MEntity_Private

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MEntity.h"
#include "List.h"
#include "MSTK_malloc.h"
#include "MSTK.h"

#ifdef __cplusplus
extern "C" {
#endif

  int MEnt_Dim(MEntity_ptr ent) {
    return (ent->dim);
  }

  int MEnt_ID(MEntity_ptr ent) {
    return (ent->id);
  }

  int MEnt_IsMarked(MEntity_ptr ent, int markerID) {
    return (ent->marker & 1<<(markerID-1));
  }

  void MEnt_Mark(MEntity_ptr ent, int markerID) {
    ent->marker = ent->marker | 1<<(markerID-1);
  }

  void MEnt_Unmark(MEntity_ptr ent, int markerID) {
    ent->marker = ent->marker & ~(1<<(markerID-1));
  }

  void List_Mark(List_ptr list, int markerID) {
    MEntity_ptr ent;
    int i, n = List_Num_Entries(list);
    
    for (i = 0; i < n; i++) {
      ent = (MEntity_ptr) List_Entry(list,i);
      MEnt_Mark(ent,markerID);
    }
  }

  void List_Unmark(List_ptr list, int markerID) {
    MEntity_ptr ent;
    int i, n = List_Num_Entries(list);
    
    for (i = 0; i < n; i++) {
      ent = (MEntity_ptr) List_Entry(list,i);
      MEnt_Unmark(ent,markerID);
    }
  }


#ifdef __cplusplus
}
#endif
