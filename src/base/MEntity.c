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
    if (ent->dim < 0)
      return MDELETED;
    else
      return ent->dim;
  }

  int MEnt_OrigDim(MEntity_ptr ent) {
    if (ent->dim < 0)
      return -(ent->dim/10)-1;
    else {
      MSTK_Report("MEnt_OrigDim","This is not a deleted entity",WARN);
      return (ent->dim);
    }
  }	

  int MEnt_ID(MEntity_ptr ent) {
    return (ent->id);
  }

  void MEnt_Delete(MEntity_ptr ent, int keep) {
    int dim;

    dim = MEnt_Dim(ent);
    if (dim < 0)
      dim = MEnt_OrigDim(ent);

    switch (dim) {
    case MVERTEX:
      MV_Delete(ent,keep);
      break;
    case MEDGE:
      ME_Delete(ent,keep);
      break;
    case MFACE:
      MF_Delete(ent,keep);
      break;
    case MREGION:
      MR_Delete(ent,keep);
      break;
    }
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
