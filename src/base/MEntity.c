#define _H_MEntity_Private

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MEntity.h"
#include "List.h"
#include "MSTK_malloc.h"
#include "MAttrib.h"
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
#ifdef DEBUG
    if (ent->dim < MVERTEX && ent->dim > MREGION)
      MSTK_Report("MEnt_Unmark","Not a valid topological entity",ERROR);
#endif

    ent->marker = ent->marker | 1<<(markerID-1);
  }

  void MEnt_Unmark(MEntity_ptr ent, int markerID) {
#ifdef DEBUG
    if (ent->dim < MVERTEX && ent->dim > MREGION)
      MSTK_Report("MEnt_Unmark","Not a valid topological entity",ERROR);
#endif

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


  /* Set value of attribute */
  void MEnt_Set_AttVal(MEntity_ptr ent, MAttrib_ptr attrib, int ival, double lval, void *pval) {
    int idx, found;
    MAttIns_ptr attins;
    List_ptr attinslist;
    
    if (!ent->AttInsList)
      ent->AttInsList = List_New(3);
      
    attinslist = ent->AttInsList;
    idx = 0; found = 0;
    while ((attins = List_Next_Entry(attinslist,&idx))) {
      if (MAttIns_Attrib(attins) == attrib) {
	found = 1;
	break;
      }
    }

    if (!found) {
      attins = MAttIns_New(attrib);
      List_Add(attinslist,attins);
    }
 
    MAttIns_Set_Value(attins, ival, lval, pval);
  }

  /* Clear value of attribute */
  void MEnt_Rem_AttVal(MEntity_ptr ent, MAttrib_ptr attrib) {
    int i, idx, found;
    MAttIns_ptr attins;
    List_ptr attinslist;
    
    if (!ent->AttInsList)
      return;
    
    attinslist = ent->AttInsList;
    idx = 0; i = 0; found = 0;
    while ((attins = List_Next_Entry(attinslist,&idx))) {
      if (MAttIns_Attrib(attins) == attrib) {
	found = 1;
	break;
      }
      else
	i++;
    }

    if (!found)
      return;

    List_Remi(attinslist,i);
    MAttIns_Delete(attins);
  }


  /* Query the value of the attribute */
  int MEnt_Get_AttVal(MEntity_ptr ent, MAttrib_ptr attrib, int *ival, double *lval, void **pval) {
    int idx, found;
    MAttIns_ptr attins;
    List_ptr attinslist;
    
    if (ival) *ival = 0;
    if (lval) *lval = 0;
    if (pval) *pval = NULL;
    
    if (!ent->AttInsList)
      return 0;
      
    attinslist = ent->AttInsList;
    idx = 0; found = 0;
    while ((attins = List_Next_Entry(attinslist,&idx))) {
      if (MAttIns_Attrib(attins) == attrib) {
	found = 1;
	break;
      }
    }

    if (!found)
      return 0;
 
    MAttIns_Get_Value(attins, ival, lval, pval);

    return 1;
  }
  
    

#ifdef __cplusplus
}
#endif
