/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#define _H_MEdge_Private

#include "MEdge.h"
#include "MEdge_jmp.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  int ME_Num_Regions_RN(MEdge_ptr e) {
    int ne, idx;
    MRegion_ptr reg;
    List_ptr vregs;

    ne = 0;

    vregs = MV_Regions(e->vertex[0]);
    idx = 0;
    while ((reg = List_Next_Entry(vregs,&idx))) {
      if (MR_UsesEntity(reg,e->vertex[1],MVERTEX))
	ne++;
    }

    return ne;
  }

  List_ptr ME_Regions_RN(MEdge_ptr e) {
    int idx;
    MRegion_ptr reg;
    List_ptr eregs, vregs;

    eregs = List_New(0);

    vregs = MV_Regions(e->vertex[0]);
    idx = 0;
    while ((reg = List_Next_Entry(vregs,&idx))) {
      if (MR_UsesEntity(reg,e->vertex[1],MVERTEX))
	List_Add(eregs,reg);
    }

    if (List_Num_Entries(eregs))
      return eregs;    
    else {
      List_Delete(eregs);
      return NULL;
    }
  }

  void ME_RegionIDs_RN(MEdge_ptr e, int *ner, int *eregids) {
    int i;
    List_ptr eregs = ME_Regions_RN(e);
    if (eregs) {
      *ner = List_Num_Entries(eregs);
      for (i = 0; i < *ner; i++)
        eregids[i] = MEnt_ID(List_Entry(eregs,i));
    }
    else
      *ner = 0;
  }

#ifdef __cplusplus
}
#endif

