/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#define _H_MEdge_Private

#include <stdlib.h>
#include "MEdge.h"
#include "MEdge_jmp.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  void ME_Set_RepType_R1(MEdge_ptr e) {
    MEdge_Adj_RN *adj;

    adj = e->adj = (MEdge_Adj_RN *) malloc(sizeof(MEdge_Adj_RN));
    adj->hnext = NULL;
    adj->lock = 0;
  }

  void ME_Delete_R1(MEdge_ptr e, int keep) {
    MEdge_Adj_RN *adj;

    if (!keep) {
      adj = (MEdge_Adj_RN *) e->adj;
      if (adj) {
	free(adj);
      }
      MEnt_Set_DelFlag((MEntity_ptr) e);
    }
    return;
  }

  void ME_Destroy_For_MESH_Delete_R1(MEdge_ptr e) {
    return;
  }

  int ME_Num_Faces_R1(MEdge_ptr e) {
    return ME_Num_Faces_R1R2(e);
  }

  List_ptr ME_Faces_R1(MEdge_ptr e) {
    return ME_Faces_R1R2(e);
  }

  int ME_Num_Regions_R1(MEdge_ptr e) {
    return ME_Num_Regions_RN(e);
  }

  List_ptr ME_Regions_R1(MEdge_ptr e) {
    return ME_Regions_RN(e);
  }

  void ME_Add_Face_R1(MEdge_ptr e, MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("ME_Add_Face",
		"Function call not suitable for this representation",MSTK_WARN);
#endif
  }
  
  void ME_Rem_Face_R1(MEdge_ptr e, MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("ME_Rem_Face",
		"Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void ME_Add_Region_R1(MEdge_ptr e, MRegion_ptr r) {
#ifdef DEBUG
    MSTK_Report("ME_Add_Region",
		"Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void ME_Rem_Region_R1(MEdge_ptr e, MRegion_ptr r) {
#ifdef DEBUG
    MSTK_Report("ME_Rem_Region",
		"Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  MEdge_ptr ME_NextInHash_R1(MEdge_ptr e) {
    MEdge_ptr hnext = ((MEdge_Adj_RN *)e->adj)->hnext; 
    return hnext;
  }

  void ME_Set_NextInHash_R1(MEdge_ptr e, MEdge_ptr next) {
    MEdge_Adj_RN *adj = (MEdge_Adj_RN *)e->adj;
    adj->hnext = next;
  }

  void ME_Lock_R1(MEdge_ptr e) {
    MEdge_Adj_RN *adj = (MEdge_Adj_RN *)e->adj;
    Hash_Lock(&adj->lock);
  }

  void ME_UnLock_R1(MEdge_ptr e) {
    MEdge_Adj_RN *adj = (MEdge_Adj_RN *)e->adj;
    Hash_UnLock(&adj->lock);
  }

  int ME_IsLocked_R1(MEdge_ptr e) {
    MEdge_Adj_RN *adj = (MEdge_Adj_RN *)e->adj;
    return Hash_IsLocked(adj->lock);
  }

#ifdef __cplusplus
}
#endif

