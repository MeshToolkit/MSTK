/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#define _H_MFace_Private

#include <stdlib.h>
#include "MFace.h"
#include "MFace_jmp.h"
#include "MSTK_private.h"
 
#ifdef __cplusplus
extern "C" {
#endif


  void MF_Set_RepType_R4(MFace_ptr f) {
    MFace_Adj_R4 *adj;

    adj = f->adj = (MFace_Adj_R4 *) malloc(sizeof(MFace_Adj_R4));
    adj->fregions[0] = (MRegion_ptr) NULL;
    adj->fregions[1] = (MRegion_ptr) NULL;
    adj->fvertices = NULL;
  }

  void MF_Delete_R4(MFace_ptr f, int keep) {
    MFace_Adj_R4 *adj;
    int idx;
    MVertex_ptr v;

    adj = (MFace_Adj_R4 *) f->adj;

    if (MEnt_Dim((MEntity_ptr) f) != MDELETED) { /* if face hasnt been temporarily deleted */
      if (adj) {
	idx = 0;
	while ((v = List_Next_Entry(adj->fvertices,&idx))) 
	  MV_Rem_Face(v,f);
      }
    }

    if (!keep) {
      if (adj) {
	if (adj->fvertices)
	  List_Delete(adj->fvertices);
	free(adj);
      }
    }
  }

  void MF_Restore_R4(MFace_ptr f) {
    MFace_Adj_R4 *adj;
    int idx;
    MVertex_ptr v;

    MEnt_Set_Dim((MEntity_ptr) f,MFACE);

    adj = (MFace_Adj_R4 *) f->adj;
    idx = 0;
    while ((v = List_Next_Entry(adj->fvertices,&idx))) 
      MV_Add_Face(v,f);
  }
  
  void MF_Destroy_For_MESH_Delete_R4(MFace_ptr f) {
    MFace_Adj_R4 *adj;

    adj = (MFace_Adj_R4 *) f->adj;
    if (adj) {
      if (adj->fvertices)
	List_Delete(adj->fvertices);
      free(adj);
    }
  }

  void MF_Set_Edges_R4(MFace_ptr f, int n, MEdge_ptr *e, int *dir) {
#ifdef DEBUG
    MSTK_Report("MF_Set_Edges",
		"Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void MF_Rem_Edge_R4(MFace_ptr f, MEdge_ptr e) {
#ifdef DEBUG
    MSTK_Report("MF_Rem_Edge_R4",
		"Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void MF_Replace_Edges_i_R4(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges) {
#ifdef DEBUG
    MSTK_Report("MF_Replace_Edges",
		"Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void MF_Replace_Edges_R4(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges) {
#ifdef DEBUG
    MSTK_Report("MF_Replace_Edge",
		"Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  int MFs_AreSame_R4(MFace_ptr f1, MFace_ptr f2) {
    return (f1 == f2);
  }

  int MF_Num_AdjFaces_R4(MFace_ptr f) {
    MSTK_Report("MF_Num_AdjFaces_R4",
		"Not yet implemented for this representation",MSTK_WARN);
    return 0;
  }

  List_ptr MF_AdjFaces_R4(MFace_ptr f) {
    MSTK_Report("MF_AdjFaces_R4",
		"Not yet implemented for this representation",MSTK_WARN);
    return 0;
  }

  void MF_Add_AdjFace_R4(MFace_ptr f, int side, MFace_ptr af) {
#ifdef DEBUG
    MSTK_Report("MF_Add_AdjFace_R4",
		"Not yet implemented for this representation",MSTK_WARN);
#endif
  }

  void MF_Rem_AdjFace_R4(MFace_ptr f, int side, MFace_ptr af) {
#ifdef DEBUG
    MSTK_Report("MF_Rem_AdjFace_R4",
		"Not yet implemented for this representation",MSTK_WARN);
#endif
  }

  MFace_ptr MF_NextInHash_R4(MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("MF_NextInHash", "Function call not suitable for this representation", MSTK_WARN);
#endif
    return NULL;
  }

  void MF_Set_NextInHash_R4(MFace_ptr f, MFace_ptr next) {
#ifdef DEBUG
    MSTK_Report("MF_NextInHash", "Function call not suitable for this representation", MSTK_WARN);
#endif
  }

  void MF_HashKey_R4(MFace_ptr f, unsigned int *pn, void* **pp) {
#ifdef DEBUG
    MSTK_Report("MF_NextInHash", "Function call not suitable for this representation", MSTK_WARN);
#endif
  }

  int MF_IsLocked_R4(MFace_ptr f) {
    return 0;
  }

#ifdef __cplusplus
}
#endif
