#define _H_MVertex_Private

#include "MVertex.h"
#include "MVertex_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MV_Set_RepType_R1(MVertex_ptr v) {
    MVertex_UpAdj_R1R2 *upadj;

    upadj = v->upadj = (MVertex_UpAdj_R1R2 *) MSTK_malloc(sizeof(MVertex_UpAdj_R1R2));
    upadj->nel = (unsigned int) 0;
    upadj->velements = List_New(10);
  }

  void MV_Delete_R1(MVertex_ptr v) {
    MVertex_UpAdj_R1R2 *upadj;

    upadj = (MVertex_UpAdj_R1R2 *) v->upadj;
    List_Delete(upadj->velements);
    MSTK_free(upadj);
  }

  int MV_Num_AdjVertices_R1(MVertex_ptr v) {
#ifdef DEBUG
    MSTK_Report("MV_Num_AdjVertices",
		"Inefficient to call this routine with this representation",
		MESG);
#endif
    return MV_Num_Edges_R1(v);
  }

  int MV_Num_Edges_R1(MVertex_ptr v) {
    int ne;
    List_ptr vedges;

#ifdef DEBUG
    MSTK_Report("MV_Num_Edges",
		"Inefficient to call this routine with this representation",
		MESG);
#endif
    
    vedges = MV_Edges_R1(v);
    ne = List_Num_Entries(vedges);
    List_Delete(vedges);
    return ne;
  }

  int MV_Num_Faces_R1(MVertex_ptr v) {
    int nf;
    List_ptr vfaces;

#ifdef DEBUG
    MSTK_Report("MV_Num_Faces",
		"Inefficient to call this routine with this representation",
		MESG);
#endif
    
    /* Have to account for the fact that this may be a surface mesh */
    /* Then v->upadj->velements contains faces of vertex which can be
       retrieved fairly efficiently */

    vfaces = MV_Faces_R1(v);
    nf = List_Num_Entries(vfaces);
    List_Delete(vfaces);
    return nf;
  }
  
  int MV_Num_Regions_R1(MVertex_ptr v) {
    MVertex_UpAdj_R1R2 *upadj;
    int i, nr = 0;
    MEntity_ptr ent;

    upadj = (MVertex_UpAdj_R1R2 *) v->upadj;
    for (i = 0; i < upadj->nel; i++) {
      ent = (MEntity_ptr) List_Entry(upadj->velements,i);
      if (MEnt_Dim(ent) == 3)
	nr++;
    }
    return nr;
  }

  List_ptr MV_AdjVertices_R1(MVertex_ptr v) {
    MSTK_Report("MV_AdjVertices_R1",
		"Not yet implemented for this representation",MESG);
    return 0;
  }

  List_ptr MV_Edges_R1(MVertex_ptr v) {
    MSTK_Report("MV_Edges_R1",
		"Not yet implemented for this representation",MESG);
    return 0;
  }

  List_ptr MV_Faces_R1(MVertex_ptr v) {
    MSTK_Report("MV_Faces_R1",
		"Not yet implemented for this representation",MESG);
    return 0;
  }

  List_ptr MV_Regions_R1(MVertex_ptr v) {
    MVertex_UpAdj_R1R2 *upadj;
    int i, nr = 0;
    MEntity_ptr ent;
    List_ptr vregions;

    upadj = (MVertex_UpAdj_R1R2*) v->upadj;
    vregions = List_New(10);

    for (i = 0; i < upadj->nel; i++) {
      ent = (MEntity_ptr) List_Entry(upadj->velements,i);
      if (MEnt_Dim(ent) == MREGION) 
	List_Add(vregions,ent);
    }
    if (nr)
      return vregions;
    else {
      List_Delete(vregions);
      return 0;
    }      
  }

  void MV_Add_Region_R1(MVertex_ptr v, MRegion_ptr mregion) {
    MVertex_UpAdj_R1R2 *upadj;

    upadj = (MVertex_UpAdj_R1R2 *) v->upadj;
    List_Add(upadj->velements,mregion);
    (upadj->nel)++;
  }

  void MV_Rem_Region_R1(MVertex_ptr v, MRegion_ptr mregion) {
   MVertex_UpAdj_R1R2 *upadj;

    upadj = (MVertex_UpAdj_R1R2 *) v->upadj;
    if (List_Rem(upadj->velements,mregion))
      (upadj->nel)--;
  }

  void MV_Add_Face_R1(MVertex_ptr v, MFace_ptr mface) {
    MVertex_UpAdj_R1R2 *upadj;

    if (MF_Region(mface,0) || MF_Region(mface,1)) {
      MSTK_Report("MV_Add_Face_R1",
		 "Can only add faces with no regions in this representation",
		 ERROR);
      return;
    }

    upadj = v->upadj;
    List_Add(upadj->velements,mface);
    (upadj->nel)++;
  }

  void MV_Rem_Face_R1(MVertex_ptr v, MFace_ptr mface) {
    MVertex_UpAdj_R1R2 *upadj;

    upadj = (MVertex_UpAdj_R1R2 *)v->upadj;
    if (List_Rem(upadj->velements,mface))
      (upadj->nel)--;
  }

  void MV_Add_AdjVertex_R1(MVertex_ptr v, MVertex_ptr av) {
#ifdef DEBUG
    MSTK_Report("MV_Add_AdjVertex","Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Rem_AdjVertex_R1(MVertex_ptr v, MVertex_ptr av) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_AdjVertex","Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Add_Edge_R1(MVertex_ptr v, MEdge_ptr e) {
#ifdef DEBUG
    MSTK_Report("MV_Add_Edge","Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Rem_Edge_R1(MVertex_ptr v, MEdge_ptr e) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_Edge","Function call not suitable for this representation",WARN);
#endif
  }


#ifdef __cplusplus
}
#endif
