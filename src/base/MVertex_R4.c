#define _H_MVertex_Private

#include "MVertex.h"
#include "MVertex_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Functions */

  void MV_Set_RepType_R4(MVertex_ptr v) {
    MVertex_UpAdj_R3R4 *upadj;
    MVertex_SameAdj_R2R4 *sameadj;

    upadj = v->upadj = (MVertex_UpAdj_R3R4 *) MSTK_malloc(sizeof(MVertex_UpAdj_R3R4));
    upadj->nf = (unsigned int) 0;
    upadj->vfaces = List_New(10);
    sameadj = v->sameadj = (MVertex_SameAdj_R2R4 *) MSTK_malloc(sizeof(MVertex_SameAdj_R2R4));
    sameadj->nvadj = (unsigned char) 0;
    sameadj->adjverts = List_New(10);
  }

  void MV_Delete_R4(MVertex_ptr v) {
    MVertex_UpAdj_R3R4 *upadj;
    MVertex_SameAdj_R2R4 *sameadj;
   
    upadj = (MVertex_UpAdj_R3R4 *) v->upadj;
    List_Delete(upadj->vfaces);
    MSTK_free(upadj);

    sameadj = (MVertex_SameAdj_R2R4 *) v->sameadj;
    List_Delete(sameadj->adjverts);
    MSTK_free(sameadj);
  }

  int MV_Num_AdjVertices_R4(MVertex_ptr v) {
    return ((MVertex_SameAdj_R2R4 *)v->sameadj)->nvadj;
  }

  int MV_Num_Edges_R4(MVertex_ptr v) {
    return ((MVertex_SameAdj_R2R4 *)v->sameadj)->nvadj;
  }

  int MV_Num_Faces_R4(MVertex_ptr v) {
    return ((MVertex_UpAdj_R3R4 *)v->upadj)->nf;
  }
  
  int MV_Num_Regions_R4(MVertex_ptr v) {
    return MV_Num_Regions_R3R4(v);
  }

  List_ptr MV_AdjVertices_R4(MVertex_ptr v) {
    MVertex_SameAdj_R2R4 *sameadj;

    sameadj = (MVertex_SameAdj_R2R4 *)v->sameadj;
    return List_Copy(sameadj->adjverts);
  }

  List_ptr MV_Edges_R4(MVertex_ptr v) {
    /* Have to create volatile edges */
    return MV_Edges_R2R4(v);
  }

  List_ptr MV_Faces_R4(MVertex_ptr v) {
    MVertex_UpAdj_R3R4 *upadj;
    upadj = (MVertex_UpAdj_R3R4 *)v->upadj;

    return List_Copy(upadj->vfaces);
  }

  List_ptr MV_Regions_R4(MVertex_ptr v) {
    return MV_Regions_R3R4(v);
  }

  void MV_Add_AdjVertex_R4(MVertex_ptr v, MVertex_ptr adjv) {
    MVertex_SameAdj_R2R4 *sameadj;

    sameadj = (MVertex_SameAdj_R2R4 *) v->sameadj;
    List_Add(sameadj->adjverts,adjv);
    (sameadj->nvadj)++;
  }

  void MV_Rem_AdjVertex_R4(MVertex_ptr v, MVertex_ptr adjv) {
    MVertex_SameAdj_R2R4 *sameadj;

    sameadj = (MVertex_SameAdj_R2R4 *) v->sameadj;
    if (List_Rem(sameadj->adjverts,adjv))
      (sameadj->nvadj)--;
  }

  void MV_Add_Edge_R4(MVertex_ptr v, MEdge_ptr medge) {
#ifdef DEBUG
    MSTK_Report("MV_Add_Edge","Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Rem_Edge_R4(MVertex_ptr v, MEdge_ptr medge) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_Edge","Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Add_Face_R4(MVertex_ptr v, MFace_ptr mface) {
    MVertex_UpAdj_R3R4 *upadj;

    upadj = (MVertex_UpAdj_R3R4 *) v->upadj;
    List_Add(upadj->vfaces,mface);
    (upadj->nf)++;
  }

  void MV_Rem_Face_R4(MVertex_ptr v, MFace_ptr mface) {
    MVertex_UpAdj_R3R4 *upadj;

    upadj = (MVertex_UpAdj_R3R4 *) v->upadj;

    if (List_Rem(upadj->vfaces,mface))
      (upadj->nf)--;
  }

  void MV_Add_Region_R4(MVertex_ptr v, MRegion_ptr mregion) {
#ifdef DEBUG
    MSTK_Report("MV_Add_Region","Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Rem_Region_R4(MVertex_ptr v, MRegion_ptr mregion) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_Region","Function call not suitable for this representation",WARN);
#endif
  }

#ifdef __cplusplus
}
#endif
