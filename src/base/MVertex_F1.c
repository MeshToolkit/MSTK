#define _H_MVertex_Private

#include "MVertex.h"
#include "MVertex_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif


  void MV_Set_RepType_F1(MVertex_ptr v) {
    MVertex_UpAdj_F1F4 *upadj;

    upadj = v->upadj = (MVertex_UpAdj_F1F4 *) MSTK_malloc(sizeof(MVertex_UpAdj_F1F4));
    upadj->ne = (unsigned int) 0;
    upadj->vedges = Set_New(10);
  }

  void MV_Delete_F1(MVertex_ptr v) {
    MVertex_UpAdj_F1F4 *upadj;

    upadj = (MVertex_UpAdj_F1F4 *) v->upadj;
    Set_Delete(upadj->vedges);
    MSTK_free(upadj);

    MSTK_free(v);
  }

  int MV_Num_AdjVertices_F1(MVertex_ptr v) {
    return ((MVertex_UpAdj_F1F4 *) v->upadj)->ne; 
  }

  int MV_Num_Edges_F1(MVertex_ptr v) {
    return ((MVertex_UpAdj_F1F4 *) v->upadj)->ne;
  }

  int MV_Num_Faces_F1(MVertex_ptr v) {
    Set_ptr vfaces;
    int nf;

#ifdef DEBUG
    MSTK_Report("MV_Num_Faces",
		"Ineficient to call this routine with this representation",
		MESG);
#endif
    
    vfaces = MV_Faces_F1(v);
    nf = Set_Num_Entries(vfaces);
    Set_Delete(vfaces);
    return nf;
  }

  int MV_Num_Regions_F1(MVertex_ptr v) {
    int nr;
    Set_ptr vregions;

#ifdef DEBUG
    MSTK_Report("MV_Num_Regions",
		"Inefficient to call this routine with this representation",
		MESG);
#endif
	
    vregions = MV_Regions_F1(v);
    nr = Set_Num_Entries(vregions);
    Set_Delete(vregions);

    return nr;
  }

  Set_ptr MV_AdjVertices_F1(MVertex_ptr v) {
    MVertex_UpAdj_F1F4 *upadj;
    Set_ptr vedges, adjv;
    int ne, i;
    MEdge_ptr vedge;
    MVertex_ptr ov;

    upadj = (MVertex_UpAdj_F1F4 *) v->upadj;
    vedges = upadj->vedges;
    if (vedges == 0)
      return 0;

    ne = upadj->ne;
    adjv = Set_New(ne);
    for (i = 0; i < ne; i++) {
      vedge = Set_Entry(vedges,i);
      ov = ME_OppVertex(vedge,v);
      Set_Add(adjv,ov);
    }

    return adjv;
  }
    

  Set_ptr MV_Edges_F1(MVertex_ptr v) {
    Set_ptr vedges;
    MVertex_UpAdj_F1F4 *upadj;

    upadj = (MVertex_UpAdj_F1F4 *) v->upadj;
    vedges = Set_Copy(upadj->vedges);
    return vedges;
  }

  Set_ptr MV_Faces_F1(MVertex_ptr v) {
    MVertex_UpAdj_F1F4 *upadj;
    int i, j, ne, nf, n, mkr;
    Set_ptr vedges, efaces, vfaces;
    MEdge_ptr edge;
    MFace_ptr face;

    upadj = (MVertex_UpAdj_F1F4 *) v->upadj;
    ne = upadj->ne;
    vedges = upadj->vedges;

    n = 0;
    vfaces = Set_New(10);
    mkr = MSTK_GetMarker();

    for (i = 0; i < ne; i++) {
      edge = Set_Entry(vedges,i);

      efaces = ME_Faces(edge);
      nf = Set_Num_Entries(efaces);
	
      for (j = 0; j < nf; j++) {
	face = Set_Entry(efaces,j);
	if (!MEnt_IsMarked(face,mkr)) {
	  MEnt_Mark(face,mkr);
	  Set_Add(vfaces,face);
	  n++;
	}
      }
      Set_Delete(efaces);
    }
    Set_Unmark(vfaces,mkr);
    MSTK_FreeMarker(mkr);
    if (n > 0)
      return vfaces;
    else {
      Set_Delete(vfaces);
      return 0;
    }
  }

  Set_ptr MV_Regions_F1(MVertex_ptr v) {
    MVertex_UpAdj_F1F4 *upadj;
    int i, j, k, ne, nf, n, mkr;
    Set_ptr vedges, efaces, vregions;
    MEdge_ptr edge;
    MFace_ptr eface;
    MRegion_ptr region;

    upadj = (MVertex_UpAdj_F1F4 *) v->upadj;
    ne = upadj->ne;
    vedges = upadj->vedges;

    n = 0;
    vregions = Set_New(10);
    mkr = MSTK_GetMarker();

    for (i = 0; i < ne; i++) {
      edge = Set_Entry(vedges,i);

      efaces = ME_Faces(edge);
      nf = Set_Num_Entries(efaces);
	
      for (j = 0; j < nf; j++) {
	eface = Set_Entry(efaces,j);
	for (k = 0; k < 2; k++) {
	  region = MF_Region(eface,k);
	  if (region && !MEnt_IsMarked(region,mkr)) {
	    MEnt_Mark(region,mkr);
	    Set_Add(vregions,region);
	    n++;
	  }
	}
      }
      Set_Delete(efaces);
    }
    Set_Unmark(vregions,mkr);
    MSTK_FreeMarker(mkr);

    if (n > 0)
      return vregions;
    else {
      Set_Delete(vregions);
      return 0;
    }
  }

  void MV_Add_AdjVertex_F1(MVertex_ptr v, MVertex_ptr av) {
#ifdef DEBUG
    MSTK_Report("MV_Add_AdjVertex","Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Rem_AdjVertex_F1(MVertex_ptr v, MVertex_ptr av) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_AdjVertex","Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Add_Edge_F1(MVertex_ptr v, MEdge_ptr e) {
    MVertex_UpAdj_F1F4 *upadj;
    
    upadj = v->upadj;
    if (upadj->vedges == NULL)
      upadj->vedges = Set_New(10);
    Set_Add(upadj->vedges,e);
    (upadj->ne)++;
  }

  void MV_Rem_Edge_F1(MVertex_ptr v, MEdge_ptr e) {
    MVertex_UpAdj_F1F4 *upadj;

    upadj = (MVertex_UpAdj_F1F4 *)v->upadj;
    if (upadj->vedges == NULL)
      return;
    if (Set_Rem(upadj->vedges,e))
      (upadj->ne)--;
  }

  void MV_Add_Face_F1(MVertex_ptr v, MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("MV_Add_Face","Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Rem_Face_F1(MVertex_ptr v, MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_Face","Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Add_Region_F1(MVertex_ptr v, MRegion_ptr f) {
#ifdef DEBUG
    MSTK_Report("MV_Add_Region","Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Rem_Region_F1(MVertex_ptr v, MRegion_ptr f) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_Region","Function call not suitable for this representation",WARN);
#endif
  }

#ifdef __cplusplus
}
#endif
