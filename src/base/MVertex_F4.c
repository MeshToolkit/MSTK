#define _H_MVertex_Private

#include "MVertex.h"
#include "MVertex_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif


  void MV_Set_RepType_F4(MVertex_ptr v) {
    MVertex_UpAdj_F1F4 *upadj;

    upadj = v->upadj = (MVertex_UpAdj_F1F4 *) MSTK_malloc(sizeof(MVertex_UpAdj_F1F4));
    upadj->ne = (unsigned int) 0;
    upadj->vedges = Set_New(10);
  }

  void MV_Delete_F4(MVertex_ptr v) {
    MVertex_UpAdj_F1F4 *upadj;

    upadj = (MVertex_UpAdj_F1F4 *) v->upadj;
    Set_Delete(upadj->vedges);
    MSTK_free(upadj);
  }

  int MV_Num_AdjVertices_F4(MVertex_ptr v) {
    return ((MVertex_UpAdj_F1F4 *) v->upadj)->ne;
  }

  int MV_Num_Edges_F4(MVertex_ptr v) {
    return ((MVertex_UpAdj_F1F4 *) v->upadj)->ne;
  }

  int MV_Num_Faces_F4(MVertex_ptr v) {
    Set_ptr vfaces;
    int nf;

#ifdef DEBUG
    MSTK_Report("MV_Num_Faces",
		"Ineficient to call this routine with this representation",
		MESG);
#endif
    
    vfaces = MV_Faces_F4(v);
    nf = Set_Num_Entries(vfaces);
    Set_Delete(vfaces);
    return nf;
  }

  int MV_Num_Regions_F4(MVertex_ptr v) {
    Set_ptr vregions;
    int nr;

#ifdef DEBUG
    MSTK_Report("MV_Num_Regions",
		"Inefficient to call this routine with this representation",
		MESG);
#endif
	
    vregions = MV_Regions_F4(v);
    nr = Set_Num_Entries(vregions);
    Set_Delete(vregions);
    return nr;
  }

  Set_ptr MV_AdjVertices_F4(MVertex_ptr v) {
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
    

  Set_ptr MV_Edges_F4(MVertex_ptr v) {
    MVertex_UpAdj_F1F4 *upadj;
    Set_ptr vedges;

    upadj = (MVertex_UpAdj_F1F4 *) v->upadj;
    vedges = Set_Copy(upadj->vedges);
    return vedges;
  }

  Set_ptr MV_Faces_F4(MVertex_ptr v) {
    MVertex_UpAdj_F1F4 *upadj;
    int i, j, k, ne, nf, nr, n, mkr;
    Set_ptr vedges, eregions, rfaces, efaces, vfaces;
    MEdge_ptr edge;
    MFace_ptr face;
    MRegion_ptr region;

    upadj = (MVertex_UpAdj_F1F4 *) v->upadj;
    ne = upadj->ne;
    vedges = upadj->vedges;

    n = 0;
    vfaces = Set_New(10);
    mkr = MSTK_GetMarker();

    for (i = 0; i < ne; i++) {
      edge = Set_Entry(vedges,i);

      eregions = ME_Regions(edge);
      if (eregions) {
	nr = Set_Num_Entries(eregions);
	
	for (j = 0; j < nr; j++) {
	  region = Set_Entry(eregions,j);
	  
	  rfaces = MR_Faces(region);
	  nf = Set_Num_Entries(rfaces);
	  
	  for (k = 0; k < nf; k++) {
	    face = Set_Entry(rfaces,k);
	    if (!MEnt_IsMarked(face,mkr)) {
	      if (MF_UsesEntity(face,v,0)) {
		MEnt_Mark(face,mkr);
		Set_Add(vfaces,face);
		n++;
	      }
	    }
	  }
	  Set_Delete(rfaces);
	}
	Set_Delete(eregions);
      }
      else {
	/* perhaps the edge has boundary faces (not connected to regions) */
	efaces = ME_Faces(edge);
	if (efaces) {
	  nf = Set_Num_Entries(efaces);
	  
	  for (k = 0; k < nf; k++) {
	    face = Set_Entry(efaces,k);
	    if (!MEnt_IsMarked(face,mkr)) {
	      if (MF_UsesEntity(face,v,0)) {
		MEnt_Mark(face,mkr);
		Set_Add(vfaces,face);
		n++;
	      }
	    }
	  }
	  Set_Delete(efaces);
	}
      }
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

  Set_ptr MV_Regions_F4(MVertex_ptr v) {
    MVertex_UpAdj_F1F4 *upadj;
    int i, j, ne, nr, n, mkr;
    Set_ptr vedges, eregions, vregions;
    MEdge_ptr edge;
    MRegion_ptr region;

    upadj = (MVertex_UpAdj_F1F4 *) v->upadj;
    ne = upadj->ne;
    vedges = upadj->vedges;

    n = 0;
    vregions = Set_New(10);
    mkr = MSTK_GetMarker();

    for (i = 0; i < ne; i++) {
      edge = Set_Entry(vedges,i);

      eregions = ME_Regions(edge);
      if (eregions) {
	nr = Set_Num_Entries(eregions);
	
	for (j = 0; j < nr; j++) {
	  region = Set_Entry(eregions,j);
	  if (!MEnt_IsMarked(region,mkr)) {
	    MEnt_Mark(region,mkr);
	    Set_Add(vregions,region);
	    n++;
	  }
	}
	
	Set_Delete(eregions);
      }
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

  void MV_Add_AdjVertex_F4(MVertex_ptr v, MVertex_ptr av) {
#ifdef DEBUG
    MSTK_Report("MV_Add_AdjVertex","Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Rem_AdjVertex_F4(MVertex_ptr v, MVertex_ptr av) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_AdjVertex","Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Add_Edge_F4(MVertex_ptr v, MEdge_ptr e) {
    MVertex_UpAdj_F1F4 *upadj;

    upadj = (MVertex_UpAdj_F1F4 *) v->upadj;

    if (upadj->vedges == NULL)
      upadj->vedges = Set_New(10);
    Set_Add(upadj->vedges,e);
    (upadj->ne)++;
  }

  void MV_Rem_Edge_F4(MVertex_ptr v, MEdge_ptr e) {
    MVertex_UpAdj_F1F4 *upadj;

    upadj = (MVertex_UpAdj_F1F4 *) v->upadj;

    if (upadj->vedges == NULL)
      return;
    if (Set_Rem(upadj->vedges,e))
      (upadj->ne)--;
  }

  void MV_Add_Face_F4(MVertex_ptr v, MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("MV_Add_Face","Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Rem_Face_F4(MVertex_ptr v, MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_Face","Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Add_Region_F4(MVertex_ptr v, MRegion_ptr f) {
#ifdef DEBUG
    MSTK_Report("MV_Add_Region","Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Rem_Region_F4(MVertex_ptr v, MRegion_ptr f) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_Region","Function call not suitable for this representation",WARN);
#endif
  }

#ifdef __cplusplus
}
#endif
