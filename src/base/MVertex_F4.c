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
    upadj->vedges = List_New(10);
  }

  void MV_Delete_F4(MVertex_ptr v, int keep) {
    MVertex_UpAdj_F1F4 *upadj;

    if (keep) {
      MSTK_KEEP_DELETED = 1;
      v->dim = MDELVERTEX;
    }
    else {
#ifdef DEBUG
      v->dim = MDELVERTEX;
#endif

      upadj = (MVertex_UpAdj_F1F4 *) v->upadj;
      List_Delete(upadj->vedges);
      MSTK_free(upadj);
    }
  }

  void MV_Restore_F4(MVertex_ptr v) {
    if (v->dim != MDELVERTEX)
      return;
    v->dim = MVERTEX;
  }

  int MV_Num_AdjVertices_F4(MVertex_ptr v) {
    return ((MVertex_UpAdj_F1F4 *) v->upadj)->ne;
  }

  int MV_Num_Edges_F4(MVertex_ptr v) {
    return ((MVertex_UpAdj_F1F4 *) v->upadj)->ne;
  }

  int MV_Num_Faces_F4(MVertex_ptr v) {
    List_ptr vfaces;
    int nf;

#ifdef DEBUG
    MSTK_Report("MV_Num_Faces",
		"Ineficient to call this routine with this representation",
		MESG);
#endif
    
    vfaces = MV_Faces_F4(v);
    nf = List_Num_Entries(vfaces);
    List_Delete(vfaces);
    return nf;
  }

  int MV_Num_Regions_F4(MVertex_ptr v) {
    List_ptr vregions;
    int nr;

#ifdef DEBUG
    MSTK_Report("MV_Num_Regions",
		"Inefficient to call this routine with this representation",
		MESG);
#endif
	
    vregions = MV_Regions_F4(v);
    nr = List_Num_Entries(vregions);
    List_Delete(vregions);
    return nr;
  }

  List_ptr MV_AdjVertices_F4(MVertex_ptr v) {
    MVertex_UpAdj_F1F4 *upadj;
    List_ptr vedges, adjv;
    int ne, i;
    MEdge_ptr vedge;
    MVertex_ptr ov;

    upadj = (MVertex_UpAdj_F1F4 *) v->upadj;
    vedges = upadj->vedges;
    if (vedges == 0)
      return 0;

    ne = upadj->ne;
    adjv = List_New(ne);
    for (i = 0; i < ne; i++) {
      vedge = List_Entry(vedges,i);
      ov = ME_OppVertex(vedge,v);
      List_Add(adjv,ov);
    }

    return adjv;
  }
    

  List_ptr MV_Edges_F4(MVertex_ptr v) {
    MVertex_UpAdj_F1F4 *upadj;
    List_ptr vedges;

    upadj = (MVertex_UpAdj_F1F4 *) v->upadj;
    vedges = List_Copy(upadj->vedges);
    return vedges;
  }

  List_ptr MV_Faces_F4(MVertex_ptr v) {
    MVertex_UpAdj_F1F4 *upadj;
    int i, j, k, ne, nf, nr, n, mkr;
    List_ptr vedges, eregions, rfaces, efaces, vfaces;
    MEdge_ptr edge;
    MFace_ptr face;
    MRegion_ptr region;

    upadj = (MVertex_UpAdj_F1F4 *) v->upadj;
    ne = upadj->ne;
    vedges = upadj->vedges;

    n = 0;
    vfaces = List_New(10);
    mkr = MSTK_GetMarker();

    for (i = 0; i < ne; i++) {
      edge = List_Entry(vedges,i);
      eregions = ME_Regions(edge);
      if (eregions) {
	nr = List_Num_Entries(eregions);
	
	for (j = 0; j < nr; j++) {
	  region = List_Entry(eregions,j);
	  
	  rfaces = MR_Faces(region);
	  nf = List_Num_Entries(rfaces);
	  
	  for (k = 0; k < nf; k++) {
	    face = List_Entry(rfaces,k);

	    if (!MEnt_IsMarked(face,mkr)) {
	      if (MF_UsesEntity(face,v,0)) {		
		MEnt_Mark(face,mkr);
		List_Add(vfaces,face);
		n++;
	      }
	    }
	  }
	  List_Delete(rfaces);
	}
	List_Delete(eregions);
      }
      else {
	/* perhaps the edge has boundary faces (not connected to regions) */
	efaces = ME_Faces(edge);
	if (efaces) {
	  nf = List_Num_Entries(efaces);
	  
	  for (k = 0; k < nf; k++) {
	    face = List_Entry(efaces,k);
	    if (!MEnt_IsMarked(face,mkr)) {
	      if (MF_UsesEntity(face,v,0)) {
		MEnt_Mark(face,mkr);
		List_Add(vfaces,face);
		n++;
	      }
	    }
	  }
	  List_Delete(efaces);
	}
      }
    }
    List_Unmark(vfaces,mkr);
    MSTK_FreeMarker(mkr);
    if (n > 0)
      return vfaces;
    else {
      List_Delete(vfaces);
      return 0;
    }
  }

  List_ptr MV_Regions_F4(MVertex_ptr v) {
    MVertex_UpAdj_F1F4 *upadj;
    int i, j, ne, nr, n, mkr;
    List_ptr vedges, eregions, vregions;
    MEdge_ptr edge;
    MRegion_ptr region;

    upadj = (MVertex_UpAdj_F1F4 *) v->upadj;
    ne = upadj->ne;
    vedges = upadj->vedges;

    n = 0;
    vregions = List_New(10);
    mkr = MSTK_GetMarker();

    for (i = 0; i < ne; i++) {
      edge = List_Entry(vedges,i);

      eregions = ME_Regions(edge);
      if (eregions) {
	nr = List_Num_Entries(eregions);
	
	for (j = 0; j < nr; j++) {
	  region = List_Entry(eregions,j);

	  if (!MEnt_IsMarked(region,mkr)) {
	    MEnt_Mark(region,mkr);
	    List_Add(vregions,region);
	    n++;
	  }
	}
	
	List_Delete(eregions);
      }
    }
    List_Unmark(vregions,mkr);
    MSTK_FreeMarker(mkr);

    if (n > 0)
      return vregions;
    else {
      List_Delete(vregions);
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
      upadj->vedges = List_New(10);
    List_Add(upadj->vedges,e);
    (upadj->ne)++;
  }

  void MV_Rem_Edge_F4(MVertex_ptr v, MEdge_ptr e) {
    MVertex_UpAdj_F1F4 *upadj;

    upadj = (MVertex_UpAdj_F1F4 *) v->upadj;

    if (upadj->vedges == NULL)
      return;
    if (List_Rem(upadj->vedges,e))
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
