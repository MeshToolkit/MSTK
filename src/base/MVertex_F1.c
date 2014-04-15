#define _H_MVertex_Private

#include "MVertex.h"
#include "MVertex_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif


  void MV_Set_RepType_F1(MVertex_ptr v) {
    MVertex_Adj_F1F4 *adj;

    adj = v->adj = (MVertex_Adj_F1F4 *) MSTK_malloc(sizeof(MVertex_Adj_F1F4));
    adj->vedges = List_New(5);
  }

  void MV_Delete_F1(MVertex_ptr v, int keep) {
    MVertex_Adj_F1F4 *adj;

    if (!keep) {
      adj = (MVertex_Adj_F1F4 *) v->adj;
      if (adj) {
	if (adj->vedges)
	  List_Delete(adj->vedges);
	MSTK_free(adj);
      }
    }
  }

  void MV_Restore_F1(MVertex_ptr v) {
    MEnt_Set_Dim((MEntity_ptr) v,MVERTEX);
  }

  void MV_Destroy_For_MESH_Delete_F1(MVertex_ptr v) {
    MVertex_Adj_F1F4 *adj;

    adj = (MVertex_Adj_F1F4 *) v->adj;
    if (adj) {
      if (adj->vedges)
	List_Delete(adj->vedges);
      MSTK_free(adj);
    }
  }

  int MV_Num_AdjVertices_F1(MVertex_ptr v) {
    List_ptr vedges = ((MVertex_Adj_F1F4 *) v->adj)->vedges;
    return List_Num_Entries(vedges);
  }

  int MV_Num_Edges_F1(MVertex_ptr v) {
    List_ptr vedges = ((MVertex_Adj_F1F4 *) v->adj)->vedges;
    return List_Num_Entries(vedges);
  }

  int MV_Num_Faces_F1(MVertex_ptr v) {
    List_ptr vfaces;
    int nf;

#ifdef DEBUG
    MSTK_Report("MV_Num_Faces",
		"Ineficient to call this routine with this representation",
		MSTK_MESG);
#endif
    
    vfaces = MV_Faces_F1(v);
    if (vfaces) {
      nf = List_Num_Entries(vfaces);
      List_Delete(vfaces);
    }
    else 
      nf = 0;
    return nf;
  }

  int MV_Num_Regions_F1(MVertex_ptr v) {
    int nr;
    List_ptr vregions;

#ifdef DEBUG
    MSTK_Report("MV_Num_Regions",
		"Inefficient to call this routine with this representation",
		MSTK_MESG);
#endif
	
    vregions = MV_Regions_F1(v);
    nr = List_Num_Entries(vregions);
    List_Delete(vregions);

    return nr;
  }

  List_ptr MV_AdjVertices_F1(MVertex_ptr v) {
    MVertex_Adj_F1F4 *adj;
    List_ptr vedges, adjv;
    int ne, i;
    MEdge_ptr vedge;
    MVertex_ptr ov;

    adj = (MVertex_Adj_F1F4 *) v->adj;
    vedges = adj->vedges;
    if (vedges == 0)
      return 0;

    ne = List_Num_Entries(vedges);
    adjv = List_New(ne);
    for (i = 0; i < ne; i++) {
      vedge = List_Entry(vedges,i);
      ov = ME_OppVertex(vedge,v);
      List_Add(adjv,ov);
    }

    return adjv;
  }

  void MV_AdjVertexIDs_F1(MVertex_ptr v, int *nvadj, int *adjvids) {
    MVertex_Adj_F1F4 *adj;
    List_ptr vedges, adjv;
    int ne, i;
    MEdge_ptr vedge;
    MVertex_ptr ov;

    adj = (MVertex_Adj_F1F4 *) v->adj;

    *nvadj = 0;
    vedges = adj->vedges;
    if (vedges) {
      *nvadj = List_Num_Entries(vedges);
      for (i = 0; i < *nvadj; i++) {
        vedge = List_Entry(vedges,i);
        ov = ME_OppVertex(vedge,v);
        adjvids[i] = MEnt_ID((MEntity_ptr)ov);
      }      
    }
  }
    

  List_ptr MV_Edges_F1(MVertex_ptr v) {
    MVertex_Adj_F1F4 *adj;
    List_ptr vedges=NULL;

    adj = (MVertex_Adj_F1F4 *) v->adj;
    if (List_Num_Entries(adj->vedges))
      vedges = List_Copy(adj->vedges);
    return vedges;
  }

  void MV_EdgeIDs_F1(MVertex_ptr v, int *nve, int *vedgeids) {
    MVertex_Adj_F1F4 *adj;
    adj = (MVertex_Adj_F1F4 *) v->adj;

    if (adj->vedges) {
      int i;
      *nve = List_Num_Entries(adj->vedges);
      for (i = 0; i < *nve; i++) 
        vedgeids[i] = MEnt_ID(List_Entry(adj->vedges,i));
    }
    else
      *nve = 0;
  }

  List_ptr MV_Faces_F1(MVertex_ptr v) {
    MVertex_Adj_F1F4 *adj;
    int i, j, ne, nf, n, mkr;
    List_ptr vedges, efaces, vfaces;
    MEdge_ptr edge;
    MFace_ptr face;

    adj = (MVertex_Adj_F1F4 *) v->adj;
    vedges = adj->vedges;
    if (!adj->vedges)
      return NULL;
    ne = List_Num_Entries(vedges);

    n = 0;
    vfaces = List_New(ne);
    mkr = MSTK_GetMarker();

    for (i = 0; i < ne; i++) {
      edge = List_Entry(vedges,i);

      efaces = ME_Faces(edge);
      if (!efaces)
        continue;

      nf = List_Num_Entries(efaces);
	
      for (j = 0; j < nf; j++) {
	face = List_Entry(efaces,j);
	if (!MEnt_IsMarked(face,mkr)) {
	  MEnt_Mark(face,mkr);
	  List_Add(vfaces,face);
	  n++;
	}
      }
      List_Delete(efaces);
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

  void MV_FaceIDs_F1(MVertex_ptr v, int *nvf, int *vfaceids) {
    List_ptr vfaces = MV_Faces_F1(v);
    if (vfaces) {
      int i;
      *nvf = List_Num_Entries(vfaces);
      for (i = 0; i < *nvf; i++)
        vfaceids[i] = MEnt_ID(List_Entry(vfaces,i));
      List_Delete(vfaces);
    }
    else
      *nvf = 0;
  }

  List_ptr MV_Regions_F1(MVertex_ptr v) {
    MVertex_Adj_F1F4 *adj;
    int i, j, k, ne, nf, n, mkr;
    List_ptr vedges, efaces, vregions;
    MEdge_ptr edge;
    MFace_ptr eface;
    MRegion_ptr region;

    adj = (MVertex_Adj_F1F4 *) v->adj;
    vedges = adj->vedges;
    if (!adj->vedges)
      return NULL;
    ne = List_Num_Entries(adj->vedges);

    n = 0;
    vregions = List_New(ne);
    mkr = MSTK_GetMarker();

    for (i = 0; i < ne; i++) {
      edge = List_Entry(vedges,i);

      efaces = ME_Faces(edge);
      if (!efaces)
        continue;
      nf = List_Num_Entries(efaces);
	
      for (j = 0; j < nf; j++) {
	eface = List_Entry(efaces,j);
	for (k = 0; k < 2; k++) {
	  region = MF_Region(eface,k);
	  if (region && !MEnt_IsMarked(region,mkr)) {
	    MEnt_Mark(region,mkr);
	    List_Add(vregions,region);
	    n++;
	  }
	}
      }
      List_Delete(efaces);
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

  void MV_RegionIDs_F1(MVertex_ptr v, int *nvr, int *vregionids) {
    List_ptr vregions = MV_Regions_F1(v);
    if (vregions) {
      int i;
      *nvr = List_Num_Entries(vregions);
      for (i = 0; i < *nvr; i++)
        vregionids[i] = MEnt_ID(List_Entry(vregions,i));
      List_Delete(vregions);
    }
    else
      *nvr = 0;
  }

  void MV_Add_AdjVertex_F1(MVertex_ptr v, MVertex_ptr av) {
#ifdef DEBUG
    MSTK_Report("MV_Add_AdjVertex","Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void MV_Rem_AdjVertex_F1(MVertex_ptr v, MVertex_ptr av) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_AdjVertex","Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void MV_Add_Edge_F1(MVertex_ptr v, MEdge_ptr e) {
    MVertex_Adj_F1F4 *adj;
    
    adj = v->adj;
    if (adj->vedges == NULL)
      adj->vedges = List_New(10);
    List_ChknAdd(adj->vedges,e);
  }

  void MV_Rem_Edge_F1(MVertex_ptr v, MEdge_ptr e) {
    MVertex_Adj_F1F4 *adj;

    adj = (MVertex_Adj_F1F4 *)v->adj;
    if (adj->vedges == NULL)
      return;
    List_Rem(adj->vedges,e);
  }

  void MV_Add_Face_F1(MVertex_ptr v, MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("MV_Add_Face","Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void MV_Rem_Face_F1(MVertex_ptr v, MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_Face","Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void MV_Add_Region_F1(MVertex_ptr v, MRegion_ptr f) {
#ifdef DEBUG
    MSTK_Report("MV_Add_Region","Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void MV_Rem_Region_F1(MVertex_ptr v, MRegion_ptr f) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_Region","Function call not suitable for this representation",MSTK_WARN);
#endif
  }

#ifdef __cplusplus
}
#endif
