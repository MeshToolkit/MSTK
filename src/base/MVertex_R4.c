#define _H_MVertex_Private

#include <stdlib.h>

#include "MVertex.h"
#include "MVertex_jmp.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Functions */

  void MV_Set_RepType_R4(MVertex_ptr v) {
    MVertex_Adj_R4 *adj;

    adj = v->adj = (MVertex_Adj_R4 *) malloc(sizeof(MVertex_Adj_R4));
    adj->vfaces = List_New(10);
    adj->adjverts = List_New(10);
  }

  void MV_Delete_R4(MVertex_ptr v, int keep) {
    MVertex_Adj_R4 *adj;
    int idx;
    MVertex_ptr adjv;

    if (MEnt_Dim((MEntity_ptr) v) != MDELETED) { /* if vtx has not been temporarily deleted*/
      adj = (MVertex_Adj_R4 *) v->adj;
      if (adj) {
	idx = 0;
	while ((adjv = List_Next_Entry(adj->adjverts,&idx)))
	  MV_Rem_AdjVertex_R4(adjv,v);
      }
    }

    if (!keep) {
      adj = (MVertex_Adj_R4 *) v->adj;
      if (adj) {
	if (adj->vfaces)
	  List_Delete(adj->vfaces);
	if (adj->adjverts)
	  List_Delete(adj->adjverts);
	free(adj);
      }
    }
  }

  void MV_Restore_R4(MVertex_ptr v) {
    MVertex_Adj_R4 *adj;
    int idx;
    MVertex_ptr adjv;

    MEnt_Set_Dim((MEntity_ptr) v,MVERTEX);

    adj = (MVertex_Adj_R4 *) v->adj;
    idx = 0;
    while ((adjv = List_Next_Entry(adj->adjverts,&idx)))
      MV_Add_AdjVertex_R4(adjv,v);
  }

  void MV_Destroy_For_MESH_Delete_R4(MVertex_ptr v) {
    MVertex_Adj_R4 *adj;

    adj = (MVertex_Adj_R4 *) v->adj;
    if (adj) {
      if (adj->vfaces)
	List_Delete(adj->vfaces);
      if (adj->adjverts)
	List_Delete(adj->adjverts);
      free(adj);
    }
  }

  int MV_Num_AdjVertices_R4(MVertex_ptr v) {
    List_ptr adjverts = ((MVertex_Adj_R4 *)v->adj)->adjverts;
    return List_Num_Entries(adjverts);
  }

  int MV_Num_Edges_R4(MVertex_ptr v) {
    List_ptr adjverts = ((MVertex_Adj_R4 *)v->adj)->adjverts;
    return List_Num_Entries(adjverts);
  }

  int MV_Num_Faces_R4(MVertex_ptr v) {
    List_ptr vfaces = ((MVertex_Adj_R4 *)v->adj)->vfaces;
    return List_Num_Entries(vfaces);
  }
  
  int MV_Num_Regions_R4(MVertex_ptr v) {
    List_ptr vregions;
    int nr;

#ifdef DEBUG
    MSTK_Report("MV_Num_Regions_R4",
		"Inefficient to call this routine with this representation",
		MSTK_MESG);
#endif
    vregions = MV_Regions_R4(v);
    nr = List_Num_Entries(vregions);
    List_Delete(vregions);

    return nr;
  }

  List_ptr MV_AdjVertices_R4(MVertex_ptr v) {
    MVertex_Adj_R4 *adj;

    adj = (MVertex_Adj_R4 *)v->adj;
    return List_Copy(adj->adjverts);
  }

  void MV_AdjVertexIDs_R4(MVertex_ptr v, int *nvadj, int *adjvids) {
    int i;
    MVertex_Adj_R4 *adj;
    adj = (MVertex_Adj_R4 *)v->adj;
    *nvadj = List_Num_Entries(adj->adjverts);
    for (i = 0; i < *nvadj; i++)
      adjvids[i] = MEnt_ID(List_Entry(adj->adjverts,i));
  }

  List_ptr MV_Edges_R4(MVertex_ptr v) {
    int i, nvadj;
    MVertex_ptr adjv;
    List_ptr vedges;
    MVertex_Adj_R4 *adj;

    MSTK_Report("MV_Edges_R4","Incomplete implementation",MSTK_ERROR);

    adj = (MVertex_Adj_R4 *) v->adj;
    nvadj = List_Num_Entries(adj->adjverts);

    /* Have to create volatile edges - SEEMS INCOMPLETE */

    vedges = List_New(nvadj);
    for (i = 0; i < nvadj; i++) {
      /* create a volatile edge from v and the adjacent vertex */
      adjv = List_Entry(adj->adjverts,i);
    }

    return vedges;
  }

  void MV_EdgeIDs_R4(MVertex_ptr v, int *nve, int *vedgeids) {
    List_ptr vedges = MV_Edges_R4(v);

    if (vedges) {
      int i;
      *nve = List_Num_Entries(vedges);
      for (i = 0; i < *nve; i++) 
        vedgeids[i] = MEnt_ID(List_Entry(vedges,i));
      List_Delete(vedges);
    }
    else
      *nve = 0;
  }

  List_ptr MV_Faces_R4(MVertex_ptr v) {
    MVertex_Adj_R4 *adj;
    adj = (MVertex_Adj_R4 *)v->adj;

    return List_Copy(adj->vfaces);
  }


  void MV_FaceIDs_R4(MVertex_ptr v, int *nvf, int *vfaceids) {
    int i;
    MVertex_Adj_R4 *adj;
    adj = (MVertex_Adj_R4 *)v->adj;
    
    *nvf = List_Num_Entries(adj->vfaces);
    for (i = 0; i < *nvf; i++)
      vfaceids[i] = MEnt_ID(List_Entry(vfaceids,i));
  }

  List_ptr MV_Regions_R4(MVertex_ptr v) {
    int i, j, nr, nf, mkr;
    MFace_ptr face;
    MRegion_ptr reg;
    List_ptr vregions, vfaces;
 
    vregions = List_New(10);
    nr = 0;
    mkr = MSTK_GetMarker();

    vfaces = MV_Faces_R4(v);
    nf = List_Num_Entries(vfaces);
    for (i = 0; i < nf; i++) {
      face = List_Entry(vfaces,i);
      for (j = 0; j < 2; j++) {
	reg = MF_Region(face,j);
	if (reg && !MEnt_IsMarked(reg,mkr)) {
	  MEnt_Mark(reg,mkr);
	  List_Add(vregions,reg);
	  nr++;
	}
      }
    }
    List_Delete(vfaces);
    List_Unmark(vregions,mkr);
    MSTK_FreeMarker(mkr);

    if (nr > 0)
      return vregions;
    else {
      List_Delete(vregions);
      return 0;
    }
  }

  void MV_RegionIDs_R4(MVertex_ptr v, int *nvr, int *vregionids) {
    int i;
    List_ptr vregions = MV_Regions_R4(v);
    if (vregions) {
      *nvr = List_Num_Entries(vregions);
      for (i = 0; i < *nvr; i++)
        vregionids[i] = MEnt_ID(List_Entry(vregions,i));
      List_Delete(vregions);
    }
    else
      *nvr = 0;
  }

  void MV_Add_AdjVertex_R4(MVertex_ptr v, MVertex_ptr adjv) {
    MVertex_Adj_R4 *adj;

    adj = (MVertex_Adj_R4 *) v->adj;
    List_ChknAdd(adj->adjverts,adjv);
  }

  void MV_Rem_AdjVertex_R4(MVertex_ptr v, MVertex_ptr adjv) {
    MVertex_Adj_R4 *adj;

    adj = (MVertex_Adj_R4 *) v->adj;
    List_Rem(adj->adjverts,adjv);
  }

  void MV_Add_Edge_R4(MVertex_ptr v, MEdge_ptr medge) {
#ifdef DEBUG
    /* This function is explicitly called from ME_Set_Vertex */
    /*MSTK_Report("MV_Add_Edge","Function call not suitable for this representation",MSTK_WARN);*/
#endif
  }

  void MV_Rem_Edge_R4(MVertex_ptr v, MEdge_ptr medge) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_Edge","Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void MV_Add_Face_R4(MVertex_ptr v, MFace_ptr mface) {
    MVertex_Adj_R4 *adj;

    adj = (MVertex_Adj_R4 *) v->adj;
    List_ChknAdd(adj->vfaces,mface);
  }

  void MV_Rem_Face_R4(MVertex_ptr v, MFace_ptr mface) {
    MVertex_Adj_R4 *adj;

    adj = (MVertex_Adj_R4 *) v->adj;

    List_Rem(adj->vfaces,mface);
  }

  void MV_Add_Region_R4(MVertex_ptr v, MRegion_ptr mregion) {
#ifdef DEBUG
    MSTK_Report("MV_Add_Region","Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void MV_Rem_Region_R4(MVertex_ptr v, MRegion_ptr mregion) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_Region","Function call not suitable for this representation",MSTK_WARN);
#endif
  }

#ifdef __cplusplus
}
#endif
