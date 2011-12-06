#define _H_MVertex_Private

#include "MVertex.h"
#include "MVertex_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MV_Set_RepType_R2(MVertex_ptr v) {
    MVertex_Adj_R2 *adj;

    adj = v->adj = (MVertex_Adj_R2 *) MSTK_malloc(sizeof(MVertex_Adj_R2));
    adj->velements = List_New(10);
    adj->adjverts = List_New(10);
  }

  void MV_Delete_R2(MVertex_ptr v, int keep) {
    MVertex_Adj_R2 *adj;
    int idx;
    MVertex_ptr adjv;

    adj = (MVertex_Adj_R2 *) v->adj;

    if (MEnt_Dim((MEntity_ptr) v) != MDELETED) { /* if vtx has not been temporarily deleted*/
      if (adj) {
	idx = 0;
	while ((adjv = List_Next_Entry(adj->adjverts,&idx))) 
	  MV_Rem_AdjVertex_R2(adjv,v);
      }
    }

    if (!keep) {
      adj = (MVertex_Adj_R2 *) v->adj;
      if (adj) {
	if (adj->velements) 
	  List_Delete(adj->velements);
        if (adj->adjverts)
	  List_Delete(adj->adjverts); 
	MSTK_free(adj);
      }
    }
  }
    
  void MV_Restore_R2(MVertex_ptr v) {
    MVertex_Adj_R2 *adj;
    int idx;
    MVertex_ptr adjv;

    MEnt_Set_Dim((MEntity_ptr) v,MVERTEX);

    adj = (MVertex_Adj_R2 *) v->adj;
    if (adj) {
      idx = 0;
      while ((adjv = List_Next_Entry(adj->adjverts,&idx)))
	MV_Add_AdjVertex_R2(adjv,v);
    }
  }    

  void MV_Destroy_For_MESH_Delete_R2(MVertex_ptr v) {
    MVertex_Adj_R2 *adj;

    adj = (MVertex_Adj_R2 *) v->adj;
    if (adj) {
      if (adj->velements)
	List_Delete(adj->velements);
      if (adj->adjverts)
        List_Delete(adj->adjverts);
      MSTK_free(adj);
    }
  }
    
  int MV_Num_AdjVertices_R2(MVertex_ptr v) {
    List_ptr adjverts = ((MVertex_Adj_R2 *) v->adj)->adjverts;
    return List_Num_Entries(adjverts);
  }

  int MV_Num_Edges_R2(MVertex_ptr v) {
    List_ptr adjverts = ((MVertex_Adj_R2 *) v->adj)->adjverts;
    return List_Num_Entries(adjverts);
  }

  int MV_Num_Faces_R2(MVertex_ptr v) {
    int i, nf;
    List_ptr vfaces;
    MFace_ptr vface;

#ifdef DEBUG
    MSTK_Report("MV_Num_Faces_R2",
		"Inefficient to call this routine with this representation",
		MSTK_MESG);
#endif
    vfaces = MV_Faces_R2(v);
    nf = List_Num_Entries(vfaces);
    for (i = 0; i < nf; i++) {
      vface = List_Entry(vfaces,i);
      if (MEnt_IsVolatile(vface))
	MF_Delete(vface,0);
    }
    List_Delete(vfaces);

    return nf;
  }
  
  int MV_Num_Regions_R2(MVertex_ptr v) {
    MVertex_Adj_R2 *adj;
    int idx, nr = 0;
    MEntity_ptr ent;

    adj = (MVertex_Adj_R2 *) v->adj;
    idx = 0;
    while ((ent = (MEntity_ptr) List_Next_Entry(adj->velements,&idx))) {
      if (MEnt_Dim(ent) == MREGION)
	nr++;
    }
    return nr;
  }

  List_ptr MV_AdjVertices_R2(MVertex_ptr v) {
    MVertex_Adj_R2 *adj;
    adj = (MVertex_Adj_R2 *) v->adj;

    return List_Copy(adj->adjverts);
  }

  List_ptr MV_Edges_R2(MVertex_ptr v) {
    int idx, ne;
    MVertex_ptr adjv[2], vtmp;
    MEdge_ptr e;
    List_ptr vedges;
    MVertex_Adj_R2 *adj = (MVertex_Adj_R2 *) v->adj;
    Mesh_ptr mesh = MEnt_Mesh((MEntity_ptr) v);

    if (!adj->adjverts)
      return NULL;

    ne = List_Num_Entries(adj->adjverts);
    vedges = List_New(ne);
    idx = 0;
    adjv[0] = v;
    while ((adjv[1] = List_Next_Entry(adj->adjverts,&idx))) {
#ifdef HASHTABLE
      if (adjv[0]>adjv[1]) {
	vtmp = adjv[0];
	adjv[0] = adjv[1];
	adjv[1] = vtmp;
      }

      e = Hash_Entry(MESH_Hash_Edges(mesh), 2, (void**)adjv);
      if (e == NULL) {
	e = ME_New(mesh);
	MEnt_Set_Volatile(e);
	ME_Set_Vertex(e,0,adjv[0]);
	ME_Set_Vertex(e,1,adjv[1]);
	ME_Set_GInfo_Auto(e);
	Hash_Add(MESH_Hash_Edges(mesh), e, 2, (void**)adjv);
      }
#else
      e = ME_New(mesh);
      MEnt_Set_Volatile(e);
      ME_Set_Vertex(e,0,adjv[0]);
      ME_Set_Vertex(e,1,adjv[1]);
      ME_Set_GInfo_Auto(e);
#endif
      List_Add(vedges,e);
      ME_Lock(e);
    }

    if (!MESH_AutoLock(mesh)) {
       idx = 0;
       while ((e = List_Next_Entry(vedges, &idx))) {
	 ME_UnLock(e);
       }
    }
    return vedges;
  }

  List_ptr MV_Faces_R2(MVertex_ptr v) {
    MVertex_Adj_R2 *adj;
    int idx, idx1, idx2, found;
    MEntity_ptr ent;
    MFace_ptr rface, lstface;
    List_ptr rfaces, vfaces;

    vfaces = List_New(0);

    adj = (MVertex_Adj_R2 *) v->adj;
    idx = 0;
    while ((ent = (MEntity_ptr) List_Next_Entry(adj->velements,&idx))) {

      if (MEnt_Dim(ent) == MREGION) {

	rfaces = MR_Faces(ent);

	idx1 = 0;
	while ((rface = List_Next_Entry(rfaces,&idx1))) {
	  if (MF_UsesEntity(rface,(MEntity_ptr) v,MVERTEX)) {
	    
	    idx2 = 0; found = 0;
	    while ((lstface = List_Next_Entry(vfaces,&idx2))) {
	      if (MFs_AreSame(rface,lstface)) {
		found = 1;
		break;
	      }
	    }

	    if (!found)
	      List_Add(vfaces,rface);
	  }
	}
	
	List_Delete(rfaces);
      }
      else { /* Must be a face */
	List_Add(vfaces,ent);
      }
    }

    if (List_Num_Entries(vfaces))
      return vfaces;
    else {
      List_Delete(vfaces);
      return NULL;
    }
  }

  List_ptr MV_Regions_R2(MVertex_ptr v) {
    MVertex_Adj_R2 *adj;
    int idx, nel, nr = 0;
    MEntity_ptr ent;
    List_ptr vregions;

    adj = (MVertex_Adj_R2*) v->adj;
    nel = List_Num_Entries(adj->velements);
    vregions = List_New(nel);

    idx = 0;
    while ((ent = (MEntity_ptr) List_Next_Entry(adj->velements,&idx))) {
      if (MEnt_Dim(ent) == MREGION) {
	List_Add(vregions,ent);
	nr++;
      }
    }
    if (nr)
      return vregions;
    else {
      List_Delete(vregions);
      return 0;
    }      
  }

  void MV_Add_AdjVertex_R2(MVertex_ptr v, MVertex_ptr adjv) {
    MVertex_Adj_R2 *adj;
    
    adj = (MVertex_Adj_R2 *) v->adj;
    List_ChknAdd(adj->adjverts,adjv);
  }

  void MV_Rem_AdjVertex_R2(MVertex_ptr v, MVertex_ptr adjv) {
    MVertex_Adj_R2 *adj;
    
    adj = (MVertex_Adj_R2 *) v->adj;
    List_Rem(adj->adjverts,adjv);
  }


  void MV_Add_Edge_R2(MVertex_ptr v, MEdge_ptr e) {
#ifdef DEBUG
    /* This function is explicitly called from ME_Set_Vertex */
    /*MSTK_Report("MV_Add_Edge",
		"Function call not suitable for this representation",MSTK_WARN);*/
#endif
  }

  void MV_Rem_Edge_R2(MVertex_ptr v, MEdge_ptr e) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_Edge",
		"Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void MV_Add_Face_R2(MVertex_ptr v, MFace_ptr mface) {
    MVertex_Adj_R2 *adj;

    if (MF_Region(mface,0) || MF_Region(mface,1)) {
      MSTK_Report("MV_Add_Face_R2",
		 "Can only add faces with no regions in this representation",
		 MSTK_ERROR);
      return;
    }

    adj = (MVertex_Adj_R2 *)v->adj;
    List_ChknAdd(adj->velements,mface);
  }

  void MV_Rem_Face_R2(MVertex_ptr v, MFace_ptr mface) {
    MVertex_Adj_R2 *adj;

    if (MF_Region(mface,0) || MF_Region(mface,1)) {
      MSTK_Report("MV_Rem_Face_R2",
      "Set should contain only faces with no regions in this representation",
		 MSTK_ERROR);
      return;
    }

    adj = (MVertex_Adj_R2 *) v->adj;
    List_Rem(adj->velements,mface);
  }

  void MV_Add_Region_R2(MVertex_ptr v, MRegion_ptr mregion) {
    MVertex_Adj_R2 *adj;

    adj = (MVertex_Adj_R2 *) v->adj;
    List_ChknAdd(adj->velements,mregion);
  }

  void MV_Rem_Region_R2(MVertex_ptr v, MRegion_ptr mregion) {
   MVertex_Adj_R2 *adj;

    adj = (MVertex_Adj_R2 *) v->adj;
    List_Rem(adj->velements,mregion);
  }


#ifdef __cplusplus
}
#endif
