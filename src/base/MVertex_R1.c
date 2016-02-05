#define _H_MVertex_Private

#include <stdlib.h>

#include "MVertex.h"
#include "MVertex_jmp.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MV_Set_RepType_R1(MVertex_ptr v) {
    MVertex_Adj_R1 *adj;

    adj = v->adj = (MVertex_Adj_R1 *) malloc(sizeof(MVertex_Adj_R1));
    adj->velements = List_New(10);
  }

  void MV_Delete_R1(MVertex_ptr v, int keep) {
    MVertex_Adj_R1 *adj;

    if (!keep) {
      adj = (MVertex_Adj_R1 *) v->adj;
      if (adj) {
	if (adj->velements)
	  List_Delete(adj->velements);
	free(adj);
      }
    }
  }

  void MV_Restore_R1(MVertex_ptr v) {
    MEnt_Set_Dim((MEntity_ptr) v,MVERTEX);
  }

  void MV_Destroy_For_MESH_Delete_R1(MVertex_ptr v) {
    MVertex_Adj_R1 *adj;

    adj = (MVertex_Adj_R1 *) v->adj;
    if (adj) {
      if (adj->velements)
	List_Delete(adj->velements);
      free(adj);
    }
  }

  int MV_Num_AdjVertices_R1(MVertex_ptr v) {
#ifdef DEBUG
    MSTK_Report("MV_Num_AdjVertices",
		"Inefficient to call this routine with this representation",
		MSTK_MESG);
#endif
    return MV_Num_Edges_R1(v);
  }

  int MV_Num_Edges_R1(MVertex_ptr v) {
    int ne;
    List_ptr vedges;

#ifdef DEBUG
    MSTK_Report("MV_Num_Edges",
		"Inefficient to call this routine with this representation",
		MSTK_MESG);
#endif
    
    vedges = MV_Edges_R1(v);
    ne = List_Num_Entries(vedges);
    List_Delete(vedges);
    return ne;
  }

  int MV_Num_Faces_R1(MVertex_ptr v) {
    int i, nf;
    List_ptr vfaces;
    MFace_ptr vface;

#ifdef DEBUG
    MSTK_Report("MV_Num_Faces_R1",
		"Inefficient to call this routine with this representation",
		MSTK_MESG);
#endif
    vfaces = MV_Faces_R1(v);
    nf = List_Num_Entries(vfaces);
    for (i = 0; i < nf; i++) {
      vface = List_Entry(vfaces,i);
      if (MEnt_IsVolatile(vface))
	MF_Delete(vface,0);
    }
    List_Delete(vfaces);

    return nf;
  }
  
  int MV_Num_Regions_R1(MVertex_ptr v) {
    MVertex_Adj_R1 *adj;
    int idx, nr = 0;
    MEntity_ptr ent;

    adj = (MVertex_Adj_R1 *) v->adj;
    idx = 0;
    while ((ent = (MEntity_ptr) List_Next_Entry(adj->velements,&idx))) {
      if (MEnt_Dim(ent) == MREGION)
	nr++;
    }
    return nr;
  }

  List_ptr MV_AdjVertices_R1(MVertex_ptr v) {
    List_ptr vedges = MV_Edges_R1(v);

    if (vedges) {      
      int idx = 0, nve = List_Num_Entries(vedges);
      List_ptr adjv = List_New(nve);
      MEdge_ptr edge;

      while ((edge = List_Next_Entry(vedges,&idx)))
	List_Add(adjv,ME_OppVertex(edge,v));
      List_Delete(vedges);
      
      return adjv;
    }
    else
      return NULL;
  }

  void MV_AdjVertexIDs_R1(MVertex_ptr v, int *nvadj, int *adjvids) {
    List_ptr vedges, adjv;
    int ne, i;
    MEdge_ptr vedge;
    MVertex_ptr ov;

    vedges = MV_Edges_R1(v);
    if (vedges) {
      *nvadj = List_Num_Entries(vedges);
      for (i = 0; i < *nvadj; i++) {
        vedge = List_Entry(vedges,i);
        ov = ME_OppVertex(vedge,v);
        adjvids[i] = MEnt_ID((MEntity_ptr)ov);
      }
      List_Delete(vedges);
    }
    else
      *nvadj = 0;

  }
    
  List_ptr MV_Edges_R1(MVertex_ptr v) {
    MVertex_Adj_R1 *adj;
    int idx, idx1, idx2, found, nfv;
    MEntity_ptr ent;
    MEdge_ptr edge, redge, lstedge;
    MVertex_ptr ev[2], vtmp;
    List_ptr redges, vedges, fverts;
    Mesh_ptr mesh = MEnt_Mesh((MEntity_ptr) v);

    vedges = List_New(0);

    adj = (MVertex_Adj_R1 *) v->adj;
    idx = 0;
    while ((ent = (MEntity_ptr) List_Next_Entry(adj->velements,&idx))) {
      if (MEnt_Dim(ent) == MREGION) {

	redges = MR_Edges(ent);

	idx1 = 0;
	while ((redge = List_Next_Entry(redges,&idx1))) {
	  if (ME_UsesEntity(redge,(MEntity_ptr) v,MVERTEX)) {
	    
	    idx2 = 0; found = 0;
	    while ((lstedge = List_Next_Entry(vedges,&idx2))) {
	      if (MEs_AreSame(redge,lstedge)) {
		found = 1;
		break;
	      }
	    }

	    if (!found) 
	      List_Add(vedges,redge);
	  }
	}
	
	List_Delete(redges);
      }
      else { /* Must be a face */
	fverts = MF_Vertices(ent,1,v);
	nfv = List_Num_Entries(fverts);

	ev[0] = List_Entry(fverts,0);
	ev[1] = List_Entry(fverts,1);

	idx2 = 0; found = 0;
	while ((lstedge = List_Next_Entry(vedges,&idx2))) {
	  if (ME_UsesEntity(lstedge,(MEntity_ptr) ev[0],MVERTEX) &&
	      ME_UsesEntity(lstedge,(MEntity_ptr) ev[1],MVERTEX)) {
	    found = 1;
	    break;
	  }
	}

	if (!found) {
#ifdef HASTABLE
	  if (ev[0]>ev[1]) {
	    vtmp = ev[0];
	    ev[0] = ev[1];
	    ev[1] = vtmp;
	  }

	  edge = Hash_Entry(MESH_Hash_Edges(mesh), 2, (void**)ev);
	  if (edge == NULL) {
	    edge = ME_New(mesh);
	    ME_Set_Vertex(edge,0,ev[0]);
	    ME_Set_Vertex(edge,1,ev[1]);
	    ME_Set_GInfo_Auto(edge);
	    Hash_Add(MESH_Hash_Edges(mesh), edge, 2, (void**)ev);
	  }
#else
	  edge = ME_New(mesh);
	  ME_Set_Vertex(edge,0,ev[0]);
	  ME_Set_Vertex(edge,1,ev[1]);
	  ME_Set_GInfo_Auto(edge);
#endif
	  List_Add(vedges,edge);
	  ME_Lock(edge);
	}

	ev[1] = ev[0];
	ev[0] = List_Entry(fverts,nfv-1);

	idx2 = 0; found = 0;
	while ((lstedge = List_Next_Entry(vedges,&idx2))) {
	  if (ME_UsesEntity(lstedge,(MEntity_ptr) ev[0],MVERTEX) &&
	      ME_UsesEntity(lstedge,(MEntity_ptr) ev[1],MVERTEX)) {
	    found = 1;
	    break;
	  }
	}

	if (!found) {
#ifdef HASHTABLE
	  if (ev[0]>ev[1]) {
	    vtmp = ev[0];
	    ev[0] = ev[1];
	    ev[1] = vtmp;
	  }

	  edge = Hash_Entry(MESH_Hash_Edges(mesh), 2, (void**)ev);
	  if (edge == NULL) {
	    edge = ME_New(mesh);
	    ME_Set_Vertex(edge,0,ev[0]);
	    ME_Set_Vertex(edge,1,ev[1]);
	    ME_Set_GInfo_Auto(edge);
	    Hash_Add(MESH_Hash_Edges(mesh), edge, 2, (void**)ev);
	  }
#else
	  edge = ME_New(mesh);
	  ME_Set_Vertex(edge,0,ev[0]);
	  ME_Set_Vertex(edge,1,ev[1]);
	  ME_Set_GInfo_Auto(edge);
#endif
	  List_Add(vedges,edge);
	  ME_Lock(edge);
	}

	List_Delete(fverts);
      }
    }

    if (!MESH_AutoLock(mesh)) {
       idx = 0;
       while ((edge = List_Next_Entry(vedges, &idx))) {
	 ME_UnLock(edge);
       }
    }

    if (List_Num_Entries(vedges)) 
      return vedges;
    else {
      List_Delete(vedges);
      return NULL;
    }
      
  }

  void MV_EdgeIDs_R1(MVertex_ptr v, int *nve, int *vedgeids) {
    List_ptr vedges = MV_Edges_R1(v);

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


  List_ptr MV_Faces_R1(MVertex_ptr v) {
    MVertex_Adj_R1 *adj;
    int idx, idx1, idx2, found;
    MEntity_ptr ent;
    MFace_ptr rface, lstface;
    List_ptr rfaces, vfaces;

    vfaces = List_New(0);

    adj = (MVertex_Adj_R1 *) v->adj;
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

  void MV_FaceIDs_R1(MVertex_ptr v, int *nvf, int *vfaceids) {
    List_ptr vfaces = MV_Faces_R1(v);
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
  

  List_ptr MV_Regions_R1(MVertex_ptr v) {
    MVertex_Adj_R1 *adj;
    int idx, nel, nr = 0;
    MEntity_ptr ent;
    List_ptr vregions;

    adj = (MVertex_Adj_R1 *) v->adj;
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

  void MV_RegionIDs_R1(MVertex_ptr v, int *nvr, int *vregionids) {
    MVertex_Adj_R1 *adj;
    int idx;
    MEntity_ptr ent;

    adj = (MVertex_Adj_R1 *) v->adj;

    *nvr = 0;
    idx = 0;
    while ((ent = (MEntity_ptr) List_Next_Entry(adj->velements,&idx))) {
      if (MEnt_Dim(ent) == MREGION)
	vregionids[(*nvr)++] = MEnt_ID(ent);
    }
  }

  void MV_Add_Region_R1(MVertex_ptr v, MRegion_ptr mregion) {
    MVertex_Adj_R1 *adj;

    adj = (MVertex_Adj_R1 *) v->adj;
    List_ChknAdd(adj->velements,mregion);
  }

  void MV_Rem_Region_R1(MVertex_ptr v, MRegion_ptr mregion) {
   MVertex_Adj_R1 *adj;

    adj = (MVertex_Adj_R1 *) v->adj;
    List_Rem(adj->velements,mregion);
  }

  void MV_Add_Face_R1(MVertex_ptr v, MFace_ptr mface) {
    MVertex_Adj_R1 *adj;

    if (MF_Region(mface,0) || MF_Region(mface,1)) {
      MSTK_Report("MV_Add_Face_R1",
		 "Can only add faces with no regions in this representation",
		 MSTK_ERROR);
      return;
    }

    adj = (MVertex_Adj_R1 *)v->adj;
    List_ChknAdd(adj->velements,mface);
  }

  void MV_Rem_Face_R1(MVertex_ptr v, MFace_ptr mface) {
    MVertex_Adj_R1 *adj;

    if (MF_Region(mface,0) || MF_Region(mface,1)) {
      MSTK_Report("MV_Rem_Face_R1",
      "Set should contain only faces with no regions in this representation",
		 MSTK_ERROR);
      return;
    }

    adj = (MVertex_Adj_R1 *) v->adj;
    List_Rem(adj->velements,mface);
  }

  void MV_Add_AdjVertex_R1(MVertex_ptr v, MVertex_ptr av) {
#ifdef DEBUG
    MSTK_Report("MV_Add_AdjVertex",
		"Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void MV_Rem_AdjVertex_R1(MVertex_ptr v, MVertex_ptr av) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_AdjVertex",
		"Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void MV_Add_Edge_R1(MVertex_ptr v, MEdge_ptr e) {
#ifdef DEBUG
    /* This function is explicitly called from ME_Set_Vertex */
/*    MSTK_Report("MV_Add_Edge",
		"Function call not suitable for this representation",MSTK_WARN);*/
#endif
  }

  void MV_Rem_Edge_R1(MVertex_ptr v, MEdge_ptr e) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_Edge",
		"Function call not suitable for this representation",MSTK_WARN);
#endif
  }


#ifdef __cplusplus
}
#endif
