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
    upadj->velements = List_New(10);
  }

  void MV_Delete_R1(MVertex_ptr v, int keep) {
    MVertex_UpAdj_R1R2 *upadj;

    if (!keep) {
      upadj = (MVertex_UpAdj_R1R2 *) v->upadj;
      if (upadj) {
	if (upadj->velements)
	  List_Delete(upadj->velements);
	MSTK_free(upadj);
      }
    }
  }

  void MV_Restore_R1(MVertex_ptr v) {
    MEnt_Set_Dim(v,MVERTEX);
  }

  void MV_Destroy_For_MESH_Delete_R1(MVertex_ptr v) {
    MVertex_UpAdj_R1R2 *upadj;

    upadj = (MVertex_UpAdj_R1R2 *) v->upadj;
    if (upadj) {
      if (upadj->velements)
	List_Delete(upadj->velements);
      MSTK_free(upadj);
    }
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
    return MV_Num_Faces_R1R2(v);
  }
  
  int MV_Num_Regions_R1(MVertex_ptr v) {
    return MV_Num_Regions_R1R2(v);
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

  List_ptr MV_Edges_R1(MVertex_ptr v) {
    MVertex_UpAdj_R1R2 *upadj;
    int idx, idx1, idx2, found, nfv;
    MEntity_ptr ent;
    MFace_ptr lstface;
    MEdge_ptr edge, redge, lstedge;
    MVertex_ptr ev[2], vtmp;
    List_ptr redges, vedges, fverts;
    Mesh_ptr mesh = MEnt_Mesh(v);

    vedges = List_New(0);

    upadj = (MVertex_UpAdj_R1R2 *) v->upadj;
    idx = 0;
    while ((ent = (MEntity_ptr) List_Next_Entry(upadj->velements,&idx))) {
      if (MEnt_Dim(ent) == MREGION) {

	redges = MR_Edges(ent);

	idx1 = 0;
	while ((redge = List_Next_Entry(redges,&idx1))) {
	  if (ME_UsesEntity(redge,v,MVERTEX)) {
	    
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
	  if (ME_UsesEntity(lstedge,ev[0],MVERTEX) &&
	      ME_UsesEntity(lstedge,ev[1],MVERTEX)) {
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
	  if (ME_UsesEntity(lstedge,ev[0],MVERTEX) &&
	      ME_UsesEntity(lstedge,ev[1],MVERTEX)) {
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

  List_ptr MV_Faces_R1(MVertex_ptr v) {
    return MV_Faces_R1R2(v);
  }

  List_ptr MV_Regions_R1(MVertex_ptr v) {
    return MV_Regions_R1R2(v);
  }

  void MV_Add_Region_R1(MVertex_ptr v, MRegion_ptr mregion) {
    MV_Add_Region_R1R2(v,mregion);
  }

  void MV_Rem_Region_R1(MVertex_ptr v, MRegion_ptr mregion) {
    MV_Rem_Region_R1R2(v,mregion);
  }

  void MV_Add_Face_R1(MVertex_ptr v, MFace_ptr mface) {
    MV_Add_Face_R1R2(v,mface);
  }

  void MV_Rem_Face_R1(MVertex_ptr v, MFace_ptr mface) {
    MV_Rem_Face_R1R2(v,mface);
  }

  void MV_Add_AdjVertex_R1(MVertex_ptr v, MVertex_ptr av) {
#ifdef DEBUG
    MSTK_Report("MV_Add_AdjVertex",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Rem_AdjVertex_R1(MVertex_ptr v, MVertex_ptr av) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_AdjVertex",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Add_Edge_R1(MVertex_ptr v, MEdge_ptr e) {
#ifdef DEBUG
    /* This function is explicitly called from ME_Set_Vertex */
/*    MSTK_Report("MV_Add_Edge",
		"Function call not suitable for this representation",WARN);*/
#endif
  }

  void MV_Rem_Edge_R1(MVertex_ptr v, MEdge_ptr e) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_Edge",
		"Function call not suitable for this representation",WARN);
#endif
  }


#ifdef __cplusplus
}
#endif
