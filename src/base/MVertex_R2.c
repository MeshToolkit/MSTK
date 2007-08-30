#define _H_MVertex_Private

#include "MVertex.h"
#include "MVertex_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MV_Set_RepType_R2(MVertex_ptr v) {
    MVertex_UpAdj_R1R2 *upadj;
    MVertex_SameAdj_R2R4 *sameadj;

    upadj = v->upadj = (MVertex_UpAdj_R1R2 *) MSTK_malloc(sizeof(MVertex_UpAdj_R1R2));
    upadj->velements = List_New(10);
    sameadj = v->sameadj = (MVertex_SameAdj_R2R4 *) MSTK_malloc(sizeof(MVertex_SameAdj_R2R4));
    sameadj->adjverts = List_New(10);
  }

  void MV_Delete_R2(MVertex_ptr v, int keep) {
    MVertex_UpAdj_R1R2 *upadj;
    MVertex_SameAdj_R2R4 *sameadj;
    int idx;
    MVertex_ptr adjv;

    sameadj = (MVertex_SameAdj_R2R4 *) v->sameadj;

    if (MEnt_Dim(v) != MDELETED) { /* if vtx has not been temporarily deleted*/
      if (sameadj) {
	idx = 0;
	while ((adjv = List_Next_Entry(sameadj->adjverts,&idx))) 
	  MV_Rem_AdjVertex_R2(adjv,v);
      }
    }

    if (!keep) {
      upadj = (MVertex_UpAdj_R1R2 *) v->upadj;
      if (upadj) {
	if (upadj->velements) 
	  List_Delete(upadj->velements);
	MSTK_free(upadj);
      }
      
      if (sameadj) {
	List_Delete(sameadj->adjverts);
	MSTK_free(sameadj);
      }
    }
  }
    
  void MV_Restore_R2(MVertex_ptr v) {
    MVertex_SameAdj_R2R4 *sameadj;
    int idx;
    MVertex_ptr adjv;

    MEnt_Set_Dim(v,MVERTEX);

    sameadj = (MVertex_SameAdj_R2R4 *) v->sameadj;
    if (sameadj) {
      idx = 0;
      while ((adjv = List_Next_Entry(sameadj->adjverts,&idx)))
	MV_Add_AdjVertex_R2(adjv,v);
    }
  }    

  void MV_Destroy_For_MESH_Delete_R2(MVertex_ptr v) {
    MVertex_UpAdj_R1R2 *upadj;
    MVertex_SameAdj_R2R4 *sameadj;

    upadj = (MVertex_UpAdj_R1R2 *) v->upadj;
    if (upadj) {
      if (upadj->velements)
	List_Delete(upadj->velements);
      MSTK_free(upadj);
    }
      
    sameadj = (MVertex_SameAdj_R2R4 *) v->sameadj;
    if (sameadj) {
      if (sameadj->adjverts)
	List_Delete(sameadj->adjverts);
      MSTK_free(sameadj);
    }
  }
    
  int MV_Num_AdjVertices_R2(MVertex_ptr v) {
    List_ptr adjverts = ((MVertex_SameAdj_R2R4 *) v->sameadj)->adjverts;
    return List_Num_Entries(adjverts);
  }

  int MV_Num_Edges_R2(MVertex_ptr v) {
    List_ptr adjverts = ((MVertex_SameAdj_R2R4 *) v->sameadj)->adjverts;
    return List_Num_Entries(adjverts);
  }

  int MV_Num_Faces_R2(MVertex_ptr v) {
    return MV_Num_Faces_R1R2(v);
  }
  
  int MV_Num_Regions_R2(MVertex_ptr v) {
    return MV_Num_Regions_R1R2(v);
  }

  List_ptr MV_AdjVertices_R2(MVertex_ptr v) {
    MVertex_SameAdj_R2R4 *sameadj;
    sameadj = (MVertex_SameAdj_R2R4 *) v->sameadj;

    return List_Copy(sameadj->adjverts);
  }

  List_ptr MV_Edges_R2(MVertex_ptr v) {
    int idx, ne;
    MVertex_ptr adjv[2], vtmp;
    MEdge_ptr e;
    List_ptr vedges;
    MVertex_SameAdj_R2R4 *sameadj = (MVertex_SameAdj_R2R4 *) v->sameadj;
    Mesh_ptr mesh = MEnt_Mesh(v);

    if (!sameadj->adjverts)
      return NULL;

    ne = List_Num_Entries(sameadj->adjverts);
    vedges = List_New(ne);
    idx = 0;
    adjv[0] = v;
    while ((adjv[1] = List_Next_Entry(sameadj->adjverts,&idx))) {
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
    }

    return vedges;
  }

  List_ptr MV_Faces_R2(MVertex_ptr v) {
    return MV_Faces_R1R2(v);
  }

  List_ptr MV_Regions_R2(MVertex_ptr v) {
    return MV_Regions_R1R2(v);
  }

  void MV_Add_AdjVertex_R2(MVertex_ptr v, MVertex_ptr adjv) {
    MVertex_SameAdj_R2R4 *sameadj;
    
    sameadj = (MVertex_SameAdj_R2R4 *) v->sameadj;
    List_Add(sameadj->adjverts,adjv);
  }

  void MV_Rem_AdjVertex_R2(MVertex_ptr v, MVertex_ptr adjv) {
    MVertex_SameAdj_R2R4 *sameadj;
    
    sameadj = (MVertex_SameAdj_R2R4 *) v->sameadj;
    List_Rem(sameadj->adjverts,adjv);
  }


  void MV_Add_Edge_R2(MVertex_ptr v, MEdge_ptr e) {
#ifdef DEBUG
    MSTK_Report("MV_Add_Edge",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Rem_Edge_R2(MVertex_ptr v, MEdge_ptr e) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_Edge",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MV_Add_Face_R2(MVertex_ptr v, MFace_ptr mface) {
    MV_Add_Face_R1R2(v,mface);
  }

  void MV_Rem_Face_R2(MVertex_ptr v, MFace_ptr mface) {
    MV_Rem_Face_R1R2(v,mface);
   }

  void MV_Add_Region_R2(MVertex_ptr v, MRegion_ptr mregion) {
    MV_Add_Region_R1R2(v,mregion);
  }

  void MV_Rem_Region_R2(MVertex_ptr v, MRegion_ptr mregion) {
    MV_Rem_Region_R1R2(v,mregion);
  }


#ifdef __cplusplus
}
#endif
