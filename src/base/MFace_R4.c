#define _H_MFace_Private

#include "MFace.h"
#include "MFace_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"
 
#ifdef __cplusplus
extern "C" {
#endif


  void MF_Set_RepType_R4(MFace_ptr f) {
    MFace_Adj_R4 *adj;

    adj = f->adj = (MFace_Adj_R4 *) MSTK_malloc(sizeof(MFace_Adj_R4));
    adj->fregions[0] = (MRegion_ptr) NULL;
    adj->fregions[1] = (MRegion_ptr) NULL;
    adj->fvertices = NULL;
  }

  void MF_Delete_R4(MFace_ptr f, int keep) {
    MFace_Adj_R4 *adj;
    int idx;
    MVertex_ptr v;

    adj = (MFace_Adj_R4 *) f->adj;

    if (MEnt_Dim(f) != MDELETED) { /* if face hasnt been temporarily deleted */
      if (adj) {
	idx = 0;
	while ((v = List_Next_Entry(adj->fvertices,&idx))) 
	  MV_Rem_Face(v,f);
      }
    }

    if (!keep) {
      if (adj) {
	if (adj->fvertices)
	  List_Delete(adj->fvertices);
	MSTK_free(adj);
      }
    }
  }

  void MF_Restore_R4(MFace_ptr f) {
    MFace_Adj_R4 *adj;
    int idx;
    MVertex_ptr v;

    MEnt_Set_Dim(f,MFACE);

    adj = (MFace_Adj_R4 *) f->adj;
    idx = 0;
    while ((v = List_Next_Entry(adj->fvertices,&idx))) 
      MV_Add_Face(v,f);
  }
  
  void MF_Destroy_For_MESH_Delete_R4(MFace_ptr f) {
    MFace_Adj_R4 *adj;

    adj = (MFace_Adj_R4 *) f->adj;
    if (adj) {
      if (adj->fvertices)
	List_Delete(adj->fvertices);
      MSTK_free(adj);
    }
  }

  int MF_Set_GInfo_Auto_R4(MFace_ptr f) {
    return MF_Set_GInfo_Auto_RN(f);
  }

  void MF_Set_Edges_R4(MFace_ptr f, int n, MEdge_ptr *e, int *dir) {
#ifdef DEBUG
    MSTK_Report("MF_Set_Edges",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Replace_Edges_i_R4(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges) {
#ifdef DEBUG
    MSTK_Report("MF_Replace_Edges",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Replace_Edges_R4(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges) {
#ifdef DEBUG
    MSTK_Report("MF_Replace_Edge",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Set_Vertices_R4(MFace_ptr f, int n, MVertex_ptr *v) {
    MF_Set_Vertices_RN(f,n,v);
  }

  void MF_Replace_Vertex_R4(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv) {
    MF_Replace_Vertex_RN(f,v,nuv);
  }

  void MF_Replace_Vertex_i_R4(MFace_ptr f, int i, MVertex_ptr nuv) {
    MF_Replace_Vertex_i_RN(f,i,nuv);
  }

  void MF_Insert_Vertex_R4(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v) {
    MF_Insert_Vertex_RN(f,nuv,b4v);
  }

  void MF_Insert_Vertex_i_R4(MFace_ptr f, MVertex_ptr nuv, int i) {
    MF_Insert_Vertex_i_RN(f,nuv,i);
  }

  int MF_Rev_EdgeDir_R4(MFace_ptr f, MEdge_ptr e) {
    return MF_Rev_EdgeDir_RN(f,e);
  }

  int MF_Rev_EdgeDir_i_R4(MFace_ptr f, int i) {
    return MF_Rev_EdgeDir_i_RN(f,i);
  }

  void MF_Add_Region_R4(MFace_ptr f, MRegion_ptr r, int side) {
    MSTK_Report("MF_Add_Region_R4", 
		"Not yet implemented for this representation",WARN);
  }

  void MF_Rem_Region_R4(MFace_ptr f, MRegion_ptr r) {
    MSTK_Report("MF_Rem_Region_R4",
		"Not yet implemented for this representation",WARN);
  }

  int MFs_AreSame_R4(MFace_ptr f1, MFace_ptr f2) {
    return (f1 == f2);
  }

  int MF_Num_Vertices_R4(MFace_ptr f) {
    return MF_Num_Vertices_RN(f);
  }

  int MF_Num_Edges_R4(MFace_ptr f) {
    return MF_Num_Edges_RN(f);
  }

  List_ptr MF_Vertices_R4(MFace_ptr f, int dir, MVertex_ptr v0) {
    return MF_Vertices_RN(f,dir,v0);
  }
	

  List_ptr MF_Edges_R4(MFace_ptr f, int dir, MVertex_ptr v0) {
    return MF_Edges_RN(f,dir,v0); 
  }

  int MF_EdgeDir_R4(MFace_ptr f, MEdge_ptr e) {
    return MF_EdgeDir_RN(f,e);
  }

  int MF_EdgeDir_i_R4(MFace_ptr f, int i) {
    return MF_EdgeDir_i_RN(f,i);
  }

  int MF_UsesEdge_R4(MFace_ptr f, MEdge_ptr e) {
    return MF_UsesEdge_RN(f,e);
  }

  int MF_UsesVertex_R4(MFace_ptr f, MVertex_ptr v) {
    return MF_UsesVertex_RN(f,v);
  }

  List_ptr MF_Regions_R4(MFace_ptr f) {
    return MF_Regions_R3R4(f);
  }

  MRegion_ptr MF_Region_R4(MFace_ptr f, int dir) {
    return MF_Region_R3R4(f,dir);
  }


  int MF_Num_AdjFaces_R4(MFace_ptr f) {
    MSTK_Report("MF_Num_AdjFaces_R4",
		"Not yet implemented for this representation",WARN);
    return 0;
  }

  List_ptr MF_AdjFaces_R4(MFace_ptr f) {
    MSTK_Report("MF_AdjFaces_R4",
		"Not yet implemented for this representation",WARN);
    return 0;
  }

  void MF_Add_AdjFace_R4(MFace_ptr f, int side, MFace_ptr af) {
#ifdef DEBUG
    MSTK_Report("MF_Add_AdjFace_R4",
		"Not yet implemented for this representation",WARN);
#endif
  }

  void MF_Rem_AdjFace_R4(MFace_ptr f, int side, MFace_ptr af) {
#ifdef DEBUG
    MSTK_Report("MF_Rem_AdjFace_R4",
		"Not yet implemented for this representation",WARN);
#endif
  }

  MFace_ptr MF_NextInHash_R4(MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("MF_NextInHash", "Function call not suitable for this representation", WARN);
#endif
    return NULL;
  }

  void MF_Set_NextInHash_R4(MFace_ptr f, MFace_ptr next) {
#ifdef DEBUG
    MSTK_Report("MF_NextInHash", "Function call not suitable for this representation", WARN);
#endif
  }

  void MF_HashKey_R4(MFace_ptr f, unsigned int *pn, void* **pp) {
#ifdef DEBUG
    MSTK_Report("MF_NextInHash", "Function call not suitable for this representation", WARN);
#endif
  }

  int MF_IsLocked_R4(MFace_ptr f) {
    return 0;
  }

#ifdef __cplusplus
}
#endif
