#define _H_MFace_Private

#include "MFace.h"
#include "MFace_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MF_Set_RepType_R1(MFace_ptr f) {
    MFace_Adj_R1 *adj;

    adj = f->adj = (MFace_Adj_R1 *) MSTK_malloc(sizeof(MFace_Adj_R1));
    adj->fvertices = NULL;
    adj->hnext = NULL;
    adj->lock = 0;
  }


  void MF_Delete_R1(MFace_ptr f, int keep) {
    MFace_Adj_R1 *adj;

    adj = (MFace_Adj_R1 *) f->adj;

    if (adj) {
      if (adj->fvertices) List_Delete(adj->fvertices);
      MSTK_free(adj);
    }

    if (!keep) {
      MEnt_Set_DelFlag(f);
    }
    if (keep) 
      MSTK_Report("MF_Delete_R1","Deletion of faces is permanent in this representation",ERROR);
  }

  void MF_Restore_R1(MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("MF_Restore_R1",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Destroy_For_MESH_Delete_R1(MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("MF_Destroy_For_MESH_Delete_R1",
		"Function call not suitable for this representation",ERROR);
#endif
  }

  int MF_Set_GInfo_Auto_R1(MFace_ptr f) {
    return MF_Set_GInfo_Auto_RN(f);
  }

  void MF_Set_Edges_R1(MFace_ptr f, int n, MEdge_ptr *e, int *dir) {
#ifdef DEBUG
    MSTK_Report("MF_Set_Edges_R1",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Replace_Edges_R1(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges) {
#ifdef DEBUG
    MSTK_Report("MF_Replace_Edges_R1",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Replace_Edges_i_R1(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges) {
#ifdef DEBUG
    MSTK_Report("MF_Replace_Edges_i_R1",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Set_Vertices_R1(MFace_ptr f, int n, MVertex_ptr *v) {
    MF_Set_Vertices_RN(f,n,v);
  }

  void MF_Replace_Vertex_R1(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv) {
#ifdef DEBUG
    MSTK_Report("MF_Replace_Vertex_R1","Modifying a temporary entity",WARN);
#endif
    MF_Replace_Vertex_RN(f,v,nuv);
  }

  void MF_Replace_Vertex_i_R1(MFace_ptr f, int i, MVertex_ptr nuv) {
#ifdef DEBUG
    MSTK_Report("MF_Replace_Vertex_i_R1","Modifying a temporary entity",WARN);
#endif
    MF_Replace_Vertex_i_RN(f,i,nuv);
  }

  void MF_Insert_Vertex_R1(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v) {
#ifdef DEBUG
    MSTK_Report("MF_Insert_Vertex_R1","Modifying a temporary entity",WARN);
#endif
    MF_Insert_Vertex_RN(f,nuv,b4v);
  }

  void MF_Insert_Vertex_i_R1(MFace_ptr f, MVertex_ptr nuv, int i) {
#ifdef DEBUG
    MSTK_Report("MF_Insert_Vertex_i_R1","Modifying a temporary entity",WARN);
#endif
    MF_Insert_Vertex_i_RN(f,nuv,i);
  }

  void MF_Add_Region_R1(MFace_ptr f, MRegion_ptr r, int side) {
#ifdef DEBUG
    MSTK_Report("MF_Add_Region_R1",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Rem_Region_R1(MFace_ptr f, MRegion_ptr r) {
#ifdef DEBUG
    MSTK_Report("MF_Rem_Region_R1",
		"Function call not suitable for this representation",WARN);
#endif
  }

  int MFs_AreSame_R1(MFace_ptr f1, MFace_ptr f2) {
    if (f1 == f2)
      return 1;
    else
      return MFs_AreSame_R1R2(f1,f2);
  }

  int MF_Num_Vertices_R1(MFace_ptr f) {
    return MF_Num_Vertices_RN(f);
  }

  int MF_Num_Edges_R1(MFace_ptr f) {
    return MF_Num_Edges_RN(f);
  }

  List_ptr MF_Vertices_R1(MFace_ptr f, int dir, MVertex_ptr v0) {
    return MF_Vertices_RN(f,dir,v0);
  }
	

  List_ptr MF_Edges_R1(MFace_ptr f, int dir, MVertex_ptr v0) {
    return MF_Edges_RN(f,dir,v0); 
  }

  int MF_EdgeDir_R1(MFace_ptr f, MEdge_ptr e) {
    return MF_EdgeDir_RN(f,e);
  }

  int MF_EdgeDir_i_R1(MFace_ptr f, int i) {
    return MF_EdgeDir_i_RN(f,i);
  }

  int MF_UsesEdge_R1(MFace_ptr f, MEdge_ptr e) {
    return MF_UsesEdge_RN(f,e);
  }

  int MF_UsesVertex_R1(MFace_ptr f, MVertex_ptr v) {
    return MF_UsesVertex_RN(f,v);
  }

  List_ptr MF_Regions_R1(MFace_ptr f) {
    return MF_Regions_R1R2(f);
  }

  MRegion_ptr MF_Region_R1(MFace_ptr f, int dir) {
    return MF_Region_R1R2(f,dir);
  }

  int MF_Num_AdjFaces_R1(MFace_ptr f) {
    MSTK_Report("MF_Num_AdjFaces_R1",
		"Not yet implemented for this representation",WARN);
    return 0;
  }

  List_ptr MF_AdjFaces_R1(MFace_ptr f) {
    MSTK_Report("MF_AdjFaces_R1",
		"Not yet implemented for this representation",WARN);
    return 0;
  }

  void MF_Add_AdjFace_R1(MFace_ptr f, int side, MFace_ptr af) {
#ifdef DEBUG
    MSTK_Report("MF_Add_AdjFace_R1",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Rem_AdjFace_R1(MFace_ptr f, int side, MFace_ptr af) {
#ifdef DEBUG
    MSTK_Report("MF_Rem_AdjFace_R1",
		"Function call not suitable for this representation",WARN);
#endif
  }

  MFace_ptr MF_NextInHash_R1(MFace_ptr f) {
    MFace_ptr hnext = ((MFace_Adj_R1 *)f->adj)->hnext; 
    return hnext;
  }

  void MF_Set_NextInHash_R1(MFace_ptr f, MFace_ptr next) {
    MFace_Adj_R1 *adj = (MFace_Adj_R1 *)f->adj;
    adj->hnext = next;
  }

  void MF_HashKey_R1(MFace_ptr f, unsigned int *pn, void* **pp) {
    MFace_Adj_R1 *adj = (MFace_Adj_R1 *)f->adj;
    *pn = List_Num_Entries(adj->fvertices);
    *pp = List_Entries(adj->fvertices);
  }

  void MF_Lock_R1(MFace_ptr f) {
    MFace_Adj_R1 *adj = (MFace_Adj_R1 *)f->adj;
    Hash_Lock(&adj->lock);
  }

  void MF_UnLock_R1(MFace_ptr f) {
    MFace_Adj_R1 *adj = (MFace_Adj_R1 *)f->adj;
    Hash_UnLock(&adj->lock);
  }

  int MF_IsLocked_R1(MFace_ptr f) {
    MFace_Adj_R1 *adj = (MFace_Adj_R1 *)f->adj;
    return Hash_IsLocked(adj->lock);
  }

#ifdef __cplusplus
}
#endif
