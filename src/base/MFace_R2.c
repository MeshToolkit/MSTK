#define _H_MFace_Private

#include "MFace.h"
#include "MFace_jmp.h"
#include "MSTK_malloc.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MF_Set_RepType_R2(MFace_ptr f) {
    /* What info should be created for virtual entities */
    MSTK_Report("MF_Set_RepType_R2","Not Implemented",ERROR);
  }

  void MF_Delete_R2(MFace_ptr f, int keep) {
    MSTK_Report("MF_Delete_R2","Not Implemented",ERROR);
  }

  void MF_Set_Edges_R2(MFace_ptr f, int n, MEdge_ptr *e, int *dir) {
#ifdef DEBUG
    MSTK_Report("MF_Set_Edges","Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Replace_Edge_R2(MFace_ptr f, MEdge_ptr e, int nnu, MEdge_ptr *nuedges, int *nudirs) {
#ifdef DEBUG
    MSTK_Report("MF_Replace_Edge","Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Replace_Edge_i_R2(MFace_ptr f, int i, int nnu, MEdge_ptr *nuedges, int *nudirs) {
#ifdef DEBUG
    MSTK_Report("MF_Replace_Edge","Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Set_Vertices_R2(MFace_ptr f, int n, MVertex_ptr *v) {
    MSTK_Report("MF_Set_Vertices","Not implemented for this representation",WARN);
  }

  void MF_Replace_Vertex_R2(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv) {
    MSTK_Report("MF_Replace_Vertex","Not implemented for this representation",WARN);
  }

  void MF_Replace_Vertex_i_R2(MFace_ptr f, int i, MVertex_ptr v) {
    MSTK_Report("MF_Replace_Vertex","Not implemented for this representation",WARN);
  }

  void MF_Insert_Vertex_R2(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v) {
    MSTK_Report("MF_Insert_Vertex","Not implemented for this representation",WARN);
  }

  void MF_Insert_Vertex_i_R2(MFace_ptr f, MVertex_ptr nuv, int i) {
    MSTK_Report("MF_Insert_Vertex_i","Not implemented for this representation",WARN);
  }

  void MF_Add_Region_R2(MFace_ptr f, MRegion_ptr r, int side) {
    MSTK_Report("MF_Add_Region_R2", "Not Implemented for this representation",WARN);
  }

  void MF_Rem_Region_R2(MFace_ptr f, MRegion_ptr r) {
    MSTK_Report("MF_Rem_Region_R2","Not Implemented for this representation",WARN);
  }

  int MF_Num_Vertices_R2(MFace_ptr f) {
    return MF_Num_Vertices_R1R2(f);
  }

  int MF_Num_Edges_R2(MFace_ptr f) {
    return MF_Num_Edges_R1R2(f);
  }

  List_ptr MF_Vertices_R2(MFace_ptr f, int dir, MVertex_ptr v0) {
    return MF_Vertices_R1R2(f,dir,v0);
  }
	

  List_ptr MF_Edges_R2(MFace_ptr f, int dir, MVertex_ptr v0) {
    return MF_Edges_R1R2(f,dir,v0); 
  }

  int MF_EdgeDir_R2(MFace_ptr f, MEdge_ptr e) {
    return MF_EdgeDir_R1R2(f,e);
  }

  int MF_EdgeDir_i_R2(MFace_ptr f, int i) {
    return MF_EdgeDir_i_R1R2(f,i);
  }

  int MF_UsesEdge_R2(MFace_ptr f, MEdge_ptr e) {
    return MF_UsesEdge_R1R2(f,e);
  }

  int MF_UsesVertex_R2(MFace_ptr f, MVertex_ptr v) {
    return MF_UsesVertex_R1R2(f,v);
  }

  List_ptr MF_Regions_R2(MFace_ptr f) {
    return MF_Regions_R1R2(f);
  }

  MRegion_ptr MF_Region_R2(MFace_ptr f, int dir) {
    return MF_Region_R1R2(f,dir);
  }


  int MF_Num_AdjFaces_R2(MFace_ptr f) {
    MSTK_Report("MF_Num_AdjFaces_R2","Not yet implemented for this representation",WARN);
    return 0;
  }

  List_ptr MF_AdjFaces_R2(MFace_ptr f) {
    MSTK_Report("MF_AdjFaces_R2","Not yet implemented for this representation",WARN);
    return 0;
  }

  void MF_Add_AdjFace_R2(MFace_ptr f, int side, MFace_ptr af) {
#ifdef DEBUG
    MSTK_Report("MF_Add_AdjFace_R2","Not yet implemented for this representation",WARN);
#endif
  }

  void MF_Rem_AdjFace_R2(MFace_ptr f, int side, MFace_ptr af) {
#ifdef DEBUG
    MSTK_Report("MF_Rem_AdjFace_R2","Not yet implemented for this representation",WARN);
#endif
  }

#ifdef __cplusplus
}
#endif
