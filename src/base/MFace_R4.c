#define _H_MFace_Private

#include "MFace.h"
#include "MFace_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"
 
#ifdef __cplusplus
extern "C" {
#endif


  void MF_Set_RepType_R4(MFace_ptr f) {
    MFace_UpAdj_R3R4 *upadj;
    MFace_DownAdj_RN *downadj;

    upadj = f->upadj = (MFace_UpAdj_R3R4 *) MSTK_malloc(sizeof(MFace_UpAdj_R3R4));
    upadj->fregions[0] = (MRegion_ptr) NULL;
    upadj->fregions[1] = (MRegion_ptr) NULL;
    downadj = f->downadj = (MFace_DownAdj_RN *) MSTK_malloc(sizeof(MFace_DownAdj_RN));
    downadj->fvertices = NULL;
  }

  void MF_Delete_R4(MFace_ptr f, int keep) {
    MFace_UpAdj_R3R4 *upadj;
    MFace_DownAdj_RN *downadj;
    int idx;
    MVertex_ptr v;

    downadj = (MFace_DownAdj_RN *) f->downadj;

    if (MEnt_Dim(f) != MDELETED) { /* if face hasnt been temporarily deleted */
      if (downadj) {
	idx = 0;
	while ((v = List_Next_Entry(downadj->fvertices,&idx))) 
	  MV_Rem_Face(v,f);
      }
    }

    if (!keep) {
      upadj  = (MFace_UpAdj_R3R4 *) f->upadj;
      if (upadj)
	MSTK_free(upadj);

      if (downadj) {
	if (downadj->fvertices)
	  List_Delete(downadj->fvertices);
	MSTK_free(downadj);
      }
    }
  }

  void MF_Restore_R4(MFace_ptr f) {
    MFace_DownAdj_RN *downadj;
    int idx;
    MVertex_ptr v;

    MEnt_Set_Dim(f,MFACE);

    downadj = (MFace_DownAdj_RN *) f->downadj;
    idx = 0;
    while ((v = List_Next_Entry(downadj->fvertices,&idx))) 
      MV_Add_Face(v,f);
  }
  
  void MF_Destroy_For_MESH_Delete_R4(MFace_ptr f) {
    MFace_UpAdj_R3R4 *upadj;
    MFace_DownAdj_RN *downadj;

    upadj  = (MFace_UpAdj_R3R4 *) f->upadj;
    if (upadj)
      MSTK_free(upadj);

    downadj = (MFace_DownAdj_RN *) f->downadj;
    if (downadj) {
      if (downadj->fvertices)
	List_Delete(downadj->fvertices);
      MSTK_free(downadj);
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

  void MF_Add_Region_R4(MFace_ptr f, MRegion_ptr r, int side) {
    MSTK_Report("MF_Add_Region_R2", 
		"Not yet implemented for this representation",WARN);
  }

  void MF_Rem_Region_R4(MFace_ptr f, MRegion_ptr r) {
    MSTK_Report("MF_Rem_Region_R2",
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
    MSTK_Report("MF_Num_AdjFaces_R1",
		"Not yet implemented for this representation",WARN);
    return 0;
  }

  List_ptr MF_AdjFaces_R4(MFace_ptr f) {
    MSTK_Report("MF_AdjFaces_R1",
		"Not yet implemented for this representation",WARN);
    return 0;
  }

  void MF_Add_AdjFace_R4(MFace_ptr f, int side, MFace_ptr af) {
#ifdef DEBUG
    MSTK_Report("MF_Add_AdjFace_R2",
		"Not yet implemented for this representation",WARN);
#endif
  }

  void MF_Rem_AdjFace_R4(MFace_ptr f, int side, MFace_ptr af) {
#ifdef DEBUG
    MSTK_Report("MF_Rem_AdjFace_R2",
		"Not yet implemented for this representation",WARN);
#endif
  }

#ifdef __cplusplus
}
#endif
