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
    MFace_DownAdj_R3R4 *downadj;

    upadj = f->upadj =(MFace_UpAdj_R3R4 *) MSTK_malloc(sizeof(MFace_UpAdj_R3R4));
    upadj->fregions[0] = (MRegion_ptr) NULL;
    upadj->fregions[1] = (MRegion_ptr) NULL;
    downadj = (MFace_DownAdj_R3R4 *) MSTK_malloc(sizeof(MFace_DownAdj_R3R4));
    downadj->nv = 0;
    downadj->fvertices = NULL;
  }

  void MF_Delete_R4(MFace_ptr f, int keep) {
    MFace_UpAdj_R3R4 *upadj;
    MFace_DownAdj_R3R4 *downadj;
    int i, nv;
    MVertex_ptr v;

    downadj = (MFace_DownAdj_R3R4 *) f->downadj;

    if (f->dim != MDELFACE) { /* if face has not been temporarily deleted */
      nv = downadj->nv;
      for (i = 0; i < nv; i++) {
	v = List_Entry(downadj->fvertices,i);
	MV_Rem_Face(v,f);
      }
    }

    if (keep) {
      MSTK_KEEP_DELETED = 1;
      f->dim = MDELFACE;
    }
    else {
#ifdef DEBUG
      f->dim = MDELFACE;
#endif

      upadj  = (MFace_UpAdj_R3R4 *) f->upadj;
      MSTK_free(upadj);

      List_Delete(downadj->fvertices);
      MSTK_free(downadj);

      MSTK_free(f);
    }
  }

  void MF_Restore_R4(MFace_ptr f) {
    MFace_DownAdj_R3R4 *downadj;
    int i, nv;
    MVertex_ptr v;

    if (f->dim != MDELFACE)
      return;

    f->dim = MFACE;

    downadj = (MFace_DownAdj_R3R4 *) f->downadj;
    nv = downadj->nv;
    for (i = 0; i < nv; i++) {
      v = List_Entry(downadj->fvertices,i);
      MV_Add_Face(v,f);
    }
  }
  

  void MF_Set_Edges_R4(MFace_ptr f, int n, MEdge_ptr *e, int *dir) {
#ifdef DEBUG
    MSTK_Report("MF_Set_Edges",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Replace_Edge_i_R4(MFace_ptr f, int i, int nnu, MEdge_ptr *nuedges, int *nudirs) {
#ifdef DEBUG
    MSTK_Report("MF_Replace_Edge",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Replace_Edge_R4(MFace_ptr f, MEdge_ptr e, int nnu, MEdge_ptr *nuedges, int *nudirs) {
#ifdef DEBUG
    MSTK_Report("MF_Replace_Edge",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Set_Vertices_R4(MFace_ptr f, int n, MVertex_ptr *v) {
    MF_Set_Vertices_R3R4(f,n,v);
  }

  void MF_Replace_Vertex_R4(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv) {
    MF_Replace_Vertex_R3R4(f,v,nuv);
  }

  void MF_Replace_Vertex_i_R4(MFace_ptr f, int i, MVertex_ptr nuv) {
    MF_Replace_Vertex_i_R3R4(f,i,nuv);
  }

  void MF_Insert_Vertex_R4(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v) {
    MF_Insert_Vertex_R3R4(f,nuv,b4v);
  }

  void MF_Insert_Vertex_i_R4(MFace_ptr f, MVertex_ptr nuv, int i) {
    MF_Insert_Vertex_i_R3R4(f,nuv,i);
  }

  void MF_Add_Region_R4(MFace_ptr f, MRegion_ptr r, int side) {
    MSTK_Report("MF_Add_Region_R2", "Not Implemented for this representation",WARN);
  }

  void MF_Rem_Region_R4(MFace_ptr f, MRegion_ptr r) {
    MSTK_Report("MF_Rem_Region_R2","Not Implemented for this representation",WARN);
  }

  int MF_Num_Vertices_R4(MFace_ptr f) {
    return MF_Num_Vertices_R3R4(f);
  }

  int MF_Num_Edges_R4(MFace_ptr f) {
    return MF_Num_Edges_R3R4(f);
  }

  List_ptr MF_Vertices_R4(MFace_ptr f, int dir, MVertex_ptr v0) {
    return MF_Vertices_R3R4(f,dir,v0);
  }
	

  List_ptr MF_Edges_R4(MFace_ptr f, int dir, MVertex_ptr v0) {
    return MF_Edges_R3R4(f,dir,v0); 
  }

  int MF_EdgeDir_R4(MFace_ptr f, MEdge_ptr e) {
    return MF_EdgeDir_R3R4(f,e);
  }

  int MF_EdgeDir_i_R4(MFace_ptr f, int i) {
    return MF_EdgeDir_i_R3R4(f,i);
  }

  int MF_UsesEdge_R4(MFace_ptr f, MEdge_ptr e) {
    return MF_UsesEdge_R3R4(f,e);
  }

  int MF_UsesVertex_R4(MFace_ptr f, MVertex_ptr v) {
    return MF_UsesVertex_R3R4(f,v);
  }

  List_ptr MF_Regions_R4(MFace_ptr f) {
    return MF_Regions_R3R4(f);
  }

  MRegion_ptr MF_Region_R4(MFace_ptr f, int dir) {
    return MF_Region_R3R4(f,dir);
  }


  int MF_Num_AdjFaces_R4(MFace_ptr f) {
    MSTK_Report("MF_Num_AdjFaces_R1","Not yet implemented for this representation",WARN);
    return 0;
  }

  List_ptr MF_AdjFaces_R4(MFace_ptr f) {
    MSTK_Report("MF_AdjFaces_R1","Not yet implemented for this representation",WARN);
    return 0;
  }

  void MF_Add_AdjFace_R4(MFace_ptr f, int side, MFace_ptr af) {
#ifdef DEBUG
    MSTK_Report("MF_Add_AdjFace_R2","Not yet implemented for this representation",WARN);
#endif
  }

  void MF_Rem_AdjFace_R4(MFace_ptr f, int side, MFace_ptr af) {
#ifdef DEBUG
    MSTK_Report("MF_Rem_AdjFace_R2","Not yet implemented for this representation",WARN);
#endif
  }

#ifdef __cplusplus
}
#endif
