#define _H_MFace_Private

#include "MFace.h"
#include "MFace_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Most of these functions are identical to those of F4 */

  void MF_Set_RepType_F1(MFace_ptr f) {
    MFace_UpAdj_F1F3 *upadj;
    MFace_DownAdj_FN *downadj;

    upadj = f->upadj = (MFace_UpAdj_F1F3 *) MSTK_malloc(sizeof(MFace_UpAdj_F1F3));
    upadj->fregions[0] = (MRegion_ptr) NULL;
    upadj->fregions[1] = (MRegion_ptr) NULL;
    downadj = f->downadj = (MFace_DownAdj_FN *) MSTK_malloc(sizeof(MFace_DownAdj_FN));
    downadj->ne = 0;
    downadj->edirs = 0;
    downadj->fedges = NULL;
  }

  void MF_Delete_F1(MFace_ptr f, int keep) {
    MFace_UpAdj_F1F3 *upadj;
    MFace_DownAdj_FN *downadj;
    MEdge_ptr e;
    int i, ne;

    downadj = (MFace_DownAdj_FN *) f->downadj;

    if (f->dim != MDELFACE) { /* if face has not been temporarily deleted */
      ne = List_Num_Entries(downadj->fedges);
      for (i = 0; i < ne; i++) {
	e = List_Entry(downadj->fedges,i);
	ME_Rem_Face(e,f);
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

      upadj = (MFace_UpAdj_F1F3 *) f->upadj;
      MSTK_free(upadj);

      List_Delete(downadj->fedges);
      MSTK_free(downadj);
      
      MSTK_free(f);
    }
  }

  void MF_Restore_F1(MFace_ptr f) {
    MFace_DownAdj_FN *downadj;
    MEdge_ptr e;
    int i, ne;

    if (f->dim != MDELFACE)
      return;

    f->dim = MFACE;

    downadj = (MFace_DownAdj_FN *) f->downadj;

    ne = List_Num_Entries(downadj->fedges);
    for (i = 0; i < ne; i++) {
      e = List_Entry(downadj->fedges,i);
      ME_Add_Face(e,f);
    }
  }

  void MF_Set_Edges_F1(MFace_ptr f, int n, MEdge_ptr *e, int *dir) {
    MF_Set_Edges_FN(f,n,e,dir);
  }

  void MF_Replace_Edge_F1(MFace_ptr f, MEdge_ptr e, int nnu, MEdge_ptr *nuedges, int *nudirs) {
    MF_Replace_Edge_FN(f,e,nnu,nuedges,nudirs);
  }

  void MF_Replace_Edge_i_F1(MFace_ptr f, int i, int nnu, MEdge_ptr *nuedges, int *nudir) {
    MF_Replace_Edge_i_FN(f,i,nnu,nuedges,nudir);
  }

  void MF_Set_Vertices_F1(MFace_ptr f, int n, MVertex_ptr *v) {
#ifdef DEBUG
    MSTK_Report("MF_Set_Vertices","Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Replace_Vertex_F1(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv) {
#ifdef DEBUG
    MSTK_Report("MF_Replace_Vertex","Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Replace_Vertex_i_F1(MFace_ptr f, int i, MVertex_ptr nuv) {
#ifdef DEBUG
    MSTK_Report("MF_Replace_Vertex","Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Insert_Vertex_F1(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v) {
#ifdef DEBUG
    MSTK_Report("MF_Insert_Vertex","Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Insert_Vertex_i_F1(MFace_ptr f, MVertex_ptr nuv, int i) {
#ifdef DEBUG
    MSTK_Report("MF_Insert_Vertex_i","Function call not suitable for this representation",WARN);
#endif
  }

  int MF_Num_Vertices_F1(MFace_ptr f) {
    return ((MFace_DownAdj_FN *)f->downadj)->ne;
  }

  int MF_Num_Edges_F1(MFace_ptr f) {
    return ((MFace_DownAdj_FN *)f->downadj)->ne;
  }

  int MF_Num_AdjFaces_F1(MFace_ptr f) {

#ifdef DEBUG
    MSTK_Report("MF_Num_AdjFaces",
		"Function call unsuitable for this representation",
		WARN);
#endif
    
    return 0;
  }

  List_ptr MF_Vertices_F1(MFace_ptr f, int dir, MVertex_ptr v0) {
    return MF_Vertices_FN(f,dir,v0);
  }
	

  List_ptr MF_Edges_F1(MFace_ptr f, int dir, MVertex_ptr v0) {
    return MF_Edges_FN(f,dir,v0);
  }

  int MF_EdgeDir_F1(MFace_ptr f, MEdge_ptr e) {
    return MF_EdgeDir_FN(f,e);
  }

  int MF_EdgeDir_i_F1(MFace_ptr f, int i) {
    return MF_EdgeDir_i_FN(f,i);
  }
			
  int MF_UsesEdge_F1(MFace_ptr f, MEdge_ptr e) {
    return MF_UsesEdge_FN(f,e);
  }

  int MF_UsesVertex_F1(MFace_ptr f, MVertex_ptr v) {
    return MF_UsesVertex_FN(f,v);
  }

  List_ptr MF_AdjFaces_F1(MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("MF_AdjFaces",
		"Function call not suitable for this representation",WARN);
#endif
    return NULL;
  }

  List_ptr MF_Regions_F1(MFace_ptr f) {
    return MF_Regions_F1F3(f);
  }

  MRegion_ptr MF_Region_F1(MFace_ptr f, int i) {
    return MF_Region_F1F3(f,i);
  }

  void MF_Add_AdjFace_F1(MFace_ptr f, int edgnum, MFace_ptr af) {
#ifdef DEBUG
    MSTK_Report("MF_Add_AdjFace",
		"Function call unsuitable for this representation",
		WARN);
#endif
  }

  void MF_Rem_AdjFace_F1(MFace_ptr f, int edgnum, MFace_ptr af) {
#ifdef DEBUG
    MSTK_Report("MF_Rem_AdjFace",
		"Function call unsuitable for this representation",
		WARN);
#endif
  }

  void MF_Add_Region_F1(MFace_ptr f, MRegion_ptr r, int side) {
    MF_Add_Region_F1F3(f,r,side);
  }

  void MF_Rem_Region_F1(MFace_ptr f, MRegion_ptr r) {
    MF_Rem_Region_F1F3(f,r);
  }


#ifdef __cplusplus
}
#endif
