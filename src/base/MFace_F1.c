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
    downadj->edirs = 0;
    downadj->fedges = NULL;
  }

  void MF_Delete_F1(MFace_ptr f, int keep) {
    MFace_UpAdj_F1F3 *upadj;
    MFace_DownAdj_FN *downadj;
    MEdge_ptr e;
    int i, ne;

    downadj = (MFace_DownAdj_FN *) f->downadj;

    if (MEnt_Dim(f) != MDELETED) { /* if face hasnt been temporarily deleted */
      if (downadj) {
	if (downadj->fedges) {
	  ne = List_Num_Entries(downadj->fedges);
	  for (i = 0; i < ne; i++) {
	    e = List_Entry(downadj->fedges,i);
	    ME_Rem_Face(e,f);
	  }
	}
      }
    }

    if (!keep) {
      upadj = (MFace_UpAdj_F1F3 *) f->upadj;
      if (upadj)
	MSTK_free(upadj);

      if (downadj) {
	List_Delete(downadj->fedges);
	MSTK_free(downadj);
      }
    }
  }

  void MF_Restore_F1(MFace_ptr f) {
    MFace_DownAdj_FN *downadj;
    MEdge_ptr e;
    int i, ne;

    MEnt_Set_Dim(f,MFACE);

    downadj = (MFace_DownAdj_FN *) f->downadj;

    ne = List_Num_Entries(downadj->fedges);
    for (i = 0; i < ne; i++) {
      e = List_Entry(downadj->fedges,i);
      ME_Add_Face(e,f);
    }
  }

  void MF_Destroy_For_MESH_Delete_F1(MFace_ptr f) {
    MFace_UpAdj_F1F3 *upadj;
    MFace_DownAdj_FN *downadj;

    upadj = (MFace_UpAdj_F1F3 *) f->upadj;
    if (upadj)
      MSTK_free(upadj);
    
    downadj = (MFace_DownAdj_FN *) f->downadj;
    if (downadj) {
      if (downadj->fedges)
	List_Delete(downadj->fedges);
      MSTK_free(downadj);
    }
  }

  int MF_Set_GInfo_Auto_F1(MFace_ptr f) {
    return MF_Set_GInfo_Auto_FN(f);
  }

  void MF_Set_Edges_F1(MFace_ptr f, int n, MEdge_ptr *e, int *dir) {
    MF_Set_Edges_FN(f,n,e,dir);
  }

  void MF_Replace_Edges_F1(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges) {
    MF_Replace_Edges_FN(f,nold,oldedges,nnu,nuedges);
  }

  void MF_Replace_Edges_i_F1(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges) {
    MF_Replace_Edges_i_FN(f,nold,i,nnu,nuedges);
  }

  void MF_Set_Vertices_F1(MFace_ptr f, int n, MVertex_ptr *v) {
    MF_Set_Vertices_FN(f,n,v);
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

  int MFs_AreSame_F1(MFace_ptr f1, MFace_ptr f2) {
    return (f1 == f2);
  }

  int MF_Num_Vertices_F1(MFace_ptr f) {
    List_ptr fedges = ((MFace_DownAdj_FN *)f->downadj)->fedges;
    return List_Num_Entries(fedges);
  }

  int MF_Num_Edges_F1(MFace_ptr f) {
    List_ptr fedges = ((MFace_DownAdj_FN *)f->downadj)->fedges;
    return List_Num_Entries(fedges);
  }

  int MF_Num_AdjFaces_F1(MFace_ptr f) {
    List_ptr adjfaces;
    int nf=0;

#ifdef DEBUG
    MSTK_Report("MF_Num_AdjFaces","Inefficient call for this representation - Call MF_AdjFaces directly",WARN);
#endif

    adjfaces = MF_AdjFaces_F1(f);
    nf = List_Num_Entries(adjfaces);
    List_Delete(adjfaces);
    return nf;
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
    List_ptr adjfaces, efaces;
    MFace_DownAdj_FN *downadj;
    MEdge_ptr e;
    MFace_ptr ef;
    int i, j, ne, nf;

    downadj = (MFace_DownAdj_FN *) f->downadj;

    ne = List_Num_Entries(downadj->fedges);
    adjfaces = List_New(ne);

    for (i = 0; i < ne; i++) {
      e = List_Entry(downadj->fedges,i);

      efaces = ME_Faces(e);
      nf = List_Num_Entries(efaces);

      if (nf > 2)
	MSTK_Report("MF_AdjFaces_F1","Non-manifold or 3D mesh!!",ERROR);

      for (j = 0; j < nf; j++) {
	ef = List_Entry(efaces,j);
	if (ef != f)
	  List_Add(adjfaces,ef);
      }
      List_Delete(efaces);
    }

    if (List_Num_Entries(adjfaces))
      return adjfaces;
    else {
      List_Delete(adjfaces);
      return NULL;
    }
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
