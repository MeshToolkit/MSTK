#define _H_MFace_Private

#include "MFace.h"
#include "MFace_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif


  void MF_Set_RepType_F4(MFace_ptr f) {
    MFace_DownAdj_FN *downadj;

    downadj = f->downadj = (MFace_DownAdj_FN *) MSTK_malloc(sizeof(MFace_DownAdj_FN));
    downadj->ne = 0;
    downadj->edirs = 0;
    downadj->fedges = NULL;
  }

  void MF_Delete_F4(MFace_ptr f) {
    MFace_DownAdj_FN *downadj;

    downadj = (MFace_DownAdj_FN *) f->downadj;
    Set_Delete(downadj->fedges);
    MSTK_free(downadj);
  }

  void MF_Set_Edges_F4(MFace_ptr f, int n, MEdge_ptr *e, int *dir) {
    MF_Set_Edges_FN(f,n,e,dir);
  }

  void MF_Replace_Edge_F4(MFace_ptr f, MEdge_ptr e, MEdge_ptr nue, int dir) {
    MF_Replace_Edge_FN(f,e,nue,dir);
  }

  void MF_Replace_Edge_i_F4(MFace_ptr f, int i, MEdge_ptr nue, int dir) {
    MF_Replace_Edge_i_FN(f,i,nue,dir);
  }

  void MF_Set_Vertices_F4(MFace_ptr f, int n, MVertex_ptr *v) {
#ifdef DEBUG
    MSTK_Report("MF_Set_Vertices","Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Replace_Vertex_i_F4(MFace_ptr f, int i, MVertex_ptr v) {
#ifdef DEBUG
    MSTK_Report("MF_Replace_Vertex","Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Replace_Vertex_F4(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv) {
#ifdef DEBUG
    MSTK_Report("MF_Replace_Vertex","Function call not suitable for this representation",WARN);
#endif
  }

  int MF_Num_Vertices_F4(MFace_ptr f) {
    return ((MFace_DownAdj_FN *)f->downadj)->ne;
  }

  int MF_Num_Edges_F4(MFace_ptr f) {
    return ((MFace_DownAdj_FN *)f->downadj)->ne;
  }

  int MF_Num_AdjFaces_F4(MFace_ptr f) {

#ifdef DEBUG
    MSTK_Report("MF_Num_AdjFaces",
		"Function call unsuitable for this representation",
		WARN);
#endif
    
    return 0;
  }

  Set_ptr MF_Vertices_F4(MFace_ptr f, int dir, MVertex_ptr v0) {
    return MF_Vertices_FN(f,dir,v0);
  }
	

  Set_ptr MF_Edges_F4(MFace_ptr f, int dir, MVertex_ptr v0) {
    return MF_Edges_FN(f,dir,v0);
  }

  int MF_EdgeDir_F4(MFace_ptr f, MEdge_ptr e) {
    return MF_EdgeDir_FN(f,e);
  }

  int MF_EdgeDir_i_F4(MFace_ptr f, int i) {
    return MF_EdgeDir_i_FN(f,i);
  }

  int MF_UsesEdge_F4(MFace_ptr f, MEdge_ptr e) {
    return MF_UsesEdge_FN(f,e);
  }

  int MF_UsesVertex_F4(MFace_ptr f, MVertex_ptr v) {
    return MF_UsesVertex_FN(f,v);
  }

  Set_ptr MF_AdjFaces_F4(MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("MF_AdjFaces",
		"Function call not suitable for this representation",WARN);
#endif
    return NULL;
  }
  
  Set_ptr MF_Regions_F4(MFace_ptr f) {
    Set_ptr fregs, eregs;
    MRegion_ptr r;
    MEdge_ptr e;
    MFace_DownAdj_FN *downadj;
    int i, k=0, nr;
    
    downadj = (MFace_DownAdj_FN *) f->downadj;
    
    fregs = Set_New(2);
    e = Set_Entry(downadj->fedges,0);
    
    eregs = ME_Regions(e);
    if (eregs) {
      nr = Set_Num_Entries(eregs);
      
      for (i = 0; i < nr; i++) {
	r = Set_Entry(eregs,i);
	if (MR_UsesEntity(r,f,2)) {
	  Set_Add(fregs,r);
	  k++;
	}
	if (k == 2)
	  break;
      }
      
      Set_Delete(eregs);
    }
			
    
    if (k) 
      return fregs;
    else {
      Set_Delete(fregs);
      return 0;
    }
  }
  
  MRegion_ptr MF_Region_F4(MFace_ptr f, int dir) {
    Set_ptr eregs;
    MRegion_ptr r, r1;
    MEdge_ptr e;
    int nr, i, fdir;
    MFace_DownAdj_FN *downadj;
    
#ifdef DEBUG
    MSTK_Report("MF_Region_F4","More efficient to use MF_Regions",MESG);
#endif
    
    downadj = (MFace_DownAdj_FN *) f->downadj;
    e = Set_Entry(downadj->fedges,0);
    
    eregs = ME_Regions(e);
    if (eregs) {
      nr = Set_Num_Entries(eregs);
      
      r1 = 0;
      for (i = 0; i < nr; i++) {
	r = Set_Entry(eregs,i);
	fdir = MR_FaceDir(r,f);
	if (fdir == !dir) {
	  r1 = r;
	  break;
	}
      }
      Set_Delete(eregs);
	
      return r1;
    }
    
    return 0;
  }

  void MF_Add_AdjFace_F4(MFace_ptr f, int edgnum, MFace_ptr af) {
#ifdef DEBUG
    MSTK_Report("MF_Add_AdjFace",
		"Function call unsuitable for this representation",
		WARN);
#endif
  }

  void MF_Rem_AdjFace_F4(MFace_ptr f, int edgnum, MFace_ptr af) {
#ifdef DEBUG
    MSTK_Report("MF_Rem_AdjFace",
		"Function call unsuitable for this representation",
		WARN);
#endif
  }

  void MF_Add_Region_F4(MFace_ptr f, MRegion_ptr r, int side) {
#ifdef DEBUG
    MSTK_Report("MF_Add_Region",
		"Function call unsuitable for this representation",
		WARN);
#endif
  }

  void MF_Rem_Region_F4(MFace_ptr f, MRegion_ptr r) {
#ifdef DEBUG
    MSTK_Report("MF_Rem_Region",
		"Function call unsuitable for this representation",
		WARN);
#endif
  }

#ifdef __cplusplus
}
#endif
