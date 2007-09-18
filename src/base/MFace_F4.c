#define _H_MFace_Private

#include "MFace.h"
#include "MFace_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif


  void MF_Set_RepType_F4(MFace_ptr f) {
    MFace_Adj_F2F4 *adj;

    adj = f->adj = (MFace_Adj_F2F4 *) MSTK_malloc(sizeof(MFace_Adj_F2F4));
    adj->edirs = 0;
    adj->fedges = NULL;
  }

  void MF_Delete_F4(MFace_ptr f, int keep) {
    MFace_Adj_F2F4 *adj;

    if (!keep) {
      adj = (MFace_Adj_F2F4 *) f->adj;
      if (adj) {
	if (adj->fedges)
	  List_Delete(adj->fedges);
	MSTK_free(adj);
      }
    }
  }

  void MF_Restore_F4(MFace_ptr f) {
    MEnt_Set_Dim(f,MFACE);
  }

  void MF_Destroy_For_MESH_Delete_F4(MFace_ptr f) {
    MFace_Adj_F2F4 *adj;

    adj = (MFace_Adj_F2F4 *) f->adj;
    if (adj) {
      if (adj->fedges) 
	List_Delete(adj->fedges);
      MSTK_free(adj);
    }
  }

  int MF_Set_GInfo_Auto_F4(MFace_ptr f) {
    return MF_Set_GInfo_Auto_FN(f);
  }

  void MF_Set_Edges_F4(MFace_ptr f, int n, MEdge_ptr *e, int *dir) {
    MF_Set_Edges_FN(f,n,e,dir);
  }

  void MF_Replace_Edges_F4(MFace_ptr f, int nold, MEdge_ptr *oldedges,  int nnu, MEdge_ptr *nuedges) {
    MF_Replace_Edges_FN(f,nold,oldedges,nnu,nuedges);
  }

  void MF_Replace_Edges_i_F4(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges) {
    MF_Replace_Edges_i_FN(f,nold,i,nnu,nuedges);
  }

  void MF_Set_Vertices_F4(MFace_ptr f, int n, MVertex_ptr *v) {
    MF_Set_Vertices_FN(f,n,v);
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

  void MF_Insert_Vertex_F4(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v) {
#ifdef DEBUG
    MSTK_Report("MF_Insert_Vertex","Function call not suitable for this representation",WARN);
#endif
  }

  void MF_Insert_Vertex_i_F4(MFace_ptr f, MVertex_ptr nuv, int i) {
#ifdef DEBUG
    MSTK_Report("MF_Insert_Vertex_i","Function call not suitable for this representation",WARN);
#endif
  }

  int MFs_AreSame_F4(MFace_ptr f1, MFace_ptr f2) {
    return (f1 == f2);
  }

  int MF_Num_Vertices_F4(MFace_ptr f) {
    List_ptr fedges = ((MFace_Adj_F2F4 *)f->adj)->fedges;
    return List_Num_Entries(fedges);
  }

  int MF_Num_Edges_F4(MFace_ptr f) {
    List_ptr fedges = ((MFace_Adj_F2F4 *)f->adj)->fedges;
    return List_Num_Entries(fedges);
  }

  int MF_Num_AdjFaces_F4(MFace_ptr f) {

#ifdef DEBUG
    MSTK_Report("MF_Num_AdjFaces",
		"Function call unsuitable for this representation",
		WARN);
#endif
    
    return 0;
  }

  List_ptr MF_Vertices_F4(MFace_ptr f, int dir, MVertex_ptr v0) {
    return MF_Vertices_FN(f,dir,v0);
  }
	

  List_ptr MF_Edges_F4(MFace_ptr f, int dir, MVertex_ptr v0) {
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

  List_ptr MF_AdjFaces_F4(MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("MF_AdjFaces",
		"Function call not suitable for this representation",WARN);
#endif
    return NULL;
  }
  
  List_ptr MF_Regions_F4(MFace_ptr f) {
    List_ptr fregs, eregs;
    MRegion_ptr r;
    MEdge_ptr e;
    MFace_Adj_F2F4 *adj;
    int i, k=0, nr;
    
    adj = (MFace_Adj_F2F4 *) f->adj;
    
    fregs = List_New(2);
    e = List_Entry(adj->fedges,0);
    
    eregs = ME_Regions(e);
    if (eregs) {
      nr = List_Num_Entries(eregs);
      
      for (i = 0; i < nr; i++) {
	r = List_Entry(eregs,i);
	if (MR_UsesEntity(r,f,2)) {
	  List_Add(fregs,r);
	  k++;
	}
	if (k == 2)
	  break;
      }
      
      List_Delete(eregs);
    }
			
    
    if (k) 
      return fregs;
    else {
      List_Delete(fregs);
      return 0;
    }
  }
  
  MRegion_ptr MF_Region_F4(MFace_ptr f, int dir) {
    List_ptr eregs;
    MRegion_ptr r, r1;
    MEdge_ptr e;
    int nr, i, fdir;
    MFace_Adj_F2F4 *adj;
    
#ifdef DEBUG
    MSTK_Report("MF_Region_F4","More efficient to use MF_Regions",MESG);
#endif
    
    adj = (MFace_Adj_F2F4 *) f->adj;
    e = List_Entry(adj->fedges,0);
    
    eregs = ME_Regions(e);
    if (eregs) {
      nr = List_Num_Entries(eregs);
      
      r1 = 0;
      for (i = 0; i < nr; i++) {
	r = List_Entry(eregs,i);
	fdir = MR_FaceDir(r,f);
	if (fdir == !dir) {
	  r1 = r;
	  break;
	}
      }
      List_Delete(eregs);
	
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

  MFace_ptr MF_NextInHash_F4(MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("MF_NextInHash", "Function call not suitable for this representation", WARN);
#endif
    return NULL;
  }

  void MF_Set_NextInHash_F4(MFace_ptr f, MFace_ptr next) {
#ifdef DEBUG
    MSTK_Report("MF_NextInHash", "Function call not suitable for this representation", WARN);
#endif
  }

  void MF_HashKey_F4(MFace_ptr f, unsigned int *pn, void* **pp) {
#ifdef DEBUG
    MSTK_Report("MF_NextInHash", "Function call not suitable for this representation", WARN);
#endif
  }

  int MF_IsLocked_F4(MFace_ptr f) {
    return 0;
  }

#ifdef __cplusplus
}
#endif
