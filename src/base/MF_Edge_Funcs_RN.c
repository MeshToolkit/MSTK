#define _H_MFace_Private

#include "MFace.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif


  int MF_Num_Edges_RN(MFace_ptr f) {
    return List_Num_Entries(((MFace_DownAdj_RN *)f->downadj)->fvertices);
  }

  List_ptr MF_Edges_RN(MFace_ptr f, int dir, MVertex_ptr v0) {
    MFace_DownAdj_RN *downadj = (MFace_DownAdj_RN *) f->downadj;
    int i, j, k, ne, nv, vgdim0, vgdim1, egdim, vgid0, vgid1, egid;
    List_ptr fedges;
    MVertex_ptr v, evtx0, evtx1;
    MEdge_ptr e;

    k = 0;
    if (v0) {
      int fnd = 0, idx = 0;
      while (!fnd && (v = List_Next_Entry(downadj->fvertices,&idx))) {
	if (v == v0)
	  fnd = 1;
	else
	  k++;
      }

      if (!fnd) {
	MSTK_Report("MF_Edges_RN","Cannot find vertex in face",ERROR);
	return NULL;
      }
    }
    
    nv = ne = List_Num_Entries(downadj->fvertices);
    fedges = List_New(ne);
    for (i = 0; i < ne; i++) {
      e = ME_New(MEnt_Mesh(f));
      MEnt_Set_Volatile(e);
      
      j = dir ? (k+i)%ne : (k-i+ne)%ne;
      evtx0 = List_Entry(downadj->fvertices,j);

      ME_Set_Vertex(e,0,evtx0);

      j = dir ? (k+(i+1))%ne : (k-(i+1)+ne)%ne;
      evtx1 = List_Entry(downadj->fvertices,j);

      ME_Set_Vertex(e,1,evtx1);

      ME_Set_GInfo_Auto(e);

      List_Add(fedges,e);
    }

    return fedges;
  }

  int MF_EdgeDir_RN(MFace_ptr f, MEdge_ptr e) {
    MFace_DownAdj_RN *downadj = (MFace_DownAdj_RN *) f->downadj;
    int i, nv;
    MVertex_ptr v0, v1, v;

    v0 = ME_Vertex(e,0);
    v1 = ME_Vertex(e,1);
    nv = List_Num_Entries(downadj->fvertices);
    for (i = 0; i < nv; i++) {
      v = List_Entry(downadj->fvertices,i);
      if (v == v0) {
	if (List_Entry(downadj->fvertices,(i+1)%nv) == v1) {
	  return 1;
	}
      }
      else if (v == v1) {
	if (List_Entry(downadj->fvertices,(i+1)%nv) == v0) {
	  return 0;
	}
      }
    }
        
    return -1;
  }

  int MF_EdgeDir_i_RN(MFace_ptr f, int i) {
    return 1;
  }

  int MF_UsesEdge_RN(MFace_ptr f, MEdge_ptr e) {
    MFace_DownAdj_RN *downadj = (MFace_DownAdj_RN *) f->downadj;
    int i, nv;
    MVertex_ptr v0, v1, v;

    v0 = ME_Vertex(e,0);
    v1 = ME_Vertex(e,1);
    nv = List_Num_Entries(downadj->fvertices);
    for (i = 0; i < nv; i++) {
      v = List_Entry(downadj->fvertices,i);
      if (v == v0) {
	if (List_Entry(downadj->fvertices,(i+1)%nv) == v1) {
	  return 1;
	}
      }
      else if (v == v1) {
	if (List_Entry(downadj->fvertices,(i+1)%nv) == v0) {
	  return 1;
	}
      }
    }
        
    return 0;
  }

#ifdef __cplusplus
}
#endif
