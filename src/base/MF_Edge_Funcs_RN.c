#define _H_MFace_Private

#include "MFace.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif


  int MF_Num_Edges_RN(MFace_ptr f) {
    return List_Num_Entries(((MFace_Adj_R1 *)f->adj)->fvertices);
  }

  List_ptr MF_Edges_RN(MFace_ptr f, int dir, MVertex_ptr v0) {
    MFace_Adj_R1 *adj = (MFace_Adj_R1 *) f->adj;
    int i, j, k, ne, nv, vgdim0, vgdim1, egdim, vgid0, vgid1, egid;
    List_ptr fedges;
    MVertex_ptr v, evtx[2], vtmp;
    MEdge_ptr e;

    k = 0;
    if (v0) {
      int fnd = 0, idx = 0;
      while (!fnd && (v = List_Next_Entry(adj->fvertices,&idx))) {
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
    
    nv = ne = List_Num_Entries(adj->fvertices);
    fedges = List_New(ne);
    for (i = 0; i < ne; i++) {
      j = dir ? (k+i)%ne : (k-i+ne)%ne;
      evtx[0] = List_Entry(adj->fvertices,j);
      j = dir ? (k+(i+1))%ne : (k-(i+1)+ne)%ne;
      evtx[1] = List_Entry(adj->fvertices,j);

#ifdef HASHTABLE
      if (evtx[0]>evtx[1]) {
	vtmp = evtx[0];
	evtx[0] = evtx[1];
	evtx[1] = vtmp;
      }

      e = Hash_Entry(MESH_Hash_Edges(MEnt_Mesh(f)), 2, evtx);
      if (e == NULL) {
	e = ME_New(MEnt_Mesh(f));
	MEnt_Set_Volatile(e);

	ME_Set_Vertex(e,0,evtx[0]);
	ME_Set_Vertex(e,1,evtx[1]);

	ME_Set_GInfo_Auto(e);
	Hash_Add(MESH_Hash_Edges(MEnt_Mesh(f)), e, 2, evtx);
      }
#else
      e = ME_New(MEnt_Mesh(f));
      MEnt_Set_Volatile(e);

      ME_Set_Vertex(e,0,evtx[0]);
      ME_Set_Vertex(e,1,evtx[1]);

      ME_Set_GInfo_Auto(e);
#endif

      List_Add(fedges,e);
      ME_Lock(e);
    }

    if (!MESH_AutoLock(MEnt_Mesh(f))) {
       i = 0;
       while ((e = List_Next_Entry(fedges, &i))) {
	 ME_UnLock(e);
       }
    }

    return fedges;
  }

  int MF_EdgeDir_RN(MFace_ptr f, MEdge_ptr e) {
    MFace_Adj_R1 *adj = (MFace_Adj_R1 *) f->adj;
    int i, nv;
    MVertex_ptr v0, v1, v;

    v0 = ME_Vertex(e,0);
    v1 = ME_Vertex(e,1);
    nv = List_Num_Entries(adj->fvertices);
    for (i = 0; i < nv; i++) {
      v = List_Entry(adj->fvertices,i);
      if (v == v0) {
	if (List_Entry(adj->fvertices,(i+1)%nv) == v1) {
	  return 1;
	}
      }
      else if (v == v1) {
	if (List_Entry(adj->fvertices,(i+1)%nv) == v0) {
	  return 0;
	}
      }
    }
        
    return -1;
  }

  int MF_EdgeDir_i_RN(MFace_ptr f, int i) {
#ifdef HASHTABLE
    MFace_Adj_R1 *adj = (MFace_Adj_R1 *) f->adj;
    int j, tdir, ne;
    MVertex_ptr evtx[2];

    ne = List_Num_Entries(adj->fvertices);

    j = i;
    evtx[0] = List_Entry(adj->fvertices,j);
    j = (i+1)%ne;
    evtx[1] = List_Entry(adj->fvertices,j);

    return (evtx[0]<evtx[1]) ? 1 : 0;
#else
    return 1;
#endif
  }

  int MF_UsesEdge_RN(MFace_ptr f, MEdge_ptr e) {
    MFace_Adj_R1 *adj = (MFace_Adj_R1 *) f->adj;
    int i, nv;
    MVertex_ptr v0, v1, v;

    v0 = ME_Vertex(e,0);
    v1 = ME_Vertex(e,1);
    nv = List_Num_Entries(adj->fvertices);
    for (i = 0; i < nv; i++) {
      v = List_Entry(adj->fvertices,i);
      if (v == v0) {
	if (List_Entry(adj->fvertices,(i+1)%nv) == v1) {
	  return 1;
	}
      }
      else if (v == v1) {
	if (List_Entry(adj->fvertices,(i+1)%nv) == v0) {
	  return 1;
	}
      }
    }
        
    return 0;
  }

#ifdef __cplusplus
}
#endif
