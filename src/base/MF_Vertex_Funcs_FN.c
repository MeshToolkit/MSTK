#define _H_MFace_Private

#include "MFace.h"
#include "MFace_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif


  Set_ptr MF_Vertices_FN(MFace_ptr f, int dir, MVertex_ptr v0) {
    int i, k, ne, edir, fnd=0;
    Set_ptr fverts;
    MEdge_ptr e;
    MVertex_ptr v;
    MFace_DownAdj_FN *downadj;

    downadj = (MFace_DownAdj_FN *) f->downadj;
    ne = downadj->ne;
    fverts = Set_New(ne);

    if (!v0) {
      for (i = 0; i < ne; i++) {
	k = dir ? i : ne-1-i;
	e = Set_Entry(downadj->fedges,k);
	edir = ((downadj->edirs)>>k) & 1;
	v = ME_Vertex(e,edir^dir);
	Set_Add(fverts,v);
      }
    }
    else {
      fnd = 0;
      for (i = 0; i < ne; i++) {
        e = Set_Entry(downadj->fedges,i);
        edir = ((downadj->edirs)>>i) & 1;
        if (ME_Vertex(e,edir^dir) == v0) {
          fnd = 1;
          k = i;
          break;
        }
      }

      if (!fnd)
        MSTK_Report("MF_Edges_F1","Cannot find vertex in face!!",FATAL);

      for (i = 0; i < ne; i++) {
	e = dir ? Set_Entry(downadj->fedges,(k+i)%ne) :
	  Set_Entry(downadj->fedges,(k+ne-i)%ne);
	edir = dir ? ((downadj->edirs)>>(k+i)%ne) & 1 :
	  ((downadj->edirs)>>(k+ne-i)%ne) & 1;
	v = ME_Vertex(e,edir^dir);
	Set_Add(fverts,v);
      }
    }

    return fverts;
  }
	

  int MF_UsesVertex_FN(MFace_ptr f, MVertex_ptr v) {
    int ne, i;
    MEdge_ptr e;
    MFace_DownAdj_FN *downadj;

    /* We have to check only ne-1 edges since they form a loop */
    downadj = (MFace_DownAdj_FN *) f->downadj;
    ne = downadj->ne;
    for (i = 0; i < ne-1; i++) {
      e = Set_Entry(downadj->fedges,i);
      if (ME_UsesEntity(e,v,0))
	return 1;
    }

    return 0;
  }


#ifdef __cplusplus
}
#endif
