#define _H_MFace_Private

#include "MFace.h"
#include "MFace_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MF_Set_Vertices_R3R4(MFace_ptr f, int n, MVertex_ptr *v) {
    MFace_DownAdj_R3R4 *downadj;
    int i;

    downadj = f->downadj;
    downadj->nv =n;
    downadj->fvertices = Set_New(n);

    for (i = 0; i < n; i++)
      Set_Add(downadj->fvertices,v[i]);
  }

  void MF_Replace_Vertex_i_R3R4(MFace_ptr f, int i, MVertex_ptr v) {
    MFace_DownAdj_R3R4 *downadj;

    downadj = f->downadj;
    if (downadj->nv == 0)
      MSTK_Report("MF_Replace_Vertex_R3R4","No initial set of vertices for face",ERROR);

    Set_Replacei(downadj->fvertices,i,v);
  }

  void MF_Replace_Vertex_R3R4(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv) {
    MFace_DownAdj_R3R4 *downadj;

    downadj = f->downadj;
    if (downadj->nv == 0)
      MSTK_Report("MF_Replace_Vertex_R3R4","No initial set of vertices for face",ERROR);

    Set_Replace(downadj->fvertices,v,nuv);
  }

  int MF_Num_Vertices_R3R4(MFace_ptr f) {
    MFace_DownAdj_R3R4 *downadj;
    downadj = (MFace_DownAdj_R3R4 *) f->downadj;
    return downadj->nv;
  }

  Set_ptr MF_Vertices_R3R4(MFace_ptr f, int dir, MVertex_ptr v0) {
    MFace_DownAdj_R3R4 *downadj;
    Set_ptr fverts;
    int i, k, nv, fnd;

    downadj = (MFace_DownAdj_R3R4 *) f->downadj;
    nv = downadj->nv;

    if (!v0) {
      if (dir) 
	fverts = Set_Copy(downadj->fvertices);
      else {
	fverts = Set_New(nv);

	for (i = 0; i < nv; i++)
	  Set_Add(fverts,Set_Entry(downadj->fvertices,i));
      }
    }
    else {
      fverts = Set_New(nv);

      for (i = 0; i < nv; i++) {
	if (Set_Entry(downadj->fvertices,i) == v0) {
	  fnd = 1;
	  k = i;
	  break;
	}
      }

      if (!fnd) {
	MSTK_Report("MF_Num_Vertices_R3R4","Cannot find starting vertex",ERROR);
	return 0;
      }

      for (i = 0; i < nv; i++) {
	if (dir)
	  Set_Add(fverts,Set_Entry(downadj->fvertices,(k+i)%nv));
	else
	  Set_Add(fverts,Set_Entry(downadj->fvertices,(k+nv-i)%nv));
      }
    }

    return fverts;      
  }
	
  int MF_UsesVertex_R3R4(MFace_ptr f, MVertex_ptr v) {
    MFace_DownAdj_R3R4 *downadj;
    downadj = (MFace_DownAdj_R3R4 *) f->downadj;
    return Set_Contains(downadj->fvertices,v);
  }

#ifdef __cplusplus
}
#endif
