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
    downadj->fvertices = List_New(n);

    for (i = 0; i < n; i++)
      List_Add(downadj->fvertices,v[i]);
  }

  void MF_Replace_Vertex_i_R3R4(MFace_ptr f, int i, MVertex_ptr v) {
    MFace_DownAdj_R3R4 *downadj;

    downadj = f->downadj;
    if (downadj->nv == 0)
      MSTK_Report("MF_Replace_Vertex_R3R4","No initial set of vertices for face",ERROR);

    List_Replacei(downadj->fvertices,i,v);
  }

  void MF_Replace_Vertex_R3R4(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv) {
    MFace_DownAdj_R3R4 *downadj;

    downadj = f->downadj;
    if (downadj->nv == 0)
      MSTK_Report("MF_Replace_Vertex_R3R4","No initial set of vertices for face",ERROR);

    List_Replace(downadj->fvertices,v,nuv);
  }

  int MF_Num_Vertices_R3R4(MFace_ptr f) {
    MFace_DownAdj_R3R4 *downadj;
    downadj = (MFace_DownAdj_R3R4 *) f->downadj;
    return downadj->nv;
  }

  List_ptr MF_Vertices_R3R4(MFace_ptr f, int dir) {
    MFace_DownAdj_R3R4 *downadj;
    downadj = (MFace_DownAdj_R3R4 *) f->downadj;
    return List_Copy(downadj->fvertices);
  }
	
  int MF_UsesVertex_R3R4(MFace_ptr f, MVertex_ptr v) {
    MFace_DownAdj_R3R4 *downadj;
    downadj = (MFace_DownAdj_R3R4 *) f->downadj;
    return List_Contains(downadj->fvertices,v);
  }

#ifdef __cplusplus
}
#endif
