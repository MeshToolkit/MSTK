#define _H_MFace_Private

#include "MFace.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif


  int MF_Num_Edges_R3R4(MFace_ptr f) {
    return ((MFace_DownAdj_R3R4 *)f->downadj)->nv;
  }

  Set_ptr MF_Edges_R3R4(MFace_ptr f, int dir, MVertex_ptr v0) {
    /* Are we supposed to create edges and send them back ? */
    MSTK_Report("MF_Edges_R3R4","Not implemented",ERROR);
    return NULL;
  }

  int MF_EdgeDir_R3R4(MFace_ptr f, MEdge_ptr e) {
    MSTK_Report("MF_EdgeDir_R3R4","Not Implemented",ERROR);
    return 0;
  }

  int MF_EdgeDir_i_R3R4(MFace_ptr f, int i) {
    MSTK_Report("MF_EdgeDir_i_R3R4","Not Implemented",ERROR);
    return 0;
  }

  int MF_UsesEdge_R3R4(MFace_ptr f, MEdge_ptr e) {
    MSTK_Report("MF_UsesEdge_R3R4","Not Implemented",ERROR);
    return 0;
  }

#ifdef __cplusplus
}
#endif
