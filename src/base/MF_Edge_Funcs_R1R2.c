#define _H_MFace_Private

#include "MFace.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  int MF_Num_Edges_R1R2(MFace_ptr f) {
    MSTK_Report("MF_Num_Edges_R1R2","Not implemented",ERROR);
    return 0;
  }

  List_ptr MF_Edges_R1R2(MFace_ptr f, int dir, MVertex_ptr v0) {
    MSTK_Report("MF_Edges_R1R2","Not implemented",ERROR);
    return 0;    
  }

  int MF_EdgeDir_R1R2(MFace_ptr f, MEdge_ptr e) {
    MSTK_Report("MF_EdgeDir_R1R2","Not implemented",ERROR);
    return 0;
  }

  int MF_EdgeDir_i_R1R2(MFace_ptr f, int i) {
    MSTK_Report("MF_EdgeDir_i_R1R2","Not implemented",ERROR);
    return 0;
  }

  int MF_UsesEdge_R1R2(MFace_ptr f, MEdge_ptr e) {
    MSTK_Report("MF_UsesEdge_R1R2","Not implemented",ERROR);
    return 0;
  }

#ifdef __cplusplus
}
#endif
