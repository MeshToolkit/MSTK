#define _H_MEdge_Private

#include "MEdge.h"
#include "MEdge_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void ME_Set_RepType_R2(MEdge_ptr e) {
    MSTK_Report("ME_Set_RepType_R2","Not implemented",ERROR);
  }

  void ME_Delete_R2(MEdge_ptr e, int keep) {
    MSTK_Report("ME_Set_RepType_R2","Not implemented",ERROR);
  }

  int ME_Num_Faces_R2(MEdge_ptr e) {
    return ME_Num_Faces_R1R2(e);
  }

  int ME_Num_Regions_R2(MEdge_ptr e) {
    return ME_Num_Regions_R1R2(e);
  }

  List_ptr ME_Faces_R2(MEdge_ptr e) {
    return ME_Faces_R1R2(e);
  }

  List_ptr ME_Regions_R2(MEdge_ptr e) {
    return ME_Regions_R1R2(e);
  }

  void ME_Add_Face_R2(MEdge_ptr e, MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("ME_Add_Face","Function call not suitable for this representation",WARN);
#endif
  }
  
  void ME_Rem_Face_R2(MEdge_ptr e, MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("ME_Rem_Face","Function call not suitable for this representation",WARN);
#endif
  }

  void ME_Add_Region_R2(MEdge_ptr e, MRegion_ptr r) {
#ifdef DEBUG
    MSTK_Report("ME_Add_Region","Function call not suitable for this representation",WARN);
#endif
  }

  void ME_Rem_Region_R2(MEdge_ptr e, MRegion_ptr r) {
#ifdef DEBUG
    MSTK_Report("ME_Rem_Region","Function call not suitable for this representation",WARN);
#endif
  }

#ifdef __cplusplus
}
#endif

