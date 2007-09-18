#define _H_MEdge_Private

#include "MEdge.h"
#include "MEdge_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void ME_Set_RepType_R2(MEdge_ptr e) {
    MEdge_Adj_RN *adj;

    adj = e->adj = (MEdge_Adj_RN *) MSTK_malloc(sizeof(MEdge_Adj_RN));
    adj->hnext = NULL;
  }

  void ME_Delete_R2(MEdge_ptr e, int keep) {
    return;
  }

  void ME_Destroy_For_MESH_Delete_R2(MEdge_ptr e) {
    return;
  }

  int ME_Num_Faces_R2(MEdge_ptr e) {
    return ME_Num_Faces_R1R2(e);
  }

  int ME_Num_Regions_R2(MEdge_ptr e) {
    return ME_Num_Regions_RN(e);
  }

  List_ptr ME_Faces_R2(MEdge_ptr e) {
    return ME_Faces_R1R2(e);
  }

  List_ptr ME_Regions_R2(MEdge_ptr e) {
    return ME_Regions_RN(e);
  }

  void ME_Add_Face_R2(MEdge_ptr e, MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("ME_Add_Face",
		"Function call not suitable for this representation",WARN);
#endif
  }
  
  void ME_Rem_Face_R2(MEdge_ptr e, MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("ME_Rem_Face",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void ME_Add_Region_R2(MEdge_ptr e, MRegion_ptr r) {
#ifdef DEBUG
    MSTK_Report("ME_Add_Region",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void ME_Rem_Region_R2(MEdge_ptr e, MRegion_ptr r) {
#ifdef DEBUG
    MSTK_Report("ME_Rem_Region",
		"Function call not suitable for this representation",WARN);
#endif
  }

  MEdge_ptr ME_NextInHash_R2(MEdge_ptr e) {
    MEdge_ptr hnext = ((MEdge_Adj_RN *)e->adj)->hnext; 
    return hnext;
  }

  void ME_Set_NextInHash_R2(MEdge_ptr e, MEdge_ptr next) {
    MEdge_Adj_RN *adj = (MEdge_Adj_RN *)e->adj;
    adj->hnext = next;
  }

  void ME_Lock_R2(MEdge_ptr e) {
    MEdge_Adj_RN *adj = (MEdge_Adj_RN *)e->adj;
    Hash_Lock(&adj->lock);
  }

  void ME_UnLock_R2(MEdge_ptr e) {
    MEdge_Adj_RN *adj = (MEdge_Adj_RN *)e->adj;
    Hash_UnLock(&adj->lock);
  }

  int ME_IsLocked_R2(MEdge_ptr e) {
    MEdge_Adj_RN *adj = (MEdge_Adj_RN *)e->adj;
    return Hash_IsLocked(adj->lock);
  }

#ifdef __cplusplus
}
#endif

