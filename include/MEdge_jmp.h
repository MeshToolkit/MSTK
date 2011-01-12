#define _H_MEdge_Private

#include "MSTK_defines.h"
#include "MSTK_types.h"
#include "MEdge.h"

#ifdef __cplusplus
extern "C" {
#endif

  void ME_Dummy1(MEdge_ptr e);
  void ME_Dummy2(MEdge_ptr e, int i);

  void ME_Set_RepType_F1(MEdge_ptr e);
  void ME_Set_RepType_F4(MEdge_ptr e);
  void ME_Set_RepType_R1(MEdge_ptr e);
  void ME_Set_RepType_R2(MEdge_ptr e);
  void ME_Set_RepType_R4(MEdge_ptr e);
  static void (*ME_Set_RepType_jmp[MSTK_MAXREP])(MEdge_ptr e) =
  {ME_Set_RepType_F1, ME_Set_RepType_F4, ME_Set_RepType_R1, ME_Set_RepType_R2, ME_Set_RepType_R4};

  void ME_Delete_F1(MEdge_ptr e, int keep);
  void ME_Delete_F4(MEdge_ptr e, int keep);
  void ME_Delete_R1(MEdge_ptr e, int keep);
  void ME_Delete_R2(MEdge_ptr e, int keep);
  void ME_Delete_R4(MEdge_ptr e, int keep);
  static void (*ME_Delete_jmp[MSTK_MAXREP])(MEdge_ptr e, int keep) = 
  {ME_Delete_F1, ME_Delete_F4, ME_Delete_R1, ME_Delete_R2, ME_Delete_R4};

  void ME_Restore_F1(MEdge_ptr e);
  void ME_Restore_F4(MEdge_ptr e);
  static void (*ME_Restore_jmp[MSTK_MAXREP])(MEdge_ptr e) = 
  {ME_Restore_F1, ME_Restore_F4, ME_Dummy1, ME_Dummy1, ME_Dummy1};

  void ME_Destroy_For_MESH_Delete_F1(MEdge_ptr e);
  void ME_Destroy_For_MESH_Delete_F4(MEdge_ptr e);
  static void (*ME_Destroy_For_MESH_Delete_jmp[MSTK_MAXREP])(MEdge_ptr e) = 
  {ME_Destroy_For_MESH_Delete_F1, ME_Destroy_For_MESH_Delete_F4, ME_Dummy1, 
   ME_Dummy1, ME_Dummy1};

  int ME_Num_Faces_F1(MEdge_ptr e);
  int ME_Num_Faces_F4(MEdge_ptr e);
  int ME_Num_Faces_R1(MEdge_ptr e);
  int ME_Num_Faces_R2(MEdge_ptr e);
  int ME_Num_Faces_R4(MEdge_ptr e);
  int ME_Num_Faces_R1R2(MEdge_ptr e);
  int ME_Num_Faces_R3R4(MEdge_ptr e);
#ifdef DEBUG
  static int (*ME_Num_Faces_jmp[MSTK_MAXREP])(MEdge_ptr e) = 
  {ME_Num_Faces_F1, ME_Num_Faces_F4, ME_Num_Faces_R1, ME_Num_Faces_R2,
   ME_Num_Faces_R4};
#else
  static int (*ME_Num_Faces_jmp[MSTK_MAXREP])(MEdge_ptr e) = 
  {ME_Num_Faces_F1, ME_Num_Faces_F4, ME_Num_Faces_R1R2, ME_Num_Faces_R1R2,
   ME_Num_Faces_R3R4};
#endif

  int ME_Num_Regions_F1(MEdge_ptr e);
  int ME_Num_Regions_F4(MEdge_ptr e);
  int ME_Num_Regions_R1(MEdge_ptr e);
  int ME_Num_Regions_R2(MEdge_ptr e);
  int ME_Num_Regions_R4(MEdge_ptr e);
  int ME_Num_Regions_RN(MEdge_ptr e);
#ifdef DEBUG
  static int (*ME_Num_Regions_jmp[MSTK_MAXREP])(MEdge_ptr e) =
  {ME_Num_Regions_F1, ME_Num_Regions_F4, ME_Num_Regions_R1, ME_Num_Regions_R2,
   ME_Num_Regions_R4};
#else
  static int (*ME_Num_Regions_jmp[MSTK_MAXREP])(MEdge_ptr e) =
  {ME_Num_Regions_F1, ME_Num_Regions_F4, ME_Num_Regions_RN,
   ME_Num_Regions_RN, ME_Num_Regions_RN};
#endif

  List_ptr ME_Faces_F1(MEdge_ptr e);
  List_ptr ME_Faces_F4(MEdge_ptr e);
  List_ptr ME_Faces_R1(MEdge_ptr e);
  List_ptr ME_Faces_R2(MEdge_ptr e);
  List_ptr ME_Faces_R4(MEdge_ptr e);
  List_ptr ME_Faces_R1R2(MEdge_ptr e);
  List_ptr ME_Faces_R3R4(MEdge_ptr e);
#ifdef DEBUG
  static List_ptr (*ME_Faces_jmp[MSTK_MAXREP])(MEdge_ptr e) =
  {ME_Faces_F1, ME_Faces_F4, ME_Faces_R1R2, ME_Faces_R1R2, ME_Faces_R4};
#else
  static List_ptr (*ME_Faces_jmp[MSTK_MAXREP])(MEdge_ptr e) =
  {ME_Faces_F1, ME_Faces_F4, ME_Faces_R1R2, ME_Faces_R1R2, ME_Faces_R3R4};
#endif
  
  List_ptr ME_Regions_F1(MEdge_ptr e);
  List_ptr ME_Regions_F4(MEdge_ptr e);
  List_ptr ME_Regions_R1(MEdge_ptr e);
  List_ptr ME_Regions_R2(MEdge_ptr e);
  List_ptr ME_Regions_R4(MEdge_ptr e);
  List_ptr ME_Regions_RN(MEdge_ptr e);
#ifdef DEBUG
  static List_ptr (*ME_Regions_jmp[MSTK_MAXREP])(MEdge_ptr e) =
  {ME_Regions_F1, ME_Regions_F4, ME_Regions_R1, ME_Regions_R2, ME_Regions_R4};
#else
  static List_ptr (*ME_Regions_jmp[MSTK_MAXREP])(MEdge_ptr e) =
  {ME_Regions_F1, ME_Regions_F4, ME_Regions_RN, ME_Regions_RN, ME_Regions_RN};
#endif

  MEdge_ptr ME_NextInHash_F1(MEdge_ptr e);
  MEdge_ptr ME_NextInHash_F4(MEdge_ptr e);
  MEdge_ptr ME_NextInHash_R1(MEdge_ptr e);
  MEdge_ptr ME_NextInHash_R2(MEdge_ptr e);
  MEdge_ptr ME_NextInHash_R4(MEdge_ptr e);
  static MEdge_ptr (*ME_NextInHash_jmp[MSTK_MAXREP])(MEdge_ptr e) =
  {ME_NextInHash_F1, ME_NextInHash_F4, ME_NextInHash_R1, ME_NextInHash_R2, ME_NextInHash_R4};

  void ME_Set_NextInHash_F1(MEdge_ptr e, MEdge_ptr next);
  void ME_Set_NextInHash_F4(MEdge_ptr e, MEdge_ptr next);
  void ME_Set_NextInHash_R1(MEdge_ptr e, MEdge_ptr next);
  void ME_Set_NextInHash_R2(MEdge_ptr e, MEdge_ptr next);
  void ME_Set_NextInHash_R4(MEdge_ptr e, MEdge_ptr next);
  static void (*ME_Set_NextInHash_jmp[MSTK_MAXREP])(MEdge_ptr e, MEdge_ptr next) =
  {ME_Set_NextInHash_F1, ME_Set_NextInHash_F4, ME_Set_NextInHash_R1, ME_Set_NextInHash_R2, ME_Set_NextInHash_R4};


  void ME_Add_Face_F1(MEdge_ptr e, MFace_ptr f);
  void ME_Add_Face_F4(MEdge_ptr e, MFace_ptr f);
  void ME_Add_Face_R1(MEdge_ptr e, MFace_ptr f);
  void ME_Add_Face_R2(MEdge_ptr e, MFace_ptr f);
  void ME_Add_Face_R4(MEdge_ptr e, MFace_ptr f);
  static void (*ME_Add_Face_jmp[MSTK_MAXREP])(MEdge_ptr e, MFace_ptr f) =
  {ME_Add_Face_F1, ME_Add_Face_F4, ME_Add_Face_R1, ME_Add_Face_R2, 
   ME_Add_Face_R4};

  void ME_Rem_Face_F1(MEdge_ptr e, MFace_ptr f);
  void ME_Rem_Face_F4(MEdge_ptr e, MFace_ptr f);
  void ME_Rem_Face_R1(MEdge_ptr e, MFace_ptr f);
  void ME_Rem_Face_R2(MEdge_ptr e, MFace_ptr f);
  void ME_Rem_Face_R4(MEdge_ptr e, MFace_ptr f);
  static void (*ME_Rem_Face_jmp[MSTK_MAXREP])(MEdge_ptr e, MFace_ptr f) =
  {ME_Rem_Face_F1, ME_Rem_Face_F4, ME_Add_Face_R1, ME_Rem_Face_R2, 
   ME_Rem_Face_R4};

  void ME_Add_Region_F1(MEdge_ptr e, MRegion_ptr r);
  void ME_Add_Region_F4(MEdge_ptr e, MRegion_ptr r);
  void ME_Add_Region_R1(MEdge_ptr e, MRegion_ptr r);
  void ME_Add_Region_R2(MEdge_ptr e, MRegion_ptr r);
  void ME_Add_Region_R4(MEdge_ptr e, MRegion_ptr r);
  static void (*ME_Add_Region_jmp[MSTK_MAXREP])(MEdge_ptr e, MRegion_ptr r) =
  {ME_Add_Region_F1, ME_Add_Region_F4, ME_Add_Region_R1, ME_Add_Region_R2, 
   ME_Add_Region_R4};

  void ME_Rem_Region_F1(MEdge_ptr e, MRegion_ptr r);
  void ME_Rem_Region_F4(MEdge_ptr e, MRegion_ptr r);
  void ME_Rem_Region_R1(MEdge_ptr e, MRegion_ptr r);
  void ME_Rem_Region_R2(MEdge_ptr e, MRegion_ptr r);
  void ME_Rem_Region_R4(MEdge_ptr e, MRegion_ptr r);
  static void (*ME_Rem_Region_jmp[MSTK_MAXREP])(MEdge_ptr e, MRegion_ptr r) =
  {ME_Rem_Region_F1, ME_Rem_Region_F4, ME_Rem_Region_R1, ME_Rem_Region_R2, 
   ME_Rem_Region_R4};


  MEdge_ptr MEs_Merge_FN(MEdge_ptr e1, MEdge_ptr e2);
  MEdge_ptr MEs_Merge_RN(MEdge_ptr e1, MEdge_ptr e2);
  static MEdge_ptr (*MEs_Merge_jmp[MSTK_MAXREP])(MEdge_ptr e1, MEdge_ptr e2) = 
  {MEs_Merge_FN, MEs_Merge_FN, MEs_Merge_RN, MEs_Merge_RN, MEs_Merge_RN};

  void ME_Lock_R1(MEdge_ptr e);
  void ME_Lock_R2(MEdge_ptr e);
  void ME_Lock_R4(MEdge_ptr e);
  static void (*ME_Lock_jmp[MSTK_MAXREP])(MEdge_ptr e) =
  {ME_Dummy1, ME_Dummy1, ME_Lock_R1, ME_Lock_R2, ME_Lock_R4};

  void ME_UnLock_R1(MEdge_ptr e);
  void ME_UnLock_R2(MEdge_ptr e);
  void ME_UnLock_R4(MEdge_ptr e);
  static void (*ME_UnLock_jmp[MSTK_MAXREP])(MEdge_ptr e) =
  {ME_Dummy1, ME_Dummy1, ME_UnLock_R1, ME_UnLock_R2, ME_UnLock_R4};

  int ME_IsLocked_F1(MEdge_ptr e);
  int ME_IsLocked_F4(MEdge_ptr e);
  int ME_IsLocked_R1(MEdge_ptr e);
  int ME_IsLocked_R2(MEdge_ptr e);
  int ME_IsLocked_R4(MEdge_ptr e);
  static int (*ME_IsLocked_jmp[MSTK_MAXREP])(MEdge_ptr e) =
  {ME_IsLocked_F1, ME_IsLocked_F4, ME_IsLocked_R1, ME_IsLocked_R2, ME_IsLocked_R4};

#ifdef __cplusplus
}
#endif
