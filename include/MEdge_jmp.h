#define _H_MEdge_Private

#include "MSTK_types.h"
#include "MEdge.h"

#ifdef __cplusplus
extern "C" {
#endif

  void ME_Dummy1(MEdge_ptr e);

  void ME_Set_RepType_F1(MEdge_ptr e);
  void ME_Set_RepType_F4(MEdge_ptr e);
  static void (*ME_Set_RepType_jmp[MSTK_MAXREP])(MEdge_ptr e) =
  {ME_Set_RepType_F1, ME_Set_RepType_F4, ME_Dummy1, ME_Dummy1, ME_Dummy1};

  void ME_Delete_F1(MEdge_ptr e);
  void ME_Delete_F4(MEdge_ptr e);
  static void (*ME_Delete_jmp[MSTK_MAXREP])(MEdge_ptr e) = 
  {ME_Delete_F1, ME_Delete_F4, ME_Dummy1, ME_Dummy1, ME_Dummy1};

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
  int ME_Num_Regions_R1R2(MEdge_ptr e);
  int ME_Num_Regions_R3R4(MEdge_ptr e);
#ifdef DEBUG
  static int (*ME_Num_Regions_jmp[MSTK_MAXREP])(MEdge_ptr e) =
  {ME_Num_Regions_F1, ME_Num_Regions_F4, ME_Num_Regions_R1, ME_Num_Regions_R2,
   ME_Num_Regions_R4};
#else
  static int (*ME_Num_Regions_jmp[MSTK_MAXREP])(MEdge_ptr e) =
  {ME_Num_Regions_F1, ME_Num_Regions_F4, ME_Num_Regions_R1R2, 
   ME_Num_Regions_R1R2, ME_Num_Regions_R3R4};
#endif

  Set_ptr ME_Faces_F1(MEdge_ptr e);
  Set_ptr ME_Faces_F4(MEdge_ptr e);
  Set_ptr ME_Faces_R1(MEdge_ptr e);
  Set_ptr ME_Faces_R2(MEdge_ptr e);
  Set_ptr ME_Faces_R4(MEdge_ptr e);
  Set_ptr ME_Faces_R1R2(MEdge_ptr e);
  Set_ptr ME_Faces_R3R4(MEdge_ptr e);
#ifdef DEBUG
  static Set_ptr (*ME_Faces_jmp[MSTK_MAXREP])(MEdge_ptr e) =
  {ME_Faces_F1, ME_Faces_F4, ME_Faces_R1R2, ME_Faces_R1R2, ME_Faces_R4};
#else
  static Set_ptr (*ME_Faces_jmp[MSTK_MAXREP])(MEdge_ptr e) =
  {ME_Faces_F1, ME_Faces_F4, ME_Faces_R1R2, ME_Faces_R1R2, ME_Faces_R3R4};
#endif
  
  Set_ptr ME_Regions_F1(MEdge_ptr e);
  Set_ptr ME_Regions_F4(MEdge_ptr e);
  Set_ptr ME_Regions_R1(MEdge_ptr e);
  Set_ptr ME_Regions_R2(MEdge_ptr e);
  Set_ptr ME_Regions_R4(MEdge_ptr e);
  Set_ptr ME_Regions_R1R2(MEdge_ptr e);
  Set_ptr ME_Regions_R3R4(MEdge_ptr e);
#ifdef DEBUG
  static Set_ptr (*ME_Regions_jmp[MSTK_MAXREP])(MEdge_ptr e) =
  {ME_Regions_F1, ME_Regions_F4, ME_Regions_R1, ME_Regions_R2, ME_Regions_R4};
#else
  static Set_ptr (*ME_Regions_jmp[MSTK_MAXREP])(MEdge_ptr e) =
  {ME_Regions_F1, ME_Regions_F4, ME_Regions_R1R2, ME_Regions_R1R2, ME_Regions_R3R4};
#endif

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


#ifdef __cplusplus
}
#endif
