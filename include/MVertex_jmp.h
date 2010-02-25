#define _H_MVertex_Private

#include "MSTK_defines.h"
#include "MSTK_types.h"
#include "MVertex.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Jump tables */

  void MV_Set_RepType_F1(MVertex_ptr v);
  void MV_Set_RepType_F4(MVertex_ptr v);
  void MV_Set_RepType_R1(MVertex_ptr v);
  void MV_Set_RepType_R2(MVertex_ptr v);
  void MV_Set_RepType_R4(MVertex_ptr v);
  static void (*MV_Set_RepType_jmp[MSTK_MAXREP]) (MVertex_ptr v) = 
  {MV_Set_RepType_F1, MV_Set_RepType_F4, MV_Set_RepType_R1, 
   MV_Set_RepType_R2, MV_Set_RepType_R4};

  void MV_Delete_F1(MVertex_ptr v, int keep);
  void MV_Delete_F4(MVertex_ptr v, int keep);
  void MV_Delete_R1(MVertex_ptr v, int keep);
  void MV_Delete_R2(MVertex_ptr v, int keep);
  void MV_Delete_R4(MVertex_ptr v, int keep);
  static void (*MV_Delete_jmp[MSTK_MAXREP])(MVertex_ptr v, int keep) = 
  {MV_Delete_F1, MV_Delete_F4, MV_Delete_R1, MV_Delete_R2, MV_Delete_R4};

  void MV_Restore_F1(MVertex_ptr v);
  void MV_Restore_F4(MVertex_ptr v);
  void MV_Restore_R1(MVertex_ptr v);
  void MV_Restore_R2(MVertex_ptr v);
  void MV_Restore_R4(MVertex_ptr v);
  static void (*MV_Restore_jmp[MSTK_MAXREP])(MVertex_ptr v) = 
  {MV_Restore_F1, MV_Restore_F4, MV_Restore_R1, MV_Restore_R2, MV_Restore_R4};

  void MV_Destroy_For_MESH_Delete_F1(MVertex_ptr v);
  void MV_Destroy_For_MESH_Delete_F4(MVertex_ptr v);
  void MV_Destroy_For_MESH_Delete_R1(MVertex_ptr v);
  void MV_Destroy_For_MESH_Delete_R2(MVertex_ptr v);
  void MV_Destroy_For_MESH_Delete_R4(MVertex_ptr v);
  static void (*MV_Destroy_For_MESH_Delete_jmp[MSTK_MAXREP])(MVertex_ptr v) = 
  {MV_Destroy_For_MESH_Delete_F1, MV_Destroy_For_MESH_Delete_F4, 
   MV_Destroy_For_MESH_Delete_R1, MV_Destroy_For_MESH_Delete_R2, 
   MV_Destroy_For_MESH_Delete_R4};

  int MV_Num_AdjVertices_F1(MVertex_ptr v);
  int MV_Num_AdjVertices_F4(MVertex_ptr v);
  int MV_Num_AdjVertices_R1(MVertex_ptr v);
  int MV_Num_AdjVertices_R2(MVertex_ptr v);
  int MV_Num_AdjVertices_R4(MVertex_ptr v);
  static int (*MV_Num_AdjVertices_jmp[MSTK_MAXREP])(MVertex_ptr v) =
  {MV_Num_AdjVertices_F1, MV_Num_AdjVertices_F4, MV_Num_AdjVertices_R1,
   MV_Num_AdjVertices_R2, MV_Num_AdjVertices_R4};

  int MV_Num_Edges_F1(MVertex_ptr v);
  int MV_Num_Edges_F4(MVertex_ptr v);
  int MV_Num_Edges_R1(MVertex_ptr v);
  int MV_Num_Edges_R2(MVertex_ptr v);
  int MV_Num_Edges_R4(MVertex_ptr v);
  static int (*MV_Num_Edges_jmp[MSTK_MAXREP])(MVertex_ptr v) = 
  {MV_Num_Edges_F1, MV_Num_Edges_F4, MV_Num_Edges_R1, 
   MV_Num_Edges_R2, MV_Num_Edges_R4};
			
  int MV_Num_Faces_F1(MVertex_ptr v);
  int MV_Num_Faces_F4(MVertex_ptr v);
  int MV_Num_Faces_R1(MVertex_ptr v);
  int MV_Num_Faces_R2(MVertex_ptr v);
  int MV_Num_Faces_R4(MVertex_ptr v);

  static int (*MV_Num_Faces_jmp[MSTK_MAXREP])(MVertex_ptr v) = 
  {MV_Num_Faces_F1, MV_Num_Faces_F4, MV_Num_Faces_R1, 
   MV_Num_Faces_R2, MV_Num_Faces_R4};

			
  int MV_Num_Regions_F1(MVertex_ptr v);
  int MV_Num_Regions_F4(MVertex_ptr v);
  int MV_Num_Regions_R1(MVertex_ptr v);
  int MV_Num_Regions_R2(MVertex_ptr v);
  int MV_Num_Regions_R4(MVertex_ptr v);
  static int (*MV_Num_Regions_jmp[MSTK_MAXREP])(MVertex_ptr v) = 
  {MV_Num_Regions_F1, MV_Num_Regions_F4, MV_Num_Regions_R1, 
   MV_Num_Regions_R2, MV_Num_Regions_R4};


  List_ptr MV_AdjVertices_F1(MVertex_ptr v);
  List_ptr MV_AdjVertices_F4(MVertex_ptr v);
  List_ptr MV_AdjVertices_R1(MVertex_ptr v);
  List_ptr MV_AdjVertices_R2(MVertex_ptr v);
  List_ptr MV_AdjVertices_R4(MVertex_ptr v);
  static List_ptr (*MV_AdjVertices_jmp[MSTK_MAXREP])(MVertex_ptr v) =
  {MV_AdjVertices_F1, MV_AdjVertices_F4, MV_AdjVertices_R1,
   MV_AdjVertices_R2, MV_AdjVertices_R4};

  List_ptr MV_Edges_F1(MVertex_ptr v);
  List_ptr MV_Edges_F4(MVertex_ptr v);
  List_ptr MV_Edges_R1(MVertex_ptr v);
  List_ptr MV_Edges_R2(MVertex_ptr v);
  List_ptr MV_Edges_R4(MVertex_ptr v);
  static List_ptr (*MV_Edges_jmp[MSTK_MAXREP])(MVertex_ptr v) = 
  {MV_Edges_F1, MV_Edges_F4, MV_Edges_R1, MV_Edges_R2, 
   MV_Edges_R4};
			
  List_ptr MV_Faces_F1(MVertex_ptr v);
  List_ptr MV_Faces_F4(MVertex_ptr v);
  List_ptr MV_Faces_R1(MVertex_ptr v);
  List_ptr MV_Faces_R2(MVertex_ptr v);
  List_ptr MV_Faces_R4(MVertex_ptr v);
  static List_ptr (*MV_Faces_jmp[MSTK_MAXREP])(MVertex_ptr v) = 
  {MV_Faces_F1, MV_Faces_F4, MV_Faces_R1, MV_Faces_R2, 
   MV_Faces_R4};
			
  List_ptr MV_Regions_F1(MVertex_ptr v);
  List_ptr MV_Regions_F4(MVertex_ptr v);
  List_ptr MV_Regions_R1(MVertex_ptr v);
  List_ptr MV_Regions_R2(MVertex_ptr v);
  List_ptr MV_Regions_R4(MVertex_ptr v);
  static List_ptr (*MV_Regions_jmp[MSTK_MAXREP])(MVertex_ptr) = 
  {MV_Regions_F1, MV_Regions_F4, MV_Regions_R1, MV_Regions_R2, 
   MV_Regions_R4};

  void MV_Add_AdjVertex_F1(MVertex_ptr v, MVertex_ptr av);
  void MV_Add_AdjVertex_F4(MVertex_ptr v, MVertex_ptr av);
  void MV_Add_AdjVertex_R1(MVertex_ptr v, MVertex_ptr av);
  void MV_Add_AdjVertex_R2(MVertex_ptr v, MVertex_ptr av);
  void MV_Add_AdjVertex_R4(MVertex_ptr v, MVertex_ptr av);
  static void (*MV_Add_AdjVertex_jmp[MSTK_MAXREP])(MVertex_ptr v, MVertex_ptr av) = 
  {MV_Add_AdjVertex_F1, MV_Add_AdjVertex_F4, MV_Add_AdjVertex_R1, 
   MV_Add_AdjVertex_R2, MV_Add_AdjVertex_R4};
				
  void MV_Rem_AdjVertex_F1(MVertex_ptr v, MVertex_ptr av);
  void MV_Rem_AdjVertex_F4(MVertex_ptr v, MVertex_ptr av);
  void MV_Rem_AdjVertex_R1(MVertex_ptr v, MVertex_ptr av);
  void MV_Rem_AdjVertex_R2(MVertex_ptr v, MVertex_ptr av);
  void MV_Rem_AdjVertex_R4(MVertex_ptr v, MVertex_ptr av);
  static void (*MV_Rem_AdjVertex_jmp[MSTK_MAXREP])(MVertex_ptr v, MVertex_ptr av) = 
  {MV_Rem_AdjVertex_F1, MV_Rem_AdjVertex_F4, MV_Rem_AdjVertex_R1, 
   MV_Rem_AdjVertex_R2, MV_Rem_AdjVertex_R4};

  void MV_Add_Edge_F1(MVertex_ptr v, MEdge_ptr e);
  void MV_Add_Edge_F4(MVertex_ptr v, MEdge_ptr e);
  void MV_Add_Edge_R1(MVertex_ptr v, MEdge_ptr e);
  void MV_Add_Edge_R2(MVertex_ptr v, MEdge_ptr e);
  void MV_Add_Edge_R4(MVertex_ptr v, MEdge_ptr e);
  static void (*MV_Add_Edge_jmp[MSTK_MAXREP])(MVertex_ptr v, MEdge_ptr e) = 
  {MV_Add_Edge_F1, MV_Add_Edge_F4, MV_Add_Edge_R1, MV_Add_Edge_R2, 
   MV_Add_Edge_R4};

  void MV_Rem_Edge_F1(MVertex_ptr v, MEdge_ptr e);
  void MV_Rem_Edge_F4(MVertex_ptr v, MEdge_ptr e);
  void MV_Rem_Edge_R1(MVertex_ptr v, MEdge_ptr e);
  void MV_Rem_Edge_R2(MVertex_ptr v, MEdge_ptr e);
  void MV_Rem_Edge_R4(MVertex_ptr v, MEdge_ptr e);
  static void (*MV_Rem_Edge_jmp[MSTK_MAXREP])(MVertex_ptr v, MEdge_ptr e) = 
  {MV_Rem_Edge_F1, MV_Rem_Edge_F4, MV_Rem_Edge_R1, MV_Rem_Edge_R2,
   MV_Rem_Edge_R4};

  void MV_Add_Face_F1(MVertex_ptr v, MFace_ptr f);
  void MV_Add_Face_F4(MVertex_ptr v, MFace_ptr f);
  void MV_Add_Face_R1(MVertex_ptr v, MFace_ptr f);
  void MV_Add_Face_R2(MVertex_ptr v, MFace_ptr f);
  void MV_Add_Face_R4(MVertex_ptr v, MFace_ptr f);
  static void (*MV_Add_Face_jmp[MSTK_MAXREP])(MVertex_ptr v, MFace_ptr f) =
  {MV_Add_Face_F1, MV_Add_Face_F4, MV_Add_Face_R1, MV_Add_Face_R2, 
   MV_Add_Face_R4};

  void MV_Rem_Face_F1(MVertex_ptr v, MFace_ptr f);
  void MV_Rem_Face_F4(MVertex_ptr v, MFace_ptr f);
  void MV_Rem_Face_R1(MVertex_ptr v, MFace_ptr f);
  void MV_Rem_Face_R2(MVertex_ptr v, MFace_ptr f);
  void MV_Rem_Face_R4(MVertex_ptr v, MFace_ptr f);
  static void (*MV_Rem_Face_jmp[MSTK_MAXREP])(MVertex_ptr v, MFace_ptr f) =
  {MV_Rem_Face_F1, MV_Rem_Face_F4, MV_Rem_Face_R1, MV_Rem_Face_R2, 
   MV_Rem_Face_R4};


  void MV_Add_Region_F1(MVertex_ptr v, MRegion_ptr r);
  void MV_Add_Region_F4(MVertex_ptr v, MRegion_ptr r);
  void MV_Add_Region_R1(MVertex_ptr v, MRegion_ptr r);
  void MV_Add_Region_R2(MVertex_ptr v, MRegion_ptr r);
  void MV_Add_Region_R4(MVertex_ptr v, MRegion_ptr r);
  static void (*MV_Add_Region_jmp[MSTK_MAXREP])(MVertex_ptr v, MRegion_ptr r) =
  {MV_Add_Region_F1, MV_Add_Region_F4, MV_Add_Region_R1, MV_Add_Region_R2, 
   MV_Add_Region_R4};

  void MV_Rem_Region_F1(MVertex_ptr v, MRegion_ptr r);
  void MV_Rem_Region_F4(MVertex_ptr v, MRegion_ptr r);
  void MV_Rem_Region_R1(MVertex_ptr v, MRegion_ptr r);
  void MV_Rem_Region_R2(MVertex_ptr v, MRegion_ptr r);
  void MV_Rem_Region_R4(MVertex_ptr v, MRegion_ptr r);
  static void (*MV_Rem_Region_jmp[MSTK_MAXREP])(MVertex_ptr v, MRegion_ptr r) =
  {MV_Rem_Region_F1, MV_Rem_Region_F4, MV_Rem_Region_R1, MV_Rem_Region_R2, 
   MV_Rem_Region_R4};

  MVertex_ptr MVs_Merge_FN(MVertex_ptr v1, MVertex_ptr v2);
  MVertex_ptr MVs_Merge_R1R2(MVertex_ptr v1, MVertex_ptr v2);
  MVertex_ptr MVs_Merge_R3R4(MVertex_ptr v1, MVertex_ptr v2);
  static MVertex_ptr (*MVs_Merge_jmp[MSTK_MAXREP])(MVertex_ptr v1, MVertex_ptr v2) =
  {MVs_Merge_FN, MVs_Merge_FN, MVs_Merge_R1R2, MVs_Merge_R1R2, MVs_Merge_R3R4};


#ifdef __cplusplus
}
#endif
