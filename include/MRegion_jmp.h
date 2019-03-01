/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#define _H_MRegion_Private

#include "MSTK.h"
#include "MRegion.h"

#ifdef __cplusplus
extern "C" {
#endif


  int MR_Set_GInfo_Auto_FNR3R4(MRegion_ptr r);
  int MR_Set_GInfo_Auto_R1(MRegion_ptr r);
  int MR_Set_GInfo_Auto_R2(MRegion_ptr r);
  static int (*MR_Set_GInfo_Auto_jmp[MSTK_MAXREP])(MRegion_ptr) = 
  {MR_Set_GInfo_Auto_FNR3R4, MR_Set_GInfo_Auto_FNR3R4, MR_Set_GInfo_Auto_R1,
   MR_Set_GInfo_Auto_R2, MR_Set_GInfo_Auto_FNR3R4};

void MR_Set_RepType_FNR3R4(MRegion_ptr r);
void MR_Set_RepType_R1(MRegion_ptr r);
void MR_Set_RepType_R2(MRegion_ptr r);
static void (*MR_Set_RepType_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_Set_RepType_FNR3R4, MR_Set_RepType_FNR3R4, MR_Set_RepType_R1, 
 MR_Set_RepType_R2, MR_Set_RepType_FNR3R4};

void MR_Delete_F1F3R3R4(MRegion_ptr r, int keep);
void MR_Delete_F4(MRegion_ptr r, int keep);
void MR_Delete_R1(MRegion_ptr r, int keep);
void MR_Delete_R2(MRegion_ptr r, int keep);
static void (*MR_Delete_jmp[MSTK_MAXREP])(MRegion_ptr r, int keep) =
{MR_Delete_F1F3R3R4, MR_Delete_F4, MR_Delete_R1, MR_Delete_R2, 
 MR_Delete_F1F3R3R4};

void MR_Restore_F1F3R3R4(MRegion_ptr r);
void MR_Restore_F4(MRegion_ptr r);
void MR_Restore_R1(MRegion_ptr r);
void MR_Restore_R2(MRegion_ptr r);
static void (*MR_Restore_jmp[MSTK_MAXREP])(MRegion_ptr r) = 
{MR_Restore_F1F3R3R4, MR_Restore_F4, MR_Restore_R1, MR_Restore_R2, 
 MR_Restore_F1F3R3R4};

void MR_Destroy_For_MESH_Delete_R1(MRegion_ptr r);
void MR_Destroy_For_MESH_Delete_R2(MRegion_ptr r);
void MR_Destroy_For_MESH_Delete_FNR3R4(MRegion_ptr r);
static void (*MR_Destroy_For_MESH_Delete_jmp[MSTK_MAXREP])(MRegion_ptr r) = 
{MR_Destroy_For_MESH_Delete_FNR3R4, MR_Destroy_For_MESH_Delete_FNR3R4, MR_Destroy_For_MESH_Delete_R1, MR_Destroy_For_MESH_Delete_R2, MR_Destroy_For_MESH_Delete_FNR3R4};

void MR_Set_Faces_FNR3R4(MRegion_ptr r,int nf,MFace_ptr *mfaces,int *dirs);
void MR_Set_Faces_R1(MRegion_ptr r,int nf,MFace_ptr *mfaces,int *dirs);
void MR_Set_Faces_R2(MRegion_ptr r,int nf,MFace_ptr *mfaces,int *dirs);
static void 
(*MR_Set_Faces_jmp[MSTK_MAXREP])(MRegion_ptr r, int nf, MFace_ptr *mfaces, int *dirs) =
{MR_Set_Faces_FNR3R4, MR_Set_Faces_FNR3R4, MR_Set_Faces_R1, MR_Set_Faces_R2, 
 MR_Set_Faces_FNR3R4};

void MR_Set_Vertices_FNR3R4(MRegion_ptr r, int nv, MVertex_ptr *mvertices,
			    int nf, int **rfvtemplate);
void MR_Set_Vertices_R1(MRegion_ptr r, int nv, MVertex_ptr *mvertices, 
			    int nf, int **rfvtemplate);
void MR_Set_Vertices_R2(MRegion_ptr r, int nv, MVertex_ptr *mvertices, 
			  int nf, int **rfvtemplate);
static void 
(*MR_Set_Vertices_jmp[MSTK_MAXREP])(MRegion_ptr r, int nv, MVertex_ptr *mvertices,
				    int nf, int **rfvtemplate) =
{MR_Set_Vertices_FNR3R4, MR_Set_Vertices_FNR3R4, MR_Set_Vertices_R1, 
 MR_Set_Vertices_R2, MR_Set_Vertices_FNR3R4};

List_ptr MR_Vertices_FNR3R4(MRegion_ptr r);
List_ptr MR_Vertices_R1(MRegion_ptr r);
List_ptr MR_Vertices_R2(MRegion_ptr r);
static List_ptr (*MR_Vertices_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_Vertices_FNR3R4, MR_Vertices_FNR3R4, MR_Vertices_R1, MR_Vertices_R2, 
 MR_Vertices_FNR3R4};

void MR_VertexIDs_R1(MRegion_ptr r, int *nrv, int *rvertids);
void MR_VertexIDs_R2(MRegion_ptr r, int *nrv, int *rvertids);
void MR_VertexIDs_FNR3R4(MRegion_ptr r, int *nrv, int *rvertids);
static void (*MR_VertexIDs_jmp[MSTK_MAXREP])(MRegion_ptr r, int *nrv, int *rvertids) =
{MR_VertexIDs_FNR3R4, MR_VertexIDs_FNR3R4, MR_VertexIDs_R1, MR_VertexIDs_R2, MR_VertexIDs_FNR3R4};

List_ptr MR_Edges_FN(MRegion_ptr r);
List_ptr MR_Edges_R1(MRegion_ptr r);
List_ptr MR_Edges_R2(MRegion_ptr r);
List_ptr MR_Edges_R3R4(MRegion_ptr r);
static List_ptr (*MR_Edges_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_Edges_FN, MR_Edges_FN, MR_Edges_R1, MR_Edges_R2, MR_Edges_R3R4};

void MR_EdgeIDs_R1(MRegion_ptr r, int *nre, int *edgeids);
void MR_EdgeIDs_R2(MRegion_ptr r, int *nre, int *edgeids);
void MR_EdgeIDs_R3R4(MRegion_ptr r, int *nre, int *edgeids);
void MR_EdgeIDs_FN(MRegion_ptr r, int *nre, int *edgeids);
static void (*MR_EdgeIDs_jmp[MSTK_MAXREP])(MRegion_ptr r, int *nre, int *edgeids) =
{MR_EdgeIDs_FN, MR_EdgeIDs_FN, MR_EdgeIDs_R1, MR_EdgeIDs_R2, MR_EdgeIDs_R3R4};

List_ptr MR_Faces_FNR3R4(MRegion_ptr r);
List_ptr MR_Faces_R1(MRegion_ptr r);
List_ptr MR_Faces_R2(MRegion_ptr r);
static List_ptr (*MR_Faces_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_Faces_FNR3R4, MR_Faces_FNR3R4, MR_Faces_R1, MR_Faces_R2, MR_Faces_FNR3R4};

void MR_FaceIDs_R1(MRegion_ptr r, int *nrf, int *faceids);
void MR_FaceIDs_R2(MRegion_ptr r, int *nrf, int *faceids);
void MR_FaceIDs_FNR3R4(MRegion_ptr r, int *nrf, int *faceids);
static void (*MR_FaceIDs_jmp[MSTK_MAXREP])(MRegion_ptr r, int *nrf, int *faceids) =
{MR_FaceIDs_FNR3R4, MR_FaceIDs_FNR3R4, MR_FaceIDs_R1, MR_FaceIDs_R2, MR_FaceIDs_FNR3R4};

List_ptr MR_AdjRegions_FNR3R4(MRegion_ptr r);
List_ptr MR_AdjRegions_R1(MRegion_ptr r);
List_ptr MR_AdjRegions_R2(MRegion_ptr r);
static List_ptr (*MR_AdjRegions_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_AdjRegions_FNR3R4, MR_AdjRegions_FNR3R4, MR_AdjRegions_R1, 
 MR_AdjRegions_R2, MR_AdjRegions_FNR3R4};

void MR_AdjRegionIDs_R1(MRegion_ptr r, int *nradj, int *adjregids);
void MR_AdjRegionIDs_R2(MRegion_ptr r, int *nradj, int *adjregids);
void MR_AdjRegionIDs_FNR3R4(MRegion_ptr r, int *nradj, int *adjregids);
static void (*MR_AdjRegionIDs_jmp[MSTK_MAXREP])(MRegion_ptr r, int *nradj, int *adjregids) =
{MR_AdjRegionIDs_FNR3R4, MR_AdjRegionIDs_FNR3R4, MR_AdjRegionIDs_R1, MR_AdjRegionIDs_R2, MR_AdjRegionIDs_FNR3R4};

int MR_FaceDir_FNR3R4(MRegion_ptr r, MFace_ptr f);
int MR_FaceDir_R1(MRegion_ptr r, MFace_ptr f);
int MR_FaceDir_R2(MRegion_ptr r, MFace_ptr f);
static int (*MR_FaceDir_jmp[MSTK_MAXREP])(MRegion_ptr r, MFace_ptr f) =
{MR_FaceDir_FNR3R4, MR_FaceDir_FNR3R4, MR_FaceDir_R1, MR_FaceDir_R2, 
 MR_FaceDir_FNR3R4};

int MR_FaceDir_i_FNR3R4(MRegion_ptr r, int i);
int MR_FaceDir_i_R1(MRegion_ptr r, int i);
int MR_FaceDir_i_R2(MRegion_ptr r, int i);
static int (*MR_FaceDir_i_jmp[MSTK_MAXREP])(MRegion_ptr r, int i) =
{MR_FaceDir_i_FNR3R4, MR_FaceDir_i_FNR3R4, MR_FaceDir_i_R1, MR_FaceDir_i_R2, 
 MR_FaceDir_i_FNR3R4};

int MR_Rev_FaceDir_FNR3R4(MRegion_ptr r, MFace_ptr f);
int MR_Rev_FaceDir_R1(MRegion_ptr r, MFace_ptr f);
int MR_Rev_FaceDir_R2(MRegion_ptr r, MFace_ptr f);
static int (*MR_Rev_FaceDir_jmp[MSTK_MAXREP])(MRegion_ptr r, MFace_ptr f) =
{MR_Rev_FaceDir_FNR3R4, MR_Rev_FaceDir_FNR3R4, MR_Rev_FaceDir_R1, 
 MR_Rev_FaceDir_R2, MR_Rev_FaceDir_FNR3R4};

int MR_Rev_FaceDir_i_FNR3R4(MRegion_ptr r, int i);
int MR_Rev_FaceDir_i_R1(MRegion_ptr r, int i);
int MR_Rev_FaceDir_i_R2(MRegion_ptr r, int i);
static int (*MR_Rev_FaceDir_i_jmp[MSTK_MAXREP])(MRegion_ptr r, int i) =
{MR_Rev_FaceDir_i_FNR3R4, MR_Rev_FaceDir_i_FNR3R4, MR_Rev_FaceDir_i_R1, 
 MR_Rev_FaceDir_i_R2, MR_Rev_FaceDir_i_FNR3R4};

int MR_Set_FaceDir_FNR3R4(MRegion_ptr r, MFace_ptr f, int dir);
int MR_Set_FaceDir_R1(MRegion_ptr r, MFace_ptr f, int dir);
int MR_Set_FaceDir_R2(MRegion_ptr r, MFace_ptr f, int dir);
static int (*MR_Set_FaceDir_jmp[MSTK_MAXREP])(MRegion_ptr r, MFace_ptr f, int dir) =
{MR_Set_FaceDir_FNR3R4, MR_Set_FaceDir_FNR3R4, MR_Set_FaceDir_R1, 
 MR_Set_FaceDir_R2, MR_Set_FaceDir_FNR3R4};

int MR_Set_FaceDir_i_FNR3R4(MRegion_ptr r, int i, int dir);
int MR_Set_FaceDir_i_R1(MRegion_ptr r, int i, int dir);
int MR_Set_FaceDir_i_R2(MRegion_ptr r, int i, int dir);
static int (*MR_Set_FaceDir_i_jmp[MSTK_MAXREP])(MRegion_ptr r, int i, int dir) =
{MR_Set_FaceDir_i_FNR3R4, MR_Set_FaceDir_i_FNR3R4, MR_Set_FaceDir_i_R1, 
 MR_Set_FaceDir_i_R2, MR_Set_FaceDir_i_FNR3R4};

  void MR_Rem_Face_FNR3R4(MRegion_ptr r, MFace_ptr f);
  void MR_Rem_Face_R1(MRegion_ptr r, MFace_ptr f);
  void MR_Rem_Face_R2(MRegion_ptr r, MFace_ptr f);
  static void 
  (*MR_Rem_Face_jmp[MSTK_MAXREP])(MRegion_ptr r, MFace_ptr f) =
  {MR_Rem_Face_FNR3R4, MR_Rem_Face_FNR3R4, MR_Rem_Face_R1, 
   MR_Rem_Face_R2, MR_Rem_Face_FNR3R4};

void MR_Replace_Faces_FNR3R4(MRegion_ptr r, int nold, MFace_ptr *f, int nnu,
                             MFace_ptr *nuf, int *nudir);
void MR_Replace_Faces_R1(MRegion_ptr r, int nold, MFace_ptr *f, int nnu, 
                         MFace_ptr *nuf, int *nudir);
void MR_Replace_Faces_R2(MRegion_ptr r, int nold, MFace_ptr *f, int nnu, 
                         MFace_ptr *nuf, int *nudir);
static void 
(*MR_Replace_Faces_jmp[MSTK_MAXREP])(MRegion_ptr r, int nold, MFace_ptr *f,
                                    int nnu, MFace_ptr * nuf, 
                                    int *nudir)=
{MR_Replace_Faces_FNR3R4, MR_Replace_Faces_FNR3R4, MR_Replace_Faces_R1, 
 MR_Replace_Faces_R2, MR_Replace_Faces_FNR3R4};

void MR_Replace_Face_i_FNR3R4(MRegion_ptr r, int i,MFace_ptr nuf,int nudir);
void MR_Replace_Face_i_R1(MRegion_ptr r, int i,MFace_ptr nuf,int nudir);
void MR_Replace_Face_i_R2(MRegion_ptr r, int i,MFace_ptr nuf,int nudir);
static void 
(*MR_Replace_Face_i_jmp[MSTK_MAXREP])(MRegion_ptr r, int i, MFace_ptr nuf, int nudir) =
{MR_Replace_Face_i_FNR3R4, MR_Replace_Face_i_FNR3R4, MR_Replace_Face_i_R1, 
 MR_Replace_Face_i_R2, MR_Replace_Face_i_FNR3R4};

void MR_Replace_Vertex_FNR3R4(MRegion_ptr r, MVertex_ptr v,MVertex_ptr nuv);
void MR_Replace_Vertex_R1(MRegion_ptr r, MVertex_ptr v,MVertex_ptr nuv);
void MR_Replace_Vertex_R2(MRegion_ptr r, MVertex_ptr v,MVertex_ptr nuv);
static void 
(*MR_Replace_Vertex_jmp[MSTK_MAXREP])(MRegion_ptr r, MVertex_ptr v, MVertex_ptr nuv) =
{MR_Replace_Vertex_FNR3R4, MR_Replace_Vertex_FNR3R4, MR_Replace_Vertex_R1, 
 MR_Replace_Vertex_R1, MR_Replace_Vertex_FNR3R4};

void MR_Replace_Vertex_i_FNR3R4(MRegion_ptr r, int i, MVertex_ptr nuv);
void MR_Replace_Vertex_i_R1(MRegion_ptr r, int i, MVertex_ptr nuv);
void MR_Replace_Vertex_i_R2(MRegion_ptr r, int i, MVertex_ptr nuv);
static void 
(*MR_Replace_Vertex_i_jmp[MSTK_MAXREP])(MRegion_ptr r, int i, MVertex_ptr nuv) =
{MR_Replace_Vertex_i_FNR3R4, MR_Replace_Vertex_i_FNR3R4, 
 MR_Replace_Vertex_i_R1, MR_Replace_Vertex_i_R2, MR_Replace_Vertex_i_FNR3R4};

void MR_Add_AdjRegion_FNR3R4(MRegion_ptr r, int facenum, MRegion_ptr aregion);
void MR_Add_AdjRegion_R1(MRegion_ptr r, int facenum, MRegion_ptr aregion);
void MR_Add_AdjRegion_R2(MRegion_ptr r, int facenum, MRegion_ptr aregion);
static void (*MR_Add_AdjRegion_jmp[MSTK_MAXREP])(MRegion_ptr r, int facenum, MRegion_ptr aregion) =
{MR_Add_AdjRegion_FNR3R4, MR_Add_AdjRegion_FNR3R4, MR_Add_AdjRegion_R1, 
 MR_Add_AdjRegion_R2, MR_Add_AdjRegion_FNR3R4};

void MR_Rem_AdjRegion_FNR3R4(MRegion_ptr r, MRegion_ptr ar);
void MR_Rem_AdjRegion_R1(MRegion_ptr r, MRegion_ptr ar);
void MR_Rem_AdjRegion_R2(MRegion_ptr r, MRegion_ptr ar);
static void (*MR_Rem_AdjRegion_jmp[MSTK_MAXREP])(MRegion_ptr r, MRegion_ptr aregion) =
{MR_Rem_AdjRegion_FNR3R4, MR_Rem_AdjRegion_FNR3R4, MR_Rem_AdjRegion_R1, 
 MR_Rem_AdjRegion_R2, MR_Rem_AdjRegion_FNR3R4};

int MR_Num_Faces_FNR3R4(MRegion_ptr r);
int MR_Num_Faces_R1(MRegion_ptr r);
int MR_Num_Faces_R2(MRegion_ptr r);
static int (*MR_Num_Faces_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_Num_Faces_FNR3R4, MR_Num_Faces_FNR3R4, MR_Num_Faces_R1, MR_Num_Faces_R2,
 MR_Num_Faces_FNR3R4};

int MR_UsesFace_FNR3R4(MRegion_ptr r, MEdge_ptr e);
int MR_UsesFace_R1(MRegion_ptr r, MEdge_ptr e);
int MR_UsesFace_R2(MRegion_ptr r, MEdge_ptr e);
static int (*MR_UsesFace_jmp[MSTK_MAXREP])(MRegion_ptr r, MEdge_ptr e) =
{MR_UsesFace_FNR3R4, MR_UsesFace_FNR3R4, MR_UsesFace_R1, MR_UsesFace_R2, 
 MR_UsesFace_FNR3R4}; 

int MR_UsesEdge_FN(MRegion_ptr r, MEdge_ptr e);
int MR_UsesEdge_R1(MRegion_ptr r, MEdge_ptr e);
int MR_UsesEdge_R2(MRegion_ptr r, MEdge_ptr e);
int MR_UsesEdge_R3R4(MRegion_ptr r, MEdge_ptr e);
static int (*MR_UsesEdge_jmp[MSTK_MAXREP])(MRegion_ptr r, MEdge_ptr e) =
{MR_UsesEdge_FN, MR_UsesEdge_FN, MR_UsesEdge_R1, MR_UsesEdge_R2, 
 MR_UsesEdge_R3R4}; 

int MR_UsesVertex_FNR3R4(MRegion_ptr r, MVertex_ptr e);
int MR_UsesVertex_R1(MRegion_ptr r, MVertex_ptr e);
int MR_UsesVertex_R2(MRegion_ptr r, MVertex_ptr e);
static int (*MR_UsesVertex_jmp[MSTK_MAXREP])(MRegion_ptr r, MEdge_ptr e) =
{MR_UsesVertex_FNR3R4, MR_UsesVertex_FNR3R4, MR_UsesVertex_R1, 
 MR_UsesVertex_R2, MR_UsesVertex_FNR3R4}; 


void MR_Update_ElementType_R1(MRegion_ptr r);
void MR_Update_ElementType_R2(MRegion_ptr r);
void MR_Update_ElementType_FNR3R4(MRegion_ptr r);
static void (*MR_Update_ElementType_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_Update_ElementType_FNR3R4, MR_Update_ElementType_FNR3R4, MR_Update_ElementType_R1, MR_Update_ElementType_R2, MR_Update_ElementType_FNR3R4};


#ifdef __cplusplus
}
#endif
