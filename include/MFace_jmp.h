/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#define _H_MFace_Private

#include "MSTK_defines.h"
#include "MSTK_types.h"
#include "MFace.h"

#ifdef __cplusplus
extern "C" {
#endif

  int  MF_Set_GInfo_Auto_F1F3(MFace_ptr f);
  int  MF_Set_GInfo_Auto_F2F4(MFace_ptr f);
  int  MF_Set_GInfo_Auto_R1(MFace_ptr f);
  int  MF_Set_GInfo_Auto_R2(MFace_ptr f);
  int  MF_Set_GInfo_Auto_R3(MFace_ptr f);
  int  MF_Set_GInfo_Auto_R4(MFace_ptr f);
  static int (*MF_Set_GInfo_Auto_jmp[MSTK_MAXREP])(MFace_ptr f) =
  {MF_Set_GInfo_Auto_F1F3, MF_Set_GInfo_Auto_F2F4, MF_Set_GInfo_Auto_R1,
   MF_Set_GInfo_Auto_R2, MF_Set_GInfo_Auto_R4};

void MF_Set_Edges_F1F3(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_F2F4(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_R1(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_R2(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_R4(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
static void 
(*MF_Set_Edges_jmp[MSTK_MAXREP])(MFace_ptr f, int n, MEdge_ptr *e, int *dir) = 
  {MF_Set_Edges_F1F3, MF_Set_Edges_F2F4, MF_Set_Edges_R1, MF_Set_Edges_R2, 
   MF_Set_Edges_R4};

void MF_Set_Vertices_F1F3(MFace_ptr f, int n, MVertex_ptr *v);
void MF_Set_Vertices_F2F4(MFace_ptr f, int n, MVertex_ptr *v);
void MF_Set_Vertices_R1(MFace_ptr f, int n, MVertex_ptr *v);
void MF_Set_Vertices_R2(MFace_ptr f, int n, MVertex_ptr *v);
void MF_Set_Vertices_R4(MFace_ptr f, int n, MVertex_ptr *v);
static void 
(*MF_Set_Vertices_jmp[MSTK_MAXREP])(MFace_ptr f, int n, MVertex_ptr *v) = 
{MF_Set_Vertices_F1F3, MF_Set_Vertices_F2F4, MF_Set_Vertices_R1, 
 MF_Set_Vertices_R2, MF_Set_Vertices_R4};


void MF_Rem_Edge_F1F3(MFace_ptr f, MEdge_ptr e);
void MF_Rem_Edge_F2F4(MFace_ptr f, MEdge_ptr e);
void MF_Rem_Edge_R1(MFace_ptr f, MEdge_ptr e);
void MF_Rem_Edge_R2(MFace_ptr f, MEdge_ptr e);
void MF_Rem_Edge_R4(MFace_ptr f, MEdge_ptr e);
static void
(*MF_Rem_Edge_jmp[MSTK_MAXREP])(MFace_ptr f, MEdge_ptr e) =
{MF_Rem_Edge_F1F3, MF_Rem_Edge_F2F4, MF_Rem_Edge_R1, MF_Rem_Edge_R2, MF_Rem_Edge_R4};

void MF_Replace_Edges_F1F3(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_F2F4(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_R1(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_R2(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_R4(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
static void (*MF_Replace_Edges_jmp[MSTK_MAXREP])
     (MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges)=
{MF_Replace_Edges_F1F3, MF_Replace_Edges_F2F4, MF_Replace_Edges_R1, 
 MF_Replace_Edges_R2, MF_Replace_Edges_R4};

void MF_Replace_Edges_i_F1F3(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_i_F2F4(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_i_R1(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_i_R2(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_i_R4(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
static void 
(*MF_Replace_Edges_i_jmp[MSTK_MAXREP])(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges) = 
{MF_Replace_Edges_i_F1F3, MF_Replace_Edges_i_F2F4, MF_Replace_Edges_i_R1, 
 MF_Replace_Edges_i_R2, MF_Replace_Edges_i_R4};

void MF_Replace_Vertex_i_F1(MFace_ptr f, int i, MVertex_ptr v);
void MF_Replace_Vertex_i_F4(MFace_ptr f, int i, MVertex_ptr v);
void MF_Replace_Vertex_i_R1(MFace_ptr f, int i, MVertex_ptr v);
void MF_Replace_Vertex_i_R2(MFace_ptr f, int i, MVertex_ptr v);
void MF_Replace_Vertex_i_R4(MFace_ptr f, int i, MVertex_ptr v);
static void 
(*MF_Replace_Vertex_i_jmp[MSTK_MAXREP])(MFace_ptr f, int i, MVertex_ptr v) =
{MF_Replace_Vertex_i_F1, MF_Replace_Vertex_i_F4, MF_Replace_Vertex_i_R1, 
 MF_Replace_Vertex_i_R2, MF_Replace_Vertex_i_R4};

void MF_Replace_Vertex_F1(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
void MF_Replace_Vertex_F4(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
void MF_Replace_Vertex_R1(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
void MF_Replace_Vertex_R2(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
void MF_Replace_Vertex_R4(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
static void
(*MF_Replace_Vertex_jmp[MSTK_MAXREP])(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv) =
{MF_Replace_Vertex_F1, MF_Replace_Vertex_F4, MF_Replace_Vertex_R1, 
 MF_Replace_Vertex_R2, MF_Replace_Vertex_R4};

void MF_Insert_Vertex_F1(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
void MF_Insert_Vertex_F4(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
void MF_Insert_Vertex_R1(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
void MF_Insert_Vertex_R2(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
void MF_Insert_Vertex_R4(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
static void 
(*MF_Insert_Vertex_jmp[MSTK_MAXREP])(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v) =
{MF_Insert_Vertex_F1, MF_Insert_Vertex_F4, MF_Insert_Vertex_R1, 
 MF_Insert_Vertex_R2, MF_Insert_Vertex_R4};

void MF_Insert_Vertex_i_F1(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_F4(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_R1(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_R2(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_R4(MFace_ptr f, MVertex_ptr nuv, int i);
static void 
(*MF_Insert_Vertex_i_jmp[MSTK_MAXREP])(MFace_ptr f, MVertex_ptr nuv, int i) =
{MF_Insert_Vertex_i_F1, MF_Insert_Vertex_i_F4, MF_Insert_Vertex_i_R1, 
 MF_Insert_Vertex_i_R2, MF_Insert_Vertex_i_R4};

void MF_Rem_Vertex_F1(MFace_ptr f, MVertex_ptr v);
void MF_Rem_Vertex_F4(MFace_ptr f, MVertex_ptr v);
void MF_Rem_Vertex_R1(MFace_ptr f, MVertex_ptr v);
void MF_Rem_Vertex_R2(MFace_ptr f, MVertex_ptr v);
void MF_Rem_Vertex_R4(MFace_ptr f, MVertex_ptr v);
static void
(*MF_Rem_Vertex_jmp[MSTK_MAXREP])(MFace_ptr f, MVertex_ptr e) =
{MF_Rem_Vertex_F1, MF_Rem_Vertex_F4, MF_Rem_Vertex_R1, MF_Rem_Vertex_R2, MF_Rem_Vertex_R4};


int MF_Rev_EdgeDir_F1F3(MFace_ptr f, MEdge_ptr e);
int MF_Rev_EdgeDir_F2F4(MFace_ptr f, MEdge_ptr e);
int MF_Rev_EdgeDir_R1(MFace_ptr f, MEdge_ptr e);
int MF_Rev_EdgeDir_R2(MFace_ptr f, MEdge_ptr e);
int MF_Rev_EdgeDir_R4(MFace_ptr f, MEdge_ptr e);
static int (*MF_Rev_EdgeDir_jmp[MSTK_MAXREP])(MFace_ptr f, MEdge_ptr e) =
{MF_Rev_EdgeDir_F1F3, MF_Rev_EdgeDir_F2F4, MF_Rev_EdgeDir_R1, MF_Rev_EdgeDir_R2, MF_Rev_EdgeDir_R4};

int MF_Rev_EdgeDir_i_F1F3(MFace_ptr f, int i);
int MF_Rev_EdgeDir_i_F2F4(MFace_ptr f, int i);
int MF_Rev_EdgeDir_i_R1(MFace_ptr f, int i);
int MF_Rev_EdgeDir_i_R2(MFace_ptr f, int i);
int MF_Rev_EdgeDir_i_R4(MFace_ptr f, int i);
static int (*MF_Rev_EdgeDir_i_jmp[MSTK_MAXREP])(MFace_ptr f, int i) =
{MF_Rev_EdgeDir_i_F1F3, MF_Rev_EdgeDir_i_F2F4, MF_Rev_EdgeDir_i_R1, MF_Rev_EdgeDir_i_R2,
 MF_Rev_EdgeDir_i_R4};


int MFs_AreSame_F1(MFace_ptr f1, MFace_ptr f2);
int MFs_AreSame_F4(MFace_ptr f1, MFace_ptr f2);
int MFs_AreSame_R1(MFace_ptr f1, MFace_ptr f2);
int MFs_AreSame_R2(MFace_ptr f1, MFace_ptr f2);
int MFs_AreSame_R4(MFace_ptr f1, MFace_ptr f2);
static int 
(*MFs_AreSame_jmp[MSTK_MAXREP])(MFace_ptr f1, MFace_ptr f2) =
{MFs_AreSame_F1, MFs_AreSame_F4, MFs_AreSame_R1, MFs_AreSame_R2, 
 MFs_AreSame_R4};

int MF_Num_Vertices_F1F3(MFace_ptr f);
int MF_Num_Vertices_F2F4(MFace_ptr f);
int MF_Num_Vertices_R1(MFace_ptr f);
int MF_Num_Vertices_R2(MFace_ptr f);
int MF_Num_Vertices_R4(MFace_ptr f);
static int (*MF_Num_Vertices_jmp[MSTK_MAXREP])(MFace_ptr f) = 
{MF_Num_Vertices_F1F3, MF_Num_Vertices_F2F4, MF_Num_Vertices_R1, 
 MF_Num_Vertices_R2, MF_Num_Vertices_R4};

int MF_Num_Edges_F1F3(MFace_ptr f);
int MF_Num_Edges_F2F4(MFace_ptr f);
int MF_Num_Edges_R1(MFace_ptr f);
int MF_Num_Edges_R2(MFace_ptr f);
int MF_Num_Edges_R4(MFace_ptr f);
static int (*MF_Num_Edges_jmp[MSTK_MAXREP])(MFace_ptr f) =
{MF_Num_Edges_F1F3, MF_Num_Edges_F2F4, MF_Num_Edges_R1, MF_Num_Edges_R2, 
 MF_Num_Edges_R4};

List_ptr MF_Vertices_F1F3(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_F2F4(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_R1(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_R2(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_R4(MFace_ptr f, int dir, MVertex_ptr v);
static List_ptr (*MF_Vertices_jmp[MSTK_MAXREP])(MFace_ptr f, int dir, MVertex_ptr v) =
{MF_Vertices_F1F3, MF_Vertices_F2F4, MF_Vertices_R1, MF_Vertices_R2, 
 MF_Vertices_R4};

void MF_VertexIDs_F1F3(MFace_ptr f, int dir, int startvid, int *nfv, int *fvertids);
void MF_VertexIDs_F2F4(MFace_ptr f, int dir, int startvid, int *nfv, int *fvertids);
void MF_VertexIDs_R1(MFace_ptr f, int dir, int startvid, int *nfv, int *fvertids);
void MF_VertexIDs_R2(MFace_ptr f, int dir, int startvid, int *nfv, int *fvertids);
void MF_VertexIDs_R4(MFace_ptr f, int dir, int startvid, int *nfv, int *fvertids);
static void (*MF_VertexIDs_jmp[MSTK_MAXREP])(MFace_ptr f, int dir, int startvid, int *nfv, int *fvertids) =
{MF_VertexIDs_F1F3, MF_VertexIDs_F2F4, MF_VertexIDs_R1, MF_VertexIDs_R2, 
 MF_VertexIDs_R4};

List_ptr MF_Edges_F1F3(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_F2F4(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_R1(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_R2(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_R4(MFace_ptr f, int dir, MVertex_ptr v);
static List_ptr (*MF_Edges_jmp[MSTK_MAXREP])(MFace_ptr f, int dir, MVertex_ptr v) =
{MF_Edges_F1F3, MF_Edges_F2F4, MF_Edges_R1, MF_Edges_R2, MF_Edges_R4};

void MF_EdgeIDs_F1F3(MFace_ptr f, int dir, int startvid, int *nfv, int *fvertids);
void MF_EdgeIDs_F2F4(MFace_ptr f, int dir, int startvid, int *nfv, int *fvertids);
void MF_EdgeIDs_R1(MFace_ptr f, int dir, int startvid, int *nfe, int *fedgeids);
void MF_EdgeIDs_R2(MFace_ptr f, int dir, int startvid, int *nfe, int *fedgeids);
void MF_EdgeIDs_R4(MFace_ptr f, int dir, int startvid, int *nfe, int *fedgeids);
static void (*MF_EdgeIDs_jmp[MSTK_MAXREP])(MFace_ptr f, int dir, int startvid, int *nfe, int *fedgeids) =
{MF_EdgeIDs_F1F3, MF_EdgeIDs_F2F4, MF_EdgeIDs_R1, MF_EdgeIDs_R2, 
 MF_EdgeIDs_R4};

int MF_EdgeDir_F1F3(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_F2F4(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_R1(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_R2(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_R4(MFace_ptr f, MEdge_ptr e);
static int (*MF_EdgeDir_jmp[MSTK_MAXREP])(MFace_ptr f, MEdge_ptr e) =
{MF_EdgeDir_F1F3, MF_EdgeDir_F2F4, MF_EdgeDir_R1, MF_EdgeDir_R2, MF_EdgeDir_R4};

int MF_EdgeDir_i_F1F3(MFace_ptr f, int i);
int MF_EdgeDir_i_F2F4(MFace_ptr f, int i);
int MF_EdgeDir_i_R1(MFace_ptr f, int i);
int MF_EdgeDir_i_R2(MFace_ptr f, int i);
int MF_EdgeDir_i_R4(MFace_ptr f, int i);
static int (*MF_EdgeDir_i_jmp[MSTK_MAXREP])(MFace_ptr f, int i) =
{MF_EdgeDir_i_F1F3, MF_EdgeDir_i_F2F4, MF_EdgeDir_i_R1, MF_EdgeDir_i_R2,
 MF_EdgeDir_i_R4};

int MF_UsesEdge_F1F3(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_F2F4(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_R1(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_R2(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_R4(MFace_ptr f, MEdge_ptr e);
static int (*MF_UsesEdge_jmp[MSTK_MAXREP])(MFace_ptr f, MEdge_ptr e) =
{MF_UsesEdge_F1F3, MF_UsesEdge_F2F4, MF_UsesEdge_R1, MF_UsesEdge_R2,
 MF_UsesEdge_R4};

int MF_UsesVertex_F1F3(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_F2F4(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_R1(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_R2(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_R4(MFace_ptr f, MVertex_ptr v);
static int (*MF_UsesVertex_jmp[MSTK_MAXREP])(MFace_ptr f, MVertex_ptr v) =
{MF_UsesVertex_F1F3, MF_UsesVertex_F2F4, MF_UsesVertex_R1, MF_UsesVertex_R2, 
 MF_UsesVertex_R4};

List_ptr MF_Regions_F1F3(MFace_ptr f);
List_ptr MF_Regions_F2F4(MFace_ptr f);
List_ptr MF_Regions_R1(MFace_ptr f);
List_ptr MF_Regions_R2(MFace_ptr f);
List_ptr MF_Regions_R4(MFace_ptr f);
static List_ptr (*MF_Regions_jmp[MSTK_MAXREP])(MFace_ptr f) =
{MF_Regions_F1F3, MF_Regions_F2F4, MF_Regions_R1, MF_Regions_R2, MF_Regions_R4};

void MF_RegionIDs_F1F3(MFace_ptr f, int *nfr, int *fregids);
void MF_RegionIDs_F2F4(MFace_ptr f, int *nfe, int *fregids);
void MF_RegionIDs_R1(MFace_ptr f, int *nfe, int *fregids);
void MF_RegionIDs_R2(MFace_ptr f, int *nfe, int *fregids);
void MF_RegionIDs_R4(MFace_ptr f, int *nfe, int *fregids);
static void (*MF_RegionIDs_jmp[MSTK_MAXREP])(MFace_ptr f, int *nfe, int *fregids) =
{MF_RegionIDs_F1F3, MF_RegionIDs_F2F4, MF_RegionIDs_R1, MF_RegionIDs_R2, MF_RegionIDs_R4};

MRegion_ptr MF_Region_F1F3(MFace_ptr f, int side);
MRegion_ptr MF_Region_F2F4(MFace_ptr f, int side);
MRegion_ptr MF_Region_R1(MFace_ptr f, int side);
MRegion_ptr MF_Region_R2(MFace_ptr f, int side);
MRegion_ptr MF_Region_R4(MFace_ptr f, int side);
static MRegion_ptr (*MF_Region_jmp[MSTK_MAXREP])(MFace_ptr f, int side) =
{MF_Region_F1F3, MF_Region_F2F4, MF_Region_R1, MF_Region_R2, MF_Region_R4};

int MF_RegionID_F1F3(MFace_ptr f, int side);
int MF_RegionID_F2F4(MFace_ptr f, int side);
int MF_RegionID_R1(MFace_ptr f, int side);
int MF_RegionID_R2(MFace_ptr f, int side);
int MF_RegionID_R4(MFace_ptr f, int side);
static int (*MF_RegionID_jmp[MSTK_MAXREP])(MFace_ptr f, int side) =
{MF_RegionID_F1F3, MF_RegionID_F2F4, MF_RegionID_R1, MF_RegionID_R2, MF_RegionID_R4};


void MF_Dummy1(MFace_ptr f);
void MF_Dummy2a(MFace_ptr f, int i);


void MF_Set_RepType_F1(MFace_ptr f);
void MF_Set_RepType_F4(MFace_ptr f);
void MF_Set_RepType_R1(MFace_ptr f);
void MF_Set_RepType_R2(MFace_ptr f);
void MF_Set_RepType_R4(MFace_ptr f);
static void (*MF_Set_RepType_jmp[MSTK_MAXREP])(MFace_ptr f) =
{MF_Set_RepType_F1, MF_Set_RepType_F4, MF_Set_RepType_R1, MF_Set_RepType_R2,
 MF_Set_RepType_R4};

void MF_Delete_F1(MFace_ptr f, int keep);
void MF_Delete_F4(MFace_ptr f, int keep);
void MF_Delete_R1(MFace_ptr f, int keep);
void MF_Delete_R2(MFace_ptr f, int keep);
void MF_Delete_R4(MFace_ptr f, int keep);
static void (*MF_Delete_jmp[MSTK_MAXREP])(MFace_ptr f, int keep) = 
{MF_Delete_F1, MF_Delete_F4, MF_Delete_R1, MF_Delete_R2, MF_Delete_R4};

void MF_Restore_F1(MFace_ptr f);
void MF_Restore_F4(MFace_ptr f);
void MF_Restore_R4(MFace_ptr f);
static void (*MF_Restore_jmp[MSTK_MAXREP])(MFace_ptr f) = 
{MF_Restore_F1, MF_Restore_F4, MF_Dummy1, MF_Dummy1, MF_Restore_R4};

void MF_Destroy_For_MESH_Delete_F1(MFace_ptr f);
void MF_Destroy_For_MESH_Delete_F4(MFace_ptr f);
void MF_Destroy_For_MESH_Delete_R1(MFace_ptr f);
void MF_Destroy_For_MESH_Delete_R2(MFace_ptr f);
void MF_Destroy_For_MESH_Delete_R4(MFace_ptr f);
static void (*MF_Destroy_For_MESH_Delete_jmp[MSTK_MAXREP])(MFace_ptr f) = 
{MF_Destroy_For_MESH_Delete_F1, MF_Destroy_For_MESH_Delete_F4, 
 MF_Destroy_For_MESH_Delete_R1, MF_Destroy_For_MESH_Delete_R2, 
 MF_Destroy_For_MESH_Delete_R4};


int MF_Num_AdjFaces_F1(MFace_ptr f);
int MF_Num_AdjFaces_F4(MFace_ptr f);
int MF_Num_AdjFaces_R1(MFace_ptr f);
int MF_Num_AdjFaces_R2(MFace_ptr f);
int MF_Num_AdjFaces_R4(MFace_ptr f);
static int (*MF_Num_AdjFaces_jmp[MSTK_MAXREP])(MFace_ptr f) =
{MF_Num_AdjFaces_F1, MF_Num_AdjFaces_F4, MF_Num_AdjFaces_R1, 
 MF_Num_AdjFaces_R2, MF_Num_AdjFaces_R4};

List_ptr MF_AdjFaces_F1(MFace_ptr f);
List_ptr MF_AdjFaces_F4(MFace_ptr f);
List_ptr MF_AdjFaces_R1(MFace_ptr f);
List_ptr MF_AdjFaces_R2(MFace_ptr f);
List_ptr MF_AdjFaces_R4(MFace_ptr f);
static List_ptr (*MF_AdjFaces_jmp[MSTK_MAXREP])(MFace_ptr f) =
{MF_AdjFaces_F1, MF_AdjFaces_F4, MF_AdjFaces_R1, MF_AdjFaces_R2,
 MF_AdjFaces_R4};


void MF_Add_Region_F1(MFace_ptr f, MRegion_ptr r, int side);
void MF_Add_Region_F1F3(MFace_ptr f, MRegion_ptr r, int side);
void MF_Add_Region_F4(MFace_ptr f, MRegion_ptr r, int side);
void MF_Dummy2b(MFace_ptr f, MRegion_ptr r, int side);
void MF_Add_Region_R4(MFace_ptr f, MRegion_ptr r, int side);
static void (*MF_Add_Region_jmp[MSTK_MAXREP])(MFace_ptr f, MRegion_ptr r, int side) =
{MF_Add_Region_F1, MF_Add_Region_F4, MF_Dummy2b, MF_Dummy2b, MF_Add_Region_R4};


void MF_Rem_Region_F1(MFace_ptr f, MRegion_ptr r);
void MF_Rem_Region_F1F3(MFace_ptr f, MRegion_ptr r);
void MF_Rem_Region_F4(MFace_ptr f, MRegion_ptr r);
void MF_Rem_Region_R1(MFace_ptr f, MRegion_ptr r);
void MF_Rem_Region_R2(MFace_ptr f, MRegion_ptr r);
void MF_Rem_Region_R4(MFace_ptr f, MRegion_ptr r);
static void (*MF_Rem_Region_jmp[MSTK_MAXREP])(MFace_ptr f, MRegion_ptr r) =
{MF_Rem_Region_F1, MF_Rem_Region_F4, MF_Rem_Region_R1, MF_Rem_Region_R2, 
 MF_Rem_Region_R4};

void MF_Add_AdjFace_F1(MFace_ptr f, int edgnum, MFace_ptr af);
void MF_Add_AdjFace_F4(MFace_ptr f, int edgnum, MFace_ptr af);
void MF_Add_AdjFace_R1(MFace_ptr f, int edgnum, MFace_ptr af);
void MF_Add_AdjFace_R2(MFace_ptr f, int edgnum, MFace_ptr af);
void MF_Add_AdjFace_R4(MFace_ptr f, int edgnum, MFace_ptr af);
static void (*MF_Add_AdjFace_jmp[MSTK_MAXREP])(MFace_ptr f, int edgnum, MFace_ptr af)=
{MF_Add_AdjFace_F1, MF_Add_AdjFace_F4, MF_Add_AdjFace_R1, MF_Add_AdjFace_R2,
 MF_Add_AdjFace_R4};

void MF_Rem_AdjFace_F1(MFace_ptr f, int edgnum, MFace_ptr af);
void MF_Rem_AdjFace_F4(MFace_ptr f, int edgnum, MFace_ptr af);
void MF_Rem_AdjFace_R1(MFace_ptr f, int edgnum, MFace_ptr af);
void MF_Rem_AdjFace_R2(MFace_ptr f, int edgnum, MFace_ptr af);
void MF_Rem_AdjFace_R4(MFace_ptr f, int edgnum, MFace_ptr af);
static void (*MF_Rem_AdjFace_jmp[MSTK_MAXREP])(MFace_ptr f, int edgnum, MFace_ptr af)=
{MF_Rem_AdjFace_F1, MF_Rem_AdjFace_F4, MF_Rem_AdjFace_R1, MF_Rem_AdjFace_R2, 
 MF_Rem_AdjFace_R4};


MFace_ptr MFs_Merge_FN(MFace_ptr f1, MFace_ptr f2, int topoflag);
MFace_ptr MFs_Merge_R1R2(MFace_ptr f1, MFace_ptr f2, int topoflag);
MFace_ptr MFs_Merge_R3R4(MFace_ptr f1, MFace_ptr f2, int topoflag);
static MFace_ptr (*MFs_Merge_jmp[MSTK_MAXREP])(MFace_ptr f1, MFace_ptr f2,
                                               int topoflag) =
{MFs_Merge_FN, MFs_Merge_FN, MFs_Merge_R1R2, MFs_Merge_R1R2, MFs_Merge_R3R4};


MFace_ptr MF_NextInHash_F1(MFace_ptr f);
MFace_ptr MF_NextInHash_F4(MFace_ptr f);
MFace_ptr MF_NextInHash_R1(MFace_ptr f);
MFace_ptr MF_NextInHash_R2(MFace_ptr f);
MFace_ptr MF_NextInHash_R4(MFace_ptr f);
static MFace_ptr (*MF_NextInHash_jmp[MSTK_MAXREP])(MFace_ptr f) =
{MF_NextInHash_F1, MF_NextInHash_F4, MF_NextInHash_R1, MF_NextInHash_R2, MF_NextInHash_R4};

void MF_Set_NextInHash_F1(MFace_ptr f, MFace_ptr next);
void MF_Set_NextInHash_F4(MFace_ptr f, MFace_ptr next);
void MF_Set_NextInHash_R1(MFace_ptr f, MFace_ptr next);
void MF_Set_NextInHash_R2(MFace_ptr f, MFace_ptr next);
void MF_Set_NextInHash_R4(MFace_ptr f, MFace_ptr next);
static void (*MF_Set_NextInHash_jmp[MSTK_MAXREP])(MFace_ptr f, MFace_ptr next) =
{MF_Set_NextInHash_F1, MF_Set_NextInHash_F4, MF_Set_NextInHash_R1, MF_Set_NextInHash_R2, MF_Set_NextInHash_R4};

void MF_HashKey_F1(MFace_ptr f, unsigned int *pn, void* **pp);
void MF_HashKey_F4(MFace_ptr f, unsigned int *pn, void* **pp);
void MF_HashKey_R1(MFace_ptr f, unsigned int *pn, void* **pp);
void MF_HashKey_R2(MFace_ptr f, unsigned int *pn, void* **pp);
void MF_HashKey_R4(MFace_ptr f, unsigned int *pn, void* **pp);
static void (*MF_HashKey_jmp[MSTK_MAXREP])(MFace_ptr f, unsigned int *pn, void* **pp) =
{MF_HashKey_F1, MF_HashKey_F4, MF_HashKey_R1, MF_HashKey_R2, MF_HashKey_R4};

void MF_Lock_R1(MFace_ptr f);
void MF_Lock_R2(MFace_ptr f);
static void (*MF_Lock_jmp[MSTK_MAXREP])(MFace_ptr f) =
{MF_Dummy1, MF_Dummy1, MF_Lock_R1, MF_Lock_R2, MF_Dummy1};

void MF_UnLock_R1(MFace_ptr f);
void MF_UnLock_R2(MFace_ptr f);
static void (*MF_UnLock_jmp[MSTK_MAXREP])(MFace_ptr f) =
{MF_Dummy1, MF_Dummy1, MF_UnLock_R1, MF_UnLock_R2, MF_Dummy1};

int MF_IsLocked_F1(MFace_ptr f);
int MF_IsLocked_F4(MFace_ptr f);
int MF_IsLocked_R1(MFace_ptr f);
int MF_IsLocked_R2(MFace_ptr f);
int MF_IsLocked_R4(MFace_ptr f);
static int (*MF_IsLocked_jmp[MSTK_MAXREP])(MFace_ptr f) =
{MF_IsLocked_F1, MF_IsLocked_F4, MF_IsLocked_R1, MF_IsLocked_R2, MF_IsLocked_R4};

#ifdef __cplusplus
}
#endif
