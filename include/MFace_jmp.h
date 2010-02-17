#define _H_MFace_Private

#include "MSTK.h"
#include "MFace.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* The following functions go through the appropriate F1 or F4
     function to the common FN function. So in the non-Debug version
     call the common function directly */

#ifdef DEBUG

  int  MF_Set_GInfo_Auto_FN(MFace_ptr f);
  int  MF_Set_GInfo_Auto_F1(MFace_ptr f);
  int  MF_Set_GInfo_Auto_F4(MFace_ptr f);
  int  MF_Set_GInfo_Auto_R1(MFace_ptr f);
  int  MF_Set_GInfo_Auto_R2(MFace_ptr f);
  int  MF_Set_GInfo_Auto_R3(MFace_ptr f);
  int  MF_Set_GInfo_Auto_R4(MFace_ptr f);
  int  MF_Set_GInfo_Auto_RN(MFace_ptr f);
  static int (*MF_Set_GInfo_Auto_jmp[MSTK_MAXREP])(MFace_ptr f) =
  {MF_Set_GInfo_Auto_F1, MF_Set_GInfo_Auto_F4, MF_Set_GInfo_Auto_R1,
   MF_Set_GInfo_Auto_R2, MF_Set_GInfo_Auto_R4};

void MF_Set_Edges_FN(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_F1(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_F4(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_R1(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_R2(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_R4(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_RN(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
static void 
(*MF_Set_Edges_jmp[MSTK_MAXREP])(MFace_ptr f, int n, MEdge_ptr *e, int *dir) = 
  {MF_Set_Edges_F1, MF_Set_Edges_F4, MF_Set_Edges_R1, MF_Set_Edges_R2, 
   MF_Set_Edges_R4};

void MF_Replace_Edges_FN(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_F1(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_F4(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_R1(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_R2(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_R4(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_RN(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
static void (*MF_Replace_Edges_jmp[MSTK_MAXREP])
     (MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges)=
{MF_Replace_Edges_F1, MF_Replace_Edges_F4, MF_Replace_Edges_R1, 
 MF_Replace_Edges_R2, MF_Replace_Edges_R4};

void MF_Replace_Edges_i_FN(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_i_F1(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_i_F4(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_i_R1(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_i_R2(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_i_R4(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_i_RN(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
static void 
(*MF_Replace_Edges_i_jmp[MSTK_MAXREP])(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges) = 
{MF_Replace_Edges_i_F1, MF_Replace_Edges_i_F4, MF_Replace_Edges_i_R1, 
 MF_Replace_Edges_i_R2, MF_Replace_Edges_i_R4};

void MF_Set_Vertices_F1(MFace_ptr f, int n, MVertex_ptr *v);
void MF_Set_Vertices_F4(MFace_ptr f, int n, MVertex_ptr *v);
void MF_Set_Vertices_R1(MFace_ptr f, int n, MVertex_ptr *v);
void MF_Set_Vertices_R2(MFace_ptr f, int n, MVertex_ptr *v);
void MF_Set_Vertices_R4(MFace_ptr f, int n, MVertex_ptr *v);
void MF_Set_Vertices_R3R4(MFace_ptr f, int n, MVertex_ptr *v);
void MF_Set_Vertices_FN(MFace_ptr f, int n, MVertex_ptr *v);
void MF_Set_Vertices_RN(MFace_ptr f, int n, MVertex_ptr *v);
static void 
(*MF_Set_Vertices_jmp[MSTK_MAXREP])(MFace_ptr f, int n, MVertex_ptr *v) = 
{MF_Set_Vertices_F1, MF_Set_Vertices_F4, MF_Set_Vertices_R1, 
 MF_Set_Vertices_R2, MF_Set_Vertices_R4};

void MF_Replace_Vertex_i_F1(MFace_ptr f, int i, MVertex_ptr v);
void MF_Replace_Vertex_i_F4(MFace_ptr f, int i, MVertex_ptr v);
void MF_Replace_Vertex_i_R1(MFace_ptr f, int i, MVertex_ptr v);
void MF_Replace_Vertex_i_R2(MFace_ptr f, int i, MVertex_ptr v);
void MF_Replace_Vertex_i_R4(MFace_ptr f, int i, MVertex_ptr v);
void MF_Replace_Vertex_i_R3R4(MFace_ptr f, int i, MVertex_ptr v);
void MF_Replace_Vertex_i_RN(MFace_ptr f, int i, MVertex_ptr v);
static void 
(*MF_Replace_Vertex_i_jmp[MSTK_MAXREP])(MFace_ptr f, int i, MVertex_ptr v) =
{MF_Replace_Vertex_i_F1, MF_Replace_Vertex_i_F4, MF_Replace_Vertex_i_R1, 
 MF_Replace_Vertex_i_R2, MF_Replace_Vertex_i_R4};

void MF_Replace_Vertex_F1(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
void MF_Replace_Vertex_F4(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
void MF_Replace_Vertex_R1(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
void MF_Replace_Vertex_R2(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
void MF_Replace_Vertex_R4(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
void MF_Replace_Vertex_R3R4(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
void MF_Replace_Vertex_RN(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
static void
(*MF_Replace_Vertex_jmp[MSTK_MAXREP])(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv) =
{MF_Replace_Vertex_F1, MF_Replace_Vertex_F4, MF_Replace_Vertex_R1, 
 MF_Replace_Vertex_R2, MF_Replace_Vertex_R4};


void MF_Insert_Vertex_F1(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
void MF_Insert_Vertex_F4(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
void MF_Insert_Vertex_R1(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
void MF_Insert_Vertex_R2(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
void MF_Insert_Vertex_R4(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
void MF_Insert_Vertex_R3R4(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
void MF_Insert_Vertex_RN(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
static void 
(*MF_Insert_Vertex_jmp[MSTK_MAXREP])(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v) =
{MF_Insert_Vertex_F1, MF_Insert_Vertex_F4, MF_Insert_Vertex_R1, 
 MF_Insert_Vertex_R2, MF_Insert_Vertex_R4};

void MF_Insert_Vertex_i_F1(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_F4(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_R1(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_R2(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_R4(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_RN(MFace_ptr f, MVertex_ptr nuv, int i);
static void 
(*MF_Insert_Vertex_i_jmp[MSTK_MAXREP])(MFace_ptr f, MVertex_ptr nuv, int i) =
{MF_Insert_Vertex_i_F1, MF_Insert_Vertex_i_F4, MF_Insert_Vertex_i_R1, 
 MF_Insert_Vertex_i_R2, MF_Insert_Vertex_i_R4};

int MF_Rev_EdgeDir_FN(MFace_ptr f, MEdge_ptr e);
int MF_Rev_EdgeDir_F1(MFace_ptr f, MEdge_ptr e);
int MF_Rev_EdgeDir_F4(MFace_ptr f, MEdge_ptr e);
int MF_Rev_EdgeDir_R1(MFace_ptr f, MEdge_ptr e);
int MF_Rev_EdgeDir_R2(MFace_ptr f, MEdge_ptr e);
int MF_Rev_EdgeDir_R4(MFace_ptr f, MEdge_ptr e);
int MF_Rev_EdgeDir_RN(MFace_ptr f, MEdge_ptr e);
static int (*MF_Rev_EdgeDir_jmp[MSTK_MAXREP])(MFace_ptr f, MEdge_ptr e) =
{MF_Rev_EdgeDir_F1, MF_Rev_EdgeDir_F4, MF_Rev_EdgeDir_R1, MF_Rev_EdgeDir_R2, 
 MF_Rev_EdgeDir_R4};

int MF_Rev_EdgeDir_i_FN(MFace_ptr f, int i);
int MF_Rev_EdgeDir_i_F1(MFace_ptr f, int i);
int MF_Rev_EdgeDir_i_F4(MFace_ptr f, int i);
int MF_Rev_EdgeDir_i_R1(MFace_ptr f, int i);
int MF_Rev_EdgeDir_i_R2(MFace_ptr f, int i);
int MF_Rev_EdgeDir_i_R4(MFace_ptr f, int i);
int MF_Rev_EdgeDir_i_RN(MFace_ptr f, int i);
static int (*MF_Rev_EdgeDir_i_jmp[MSTK_MAXREP])(MFace_ptr f, int i) =
{MF_Rev_EdgeDir_i_F1, MF_Rev_EdgeDir_i_F4, MF_Rev_EdgeDir_i_R1, 
 MF_Rev_EdgeDir_i_R2, MF_Rev_EdgeDir_i_R4};


int MFs_AreSame_F1(MFace_ptr f1, MFace_ptr f2);
int MFs_AreSame_F4(MFace_ptr f1, MFace_ptr f2);
int MFs_AreSame_R1(MFace_ptr f1, MFace_ptr f2);
int MFs_AreSame_R2(MFace_ptr f1, MFace_ptr f2);
int MFs_AreSame_R4(MFace_ptr f1, MFace_ptr f2);
int MFs_AreSame_R1R2(MFace_ptr f1, MFace_ptr f2);
static int 
(*MFs_AreSame_jmp[MSTK_MAXREP])(MFace_ptr f1, MFace_ptr f2) =
{MFs_AreSame_F1, MFs_AreSame_F4, MFs_AreSame_R1, MFs_AreSame_R2, MFs_AreSame_R4};

int MF_Num_Vertices_F1(MFace_ptr f);
int MF_Num_Vertices_F4(MFace_ptr f);
int MF_Num_Vertices_R1(MFace_ptr f);
int MF_Num_Vertices_R2(MFace_ptr f);
int MF_Num_Vertices_R4(MFace_ptr f);
int MF_Num_Vertices_RN(MFace_ptr f);
static int (*MF_Num_Vertices_jmp[MSTK_MAXREP])(MFace_ptr f) = 
{MF_Num_Vertices_F1, MF_Num_Vertices_F4, MF_Num_Vertices_R1, 
 MF_Num_Vertices_R2, MF_Num_Vertices_R4};

int MF_Num_Edges_F1(MFace_ptr f);
int MF_Num_Edges_F4(MFace_ptr f);
int MF_Num_Edges_R1(MFace_ptr f);
int MF_Num_Edges_R2(MFace_ptr f);
int MF_Num_Edges_R4(MFace_ptr f);
int MF_Num_Edges_RN(MFace_ptr f);
static int (*MF_Num_Edges_jmp[MSTK_MAXREP])(MFace_ptr f) =
{MF_Num_Edges_F1, MF_Num_Edges_F4, MF_Num_Edges_R1, MF_Num_Edges_R2, 
 MF_Num_Edges_R4};

List_ptr MF_Vertices_FN(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_F1(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_F4(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_R1(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_R2(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_R4(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_RN(MFace_ptr f, int dir, MVertex_ptr v);
static List_ptr (*MF_Vertices_jmp[MSTK_MAXREP])(MFace_ptr f, int dir, MVertex_ptr v) =
{MF_Vertices_F1, MF_Vertices_F4, MF_Vertices_R1, MF_Vertices_R2, 
 MF_Vertices_R4};

List_ptr MF_Edges_FN(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_F1(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_F4(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_R1(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_R2(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_R4(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_RN(MFace_ptr f, int dir, MVertex_ptr v);
static List_ptr (*MF_Edges_jmp[MSTK_MAXREP])(MFace_ptr f, int dir, MVertex_ptr v) =
{MF_Edges_F1, MF_Edges_F4, MF_Edges_R1, MF_Edges_R2, MF_Edges_R4};

int MF_EdgeDir_FN(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_F1(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_F4(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_R1(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_R2(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_R4(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_RN(MFace_ptr f, MEdge_ptr e);
static int (*MF_EdgeDir_jmp[MSTK_MAXREP])(MFace_ptr f, MEdge_ptr e) =
{MF_EdgeDir_F1, MF_EdgeDir_F4, MF_EdgeDir_R1, MF_EdgeDir_R2, MF_EdgeDir_R4};

int MF_EdgeDir_i_FN(MFace_ptr f, int i);
int MF_EdgeDir_i_F1(MFace_ptr f, int i);
int MF_EdgeDir_i_F4(MFace_ptr f, int i);
int MF_EdgeDir_i_R1(MFace_ptr f, int i);
int MF_EdgeDir_i_R2(MFace_ptr f, int i);
int MF_EdgeDir_i_R4(MFace_ptr f, int i);
int MF_EdgeDir_i_RN(MFace_ptr f, int i);
static int (*MF_EdgeDir_i_jmp[MSTK_MAXREP])(MFace_ptr f, int i) =
{MF_EdgeDir_i_F1, MF_EdgeDir_i_F4, MF_EdgeDir_i_R1, MF_EdgeDir_i_R2,
 MF_EdgeDir_i_R4};

int MF_UsesEdge_FN(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_F1(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_F4(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_R1(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_R2(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_R4(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_RN(MFace_ptr f, MEdge_ptr e);
static int (*MF_UsesEdge_jmp[MSTK_MAXREP])(MFace_ptr f, MEdge_ptr e) =
{MF_UsesEdge_F1, MF_UsesEdge_F4, MF_UsesEdge_R1, MF_UsesEdge_R2, 
 MF_UsesEdge_R4};

int MF_UsesVertex_FN(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_F1(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_F4(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_R1(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_R2(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_R4(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_RN(MFace_ptr f, MVertex_ptr v);
static int (*MF_UsesVertex_jmp[MSTK_MAXREP])(MFace_ptr f, MVertex_ptr v) =
{MF_UsesVertex_F1, MF_UsesVertex_F4, MF_UsesVertex_R1, MF_UsesVertex_R2, 
 MF_UsesVertex_R4};

List_ptr MF_Regions_F1F3(MFace_ptr f);
List_ptr MF_Regions_F1(MFace_ptr f);
List_ptr MF_Regions_F4(MFace_ptr f);
List_ptr MF_Regions_R1(MFace_ptr f);
List_ptr MF_Regions_R2(MFace_ptr f);
List_ptr MF_Regions_R4(MFace_ptr f);
List_ptr MF_Regions_R1R2(MFace_ptr f);
List_ptr MF_Regions_R3R4(MFace_ptr f);
static List_ptr (*MF_Regions_jmp[MSTK_MAXREP])(MFace_ptr f) =
{MF_Regions_F1F3, MF_Regions_F4, MF_Regions_R1, MF_Regions_R2, MF_Regions_R4};

MRegion_ptr MF_Region_F1F3(MFace_ptr f, int side);
MRegion_ptr MF_Region_F1(MFace_ptr f, int side);
MRegion_ptr MF_Region_F4(MFace_ptr f, int side);
MRegion_ptr MF_Region_R1(MFace_ptr f, int side);
MRegion_ptr MF_Region_R2(MFace_ptr f, int side);
MRegion_ptr MF_Region_R4(MFace_ptr f, int side);
MRegion_ptr MF_Region_R1R2(MFace_ptr f, int side);
MRegion_ptr MF_Region_R3R4(MFace_ptr f, int side);
static MRegion_ptr (*MF_Region_jmp[MSTK_MAXREP])(MFace_ptr f, int side) =
{MF_Region_F1F3, MF_Region_F4, MF_Region_R1, MF_Region_R2, MF_Region_R4};




#else




  int  MF_Set_GInfo_Auto_FN(MFace_ptr f);
  int  MF_Set_GInfo_Auto_RN(MFace_ptr f);
  static int (*MF_Set_GInfo_Auto_jmp[MSTK_MAXREP])(MFace_ptr f) =
  {MF_Set_GInfo_Auto_FN, MF_Set_GInfo_Auto_FN, MF_Set_GInfo_Auto_RN,
   MF_Set_GInfo_Auto_RN, MF_Set_GInfo_Auto_RN};

void MF_Set_Edges_FN(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_R1(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_R2(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_R4(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_RN(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
static void 
(*MF_Set_Edges_jmp[MSTK_MAXREP])(MFace_ptr f, int n, MEdge_ptr *e, int *dir) = 
  {MF_Set_Edges_FN, MF_Set_Edges_FN, MF_Set_Edges_R1, MF_Set_Edges_R2, 
   MF_Set_Edges_R4};

void MF_Set_Vertices_F1(MFace_ptr f, int n, MVertex_ptr *v);
void MF_Set_Vertices_F4(MFace_ptr f, int n, MVertex_ptr *v);
void MF_Set_Vertices_FN(MFace_ptr f, int n, MVertex_ptr *v);
void MF_Set_Vertices_RN(MFace_ptr f, int n, MVertex_ptr *v);
static void 
(*MF_Set_Vertices_jmp[MSTK_MAXREP])(MFace_ptr f, int n, MVertex_ptr *v) = 
{MF_Set_Vertices_FN, MF_Set_Vertices_FN, MF_Set_Vertices_RN, 
 MF_Set_Vertices_RN, MF_Set_Vertices_RN};

void MF_Replace_Edges_FN(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_F1(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_F4(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_R1(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_R2(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_R4(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_RN(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
static void (*MF_Replace_Edges_jmp[MSTK_MAXREP])
     (MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges)=
{MF_Replace_Edges_FN, MF_Replace_Edges_FN, MF_Replace_Edges_R1, 
 MF_Replace_Edges_R2, MF_Replace_Edges_R4};

void MF_Replace_Edges_i_FN(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_i_R1(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_i_R2(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_i_R4(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_i_RN(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
static void 
(*MF_Replace_Edges_i_jmp[MSTK_MAXREP])(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges) = 
{MF_Replace_Edges_i_FN, MF_Replace_Edges_i_FN, MF_Replace_Edges_i_R1, 
 MF_Replace_Edges_i_R2, MF_Replace_Edges_i_R4};

void MF_Replace_Vertex_i_F1(MFace_ptr f, int i, MVertex_ptr v);
void MF_Replace_Vertex_i_F4(MFace_ptr f, int i, MVertex_ptr v);
void MF_Replace_Vertex_i_RN(MFace_ptr f, int i, MVertex_ptr v);
static void 
(*MF_Replace_Vertex_i_jmp[MSTK_MAXREP])(MFace_ptr f, int i, MVertex_ptr v) =
{MF_Replace_Vertex_i_F1, MF_Replace_Vertex_i_F4, MF_Replace_Vertex_i_RN, 
 MF_Replace_Vertex_i_RN, MF_Replace_Vertex_i_RN};

void MF_Replace_Vertex_F1(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
void MF_Replace_Vertex_F4(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
void MF_Replace_Vertex_RN(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
static void
(*MF_Replace_Vertex_jmp[MSTK_MAXREP])(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv) =
{MF_Replace_Vertex_F1, MF_Replace_Vertex_F4, MF_Replace_Vertex_RN, 
 MF_Replace_Vertex_RN, MF_Replace_Vertex_RN};

void MF_Insert_Vertex_F1(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
void MF_Insert_Vertex_F4(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
void MF_Insert_Vertex_RN(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
static void 
(*MF_Insert_Vertex_jmp[MSTK_MAXREP])(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v) =
{MF_Insert_Vertex_F1, MF_Insert_Vertex_F4, MF_Insert_Vertex_RN, 
 MF_Insert_Vertex_RN, MF_Insert_Vertex_RN};

void MF_Insert_Vertex_i_F1(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_F4(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_RN(MFace_ptr f, MVertex_ptr nuv, int i);
static void 
(*MF_Insert_Vertex_i_jmp[MSTK_MAXREP])(MFace_ptr f, MVertex_ptr nuv, int i) =
{MF_Insert_Vertex_i_F1, MF_Insert_Vertex_i_F4, MF_Insert_Vertex_i_RN, 
 MF_Insert_Vertex_i_RN, MF_Insert_Vertex_i_RN};

int MF_Rev_EdgeDir_FN(MFace_ptr f, MEdge_ptr e);
int MF_Rev_EdgeDir_RN(MFace_ptr f, MEdge_ptr e);
static int (*MF_Rev_EdgeDir_jmp[MSTK_MAXREP])(MFace_ptr f, MEdge_ptr e) =
{MF_Rev_EdgeDir_FN, MF_Rev_EdgeDir_FN, MF_Rev_EdgeDir_RN, MF_Rev_EdgeDir_RN, MF_Rev_EdgeDir_RN};

int MF_Rev_EdgeDir_i_FN(MFace_ptr f, int i);
int MF_Rev_EdgeDir_i_RN(MFace_ptr f, int i);
static int (*MF_Rev_EdgeDir_i_jmp[MSTK_MAXREP])(MFace_ptr f, int i) =
{MF_Rev_EdgeDir_i_FN, MF_Rev_EdgeDir_i_FN, MF_Rev_EdgeDir_i_RN, MF_Rev_EdgeDir_i_RN,
 MF_Rev_EdgeDir_i_RN};


int MFs_AreSame_F1(MFace_ptr f1, MFace_ptr f2);
int MFs_AreSame_F4(MFace_ptr f1, MFace_ptr f2);
int MFs_AreSame_R1(MFace_ptr f1, MFace_ptr f2);
int MFs_AreSame_R2(MFace_ptr f1, MFace_ptr f2);
int MFs_AreSame_R4(MFace_ptr f1, MFace_ptr f2);
int MFs_AreSame_R1R2(MFace_ptr f1, MFace_ptr f2);
static int 
(*MFs_AreSame_jmp[MSTK_MAXREP])(MFace_ptr f1, MFace_ptr f2) =
{MFs_AreSame_F1, MFs_AreSame_F4, MFs_AreSame_R1R2, MFs_AreSame_R1R2, 
 MFs_AreSame_R4};

int MF_Num_Vertices_F1(MFace_ptr f);
int MF_Num_Vertices_F4(MFace_ptr f);
int MF_Num_Vertices_RN(MFace_ptr f);
static int (*MF_Num_Vertices_jmp[MSTK_MAXREP])(MFace_ptr f) = 
{MF_Num_Vertices_F1, MF_Num_Vertices_F4, MF_Num_Vertices_RN, 
 MF_Num_Vertices_RN, MF_Num_Vertices_RN};

int MF_Num_Edges_F1(MFace_ptr f);
int MF_Num_Edges_F4(MFace_ptr f);
int MF_Num_Edges_RN(MFace_ptr f);
static int (*MF_Num_Edges_jmp[MSTK_MAXREP])(MFace_ptr f) =
{MF_Num_Edges_F1, MF_Num_Edges_F4, MF_Num_Edges_RN, MF_Num_Edges_RN, 
 MF_Num_Edges_RN};

List_ptr MF_Vertices_FN(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_RN(MFace_ptr f, int dir, MVertex_ptr v);
static List_ptr (*MF_Vertices_jmp[MSTK_MAXREP])(MFace_ptr f, int dir, MVertex_ptr v) =
{MF_Vertices_FN, MF_Vertices_FN, MF_Vertices_RN, MF_Vertices_RN, 
 MF_Vertices_RN};

List_ptr MF_Edges_FN(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_RN(MFace_ptr f, int dir, MVertex_ptr v);
static List_ptr (*MF_Edges_jmp[MSTK_MAXREP])(MFace_ptr f, int dir, MVertex_ptr v) =
{MF_Edges_FN, MF_Edges_FN, MF_Edges_RN, MF_Edges_RN, MF_Edges_RN};

int MF_EdgeDir_FN(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_RN(MFace_ptr f, MEdge_ptr e);
static int (*MF_EdgeDir_jmp[MSTK_MAXREP])(MFace_ptr f, MEdge_ptr e) =
{MF_EdgeDir_FN, MF_EdgeDir_FN, MF_EdgeDir_RN, MF_EdgeDir_RN, MF_EdgeDir_RN};

int MF_EdgeDir_i_FN(MFace_ptr f, int i);
int MF_EdgeDir_i_RN(MFace_ptr f, int i);
static int (*MF_EdgeDir_i_jmp[MSTK_MAXREP])(MFace_ptr f, int i) =
{MF_EdgeDir_i_FN, MF_EdgeDir_i_FN, MF_EdgeDir_i_RN, MF_EdgeDir_i_RN,
 MF_EdgeDir_i_RN};

int MF_UsesEdge_FN(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_RN(MFace_ptr f, MEdge_ptr e);
static int (*MF_UsesEdge_jmp[MSTK_MAXREP])(MFace_ptr f, MEdge_ptr e) =
{MF_UsesEdge_FN, MF_UsesEdge_FN, MF_UsesEdge_RN, MF_UsesEdge_RN,
 MF_UsesEdge_RN};

int MF_UsesVertex_FN(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_RN(MFace_ptr f, MVertex_ptr v);
static int (*MF_UsesVertex_jmp[MSTK_MAXREP])(MFace_ptr f, MVertex_ptr v) =
{MF_UsesVertex_FN, MF_UsesVertex_FN, MF_UsesVertex_RN, MF_UsesVertex_RN, 
 MF_UsesVertex_RN};

List_ptr MF_Regions_F1F3(MFace_ptr f);
List_ptr MF_Regions_F4(MFace_ptr f);
List_ptr MF_Regions_R1R2(MFace_ptr f);
List_ptr MF_Regions_R3R4(MFace_ptr f);
static List_ptr (*MF_Regions_jmp[MSTK_MAXREP])(MFace_ptr f) =
{MF_Regions_F1F3, MF_Regions_F4, MF_Regions_R1R2, MF_Regions_R1R2, MF_Regions_R3R4};

MRegion_ptr MF_Region_F1F3(MFace_ptr f, int side);
MRegion_ptr MF_Region_F4(MFace_ptr f, int side);
MRegion_ptr MF_Region_R1R2(MFace_ptr f, int side);
MRegion_ptr MF_Region_R3R4(MFace_ptr f, int side);
static MRegion_ptr (*MF_Region_jmp[MSTK_MAXREP])(MFace_ptr f, int side) =
{MF_Region_F1F3, MF_Region_F4, MF_Region_R1R2, MF_Region_R1R2, MF_Region_R3R4};

#endif


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
