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

void MF_Set_Edges_FN(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_F1(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_F4(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_R1(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_R2(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_R4(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
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
static void 
(*MF_Insert_Vertex_jmp[MSTK_MAXREP])(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v) =
{MF_Insert_Vertex_F1, MF_Insert_Vertex_F4, MF_Insert_Vertex_R1, 
 MF_Insert_Vertex_R2, MF_Insert_Vertex_R4};

void MF_Insert_Vertex_i_F1(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_F4(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_R1(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_R2(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_R4(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_R3R4(MFace_ptr f, MVertex_ptr nuv, int i);
static void 
(*MF_Insert_Vertex_i_jmp[MSTK_MAXREP])(MFace_ptr f, MVertex_ptr nuv, int i) =
{MF_Insert_Vertex_i_F1, MF_Insert_Vertex_i_F4, MF_Insert_Vertex_i_R1, 
 MF_Insert_Vertex_i_R2, MF_Insert_Vertex_i_R4};

int MF_Num_Vertices_F1(MFace_ptr f);
int MF_Num_Vertices_F4(MFace_ptr f);
int MF_Num_Vertices_R1(MFace_ptr f);
int MF_Num_Vertices_R2(MFace_ptr f);
int MF_Num_Vertices_R4(MFace_ptr f);
int MF_Num_Vertices_R1R2(MFace_ptr f);
int MF_Num_Vertices_R3R4(MFace_ptr f);
static int (*MF_Num_Vertices_jmp[MSTK_MAXREP])(MFace_ptr f) = 
{MF_Num_Vertices_F1, MF_Num_Vertices_F4, MF_Num_Vertices_R1, 
 MF_Num_Vertices_R2, MF_Num_Vertices_R4};

int MF_Num_Edges_F1(MFace_ptr f);
int MF_Num_Edges_F4(MFace_ptr f);
int MF_Num_Edges_R1(MFace_ptr f);
int MF_Num_Edges_R2(MFace_ptr f);
int MF_Num_Edges_R4(MFace_ptr f);
int MF_Num_Edges_R1R2(MFace_ptr f);
int MF_Num_Edges_R3R4(MFace_ptr f);
static int (*MF_Num_Edges_jmp[MSTK_MAXREP])(MFace_ptr f) =
{MF_Num_Edges_F1, MF_Num_Edges_F4, MF_Num_Edges_R1, MF_Num_Edges_R2, 
 MF_Num_Edges_R4};

List_ptr MF_Vertices_FN(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_F1(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_F4(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_R1(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_R2(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_R4(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_R1R2(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_R3R4(MFace_ptr f, int dir, MVertex_ptr v);
static List_ptr (*MF_Vertices_jmp[MSTK_MAXREP])(MFace_ptr f, int dir, MVertex_ptr v) =
{MF_Vertices_F1, MF_Vertices_F4, MF_Vertices_R1, MF_Vertices_R2, 
 MF_Vertices_R4};

List_ptr MF_Edges_FN(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_F1(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_F4(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_R1(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_R2(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_R4(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_R1R2(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_R3R4(MFace_ptr f, int dir, MVertex_ptr v);
static List_ptr (*MF_Edges_jmp[MSTK_MAXREP])(MFace_ptr f, int dir, MVertex_ptr v) =
{MF_Edges_F1, MF_Edges_F4, MF_Edges_R1, MF_Edges_R2, MF_Edges_R4};

int MF_EdgeDir_FN(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_F1(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_F4(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_R1(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_R2(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_R4(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_R1R2(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_R3R4(MFace_ptr f, MEdge_ptr e);
static int (*MF_EdgeDir_jmp[MSTK_MAXREP])(MFace_ptr f, MEdge_ptr e) =
{MF_EdgeDir_F1, MF_EdgeDir_F4, MF_EdgeDir_R1, MF_EdgeDir_R2, MF_EdgeDir_R4};

int MF_EdgeDir_i_FN(MFace_ptr f, int i);
int MF_EdgeDir_i_F1(MFace_ptr f, int i);
int MF_EdgeDir_i_F4(MFace_ptr f, int i);
int MF_EdgeDir_i_R1(MFace_ptr f, int i);
int MF_EdgeDir_i_R2(MFace_ptr f, int i);
int MF_EdgeDir_i_R4(MFace_ptr f, int i);
int MF_EdgeDir_i_R1R2(MFace_ptr f, int i);
int MF_EdgeDir_i_R3R4(MFace_ptr f, int i);
static int (*MF_EdgeDir_i_jmp[MSTK_MAXREP])(MFace_ptr f, int i) =
{MF_EdgeDir_i_F1, MF_EdgeDir_i_F4, MF_EdgeDir_i_R1, MF_EdgeDir_i_R2,
 MF_EdgeDir_i_R4};

int MF_UsesEdge_FN(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_F1(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_F4(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_R1(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_R2(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_R4(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_R1R2(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_R3R4(MFace_ptr f, MEdge_ptr e);
static int (*MF_UsesEdge_jmp[MSTK_MAXREP])(MFace_ptr f, MEdge_ptr e) =
{MF_UsesEdge_F1, MF_UsesEdge_F4, MF_UsesEdge_R1, MF_UsesEdge_R2, 
 MF_UsesEdge_R4};

int MF_UsesVertex_FN(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_F1(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_F4(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_R1(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_R2(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_R4(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_R1R2(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_R3R4(MFace_ptr f, MVertex_ptr v);
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

void MF_Set_Edges_FN(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_R1(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_R2(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
void MF_Set_Edges_R4(MFace_ptr f, int n, MEdge_ptr *e, int *dir);
static void 
(*MF_Set_Edges_jmp[MSTK_MAXREP])(MFace_ptr f, int n, MEdge_ptr *e, int *dir) = 
  {MF_Set_Edges_FN, MF_Set_Edges_FN, MF_Set_Edges_R1, MF_Set_Edges_R2, 
   MF_Set_Edges_R4};

void MF_Set_Vertices_F1(MFace_ptr f, int n, MVertex_ptr *v);
void MF_Set_Vertices_F4(MFace_ptr f, int n, MVertex_ptr *v);
void MF_Set_Vertices_R1(MFace_ptr f, int n, MVertex_ptr *v);
void MF_Set_Vertices_R2(MFace_ptr f, int n, MVertex_ptr *v);
void MF_Set_Vertices_R3R4(MFace_ptr f, int n, MVertex_ptr *v);
static void 
(*MF_Set_Vertices_jmp[MSTK_MAXREP])(MFace_ptr f, int n, MVertex_ptr *v) = 
{MF_Set_Vertices_F1, MF_Set_Vertices_F4, MF_Set_Vertices_R1, 
 MF_Set_Vertices_R2, MF_Set_Vertices_R3R4};

void MF_Replace_Edges_FN(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_F1(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_F4(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_R1(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_R2(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_R4(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
static void (*MF_Replace_Edges_jmp[MSTK_MAXREP])
     (MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges)=
{MF_Replace_Edges_FN, MF_Replace_Edges_FN, MF_Replace_Edges_R1, 
 MF_Replace_Edges_R2, MF_Replace_Edges_R4};

void MF_Replace_Edges_i_FN(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_i_R1(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_i_R2(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
void MF_Replace_Edges_i_R4(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges);
static void 
(*MF_Replace_Edges_i_jmp[MSTK_MAXREP])(MFace_ptr f, int nold, int i, int nnu, MEdge_ptr *nuedges) = 
{MF_Replace_Edges_i_FN, MF_Replace_Edges_i_FN, MF_Replace_Edges_i_R1, 
 MF_Replace_Edges_i_R2, MF_Replace_Edges_i_R4};

void MF_Replace_Vertex_i_F1(MFace_ptr f, int i, MVertex_ptr v);
void MF_Replace_Vertex_i_F4(MFace_ptr f, int i, MVertex_ptr v);
void MF_Replace_Vertex_i_R1(MFace_ptr f, int i, MVertex_ptr v);
void MF_Replace_Vertex_i_R2(MFace_ptr f, int i, MVertex_ptr v);
void MF_Replace_Vertex_i_R3R4(MFace_ptr f, int i, MVertex_ptr v);
static void 
(*MF_Replace_Vertex_i_jmp[MSTK_MAXREP])(MFace_ptr f, int i, MVertex_ptr v) =
{MF_Replace_Vertex_i_F1, MF_Replace_Vertex_i_F4, MF_Replace_Vertex_i_R1, 
 MF_Replace_Vertex_i_R2, MF_Replace_Vertex_i_R3R4};

void MF_Replace_Vertex_F1(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
void MF_Replace_Vertex_F4(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
void MF_Replace_Vertex_R1(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
void MF_Replace_Vertex_R2(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
void MF_Replace_Vertex_R3R4(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
static void
(*MF_Replace_Vertex_jmp[MSTK_MAXREP])(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv) =
{MF_Replace_Vertex_F1, MF_Replace_Vertex_F4, MF_Replace_Vertex_R1, 
 MF_Replace_Vertex_R2, MF_Replace_Vertex_R3R4};

void MF_Insert_Vertex_F1(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
void MF_Insert_Vertex_F4(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
void MF_Insert_Vertex_R1(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
void MF_Insert_Vertex_R2(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
void MF_Insert_Vertex_R3R4(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v);
static void 
(*MF_Insert_Vertex_jmp[MSTK_MAXREP])(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v) =
{MF_Insert_Vertex_F1, MF_Insert_Vertex_F4, MF_Insert_Vertex_R1, 
 MF_Insert_Vertex_R2, MF_Insert_Vertex_R3R4};

void MF_Insert_Vertex_i_F1(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_F4(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_R1(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_R2(MFace_ptr f, MVertex_ptr nuv, int i);
void MF_Insert_Vertex_i_R3R4(MFace_ptr f, MVertex_ptr nuv, int i);
static void 
(*MF_Insert_Vertex_i_jmp[MSTK_MAXREP])(MFace_ptr f, MVertex_ptr nuv, int i) =
{MF_Insert_Vertex_i_F1, MF_Insert_Vertex_i_F4, MF_Insert_Vertex_i_R1, 
 MF_Insert_Vertex_i_R2, MF_Insert_Vertex_i_R3R4};

int MF_Num_Vertices_F1(MFace_ptr f);
int MF_Num_Vertices_F4(MFace_ptr f);
int MF_Num_Vertices_R1R2(MFace_ptr f);
int MF_Num_Vertices_R3R4(MFace_ptr f);
static int (*MF_Num_Vertices_jmp[MSTK_MAXREP])(MFace_ptr f) = 
{MF_Num_Vertices_F1, MF_Num_Vertices_F4, MF_Num_Vertices_R1R2, 
 MF_Num_Vertices_R1R2, MF_Num_Vertices_R3R4};

int MF_Num_Edges_F1(MFace_ptr f);
int MF_Num_Edges_F4(MFace_ptr f);
int MF_Num_Edges_R1R2(MFace_ptr f);
int MF_Num_Edges_R3R4(MFace_ptr f);
static int (*MF_Num_Edges_jmp[MSTK_MAXREP])(MFace_ptr f) =
{MF_Num_Edges_F1, MF_Num_Edges_F4, MF_Num_Edges_R1R2, MF_Num_Edges_R1R2, 
 MF_Num_Edges_R3R4};

List_ptr MF_Vertices_FN(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_R1R2(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Vertices_R3R4(MFace_ptr f, int dir, MVertex_ptr v);
static List_ptr (*MF_Vertices_jmp[MSTK_MAXREP])(MFace_ptr f, int dir, MVertex_ptr v) =
{MF_Vertices_FN, MF_Vertices_FN, MF_Vertices_R1R2, MF_Vertices_R1R2, 
 MF_Vertices_R3R4};

List_ptr MF_Edges_FN(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_R1R2(MFace_ptr f, int dir, MVertex_ptr v);
List_ptr MF_Edges_R3R4(MFace_ptr f, int dir, MVertex_ptr v);
static List_ptr (*MF_Edges_jmp[MSTK_MAXREP])(MFace_ptr f, int dir, MVertex_ptr v) =
{MF_Edges_FN, MF_Edges_FN, MF_Edges_R1R2, MF_Edges_R1R2, MF_Edges_R3R4};

int MF_EdgeDir_FN(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_R1R2(MFace_ptr f, MEdge_ptr e);
int MF_EdgeDir_R3R4(MFace_ptr f, MEdge_ptr e);
static int (*MF_EdgeDir_jmp[MSTK_MAXREP])(MFace_ptr f, MEdge_ptr e) =
{MF_EdgeDir_FN, MF_EdgeDir_FN, MF_EdgeDir_R1R2, MF_EdgeDir_R1R2, MF_EdgeDir_R3R4};

int MF_EdgeDir_i_FN(MFace_ptr f, int i);
int MF_EdgeDir_i_R1R2(MFace_ptr f, int i);
int MF_EdgeDir_i_R3R4(MFace_ptr f, int i);
static int (*MF_EdgeDir_i_jmp[MSTK_MAXREP])(MFace_ptr f, int i) =
{MF_EdgeDir_i_FN, MF_EdgeDir_i_FN, MF_EdgeDir_i_R1R2, MF_EdgeDir_i_R1R2,
 MF_EdgeDir_i_R3R4};

int MF_UsesEdge_FN(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_R1R2(MFace_ptr f, MEdge_ptr e);
int MF_UsesEdge_R3R4(MFace_ptr f, MEdge_ptr e);
static int (*MF_UsesEdge_jmp[MSTK_MAXREP])(MFace_ptr f, MEdge_ptr e) =
{MF_UsesEdge_FN, MF_UsesEdge_FN, MF_UsesEdge_R1R2, MF_UsesEdge_R1R2,
 MF_UsesEdge_R3R4};

int MF_UsesVertex_FN(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_R1R2(MFace_ptr f, MVertex_ptr v);
int MF_UsesVertex_R3R4(MFace_ptr f, MVertex_ptr v);
static int (*MF_UsesVertex_jmp[MSTK_MAXREP])(MFace_ptr f, MVertex_ptr v) =
{MF_UsesVertex_FN, MF_UsesVertex_FN, MF_UsesVertex_R1R2, MF_UsesVertex_R1R2, 
 MF_UsesVertex_R3R4};

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
void MF_Set_RepType_R4(MFace_ptr f);
static void (*MF_Set_RepType_jmp[MSTK_MAXREP])(MFace_ptr f) =
{MF_Set_RepType_F1, MF_Set_RepType_F4, MF_Dummy1, MF_Dummy1,
 MF_Set_RepType_R4};

void MF_Delete_F1(MFace_ptr f, int keep);
void MF_Delete_F4(MFace_ptr f, int keep);
void MF_Delete_R4(MFace_ptr f, int keep);
static void (*MF_Delete_jmp[MSTK_MAXREP])(MFace_ptr f, int keep) = 
{MF_Delete_F1, MF_Delete_F4, MF_Dummy2a, MF_Dummy2a, MF_Delete_R4};

void MF_Restore_F1(MFace_ptr f);
void MF_Restore_F4(MFace_ptr f);
void MF_Restore_R4(MFace_ptr f);
static void (*MF_Restore_jmp[MSTK_MAXREP])(MFace_ptr f) = 
{MF_Restore_F1, MF_Restore_F4, MF_Dummy1, MF_Dummy1, MF_Restore_R4};

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

  
#ifdef __cplusplus
}
#endif
