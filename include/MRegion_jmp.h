#define _H_MRegion_Private

#include "MSTK.h"
#include "MRegion.h"

#ifdef __cplusplus
extern "C" {
#endif


#ifdef DEBUG

void MR_Set_RepType_F1(MRegion_ptr r);
void MR_Set_RepType_F4(MRegion_ptr r);
void MR_Set_RepType_R1(MRegion_ptr r);
void MR_Set_RepType_R2(MRegion_ptr r);
void MR_Set_RepType_R4(MRegion_ptr r);
void MR_Set_RepType_FNR3R4(MRegion_ptr r);
static void (*MR_Set_RepType_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_Set_RepType_F1, MR_Set_RepType_F4, MR_Set_RepType_R1, 
 MR_Set_RepType_R2, MR_Set_RepType_R4};

void MR_Delete_F1(MRegion_ptr r);
void MR_Delete_F4(MRegion_ptr r);
void MR_Delete_R1(MRegion_ptr r);
void MR_Delete_R2(MRegion_ptr r);
void MR_Delete_R4(MRegion_ptr r);
void MR_Delete_F1F3R3R4(MRegion_ptr r);
static void (*MR_Delete_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_Delete_F1, MR_Delete_F4, MR_Delete_R1, MR_Delete_R2, MR_Delete_R4};

void MR_Set_Faces_F1(MRegion_ptr r, int nf, MFace_ptr *mfaces, int *dirs);
void MR_Set_Faces_F4(MRegion_ptr r, int nf, MFace_ptr *mfaces, int *dirs);
void MR_Set_Faces_R1(MRegion_ptr r, int nf, MFace_ptr *mfaces, int *dirs);
void MR_Set_Faces_R2(MRegion_ptr r, int nf, MFace_ptr *mfaces, int *dirs);
void MR_Set_Faces_R4(MRegion_ptr r, int nf, MFace_ptr *mfaces, int *dirs);
void MR_Set_Faces_F1F3R3R4(MRegion_ptr r, int nf, MFace_ptr *mfaces, int *dirs);
static void 
(*MR_Set_Faces_jmp[MSTK_MAXREP])(MRegion_ptr r,int nf,MFace_ptr *mfaces,int *dirs) =
{MR_Set_Faces_F1, MR_Set_Faces_F4, MR_Set_Faces_R1, MR_Set_Faces_R2, 
 MR_Set_Faces_R4};

void MR_Set_Vertices_F1(MRegion_ptr r, int nv, MVertex_ptr *mvertices);
void MR_Set_Vertices_F4(MRegion_ptr r, int nv, MVertex_ptr *mvertices);
void MR_Set_Vertices_R1(MRegion_ptr r, int nv, MVertex_ptr *mvertices);
void MR_Set_Vertices_R2(MRegion_ptr r, int nv, MVertex_ptr *mvertices);
void MR_Set_Vertices_R4(MRegion_ptr r, int nv, MVertex_ptr *mvertices);
void MR_Set_Vertices_FNR3R4(MRegion_ptr r, int nv, MVertex_ptr *mvertices);
static void 
(*MR_Set_Vertices_jmp[MSTK_MAXREP])(MRegion_ptr r, int nv, MVertex_ptr *mvertices) =
{MR_Set_Vertices_F1, MR_Set_Vertices_F4, MR_Set_Vertices_R1, 
 MR_Set_Vertices_R2, MR_Set_Vertices_R4};

Set_ptr MR_Vertices_F1(MRegion_ptr r);
Set_ptr MR_Vertices_F4(MRegion_ptr r);
Set_ptr MR_Vertices_R1(MRegion_ptr r);
Set_ptr MR_Vertices_R2(MRegion_ptr r);
Set_ptr MR_Vertices_R4(MRegion_ptr r);
Set_ptr MR_Vertices_FNR3R4(MRegion_ptr r);
static Set_ptr (*MR_Vertices_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_Vertices_F1, MR_Vertices_F4, MR_Vertices_R1, MR_Vertices_R2, 
MR_Vertices_R4};

Set_ptr MR_Edges_F1(MRegion_ptr r);
Set_ptr MR_Edges_F4(MRegion_ptr r);
Set_ptr MR_Edges_R1(MRegion_ptr r);
Set_ptr MR_Edges_R2(MRegion_ptr r);
Set_ptr MR_Edges_R4(MRegion_ptr r);
Set_ptr MR_Edges_FNR3R4(MRegion_ptr r);
static Set_ptr (*MR_Edges_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_Edges_F1, MR_Edges_F4, MR_Edges_R1, MR_Edges_R2, MR_Edges_R4};

Set_ptr MR_Faces_F1(MRegion_ptr r);
Set_ptr MR_Faces_F4(MRegion_ptr r);
Set_ptr MR_Faces_R1(MRegion_ptr r);
Set_ptr MR_Faces_R2(MRegion_ptr r);
Set_ptr MR_Faces_R4(MRegion_ptr r);
Set_ptr MR_Faces_FNR3R4(MRegion_ptr r);
static Set_ptr (*MR_Faces_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_Faces_F1, MR_Faces_F4, MR_Faces_R1, MR_Faces_R2, MR_Faces_R4};

Set_ptr MR_AdjRegions_F1(MRegion_ptr r);
Set_ptr MR_AdjRegions_F4(MRegion_ptr r);
Set_ptr MR_AdjRegions_R1(MRegion_ptr r);
Set_ptr MR_AdjRegions_R2(MRegion_ptr r);
Set_ptr MR_AdjRegions_R4(MRegion_ptr r);
Set_ptr MR_AdjRegions_FNR3R4(MRegion_ptr r);
static Set_ptr (*MR_AdjRegions_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_AdjRegions_F1, MR_AdjRegions_F4, MR_AdjRegions_R1, MR_AdjRegions_R2, 
 MR_AdjRegions_R4};

int MR_FaceDir_F1(MRegion_ptr r, MFace_ptr f);
int MR_FaceDir_F4(MRegion_ptr r, MFace_ptr f);
int MR_FaceDir_R1(MRegion_ptr r, MFace_ptr f);
int MR_FaceDir_R2(MRegion_ptr r, MFace_ptr f);
int MR_FaceDir_R4(MRegion_ptr r, MFace_ptr f);
int MR_FaceDir_FNR3R4(MRegion_ptr r, MFace_ptr f);
static int (*MR_FaceDir_jmp[MSTK_MAXREP])(MRegion_ptr r, MFace_ptr f) =
{MR_FaceDir_F1, MR_FaceDir_F4, MR_FaceDir_R1, MR_FaceDir_R2, 
 MR_FaceDir_R4};

int MR_FaceDir_i_F1(MRegion_ptr r, int i);
int MR_FaceDir_i_F4(MRegion_ptr r, int i);
int MR_FaceDir_i_R1(MRegion_ptr r, int i);
int MR_FaceDir_i_R2(MRegion_ptr r, int i);
int MR_FaceDir_i_R4(MRegion_ptr r, int i);
int MR_FaceDir_i_FNR3R4(MRegion_ptr r, int i);
static int (*MR_FaceDir_i_jmp[MSTK_MAXREP])(MRegion_ptr r, int i) =
{MR_FaceDir_i_F1, MR_FaceDir_i_F4, MR_FaceDir_i_R1, MR_FaceDir_i_R2, 
 MR_FaceDir_i_R4};

void MR_Replace_Face_F1(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf, int nudir);
void MR_Replace_Face_F4(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf, int nudir);
void MR_Replace_Face_R1(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf, int nudir);
void MR_Replace_Face_R2(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf, int nudir);
void MR_Replace_Face_R4(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf, int nudir);
void MR_Replace_Face_FNR3R4(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf, int nudir);
static void (*MR_Replace_Face_jmp[MSTK_MAXREP])(MRegion_ptr r, MFace_ptr f, 
				      MFace_ptr nuf, int nudir) =
{MR_Replace_Face_F1, MR_Replace_Face_F4, MR_Replace_Face_R1, 
 MR_Replace_Face_R2, MR_Replace_Face_R4};

void MR_Replace_Face_i_F1(MRegion_ptr r, int i,MFace_ptr nuf,int nudir);
void MR_Replace_Face_i_F4(MRegion_ptr r, int i,MFace_ptr nuf,int nudir);
void MR_Replace_Face_i_R1(MRegion_ptr r, int i,MFace_ptr nuf,int nudir);
void MR_Replace_Face_i_R2(MRegion_ptr r, int i,MFace_ptr nuf,int nudir);
void MR_Replace_Face_i_R4(MRegion_ptr r, int i,MFace_ptr nuf,int nudir);
void MR_Replace_Face_i_FNR3R4(MRegion_ptr r, int i,MFace_ptr nuf,int nudir);
static void 
(*MR_Replace_Face_i_jmp[MSTK_MAXREP])(MRegion_ptr r, int i,MFace_ptr nuf,int nudir)=
{MR_Replace_Face_i_F1, MR_Replace_Face_i_F4, MR_Replace_Face_i_R1, 
 MR_Replace_Face_i_R2, MR_Replace_Face_i_R4};

void MR_Replace_Vertex_F1(MRegion_ptr r, MVertex_ptr v,MVertex_ptr nuv);
void MR_Replace_Vertex_F4(MRegion_ptr r, MVertex_ptr v,MVertex_ptr nuv);
void MR_Replace_Vertex_R1(MRegion_ptr r, MVertex_ptr v,MVertex_ptr nuv);
void MR_Replace_Vertex_R2(MRegion_ptr r, MVertex_ptr v,MVertex_ptr nuv);
void MR_Replace_Vertex_R4(MRegion_ptr r, MVertex_ptr v,MVertex_ptr nuv);
void MR_Replace_Vertex_FNR3R4(MRegion_ptr r, MVertex_ptr v,MVertex_ptr nuv);
static void 
(*MR_Replace_Vertex_jmp[MSTK_MAXREP])(MRegion_ptr r, MVertex_ptr v,MVertex_ptr nuv)=
{MR_Replace_Vertex_F1, MR_Replace_Vertex_F4, MR_Replace_Vertex_R1, 
 MR_Replace_Vertex_R2, MR_Replace_Vertex_R4};

void MR_Replace_Vertex_i_F1(MRegion_ptr r, int i, MVertex_ptr nuv);
void MR_Replace_Vertex_i_F4(MRegion_ptr r, int i, MVertex_ptr nuv);
void MR_Replace_Vertex_i_R1(MRegion_ptr r, int i, MVertex_ptr nuv);
void MR_Replace_Vertex_i_R2(MRegion_ptr r, int i, MVertex_ptr nuv);
void MR_Replace_Vertex_i_R4(MRegion_ptr r, int i, MVertex_ptr nuv);
void MR_Replace_Vertex_i_FNR3R4(MRegion_ptr r, int i, MVertex_ptr nuv);
static void 
(*MR_Replace_Vertex_i_jmp[MSTK_MAXREP])(MRegion_ptr r, int i, MVertex_ptr nuv) =
{MR_Replace_Vertex_i_F1, MR_Replace_Vertex_i_F4, MR_Replace_Vertex_i_R1, 
 MR_Replace_Vertex_i_R2, MR_Replace_Vertex_i_R4};

void MR_Add_AdjRegion_F1(MRegion_ptr r, int facenum, MRegion_ptr aregion);
void MR_Add_AdjRegion_F4(MRegion_ptr r, int facenum, MRegion_ptr aregion);
void MR_Add_AdjRegion_R1(MRegion_ptr r, int facenum, MRegion_ptr aregion);
void MR_Add_AdjRegion_R2(MRegion_ptr r, int facenum, MRegion_ptr aregion);
void MR_Add_AdjRegion_R4(MRegion_ptr r, int facenum, MRegion_ptr aregion);
void MR_Add_AdjRegion_FNR3R4(MRegion_ptr r, int facenum, MRegion_ptr aregion);
static void (*MR_Add_AdjRegion_jmp[MSTK_MAXREP])(MRegion_ptr r, int facenum, MRegion_ptr ar) =
{MR_Add_AdjRegion_F1, MR_Add_AdjRegion_F4, MR_Add_AdjRegion_R1, 
 MR_Add_AdjRegion_R2, MR_Add_AdjRegion_R4};

void MR_Rem_AdjRegion_F1(MRegion_ptr r, MRegion_ptr ar);
void MR_Rem_AdjRegion_F4(MRegion_ptr r, MRegion_ptr ar);
void MR_Rem_AdjRegion_R1(MRegion_ptr r, MRegion_ptr ar);
void MR_Rem_AdjRegion_R2(MRegion_ptr r, MRegion_ptr ar);
void MR_Rem_AdjRegion_R4(MRegion_ptr r, MRegion_ptr ar);
void MR_Rem_AdjRegion_FNR3R4(MRegion_ptr r, MRegion_ptr ar);
static void (*MR_Rem_AdjRegion_jmp[MSTK_MAXREP])(MRegion_ptr r, MRegion_ptr ar) =
{MR_Rem_AdjRegion_F1, MR_Rem_AdjRegion_F4, MR_Rem_AdjRegion_R1, 
 MR_Rem_AdjRegion_R2, MR_Rem_AdjRegion_R4};

int MR_Num_Faces_F1(MRegion_ptr r);
int MR_Num_Faces_F4(MRegion_ptr r);
int MR_Num_Faces_R1(MRegion_ptr r);
int MR_Num_Faces_R2(MRegion_ptr r);
int MR_Num_Faces_R4(MRegion_ptr r);
int MR_Num_Faces_FNR3R4(MRegion_ptr r);
static int (*MR_Num_Faces_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_Num_Faces_FNR3R4, MR_Num_Faces_FNR3R4, MR_Num_Faces_R1, MR_Num_Faces_R2,
 MR_Num_Faces_FNR3R4};

int MR_UsesFace_F1(MRegion_ptr r, MFace_ptr f);
int MR_UsesFace_F4(MRegion_ptr r, MFace_ptr f);
int MR_UsesFace_R1(MRegion_ptr r, MFace_ptr f);
int MR_UsesFace_R2(MRegion_ptr r, MFace_ptr f);
int MR_UsesFace_R4(MRegion_ptr r, MFace_ptr f);
int MR_UsesFace_FNR3R4(MRegion_ptr r, MFace_ptr f);
static int (*MR_UsesFace_jmp[MSTK_MAXREP])(MRegion_ptr r, MFace_ptr f) =
{MR_UsesFace_F1, MR_UsesFace_F4, MR_UsesFace_R1, MR_UsesFace_R2, 
 MR_UsesFace_R4}; 

int MR_UsesEdge_F1(MRegion_ptr r, MEdge_ptr e);
int MR_UsesEdge_F4(MRegion_ptr r, MEdge_ptr e);
int MR_UsesEdge_R1(MRegion_ptr r, MEdge_ptr e);
int MR_UsesEdge_R2(MRegion_ptr r, MEdge_ptr e);
int MR_UsesEdge_R4(MRegion_ptr r, MEdge_ptr e);
int MR_UsesEdge_FNR3R4(MRegion_ptr r, MEdge_ptr e);
static int (*MR_UsesEdge_jmp[MSTK_MAXREP])(MRegion_ptr r, MEdge_ptr e) =
{MR_UsesEdge_F1, MR_UsesEdge_F4, MR_UsesEdge_R1, MR_UsesEdge_R2, 
 MR_UsesEdge_R4}; 

int MR_UsesVertex_F1(MRegion_ptr r, MVertex_ptr e);
int MR_UsesVertex_F4(MRegion_ptr r, MVertex_ptr e);
int MR_UsesVertex_R1(MRegion_ptr r, MVertex_ptr e);
int MR_UsesVertex_R2(MRegion_ptr r, MVertex_ptr e);
int MR_UsesVertex_R4(MRegion_ptr r, MVertex_ptr e);
int MR_UsesVertex_FNR3R4(MRegion_ptr r, MVertex_ptr e);
static int (*MR_UsesVertex_jmp[MSTK_MAXREP])(MRegion_ptr r, MVertex_ptr e) =
{MR_UsesVertex_F1, MR_UsesVertex_F4, MR_UsesVertex_R1, MR_UsesVertex_R2, 
 MR_UsesVertex_R4}; 

#else

void MR_Set_RepType_FNR3R4(MRegion_ptr r);
void MR_Set_RepType_R1(MRegion_ptr r);
void MR_Set_RepType_R2(MRegion_ptr r);
static void (*MR_Set_RepType_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_Set_RepType_FNR3R4, MR_Set_RepType_FNR3R4, MR_Set_RepType_R1, 
 MR_Set_RepType_R2, MR_Set_RepType_FNR3R4};

void MR_Delete_F1F3R3R4(MRegion_ptr r);
void MR_Delete_F4(MRegion_ptr r);
void MR_Delete_R1(MRegion_ptr r);
void MR_Delete_R2(MRegion_ptr r);
static void (*MR_Delete_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_Delete_F1F3R3R4, MR_Delete_F4, MR_Delete_R1, MR_Delete_R2, 
 MR_Delete_F1F3R3R4};

void MR_Set_Faces_F1F3R3R4(MRegion_ptr r,int nf,MFace_ptr *mfaces,int *dirs);
void MR_Set_Faces_F4(MRegion_ptr r,int nf,MFace_ptr *mfaces,int *dirs);
void MR_Set_Faces_R1(MRegion_ptr r,int nf,MFace_ptr *mfaces,int *dirs);
void MR_Set_Faces_R2(MRegion_ptr r,int nf,MFace_ptr *mfaces,int *dirs);
static void 
(*MR_Set_Faces_jmp[MSTK_MAXREP])(MRegion_ptr r, int nf, MFace_ptr *mfaces, int *dirs) =
{MR_Set_Faces_F1F3R3R4, MR_Set_Faces_F4, MR_Set_Faces_R1, MR_Set_Faces_R2, 
 MR_Set_Faces_F1F3R3R4};

void MR_Set_Vertices_FNR3R4(MRegion_ptr r, int nv, MVertex_ptr *mvertices);
void MR_Set_Vertices_R1(MRegion_ptr r, int nv, MVertex_ptr *mvertices);
void MR_Set_Vertices_R2(MRegion_ptr r, int nv, MVertex_ptr *mvertices);
static void 
(*MR_Set_Vertices_jmp[MSTK_MAXREP])(MRegion_ptr r, int nv, MVertex_ptr *mvertices) =
{MR_Set_Vertices_FNR3R4, MR_Set_Vertices_FNR3R4, MR_Set_Vertices_R1, 
 MR_Set_Vertices_R2, MR_Set_Vertices_FNR3R4};

Set_ptr MR_Vertices_FNR3R4(MRegion_ptr r);
Set_ptr MR_Vertices_R1(MRegion_ptr r);
Set_ptr MR_Vertices_R2(MRegion_ptr r);
static Set_ptr (*MR_Vertices_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_Vertices_FNR3R4, MR_Vertices_FNR3R4, MR_Vertices_R1, MR_Vertices_R2, 
 MR_Vertices_FNR3R4};

Set_ptr MR_Edges_FNR3R4(MRegion_ptr r);
Set_ptr MR_Edges_R1(MRegion_ptr r);
Set_ptr MR_Edges_R2(MRegion_ptr r);
static Set_ptr (*MR_Edges_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_Edges_FNR3R4, MR_Edges_FNR3R4, MR_Edges_R1, MR_Edges_R2, MR_Edges_FNR3R4};

Set_ptr MR_Faces_FNR3R4(MRegion_ptr r);
Set_ptr MR_Faces_R1(MRegion_ptr r);
Set_ptr MR_Faces_R2(MRegion_ptr r);
static Set_ptr (*MR_Faces_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_Faces_FNR3R4, MR_Faces_FNR3R4, MR_Faces_R1, MR_Faces_R2, MR_Faces_FNR3R4};

Set_ptr MR_AdjRegions_FNR3R4(MRegion_ptr r);
Set_ptr MR_AdjRegions_R1(MRegion_ptr r);
Set_ptr MR_AdjRegions_R2(MRegion_ptr r);
static Set_ptr (*MR_AdjRegions_jmp[MSTK_MAXREP])(MRegion_ptr r) =
{MR_AdjRegions_FNR3R4, MR_AdjRegions_FNR3R4, MR_AdjRegions_R1, 
 MR_AdjRegions_R2, MR_AdjRegions_FNR3R4};

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

void MR_Replace_Face_FNR3R4(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf, int nudir);
void MR_Replace_Face_R1(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf, int nudir);
void MR_Replace_Face_R2(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf, int nudir);
static void 
(*MR_Replace_Face_jmp[MSTK_MAXREP])(MRegion_ptr r, MFace_ptr f,MFace_ptr nuf,int nudir)=
{MR_Replace_Face_FNR3R4, MR_Replace_Face_FNR3R4, MR_Replace_Face_R1, 
 MR_Replace_Face_R2, MR_Replace_Face_FNR3R4};

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
 MR_Replace_Vertex_R2, MR_Replace_Vertex_FNR3R4};

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

int MR_UsesEdge_FNR3R4(MRegion_ptr r, MEdge_ptr e);
int MR_UsesEdge_R1(MRegion_ptr r, MEdge_ptr e);
int MR_UsesEdge_R2(MRegion_ptr r, MEdge_ptr e);
static int (*MR_UsesEdge_jmp[MSTK_MAXREP])(MRegion_ptr r, MEdge_ptr e) =
{MR_UsesEdge_FNR3R4, MR_UsesEdge_FNR3R4, MR_UsesEdge_R1, MR_UsesEdge_R2, 
 MR_UsesEdge_FNR3R4}; 

int MR_UsesVertex_FNR3R4(MRegion_ptr r, MVertex_ptr e);
int MR_UsesVertex_R1(MRegion_ptr r, MVertex_ptr e);
int MR_UsesVertex_R2(MRegion_ptr r, MVertex_ptr e);
static int (*MR_UsesVertex_jmp[MSTK_MAXREP])(MRegion_ptr r, MEdge_ptr e) =
{MR_UsesVertex_FNR3R4, MR_UsesVertex_FNR3R4, MR_UsesVertex_R1, 
 MR_UsesVertex_R2, MR_UsesVertex_FNR3R4}; 

#endif


#ifdef __cplusplus
}
#endif
