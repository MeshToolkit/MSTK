#ifndef _H_MSTK
#define _H_MSTK

#ifdef   __cplusplus
extern "C" {
#endif

  /* #include "gmtk.h" */
#include "MSTK_defines.h"
#include "MSTK_types.h"
#include "MSTK_util.h"
#include "List.h"
#include "MSTK_malloc.h"


/*******************************************************************/

void        MSTK_Init(void);

/********************************************************************/
/*        MESH OBJECT OPERATORS                                     */
/********************************************************************/  

  Mesh_ptr   MESH_New(RepType type);
  void       MESH_Delete(Mesh_ptr mesh);
  int        MESH_InitFromFile(Mesh_ptr mesh, const char *filename);
  void       MESH_WriteToFile(Mesh_ptr mesh, const char *filename);

  GModel_ptr MESH_GModel(Mesh_ptr mesh);
  RepType    MESH_RepType(Mesh_ptr mesh);

  int        MESH_Num_Vertices(Mesh_ptr mesh);
  int        MESH_Num_Edges(Mesh_ptr mesh);
  int        MESH_Num_Faces(Mesh_ptr mesh);
  int        MESH_Num_Regions(Mesh_ptr mesh);

  MVertex_ptr  MESH_Vertex(Mesh_ptr mesh, int i);
  MEdge_ptr    MESH_Edge(Mesh_ptr mesh, int i);
  MFace_ptr    MESH_Face(Mesh_ptr mesh, int i);
  MRegion_ptr  MESH_Region(Mesh_ptr mesh, int i);

  MVertex_ptr  MESH_Next_Vertex(Mesh_ptr mesh, int *index);
  MEdge_ptr    MESH_Next_Edge(Mesh_ptr mesh, int *index);
  MFace_ptr    MESH_Next_Face(Mesh_ptr mesh, int *index);
  MRegion_ptr  MESH_Next_Region(Mesh_ptr mesh, int *index);


  void       MESH_Set_GModel(Mesh_ptr mesh, GModel_ptr geom);
  int        MESH_Change_RepType(Mesh_ptr mesh, int nurep);

/********************************************************************/
/*        MESH VERTEX OPERATORS                                     */
/********************************************************************/

  MVertex_ptr MV_New(Mesh_ptr mesh);
  void MV_Delete(MVertex_ptr mvertex);
  void MV_Set_RepType(MVertex_ptr mvertex, RepType reptype);
  void MV_Set_Coords(MVertex_ptr mvertex, double *xyz);
  void MV_Set_GEntity(MVertex_ptr mvertex, GEntity_ptr gent);
  void MV_Set_GEntDim(MVertex_ptr mvertex, int gdim);
  void MV_Set_GEntID(MVertex_ptr mvertex, int gid);
  void MV_Add_AdjVertex(MVertex_ptr mvertex, MVertex_ptr adjvertex);
  void MV_Rem_AdjVertex(MVertex_ptr mvertex, MVertex_ptr adjvertex);
  void MV_Set_ID(MVertex_ptr mvertex, int id);

  Mesh_ptr MV_Mesh(MVertex_ptr mv);
  int MV_ID(MVertex_ptr mvertex);
  int MV_GEntDim(MVertex_ptr mvertex);
  int MV_GEntID(MVertex_ptr mvertex);
  GEntity_ptr MV_GEntity(MVertex_ptr mvertex);

  void MV_Coords(MVertex_ptr mvertex, double *xyz);

  int MV_Num_AdjVertices(MVertex_ptr mvertex);
  int MV_Num_Edges(MVertex_ptr mvertex);
  int MV_Num_Faces(MVertex_ptr mvertex);
  int MV_Num_Regions(MVertex_ptr mvertex);
  List_ptr MV_AdjVertices(MVertex_ptr mvertex);
  List_ptr MV_Edges(MVertex_ptr mvertex);
  List_ptr MV_Faces(MVertex_ptr mvertex);
  List_ptr MV_Regions(MVertex_ptr mvertex);


/********************************************************************/
/*        MESH EDGE OPERATORS                                       */
/********************************************************************/

  MEdge_ptr ME_New(Mesh_ptr mesh);
  void ME_Delete(MEdge_ptr medge);
  void ME_Set_RepType(MEdge_ptr medge, RepType reptype);
  void ME_Set_GEntity(MEdge_ptr medge, GEntity_ptr gent);
  void ME_Set_GEntDim(MEdge_ptr medge, int gdim);
  void ME_Set_GEntID(MEdge_ptr medge, int gid);
  void ME_Set_ID(MEdge_ptr medge, int id);

  void ME_Set_Vertex(MEdge_ptr medge, int i, MVertex_ptr vertex);
  void ME_Replace_Vertex(MEdge_ptr medge, MVertex_ptr vert, MVertex_ptr nuvert);

  Mesh_ptr ME_Mesh(MEdge_ptr medge);
  int ME_ID(MEdge_ptr medge);
  int ME_GEntDim(MEdge_ptr medge);
  int ME_GEntID(MEdge_ptr medge);
  GEntity_ptr ME_GEntity(MEdge_ptr medge);
  int ME_Num_Faces(MEdge_ptr medge);
  int ME_Num_Regions(MEdge_ptr medge);
  MVertex_ptr ME_Vertex(MEdge_ptr medge, int i);
  MVertex_ptr ME_OppVertex(MEdge_ptr medge, MVertex_ptr ov);
  int ME_UsesEntity(MEdge_ptr medge, MEntity_ptr mentity, int etype);
  List_ptr ME_Faces(MEdge_ptr medge);
  List_ptr ME_Regions(MEdge_ptr medge);


  MEdge_ptr MVs_CommonEdge(MVertex_ptr v1, MVertex_ptr v2);

  double ME_Len(MEdge_ptr e);
  double ME_LenSqr(MEdge_ptr e);
  void   ME_Vec(MEdge_ptr e, double *evec);

/********************************************************************/
/*        MESH FACE OPERATORS                                       */
/********************************************************************/
  MFace_ptr MF_New(Mesh_ptr mesh);
  void MF_Delete(MFace_ptr mface);
  void MF_Set_RepType(MFace_ptr mface, RepType reptype);
  void MF_Set_GEntity(MFace_ptr mface, GEntity_ptr gent);
  void MF_Set_GEntDim(MFace_ptr mface, int gdim);
  void MF_Set_GEntID(MFace_ptr mface, int gid);
  void MF_Set_ID(MFace_ptr mface, int id);

  /* These can be called by a user/application after calling MF_New() */
  void MF_Set_Edges(MFace_ptr mface, int n, MEdge_ptr *edges, int *dirs);
  void MF_Set_Vertices(MFace_ptr mface, int n, MVertex_ptr *verts);

  /* Can be called by higher level mesh modification operators */
  void MF_Replace_Edge(MFace_ptr mface, MEdge_ptr edge, MEdge_ptr nuedge, int dir);
  void MF_Replace_Vertex(MFace_ptr mface, MVertex_ptr mvertex, MVertex_ptr nuvertex);
  void MF_Replace_Edge_i(MFace_ptr mface, int i, MEdge_ptr nuedge, int dir);
  void MF_Replace_Vertex_i(MFace_ptr mface, int i, MVertex_ptr nuvertex);

  Mesh_ptr MF_Mesh(MFace_ptr mf);
  int MF_ID(MFace_ptr mface);
  int MF_GEntDim(MFace_ptr mface);
  int MF_GEntID(MFace_ptr mface);
  GEntity_ptr MF_GEntity(MFace_ptr mface);
  int MF_Num_Vertices(MFace_ptr mface);
  int MF_Num_Edges(MFace_ptr mface);
  int MF_Num_AdjFaces(MFace_ptr mface);
  List_ptr MF_Vertices(MFace_ptr mface, int dir);
  List_ptr MF_Edges(MFace_ptr mface, int dir, MVertex_ptr mvert);
  List_ptr MF_AdjFaces(MFace_ptr mface);
  int MF_EdgeDir(MFace_ptr mface, MEdge_ptr medge);
  int MF_EdgeDir_i(MFace_ptr mface, int i);
  int MF_UsesEntity(MFace_ptr mface, MEntity_ptr mentity, int type);
  List_ptr MF_Regions(MFace_ptr mface);
  MRegion_ptr MF_Region(MFace_ptr mface, int side);



  void MF_Coords(MFace_ptr mface, int *n, double (*xyz)[3]);



/********************************************************************/
/*        MESH REGN OPERATORS                                       */
/********************************************************************/

  MRegion_ptr MR_New(Mesh_ptr mesh);
  void MR_Delete(MRegion_ptr mregion);
  void MR_Set_RepType(MRegion_ptr mregion, RepType reptype);
  void MR_Set_GEntity(MRegion_ptr mregion, GEntity_ptr gent);
  void MR_Set_GEntDim(MRegion_ptr mregion, int gdim);
  void MR_Set_GEntID(MRegion_ptr mregion, int gid);
  void MR_Set_ID(MRegion_ptr mregion, int id);

  /* Can be called by user/application after creating region by MR_New(); */
  void MR_Set_Faces(MRegion_ptr mregion, int nf, MFace_ptr *mfaces, int *dirs);
  void MR_Set_Vertices(MRegion_ptr mregion, int nv, MVertex_ptr *mvertices);

  /* Can be called by mesh modification routines */
  void MR_Replace_Face(MRegion_ptr mregion, MFace_ptr mface, MFace_ptr nuface, int dir);
  void MR_Replace_Vertex(MRegion_ptr mregion, MVertex_ptr mvertex, MVertex_ptr nuvertex);
  void MR_Replace_Face_i(MRegion_ptr mregion, int i, MFace_ptr mface, int dir);
  void MR_Replace_Vertex_i(MRegion_ptr mregion, int i, MVertex_ptr mvertex);


  Mesh_ptr MR_Mesh(MRegion_ptr mregion);
  int MR_ID(MRegion_ptr mregion);
  int MR_GEntDim(MRegion_ptr mregion);
  int MR_GEntID(MRegion_ptr mregion);
  GEntity_ptr MR_GEntity(MRegion_ptr mregion);

  int MR_Num_Vertices(MRegion_ptr mregion);
  int MR_Num_Edges(MRegion_ptr mregion);
  int MR_Num_Faces(MRegion_ptr mregion);
  int MR_Num_AdjRegions(MRegion_ptr mregion);
  List_ptr MR_Vertices(MRegion_ptr mregion);
  List_ptr MR_Edges(MRegion_ptr mregion);
  List_ptr MR_Faces(MRegion_ptr mregion);
  List_ptr MR_AdjRegions(MRegion_ptr mregion);
  int MR_FaceDir(MRegion_ptr mregion, MFace_ptr mface);
  int MR_FaceDir_i(MRegion_ptr mregion, int i);
  int MR_UsesEntity(MRegion_ptr mregion, MEntity_ptr ment, int type);


  void MR_Coords(MRegion_ptr mregion, int *n, double (*xyz)[3]);


  /************************************************************************/
  /* GENERIC ENTITY OPERATORS                                             */
  /************************************************************************/

  void MEnt_Set_GEntity(MEntity_ptr mentity, GEntity_ptr gent);
  void MEnt_Set_GEntDim(MEntity_ptr mentity, int gdim);
  void MEnt_Set_ID(MEntity_ptr mentity, int id);

  int MEnt_ID(MEntity_ptr mentity);
  int MEnt_Dim(MEntity_ptr mentity);
  Mesh_ptr MEnt_Mesh(MEntity_ptr mentity);
  int MEnt_GEntDim(MEntity_ptr mentity);
  GEntity_ptr MEnt_GEntity(MEntity_ptr mentity);
  

  /************************************************************************/
  /* ENTITY MARKING                                                       */
  /************************************************************************/

  int MSTK_GetMarker();
  void MSTK_FreeMarker(int mkr);
  void MEnt_Mark(MEntity_ptr ent, int mkr);
  int MEnt_IsMarked(MEntity_ptr ent, int mkr);
  void MEnt_Unmark(MEntity_ptr ent, int mkr);
  void List_Mark(List_ptr list, int mkr);
  void List_Unmark(List_ptr list, int mkr);


/*************************************************************************/
/*  SOME GENERIC OPERATORS FOR MESH REGIONS BASED ON TYPE AND CONVENTION */
/*************************************************************************/

int         RType_NumVerts(MRType type);
int         RType_NumFaces(MRType type);
int         RType_NumFVerts(MRType type, int locfnum);
int         RType_LocFVNums(MRType type, int locfnum, MVertex_ptr *lverts);


/***********************************************************************/
/* Mesh Modification Operators                                         */
/***********************************************************************/

int ME_Swap2D(MEdge_ptr e, MEdge_ptr *enew, MFace_ptr fnew[2]);
  MFace_ptr MFs_Join(MFace_ptr f1, MFace_ptr f2, MEdge_ptr e);


/**********************************************************************/
/* Higher level Mesh Query Operators                                  */
/**********************************************************************/



/**********************************************************************/
/* Element quality evaluation                                         */
/**********************************************************************/

  void MF_CondNums(MFace_ptr f, double *condnums);

/*******************************************/
/* UTILITIES                               */
/*******************************************/

void        MSTK_Report(char *module, char *message, ErrType severity);

#ifdef DEBUG
  void        List_PrintID(List_ptr l);
#endif


#ifdef __cplusplus
	   }
#endif

#endif
