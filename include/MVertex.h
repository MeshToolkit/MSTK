#ifndef _H_MVertex
#define _H_MVertex

#include "MSTK_defines.h"
#include "MSTK_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _H_MVertex_Private
#define _H_MEntity_Private

#include "MEntity.h"

  typedef struct MVertex {

    /*------------------------------------------------------------*/
    /* Common data for all mesh entities */

    Mesh_ptr mesh;
    List_ptr AttInsList;

    unsigned int dim_id;
    unsigned int rtype_gdim_gid;
    unsigned int marker;

#ifdef MSTK_HAVE_MPI
    unsigned int ptype_masterparid;
    unsigned int globalid;  /* if -ve, it represents local id on master proc */
#endif
    /*------------------------------------------------------------*/

    void *adj; /* pointer to entity adjacency structure */

    /* Vertex coordinates */
    double xyz[3];

  } MVertex, *MVertex_ptr;

    /*----- Upward adjacency definitions --------*/

  typedef struct MVertex_Adj_F1F4 {
    List_ptr vedges;
  } MVertex_Adj_F1F4;

  typedef struct MVertex_Adj_F2 {
    List_ptr vregions;
  } MVertex_Adj_F2;

  typedef struct MVertex_Adj_F3 {
    List_ptr vfaces;
  } MVertex_Adj_F3;

  typedef struct MVertex_Adj_R1 {
    List_ptr velements;
  } MVertex_Adj_R1;

  typedef struct MVertex_Adj_R2 {
    List_ptr velements;
    List_ptr adjverts;
  } MVertex_Adj_R2;

  typedef struct MVertex_Adj_R3 {
    List_ptr vfaces;
  } MVertex_Adj_R3;

  typedef struct MVertex_Adj_R4 {
    List_ptr vfaces;
    List_ptr adjverts;
  } MVertex_Adj_R4;

#else
  typedef void *MVertex_ptr;
#endif


  /*-------- Interface Declarations -----------*/

  MVertex_ptr MV_New(Mesh_ptr mesh);
  void MV_Delete(MVertex_ptr mvertex, int keep);
  void MV_Set_Coords(MVertex_ptr mvertex, double *xyz);
  void MV_Set_GEntity(MVertex_ptr mvertex, GEntity_ptr gent);
  void MV_Set_GEntID(MVertex_ptr mvertex, int gid);
  void MV_Set_GEntDim(MVertex_ptr mvertex, int gdim);
  void MV_Set_ID(MVertex_ptr mvertex, int id);

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

  void MV_Add_AdjVertex(MVertex_ptr mvertex, MVertex_ptr adjvertex);
  void MV_Rem_AdjVertex(MVertex_ptr mvertex, MVertex_ptr adjvertex);


  /* Calling applications can only set the representation type for the
     entire mesh, not for individual mesh entities */
  void MV_Set_RepType(MVertex_ptr v, RepType rtype);

  void MV_Add_Edge(MVertex_ptr mvertex, MEdge_ptr medge);
  void MV_Add_Face(MVertex_ptr mvertex, MFace_ptr mface);
  void MV_Add_Region(MVertex_ptr mvertex, MRegion_ptr mregion);
  void MV_Rem_Edge(MVertex_ptr mvertex, MEdge_ptr medge);
  void MV_Rem_Face(MVertex_ptr mvertex, MFace_ptr mface);
  void MV_Rem_Region(MVertex_ptr mvertex, MRegion_ptr mregion);

#ifdef MSTK_HAVE_MPI
  PType MV_PType(MVertex_ptr mvertex);  
  void  MV_Set_PType(MVertex_ptr mvertex, PType ptype);

  /* Is the entity on the partition boundary? */
  int MV_OnParBoundary(MVertex_ptr mvertex);

  /* Mark/Unmark the entity as being on the partition boundary */
  void MV_Flag_OnParBoundary(MVertex_ptr mvertex);
  void MV_Unflag_OnParBoundary(MVertex_ptr mvertex);

  int   MV_MasterParID(MVertex_ptr mvertex);
  void  MV_Set_MasterParID(MVertex_ptr mvertex, int masterparid);
  int   MV_GlobalID(MVertex_ptr mvertex);
  void  MV_Set_GlobalID(MVertex_ptr mvertex, int globalid);
  MVertex_ptr MV_GhostNew(Mesh_ptr mesh);
#endif

#ifdef __cplusplus
}
#endif

#endif
