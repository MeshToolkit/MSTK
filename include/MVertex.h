#ifndef _H_MVertex
#define _H_MVertex

#ifdef __cplusplus
extern "C" {
#endif

#include "MSTK_defines.h"
#include "MSTK_types.h"
#include "List.h"
#include "Hash.h"

#ifdef _H_MVertex_Private

  typedef struct MVertex {

    /* Common data structure for all mesh entities */

    MEntity_Data_ptr entdat;

    /* specific to mesh vertices */

    void *upadj;
    void *sameadj;
    double xyz[3];

  } MVertex, *MVertex_ptr;

    /*----- Upward adjacency definitions --------*/

  typedef struct MVertex_UpAdj_F1F4 {
    List_ptr vedges;
  } MVertex_UpAdj_F1F4;

  typedef struct MVertex_UpAdj_F2 {
    List_ptr vregions;
  } MVertex_UpAdj_F2;

  typedef struct MVertex_UpAdj_F3 {
    List_ptr vfaces;
  } MVertex_UpAdj_F3;

  typedef struct MVertex_UpAdj_F5 {
    List_ptr vedges;
    List_ptr vfaces;
    List_ptr vregions;
  } MVertex_UpAdj_F5;

  typedef struct MVertex_UpAdj_F6 {
    List_ptr vfaces;
    List_ptr vregions;
  } MVertex_UpAdj_F6;

  typedef struct MVertex_UpAdj_R1R2 {
    List_ptr velements;
  } MVertex_UpAdj_R1R2;

  typedef struct MVertex_UpAdj_R3R4 {
    List_ptr vfaces;
  } MVertex_UpAdj_R3R4;

  /*-------  Same Level adjacency definitions --------*/

  typedef struct MVertex_SameAdj_R2R4 {
    List_ptr adjverts;
  } MVertex_SameAdj_R2R4;

  /*-------- Downward adjacency definitions ----------*/

  /*    NONE     */

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

#ifdef __cplusplus
}
#endif

#endif
