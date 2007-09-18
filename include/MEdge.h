#ifndef _H_MEdge
#define _H_MEdge

#ifdef __cplusplus
extern "C" {
#endif

#include "MSTK_defines.h"
#include "MSTK_types.h"
#include "List.h"

#ifdef _H_MEdge_Private
  typedef struct MEdge {

    /* Common data structure for all mesh entities */

    MEntity_Data entdat;

    /* Specific to mesh edges */

    void *adj;
    MVertex_ptr vertex[2];

  } MEdge, *MEdge_ptr;

  /*----- Adjacency definitions --------*/

  typedef struct MEdge_Adj_F1 {
    List_ptr efaces;
  } MEdge_Adj_F1;

  typedef struct MEdge_Adj_F4 {
    List_ptr elements;
  } MEdge_Adj_F4;

  typedef struct MEdge_Adj_RN {
    MEdge_ptr hnext;
    int lock;
  } MEdge_Adj_RN;
#else
  typedef void *MEdge_ptr;
#endif

  MEdge_ptr ME_New(Mesh_ptr mesh);
  void ME_Delete(MEdge_ptr medge, int keep);
  void ME_Set_GEntity(MEdge_ptr medge, GEntity_ptr gent);
  void ME_Set_GEntDim(MEdge_ptr medge, int gdim);
  void ME_Set_GEntID(MEdge_ptr medge, int gid);
  int  ME_Derive_GInfo(MEdge_ptr medge);
  void ME_Set_ID(MEdge_ptr medge, int id);

  void ME_Set_Vertex(MEdge_ptr medge, int i, MVertex_ptr vertex);
  
  int ME_ID(MEdge_ptr medge);
  int ME_GEntDim(MEdge_ptr medge);
  int ME_GEntID(MEdge_ptr medge);
  GEntity_ptr ME_GEntity(MEdge_ptr medge);
  MVertex_ptr ME_Vertex(MEdge_ptr medge, int i);
  MVertex_ptr ME_OppVertex(MEdge_ptr e, MVertex_ptr v);

  int ME_Num_Faces(MEdge_ptr medge);
  int ME_Num_Regions(MEdge_ptr medge);
  List_ptr ME_Faces(MEdge_ptr medge);
  List_ptr ME_Regions(MEdge_ptr medge);
  int ME_UsesEntity(MEdge_ptr medge, MEntity_ptr mentity, int etype);

  /* Calling applications can only set the representation type for the
     entire mesh, not for individual mesh entities */
  void ME_Set_RepType(MEdge_ptr medge, RepType rtype);

  void ME_Add_Face(MEdge_ptr medge, MFace_ptr mface);
  void ME_Add_Region(MEdge_ptr medge, MRegion_ptr mregion);
  void ME_Rem_Face(MEdge_ptr medge, MFace_ptr mface);
  void ME_Rem_Region(MEdge_ptr medge, MRegion_ptr mregion);

  void ME_Lock(MEdge_ptr e);
  void ME_UnLock(MEdge_ptr e);
  int ME_IsLocked(MEdge_ptr e);

#ifdef __cplusplus
}
#endif

#endif
