#ifndef _H_MEdge
#define _H_MEdge

#ifdef __cplusplus
extern "C" {
#endif

#include "MSTK_types.h"
#include "List.h"

#ifdef _H_MEdge_Private
  typedef struct MEdge {
    int id;
    int marker;
    Mesh_ptr mesh;
    char dim;
    char gdim;
    int gid;
    GEntity_ptr gent;
    RepType repType;
    void *upadj;
    MVertex_ptr vertex[2];
  } MEdge, *MEdge_ptr;

  /*----- Upward adjacency definitions --------*/

  typedef struct MEdge_UpAdj_F1 {
    unsigned int nf;
    List_ptr efaces;
  } MEdge_UpAdj_F1;

  typedef struct MEdge_UpAdj_F4 {
    unsigned int nel;
    List_ptr elements;
  } MEdge_UpAdj_F4;
#else
  typedef void *MEdge_ptr;
#endif

  MEdge_ptr ME_New(Mesh_ptr mesh);
  void ME_Set_GEntity(MEdge_ptr medge, GEntity_ptr gent);
  void ME_Set_GEntDim(MEdge_ptr medge, int gdim);
  void ME_Set_GEntID(MEdge_ptr medge, int gid);
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

#ifdef __cplusplus
}
#endif

#endif
