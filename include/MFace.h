#ifndef _H_MFace
#define _H_MFace

#ifdef __cplusplus
extern "C" {
#endif

#include "MSTK_types.h"
#include "Set.h"

#ifdef _H_MFace_Private
  typedef struct MFace {
    int id;
    int marker;
    Mesh_ptr mesh;
    char dim;
    char gdim;
    int gid;
    GEntity_ptr gent;
    RepType repType;
    void *upadj;
    void *sameadj;
    void *downadj;
  } MFace, *MFace_ptr;

  /*----- Upward adjacency definitions --------*/

  typedef struct MFace_UpAdj_F1F3 {
    MRegion_ptr fregions[2];
  } MFace_UpAdj_F1F3;

  typedef struct MFace_UpAdj_R3R4 {
    MRegion_ptr fregions[2];
  } MFace_UpAdj_R3R4;

  /*----- Same Level adjacency definitions ------*/

  typedef struct MFace_SameAdj_R2R4 {
    unsigned char nadj;
    Set_ptr adjfaces;
  } MFace_SameAdj_R2R4;

  /*----- Downward adjacency definitions --------*/

  typedef struct MFace_DownAdj_FN {
    unsigned char ne;
    int edirs;
    Set_ptr fedges;
  } MFace_DownAdj_FN;

  typedef struct MFace_DownAdj_R3R4 {
    unsigned char nv;
    Set_ptr *fvertices;
  } MFace_DownAdj_R3R4;
  
#else
  typedef void *MFace_ptr;
#endif  

  /*--------- Interface declaration -----------*/

  MFace_ptr MF_New(Mesh_ptr mesh);
  void MF_Set_GEntity(MFace_ptr f, GEntity_ptr gent);
  void MF_Set_GEntDim(MFace_ptr f, int gdim);
  void MF_Set_GEntID(MFace_ptr f, int gid);
  void MF_Set_ID(MFace_ptr f, int id);

  /* These can be called by a user/application after calling MF_New() */
  void MF_Set_Edges(MFace_ptr f, int n, MEdge_ptr *edges, int *dirs);
  void MF_Set_Vertices(MFace_ptr f, int n, MVertex_ptr *verts);

  /* Can be called by higher level mesh modification operators */
  void MF_Replace_Edge(MFace_ptr f, MEdge_ptr e, MEdge_ptr nue, int dir);
  void MF_Replace_Vertex(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
  void MF_Replace_Edge_i(MFace_ptr f, int i, MEdge_ptr e, int dir);
  void MF_Replace_Vertex_i(MFace_ptr f, int i, MVertex_ptr v);

  int MF_ID(MFace_ptr f);
  int MF_GEntDim(MFace_ptr f);
  GEntity_ptr MF_GEntity(MFace_ptr f);
  int MF_GEntID(MFace_ptr f);
  int MF_UsesEntity(MFace_ptr f, MEntity_ptr e, int etype);
  int MF_Num_Vertices(MFace_ptr f);
  int MF_Num_Edges(MFace_ptr f);
  int MF_Num_AdjFaces(MFace_ptr f);
  Set_ptr MF_Vertices(MFace_ptr f, int dir, MVertex_ptr v);
  Set_ptr MF_Edges(MFace_ptr f, int dir, MVertex_ptr v);
  Set_ptr MF_AdjFaces(MFace_ptr f);
  int MF_EdgeDir(MFace_ptr f, MEdge_ptr e);
  int MF_EdgeDir_i(MFace_ptr f, int i);
  int MF_UsesEntity(MFace_ptr f, MEntity_ptr e, int type);
  Set_ptr MF_Regions(MFace_ptr f);
  MRegion_ptr MF_Region(MFace_ptr f, int side);

  /* Calling applications can only set the representation type for the
     entire mesh, not for individual mesh entities */
  void MF_Set_RepType(MFace_ptr f, RepType rtype);

  /* Add, Remove AdjFace can be automatically called when faces are
     being created or deleted. They do not need user invocation */
  void MF_Add_AdjFace(MFace_ptr f, int enbr, MFace_ptr af);
  void MF_Rem_AdjFace(MFace_ptr f, int enbr, MFace_ptr af);

  /* Add, Remove Region can be called internally when a region is
     being created from a list of faces */
  void MF_Add_Region(MFace_ptr f, MRegion_ptr r, int side);
  void MF_Rem_Region(MFace_ptr f, MRegion_ptr r);

#ifdef __cplusplus
}
#endif

#endif
