#ifndef _H_MFace
#define _H_MFace

#include "MSTK_types.h"
#include "MSTK_defines.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _H_MFace_Private
#define _H_MEntity_Private

#include "MEntity.h"

  typedef struct MFace {

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

  } MFace, *MFace_ptr;

  /*----- Adjacency definitions --------*/
  /* Put in first place common fields, like downward adjacencies */
  /* The most individual fields should go to the end of structure */
  /* If there are functions which work with several different
   * representations, we must be sure, that used fields and their
   * position is the same for this representations. */
  /* Hash-related data is Ok to be in the end, since we use different
   * functions for each R-representation */

  typedef struct MFace_Adj_F1F3 {
    unsigned long long edirs;
    List_ptr fedges;
    MRegion_ptr fregions[2];
  } MFace_Adj_F1F3;

  typedef struct MFace_Adj_F2F4 {
    unsigned long long edirs;
    List_ptr fedges;
  } MFace_Adj_F2F4;

  typedef struct MFace_Adj_R1 {
    List_ptr fvertices;
    MFace_ptr hnext;
    int lock;
  } MFace_Adj_R1;

  typedef struct MFace_Adj_R2 {
    List_ptr fvertices;
    List_ptr adjfaces;
    MFace_ptr hnext;
    int lock;
  } MFace_Adj_R2;

  typedef struct MFace_Adj_R3 {
    List_ptr fvertices;
    MRegion_ptr fregions[2];
  } MFace_Adj_R3;

  typedef struct MFace_Adj_R4 {
    List_ptr fvertices;
    MRegion_ptr fregions[2];
    List_ptr adjfaces;
  } MFace_Adj_R4;

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
  void MF_Rem_Edge(MFace_ptr f, MEdge_ptr e);
  void MF_Replace_Edge(MFace_ptr f, MEdge_ptr e, int nnu, MEdge_ptr *nuedges, int *nudirs);
  void MF_Replace_Vertex(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv);
  void MF_Replace_Edge_i(MFace_ptr f, int i, int nnu, MEdge_ptr *nuedges, int *nudirs);
  void MF_Replace_Vertex_i(MFace_ptr f, int i, MVertex_ptr v);
  void MF_Insert_Vertex(MFace_ptr mface, MVertex_ptr nuv, MVertex_ptr b4v);
  void MF_Insert_Vertex_i(MFace_ptr mface, MVertex_ptr nuv, int i);
  int MF_Rev_EdgeDir(MFace_ptr f, MEdge_ptr e);
  int MF_Rev_EdgeDir_i(MFace_ptr f, int i);

  int MF_ID(MFace_ptr f);
  int MF_GEntDim(MFace_ptr f);
  GEntity_ptr MF_GEntity(MFace_ptr f);
  int MF_GEntID(MFace_ptr f);
  int MF_UsesEntity(MFace_ptr f, MEntity_ptr e, int etype);
  int MF_Num_Vertices(MFace_ptr f);
  int MF_Num_Edges(MFace_ptr f);
  int MF_Num_AdjFaces(MFace_ptr f);
  List_ptr MF_Vertices(MFace_ptr f, int dir, MVertex_ptr v);
  List_ptr MF_Edges(MFace_ptr f, int dir, MVertex_ptr v);
  List_ptr MF_AdjFaces(MFace_ptr f);
  int MF_EdgeDir(MFace_ptr f, MEdge_ptr e);
  int MF_EdgeDir_i(MFace_ptr f, int i);
  int MF_UsesEntity(MFace_ptr f, MEntity_ptr e, int type);
  List_ptr MF_Regions(MFace_ptr f);
  MRegion_ptr MF_Region(MFace_ptr f, int side);

  MFace_ptr MVs_CommonFace(int nv, MVertex_ptr *fverts);
  MFace_ptr MEs_CommonFace(int ne, MEdge_ptr *fedges);

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

  void MF_Lock(MFace_ptr f);
  void MF_UnLock(MFace_ptr f);
  int MF_IsLocked(MFace_ptr f);

#ifdef MSTK_HAVE_MPI
  PType MF_PType(MFace_ptr f);  
  void  MF_Set_PType(MFace_ptr f, PType ptype);
  int   MF_MasterParID(MFace_ptr f);
  void  MF_Set_MasterParID(MFace_ptr f, int masterparid);
  int   MF_GlobalID(MFace_ptr f);
  void  MF_Set_GlobalID(MFace_ptr f, int globalid);
  MFace_ptr MF_GhostNew(Mesh_ptr mesh);
#endif

#ifdef __cplusplus
}
#endif

#endif
