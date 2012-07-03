#ifndef _H_MRegion
#define _H_MRegion

#include "MSTK_defines.h"
#include "MSTK_types.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int mrtype_nv[6], mrtype_ne[6], mrtype_nf[6];

#ifdef _H_MRegion_Private
#define _H_MEntity_Private

#include "MEntity.h"

  typedef struct MRegion {

    /* Common data structure for all mesh entities */

    MEntity_Data entdat;

    /* Specific to mesh regions */
    void *adj;

  } MRegion, *MRegion_ptr;

  /*----- Adjacency definitions --------*/

  typedef struct MRegion_Adj_FN {
    unsigned int *fdirs;
    MFace_ptr *rfaces;
  } MRegion_Adj_FN;

  typedef struct MRegion_Adj_R1 {
    MVertex_ptr *rvertices;
    int **fvtemplate;
  } MRegion_Adj_R1;

  typedef struct MRegion_Adj_R2 {
    MVertex_ptr *rvertices;
    int **fvtemplate;
    MRegion_ptr *aregions;
  } MRegion_Adj_R2;

  typedef struct MRegion_Adj_R3R4 {
    unsigned int *fdirs;
    MFace_ptr *rfaces;
  } MRegion_Adj_R3R4;

#else
  typedef void *MRegion_ptr;
#endif


  MRegion_ptr MR_New(Mesh_ptr mesh);
  void MR_Delete(MRegion_ptr r, int keep);
  void MR_Set_GEntity(MRegion_ptr r, GEntity_ptr gent);
  void MR_Set_GEntID(MRegion_ptr r, int gid);
  void MR_Set_ID(MRegion_ptr r, int id);

  /* Can be called by user/application after creating region by MR_New(); */
  void MR_Set_Faces(MRegion_ptr r, int nf, MFace_ptr *mfaces, int *dirs);
  void MR_Set_Vertices(MRegion_ptr r, int nv, MVertex_ptr *mvertices, int nf, int **rfvtemplate);

  /* Can be called by mesh modification routines */
  void MR_Rem_Face(MRegion_ptr r, MFace_ptr f);
  void MR_Replace_Face(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf, int dir);
  void MR_Replace_Vertex(MRegion_ptr r, MVertex_ptr v, MVertex_ptr nuv);
  void MR_Replace_Face_i(MRegion_ptr r, int i, MFace_ptr f, int dir);
  void MR_Replace_Vertex_i(MRegion_ptr r, int i, MVertex_ptr v);
  int MR_Rev_FaceDir(MRegion_ptr r, MFace_ptr f);
  int MR_Rev_FaceDir_i(MRegion_ptr r, int i);

  int MR_ID(MRegion_ptr r);
  GEntity_ptr MR_GEntity(MRegion_ptr r);
  int MR_GEntID(MRegion_ptr r);
  MRType MR_ElementType(MRegion_ptr r);

  int MR_Num_Vertices(MRegion_ptr r);
  int MR_Num_Edges(MRegion_ptr r);
  int MR_Num_Faces(MRegion_ptr r);
  int MR_Num_AdjRegions(MRegion_ptr r);
  List_ptr MR_Vertices(MRegion_ptr r);
  List_ptr MR_Edges(MRegion_ptr r);
  List_ptr MR_Faces(MRegion_ptr r);
  List_ptr MR_AdjRegions(MRegion_ptr r);
  int MR_FaceDir(MRegion_ptr r, MFace_ptr f);
  int MR_FaceDir_i(MRegion_ptr r, int i);
  int MR_UsesEntity(MRegion_ptr r, MEntity_ptr e, int type);

  /* Users can only set the representation type for the entire mesh,
     not for individual mesh entities */
  void MR_Set_RepType(MRegion_ptr r, RepType rtype);

  /* Adjacent region info will be updated by private functions when
     regions are created or deleted. There is no need for invocation
     by users/applications */
  void MR_Add_AdjRegion(MRegion_ptr r, int facenum, MRegion_ptr ar);
  void MR_Rem_AdjRegion(MRegion_ptr r, MRegion_ptr ar);

  PType MR_PType(MRegion_ptr r);  
  void  MR_Set_PType(MRegion_ptr r, PType ptype);
  int   MR_MasterParID(MRegion_ptr r);
  void  MR_Set_MasterParID(MRegion_ptr r, int masterparid);
  int   MR_GlobalID(MRegion_ptr r);
  void  MR_Set_GlobalID(MRegion_ptr r, int globalid);
  MRegion_ptr MR_GhostNew(Mesh_ptr mesh);

#ifdef __cplusplus
}
#endif

#endif
