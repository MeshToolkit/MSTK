#define _H_MRegion_Private

#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

 int mrtype_nv[7]={0,0,4,5,6,8,-1};
 int mrtype_ne[7]={0,0,6,8,9,12,-1};
 int mrtype_nf[7]={0,0,4,5,5,6,-1};

  MRegion_ptr MR_New(Mesh_ptr mesh) {
    MRegion_ptr r;
    RepType rtype;

    r = (MRegion_ptr) MSTK_malloc(sizeof(MRegion));

    r->id = 0;
    r->marker = 0;
    r->mesh = mesh;
    r->dim = 3;
    r->gdim = 3;
    r->gid = 0;
    r->gent = (GEntity_ptr) NULL;
    r->sameadj = (void *) NULL;
    r->downadj = (void *) NULL;

    rtype = mesh ? MESH_RepType(mesh) : F1;
    MR_Set_RepType(r,rtype);

    if (mesh) MESH_Add_Region(mesh,r);

    return r;
  }

  void MR_Delete(MRegion_ptr r) {
    MESH_Rem_Region(r->mesh,r);
    (*MR_Delete_jmp[r->repType])(r);
  }

  void MR_Set_RepType(MRegion_ptr r, RepType rtype) {
    r->repType = rtype;
    (*MR_Set_RepType_jmp[r->repType])(r);
  }

  void MR_Set_GEntity(MRegion_ptr r, GEntity_ptr gent) {
    r->gent = gent;
  }

  void MR_Set_GEntDim(MRegion_ptr r, int gdim) {
    r->gdim = gdim;
  }

  void MR_Set_GEntID(MRegion_ptr r, int gid) {
    r->gid = gid;
  }

  void MR_Set_ID(MRegion_ptr r, int id) {
    r->id = id;
  }

  void MR_Set_Faces(MRegion_ptr r, int nf, MFace_ptr *mfaces, int *dirs) {
    (*MR_Set_Faces_jmp[r->repType])(r,nf,mfaces,dirs);
  }

  void MR_Set_Vertices(MRegion_ptr r, int nv, MFace_ptr *mvertices) {
    (*MR_Set_Vertices_jmp[r->repType])(r,nv,mvertices);
  }

  Mesh_ptr MR_Mesh(MRegion_ptr r) {
    r->mesh;
  }

  int MR_ID(MRegion_ptr r) {
    return r->id;
  }

  int MR_GEntDim(MRegion_ptr r) {
    return r->gdim;
  }

  GEntity_ptr MR_GEntity(MRegion_ptr r) {
    return r->gent;
  }

  MRType MR_ElementType(MRegion_ptr r) {
    int nt, i;
    MFace_ptr face;
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;
    switch (downadj->nf) {
    case 4:
      return TET;
    case 5:
      /* it could be a pyramid or a triangular prism */
      nt = 0;
      for (i = 0; i < downadj->nf; i++) {
	face = List_Entry(downadj->rfaces,i);
	if (MF_Num_Edges(face) == 3)
	  nt++;
      }
      if (nt == 2)
	return PRISM;
      else 
	return PYRAMID; 
    case 6:
      return HEX;
    default:
      if (downadj->nf > 6)
	return POLYHED;
      else
	return RUNKNOWN;
    }
  }

  int MR_Num_Vertices(MRegion_ptr r) {
    return mrtype_nv[MR_ElementType(r)];
  }

  int MR_Num_Edges(MRegion_ptr r) {
    return mrtype_ne[MR_ElementType(r)];
  }

  int MR_Num_Faces(MRegion_ptr r) {
    return (*MR_Num_Faces_jmp[r->repType])(r);
  }

  /*
  int MR_Num_AdjRegions(MRegion_ptr r) {
    return (*MR_Num_AdjRegions_jmp[r->repType])(r);
  }
  */

  List_ptr MR_Vertices(MRegion_ptr r) {
    return (*MR_Vertices_jmp[r->repType])(r);
  }

  List_ptr MR_Edges(MRegion_ptr r) {
    return (*MR_Edges_jmp[r->repType])(r);
  }

  List_ptr MR_Faces(MRegion_ptr r) {
    return (*MR_Faces_jmp[r->repType])(r);
  }

  List_ptr MR_AdjRegions(MRegion_ptr r) {
    return (*MR_AdjRegions_jmp[r->repType])(r);
  }

  int MR_FaceDir(MRegion_ptr r, MFace_ptr f) {
    return (*MR_FaceDir_jmp[r->repType])(r,f);
  }

  int MR_FaceDir_i(MRegion_ptr r, int i) {
    return (*MR_FaceDir_i_jmp[r->repType])(r,i);
  }

  void MR_Replace_Face(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf, int nudir) {
    (*MR_Replace_Face_jmp[r->repType])(r,f,nuf,nudir);
  }

  void MR_Replace_Vertex(MRegion_ptr r, MVertex_ptr v, MVertex_ptr nuv) {
    (*MR_Replace_Vertex_jmp[r->repType])(r,v,nuv);
  }

  void MR_Replace_Face_i(MRegion_ptr r, int i, MFace_ptr nuf, int nudir) {
    (*MR_Replace_Face_i_jmp[r->repType])(r,i,nuf,nudir);
  }

  void MR_Replace_Vertex_i(MRegion_ptr r, int i, MVertex_ptr nuv) {
    (*MR_Replace_Vertex_i_jmp[r->repType])(r,i,nuv);
  }

  void MR_Add_AdjRegion(MRegion_ptr r, int facenum, MRegion_ptr aregion) {
    (*MR_Add_AdjRegion_jmp[r->repType])(r,facenum,aregion);
  }

  void MR_Rem_AdjRegion(MRegion_ptr r, MRegion_ptr aregion) {
    (*MR_Rem_AdjRegion_jmp[r->repType])(r,aregion);
  }

  int MR_UsesEntity(MRegion_ptr r, MEntity_ptr e, int etype) {
    
    switch (etype) {
    case 3:
      return (r == (MRegion_ptr) e);
    case 2:
      return (*MR_UsesFace_jmp[r->repType])(r,(MFace_ptr)e);
    case 1:
      return (*MR_UsesEdge_jmp[r->repType])(r,(MEdge_ptr)e);
    case 0:
      return (*MR_UsesVertex_jmp[r->repType])(r,(MVertex_ptr)e);
    default:
      MSTK_Report("MR_UsesEntity","Invalid entity type",ERROR);
    }
    return 0;
  }

  
#ifdef __cplusplus
}
#endif
