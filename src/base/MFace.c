#define _H_MFace_Private

#include "MFace.h"
#include "MFace_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"


#ifdef __cplusplus
extern "C" {
#endif

  MFace_ptr MF_New(Mesh_ptr mesh) {
    MFace_ptr f;
    RepType rtype;

    f = (MFace_ptr) MSTK_malloc(sizeof(MFace));

    f->id = 0;
    f->marker = 0;
    f->mesh = mesh;
    f->dim = 2;
    f->gdim = 4; /* nonsensical value since we don;t know what it is */
    f->gid = 0;
    f->gent = (GEntity_ptr) NULL;
    f->upadj = (void *) NULL;
    f->sameadj = (void *) NULL;
    f->downadj = (void *) NULL;

    rtype = mesh ? MESH_RepType(mesh) : F1;
    MF_Set_RepType(f,rtype);
    
    if (mesh) MESH_Add_Face(mesh,f);

    return f;
  }

  void MF_Delete(MFace_ptr f) {
    MESH_Rem_Face(f->mesh,f);
    (*MF_Delete_jmp[f->repType])(f);
  }

  void MF_Set_RepType(MFace_ptr f, RepType rtype) {
    f->repType = rtype;
    (*MF_Set_RepType_jmp[f->repType])(f);
  }

  void MF_Set_GEntity(MFace_ptr f, GEntity_ptr gent) {
    f->gent = gent;
  }

  void MF_Set_GEntDim(MFace_ptr f, int gdim) {
    f->gdim = gdim;
  }

  void MF_Set_GEntID(MFace_ptr f, int gid) {
    f->gid = gid;
  }

  void MF_Set_ID(MFace_ptr f, int id) {
    f->id = id;
  }

  void MF_Set_Edges(MFace_ptr f, int n, MEdge_ptr *edges, int *dir) {
    (*MF_Set_Edges_jmp[f->repType])(f,n,edges,dir);
  }

  void MF_Replace_Edge(MFace_ptr f, MEdge_ptr e, MEdge_ptr nue, int nudir) {
    (*MF_Replace_Edge_jmp[f->repType])(f,e,nue,nudir);
  }

  void MF_Replace_Edge_i(MFace_ptr f, int i, MEdge_ptr nue, int nudir) {
    (*MF_Replace_Edge_i_jmp[f->repType])(f,i,nue,nudir);
  }

  void MF_Set_Vertices(MFace_ptr f, int n, MVertex_ptr *verts) {
    (*MF_Set_Vertices_jmp[f->repType])(f,n,verts);
  }

  void MF_Replace_Vertex(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv) {
    (*MF_Replace_Vertex_jmp[f->repType])(f,v,nuv);
  }

  void MF_Replace_Vertex_i(MFace_ptr f, int i, MVertex_ptr nuv) {
    (*MF_Replace_Vertex_i_jmp[f->repType])(f,i,nuv);
  }

  Mesh_ptr MF_Mesh(MFace_ptr f) {
    return f->mesh;
  }

  int MF_ID(MFace_ptr f) {
    return f->id;
  }

  int MF_GEntDim(MFace_ptr f) {
    return f->gdim;
  }

  GEntity_ptr MF_GEntity(MFace_ptr f) {
    return f->gent;
  }

  int MF_UsesEntity(MFace_ptr f, MEntity_ptr e, int etype) {

    switch (etype) {
    case 3:
      return 0;
    case 2:
      return (f == (MFace_ptr) e);
    case 1:
      return (*MF_UsesEdge_jmp[f->repType])(f,(MEdge_ptr)e);
    case 0:
      return (*MF_UsesVertex_jmp[f->repType])(f,(MVertex_ptr)e);
    default:
      MSTK_Report("MF_UsesEntity","Invalid entity type",ERROR);
    }
    return 0;
  }

  int MF_Num_Vertices(MFace_ptr f) {
    return (*MF_Num_Vertices_jmp[f->repType])(f);
  }

  int MF_Num_Edges(MFace_ptr f) {
    return (*MF_Num_Edges_jmp[f->repType])(f);
  }

  List_ptr MF_Vertices(MFace_ptr f, int dir) {
    return (*MF_Vertices_jmp[f->repType])(f,dir);
  }

  List_ptr MF_Edges(MFace_ptr f, int dir, MVertex_ptr v) {
    return (*MF_Edges_jmp[f->repType])(f,dir,v);
  }

  int MF_Num_AdjFaces(MFace_ptr f) {
    return (*MF_Num_AdjFaces_jmp[f->repType])(f);
  }

  int MF_EdgeDir(MFace_ptr f, MEdge_ptr e) {
    return (*MF_EdgeDir_jmp[f->repType])(f,e);
  }

  int MF_EdgeDir_i(MFace_ptr f, int i) {
    return (*MF_EdgeDir_i_jmp[f->repType])(f,i);
  }

  List_ptr MF_AdjFaces(MFace_ptr f) {
    return (*MF_AdjFaces_jmp[f->repType])(f);
  }

  List_ptr MF_Regions(MFace_ptr f) {
    return (*MF_Regions_jmp[f->repType])(f);
  }

  MRegion_ptr MF_Region(MFace_ptr f, int side) {
    return (*MF_Region_jmp[f->repType])(f,side);
  }

  void MF_Add_AdjFace(MFace_ptr f, int edgnum, MFace_ptr af) {
    (*MF_Add_AdjFace_jmp[f->repType])(f,edgnum,af);
  }

  void MF_Rem_AdjFace(MFace_ptr f, int edgnum, MFace_ptr af) {
    (*MF_Rem_AdjFace_jmp[f->repType])(f,edgnum,af);
  }

  void MF_Add_Region(MFace_ptr f, MRegion_ptr r, int side) {
    (*MF_Add_Region_jmp[f->repType])(f,r,side);
  }

  void MF_Rem_Region(MFace_ptr f, MRegion_ptr r) {
    (*MF_Rem_Region_jmp[f->repType])(f,r);
  }

  void MF_Dummy1(MFace_ptr f) {
    return;
  }

  void MF_Dummy2b(MFace_ptr f, MRegion_ptr r, int side) {
    return;
  }

#ifdef __cplusplus
}
#endif

