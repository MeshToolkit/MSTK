#define _H_MVertex_Private

#include <string.h>
#include "MVertex.h"
#include "MVertex_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  MVertex_ptr MV_New(Mesh_ptr mesh) {
    MVertex_ptr v;
    RepType rtype;
    
    v = (MVertex_ptr) MSTK_malloc(sizeof(MVertex));

    v->id = 0;
    v->marker = 0;
    v->mesh = mesh;
    v->dim  = 0;
    v->gdim = 4; /* Nonsensical value since we don't know what it is */
    v->gid  = 0;
    v->gent = (GEntity_ptr) NULL;
    v->xyz[0] = v->xyz[1] = v->xyz[2] = 0.0;
    v->upadj = NULL;
    v->sameadj = NULL;

    rtype = mesh ? MESH_RepType(mesh) : F1;
    MV_Set_RepType(v,rtype);

    if (mesh) MESH_Add_Vertex(mesh,v);

    return v;
  } 

  void MV_Delete(MVertex_ptr v) {
    MESH_Rem_Vertex(v->mesh,v);
    (*MV_Delete_jmp[v->repType])(v);
  }

  void MV_Set_RepType(MVertex_ptr v, RepType rtype) {
    v->repType = rtype;
    (*MV_Set_RepType_jmp[rtype])(v);
  }

  void MV_Set_Coords(MVertex_ptr v, double *xyz) {
    v->xyz[0] = xyz[0];
    v->xyz[1] = xyz[1];
    v->xyz[2] = xyz[2];
  }

  void MV_Set_GEntity(MVertex_ptr v, GEntity_ptr gent) {
    v->gent = gent;
  }

  void MV_Set_GEntDim(MVertex_ptr v, int gdim) {
    v->gdim = gdim;
  }

  void MV_Set_GEntID(MVertex_ptr v, int gid) {
    v->gid = gid;
  }

  void MV_Set_ID(MVertex_ptr v, int id) {
    v->id = id;
  }

  Mesh_ptr MV_Mesh(MVertex_ptr v) {
    return v->mesh;
  }

  int MV_ID(MVertex_ptr v) {
    return v->id;
  }

  int MV_GEntDim(MVertex_ptr v) {
    return v->gdim;
  }

  GEntity_ptr MV_GEntity(MVertex_ptr v) {
    return v->gent;
  }

  void MV_Coords(MVertex_ptr v, double *xyz) {
    xyz[0] = v->xyz[0]; xyz[1] = v->xyz[1]; xyz[2] = v->xyz[2];
  }   

  int MV_Num_AdjVertices(MVertex_ptr v) {
    return (*MV_Num_Edges_jmp[v->repType])(v);
  }

  int MV_Num_Edges(MVertex_ptr v) {
    return (*MV_Num_Edges_jmp[v->repType])(v);
  }

  int MV_Num_Faces(MVertex_ptr v) {
    return (*MV_Num_Faces_jmp[v->repType])(v);
  }
  
  int MV_Num_Regions(MVertex_ptr v) {
    return (*MV_Num_Regions_jmp[v->repType])(v);
  }

  List_ptr MV_AdjVertices(MVertex_ptr v) {
    return (*MV_AdjVertices_jmp[v->repType])(v);
  }

  List_ptr MV_Edges(MVertex_ptr v) {
    return (*MV_Edges_jmp[v->repType])(v);
  }

  List_ptr MV_Faces(MVertex_ptr v) {
    return (*MV_Faces_jmp[v->repType])(v);
  }

  List_ptr MV_Regions(MVertex_ptr v) {
    return (*MV_Regions_jmp[v->repType])(v);
  }

  void MV_Add_AdjVertex(MVertex_ptr v, MVertex_ptr adjv) {
    (*MV_Add_AdjVertex_jmp[v->repType])(v,adjv);
  }

  void MV_Rem_AdjVertex(MVertex_ptr v, MVertex_ptr adjv) {
    (*MV_Rem_AdjVertex_jmp[v->repType])(v,adjv);
  }

  void MV_Add_Edge(MVertex_ptr v, MEdge_ptr e) {
    (*MV_Add_Edge_jmp[v->repType])(v,e);
  }

  void MV_Add_Face(MVertex_ptr v, MFace_ptr f) {
    (*MV_Add_Face_jmp[v->repType])(v,f);
  }

  void MV_Add_Region(MVertex_ptr v, MRegion_ptr r) {
    (*MV_Add_Region_jmp[v->repType])(v,r);
  }

  void MV_Rem_Edge(MVertex_ptr v, MEdge_ptr e) {
    (*MV_Rem_Edge_jmp[v->repType])(v,e);
  }

  void MV_Rem_Face(MVertex_ptr v, MFace_ptr f) {
    (*MV_Rem_Face_jmp[v->repType])(v,f);
  }

  void MV_Rem_Region(MVertex_ptr v, MRegion_ptr r) {
    (*MV_Rem_Region_jmp[v->repType])(v,r);
  }

#ifdef __cplusplus
}
#endif
