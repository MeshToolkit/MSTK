#define _H_MEdge_Private

#include "MEdge.h"
#include "MEdge_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  MEdge_ptr ME_New(Mesh_ptr mesh) {
    MEdge_ptr e;
    RepType rtype;

    e = (MEdge_ptr) MSTK_malloc(sizeof(MEdge));

    e->id = 0;
    e->marker = 0;
    e->mesh = mesh;
    e->dim = 1;
    e->gdim = 4; /* Nonsensical value since we don't know what it is */
    e->gid = 0;
    e->gent = (GEntity_ptr) NULL;
    e->AttInsList = NULL;
    e->upadj = (void *) NULL;
    e->vertex[0] = e->vertex[1] = (MVertex_ptr) NULL;

    rtype = mesh ? MESH_RepType(mesh) : F1;
    ME_Set_RepType(e,rtype);

    if (mesh) MESH_Add_Edge(mesh,e);

    return e;
  }

  void ME_Delete(MEdge_ptr e, int keep) {
    int idx;
    MAttIns_ptr attins;

    if (e->dim != MDELEDGE)
      MESH_Rem_Edge(e->mesh,e);

    if (!keep) {
      if (e->AttInsList) {
	idx = 0;
	while ((attins = List_Next_Entry(e->AttInsList,&idx)))
	  MAttIns_Delete(attins);
	List_Delete(e->AttInsList);
      }
    }

    (*ME_Delete_jmp[e->repType])(e,keep);
  }

  void ME_Restore(MEdge_ptr e) {
    if (e->dim == MDELEDGE)
      MESH_Add_Edge(e->mesh,e);
    (*ME_Restore_jmp[e->repType])(e);
  }

  void ME_Set_RepType(MEdge_ptr e, RepType rtype) {
    e->repType = rtype;
    (*ME_Set_RepType_jmp[e->repType])(e);
  }

  void ME_Set_GEntity(MEdge_ptr e, GEntity_ptr gent) {
    e->gent = gent;
  }

  void ME_Set_GEntDim(MEdge_ptr e, int gdim) {
    e->gdim = gdim;
  }

  void ME_Set_GEntID(MEdge_ptr e, int gid) {
    e->gid = gid;
  }

  void ME_Set_ID(MEdge_ptr e, int id) {
    e->id = id;
  }

  void ME_Set_Vertex(MEdge_ptr e, int i, MVertex_ptr v) {
    e->vertex[i] = v;
    MV_Add_Edge(v,e);
  }

  void ME_Replace_Vertex(MEdge_ptr e, MVertex_ptr oldv, MVertex_ptr nuv) {
    if (e->vertex[0] == oldv)
      e->vertex[0] = nuv;
    else if (e->vertex[1] == oldv)
      e->vertex[1] = nuv;
    else
      MSTK_Report("ME_Replace_Vertex","Cannot find vertex in edge",ERROR);
  }

  Mesh_ptr ME_Mesh(MEdge_ptr e) {
    return e->mesh;
  }

  int ME_ID(MEdge_ptr e) {
    return e->id;
  }

  int ME_GEntDim(MEdge_ptr e) {
    return e->gdim;
  }

  int ME_GEntID(MEdge_ptr e) {
    return e->gid;
  }

  GEntity_ptr ME_GEntity(MEdge_ptr e) {
    return e->gent;
  }

  MVertex_ptr ME_Vertex(MEdge_ptr e, int i) {
    return e->vertex[i];
  }

  MVertex_ptr ME_OppVertex(MEdge_ptr e, MVertex_ptr v) {
    return (e->vertex[0] == v ? e->vertex[1] : e->vertex[0]);
  }

  int ME_Num_Faces(MEdge_ptr e) {
    return (*ME_Num_Faces_jmp[e->repType])(e);
  }

  int ME_Num_Regions(MEdge_ptr e) {
    return (*ME_Num_Regions_jmp[e->repType])(e);
  }

  List_ptr ME_Faces(MEdge_ptr e) {
    return (*ME_Faces_jmp[e->repType])(e);
  }

  List_ptr ME_Regions(MEdge_ptr e) {
    return (*ME_Regions_jmp[e->repType])(e);
  }

  int ME_UsesEntity(MEdge_ptr e, MEntity_ptr ent, int etype) {
    switch (etype) {
    case 0:
      if (e->vertex[0] == (MVertex_ptr) ent || 
	  e->vertex[1] == (MVertex_ptr) ent)
	return 1;
      else
	return 0;
    case 1:
      return (e == (MEdge_ptr) ent);
    case 2: case 3:
      return 0;
    default:
      MSTK_Report("ME_UsesEntity","Invalid entity type",ERROR);
    }
    return 0;
  }

  void ME_Add_Face(MEdge_ptr e, MFace_ptr f) {
    (*ME_Add_Face_jmp[e->repType])(e,f);
  }

  void ME_Rem_Face(MEdge_ptr e, MFace_ptr f) {
    (*ME_Rem_Face_jmp[e->repType])(e,f);
  }

  void ME_Add_Region(MEdge_ptr e, MRegion_ptr r) {
    (*ME_Add_Region_jmp[e->repType])(e,r);
  }

  void ME_Rem_Region(MEdge_ptr e, MRegion_ptr r) {
    (*ME_Rem_Region_jmp[e->repType])(e,r);
  }



  MEdge_ptr MVs_CommonEdge(MVertex_ptr v1, MVertex_ptr v2) {
    List_ptr vedges;
    MEdge_ptr edge=0;
    MVertex_ptr ev;
    int i, found, nve;

    vedges = MV_Edges(v1);
    if (!vedges)
      return 0;

    nve = List_Num_Entries(vedges);

    for (i = 0, found = 0; i < nve; i++) {
      edge = List_Entry(vedges,i);
      ev = ME_Vertex(edge,0);
      if (ev == v2) {
        found = 1;
        break;
      }
      else {
        ev = ME_Vertex(edge,1);
        if (ev == v2) {
          found = 1;
          break;
        }
      }
    }
    List_Delete(vedges);

    return (found ? edge : 0);
  }

  void ME_Dummy1(MEdge_ptr e) {
    return;
  }

  void ME_Dummy2(MEdge_ptr e, int i) {
    return;
  }

#ifdef __cplusplus
}
#endif

