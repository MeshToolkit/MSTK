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
    RepType RTYPE;

    e = (MEdge_ptr) MSTK_malloc(sizeof(MEdge));

    MEnt_Init_CmnData(e);
    MEnt_Set_Mesh(e,mesh);
    MEnt_Set_Dim(e,1);
    MEnt_Set_GEntDim(e,4); /* Nonsensical value as we don't know what it is */
    MEnt_Set_GEntID(e,0);

    e->upadj = (void *) NULL;
    e->vertex[0] = e->vertex[1] = (MVertex_ptr) NULL;

    RTYPE = mesh ? MESH_RepType(mesh) : F1;
    ME_Set_RepType(e,RTYPE);

    if (mesh) MESH_Add_Edge(mesh,e);

    return e;
  }

  void ME_Delete(MEdge_ptr e, int keep) {
    RepType RTYPE = MEnt_RepType(e);
    Mesh_ptr mesh;

    (*ME_Delete_jmp[RTYPE])(e,keep);

    if (MEnt_Dim(e) != MDELETED) {
      mesh = MEnt_Mesh(e);
      MESH_Rem_Edge(mesh,e);
      MEnt_Set_DelFlag(e);
    }

    if (!keep) {
      MEnt_Free_CmnData(e);
      MSTK_free(e);
    }
  }

  void ME_Restore(MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType(e);
    Mesh_ptr mesh = MEnt_Mesh(e);

#ifdef DEBUG
    if (MEnt_Dim(e) != MDELETED) {
      MSTK_Report("ME_Restore",
		  "Trying to restore edge that is not deleted",WARN);
      return;
    }
#endif

    MEnt_Rem_DelFlag(e);

    MESH_Add_Edge(mesh,e);

    (*ME_Restore_jmp[RTYPE])(e);
  }

  void ME_Destroy_For_MESH_Delete(MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType(e);

    (*ME_Destroy_For_MESH_Delete_jmp[RTYPE])(e);

    MEnt_Free_CmnData(e);

    MSTK_free(e);
  }

  void ME_Set_RepType(MEdge_ptr e, RepType RTYPE) {
    MEnt_Set_RepType_Data(e,RTYPE);
    (*ME_Set_RepType_jmp[RTYPE])(e);
  }

  void ME_Set_GEntity(MEdge_ptr e, GEntity_ptr gent) {
  }

  void ME_Set_GEntDim(MEdge_ptr e, int gdim) {
    MEnt_Set_GEntDim(e,gdim);
  }

  void ME_Set_GEntID(MEdge_ptr e, int gid) {
    MEnt_Set_GEntID(e,gid);
  }

  void ME_Set_ID(MEdge_ptr e, int id) {
    MEnt_Set_ID(e,id);
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

  int ME_Set_GInfo_Auto(MEdge_ptr e) {
    MVertex_ptr v0, v1;
    int gdim0, gdim1, gid0, gid1;

    v0 = e->vertex[0];
    gdim0 = MV_GEntDim(v0); gid0 = MV_GEntID(v0);
    v1 = e->vertex[1];
    gdim1 = MV_GEntDim(v1); gid1 = MV_GEntID(v1);

    if (gdim0 == gdim1) {
      if (gid0 == gid1) {
	MEnt_Set_GEntDim(e,gdim0);
	MEnt_Set_GEntID(e,gid0);
	return 1;
      }
      else { /* Unknown classification */
	MEnt_Set_GEntDim(e,4);
	MEnt_Set_GEntID(e,-1);
	return 0;
      }
    }
    else {
      if (gdim0 > gdim1) {
	MEnt_Set_GEntDim(e,gdim0);
	MEnt_Set_GEntID(e,gid0);
      }
      else {
	MEnt_Set_GEntDim(e,gdim1);
	MEnt_Set_GEntID(e,gid1);
      }
      return 1;
    }
  }

  Mesh_ptr ME_Mesh(MEdge_ptr e) {
    return MEnt_Mesh(e);
  }

  int ME_ID(MEdge_ptr e) {
    return MEnt_ID(e);
  }

  int ME_GEntDim(MEdge_ptr e) {
    return MEnt_GEntDim(e);
  }

  int ME_GEntID(MEdge_ptr e) {
    return MEnt_GEntID(e);
  }

  GEntity_ptr ME_GEntity(MEdge_ptr e) {
    return MEnt_GEntity(e);
  }

  int MEs_AreSame(MEdge_ptr e1, MEdge_ptr e2) {
    if (e1 == e2)
      return 1;
    else {
      if ((e1->vertex[0] == e2->vertex[0] && e1->vertex[1] == e2->vertex[1]) ||
	  (e1->vertex[1] == e2->vertex[0] && e1->vertex[0] == e2->vertex[1]))
	return 1;
      else
	return 0;
    }
  }

  MVertex_ptr ME_Vertex(MEdge_ptr e, int i) {
    return e->vertex[i];
  }

  MVertex_ptr ME_OppVertex(MEdge_ptr e, MVertex_ptr v) {
    return (e->vertex[0] == v ? e->vertex[1] : e->vertex[0]);
  }

  int ME_Num_Faces(MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType(e);
    return (*ME_Num_Faces_jmp[RTYPE])(e);
  }

  int ME_Num_Regions(MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType(e);
    return (*ME_Num_Regions_jmp[RTYPE])(e);
  }

  List_ptr ME_Faces(MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType(e);
    return (*ME_Faces_jmp[RTYPE])(e);
  }

  List_ptr ME_Regions(MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType(e);
    return (*ME_Regions_jmp[RTYPE])(e);
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
    RepType RTYPE = MEnt_RepType(e);
    (*ME_Add_Face_jmp[RTYPE])(e,f);
  }

  void ME_Rem_Face(MEdge_ptr e, MFace_ptr f) {
    RepType RTYPE = MEnt_RepType(e);
    (*ME_Rem_Face_jmp[RTYPE])(e,f);
  }

  void ME_Add_Region(MEdge_ptr e, MRegion_ptr r) {
    RepType RTYPE = MEnt_RepType(e);
    (*ME_Add_Region_jmp[RTYPE])(e,r);
  }

  void ME_Rem_Region(MEdge_ptr e, MRegion_ptr r) {
    RepType RTYPE = MEnt_RepType(e);
    (*ME_Rem_Region_jmp[RTYPE])(e,r);
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

