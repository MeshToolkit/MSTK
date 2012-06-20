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
    RepType RTYPE;
    
    v = (MVertex_ptr) MSTK_malloc(sizeof(MVertex));

    MEnt_Init_CmnData((MEntity_ptr) v);
    MEnt_Set_Mesh((MEntity_ptr) v,mesh);
    MEnt_Set_Dim((MEntity_ptr) v,MVERTEX);
    MEnt_Set_GEntDim((MEntity_ptr) v,4);
    MEnt_Set_GEntID((MEntity_ptr) v,0);

    v->xyz[0] = v->xyz[1] = v->xyz[2] = 0.0;
    v->adj = NULL;

    RTYPE = mesh ? MESH_RepType(mesh) : F1;
    MV_Set_RepType(v,RTYPE);

    if (mesh) MESH_Add_Vertex(mesh,v);

    return v;
  } 

  void MV_Delete(MVertex_ptr v, int keep) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) v);
    Mesh_ptr mesh;

    (*MV_Delete_jmp[RTYPE])(v,keep);

    if (MEnt_Dim((MEntity_ptr) v) != MDELETED) {
      mesh = MEnt_Mesh((MEntity_ptr) v);

#ifdef MSTK_HAVE_MPI
      if (MV_PType(v) == PGHOST)
	MESH_Rem_GhostVertex(mesh,v);
      else
#endif
	MESH_Rem_Vertex(mesh,v);

      MEnt_Set_DelFlag((MEntity_ptr) v);
    }

    if (!keep) {
      MEnt_Free_CmnData((MEntity_ptr) v);
      MSTK_free(v);
    }
  }

  void MV_Restore(MVertex_ptr v) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) v);
    Mesh_ptr mesh = MEnt_Mesh((MEntity_ptr) v);
    
#ifdef DEBUG
    if (MEnt_Dim((MEntity_ptr) v) != MDELETED) {
      MSTK_Report("MV_Restore",
		  "Trying to restore vertex that is not deleted",MSTK_WARN);
      return;
    }
#endif

    MEnt_Rem_DelFlag((MEntity_ptr) v);

#ifdef MSTK_HAVE_MPI
    if (ME_PType(v) == PGHOST)
      MESH_Add_GhostVertex(mesh,v);
    else
#endif
      MESH_Add_Vertex(mesh,v);

    (*MV_Restore_jmp[RTYPE])(v);
  }

  void MV_Destroy_For_MESH_Delete(MVertex_ptr v) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) v);

    (*MV_Destroy_For_MESH_Delete_jmp[RTYPE])(v);

    MEnt_Free_CmnData((MEntity_ptr) v);

    MSTK_free(v);
  }

  void MV_Set_RepType(MVertex_ptr v, RepType RTYPE) {
    MEnt_Set_RepType_Data((MEntity_ptr) v,RTYPE);
    (*MV_Set_RepType_jmp[RTYPE])(v);
  }

  void MV_Set_Coords(MVertex_ptr v, double *xyz) {
    v->xyz[0] = xyz[0];
    v->xyz[1] = xyz[1];
    v->xyz[2] = xyz[2];
  }

  void MV_Set_GEntity(MVertex_ptr v, GEntity_ptr gent) {    
  }

  void MV_Set_GEntDim(MVertex_ptr v, int gdim) {
    MEnt_Set_GEntDim((MEntity_ptr) v,gdim);
  }

  void MV_Set_GEntID(MVertex_ptr v, int gid) {
    MEnt_Set_GEntID((MEntity_ptr) v,gid);
  }

  void MV_Set_ID(MVertex_ptr v, int id) {
    MEnt_Set_ID((MEntity_ptr) v,id);
  }

  Mesh_ptr MV_Mesh(MVertex_ptr v) {
    return MEnt_Mesh((MEntity_ptr) v);
  }

  int MV_ID(MVertex_ptr v) {
    return MEnt_ID((MEntity_ptr) v);
  }

  int MV_GEntDim(MVertex_ptr v) {
    return MEnt_GEntDim((MEntity_ptr) v);
  }

  int MV_GEntID(MVertex_ptr v) {
    return MEnt_GEntID((MEntity_ptr) v);
  }

  GEntity_ptr MV_GEntity(MVertex_ptr v) {
    return MEnt_GEntity((MEntity_ptr) v);
  }

  void MV_Coords(MVertex_ptr v, double *xyz) {
    xyz[0] = v->xyz[0]; xyz[1] = v->xyz[1]; xyz[2] = v->xyz[2];
  }   

  int MV_Num_AdjVertices(MVertex_ptr v) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) v);
    return (*MV_Num_Edges_jmp[RTYPE])(v);
  }

  int MV_Num_Edges(MVertex_ptr v) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) v);
    return (*MV_Num_Edges_jmp[RTYPE])(v);
  }

  int MV_Num_Faces(MVertex_ptr v) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) v);
    return (*MV_Num_Faces_jmp[RTYPE])(v);
  }
  
  int MV_Num_Regions(MVertex_ptr v) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) v);
    return (*MV_Num_Regions_jmp[RTYPE])(v);
  }

  List_ptr MV_AdjVertices(MVertex_ptr v) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) v);
    return (*MV_AdjVertices_jmp[RTYPE])(v);
  }

  List_ptr MV_Edges(MVertex_ptr v) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) v);
    return (*MV_Edges_jmp[RTYPE])(v);
  }

  List_ptr MV_Faces(MVertex_ptr v) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) v);
    return (*MV_Faces_jmp[RTYPE])(v);
  }

  List_ptr MV_Regions(MVertex_ptr v) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) v);
    return (*MV_Regions_jmp[RTYPE])(v);
  }

  void MV_Add_AdjVertex(MVertex_ptr v, MVertex_ptr adjv) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) v);
    (*MV_Add_AdjVertex_jmp[RTYPE])(v,adjv);
  }

  void MV_Rem_AdjVertex(MVertex_ptr v, MVertex_ptr adjv) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) v);
    (*MV_Rem_AdjVertex_jmp[RTYPE])(v,adjv);
  }

  void MV_Add_Edge(MVertex_ptr v, MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) v);
    (*MV_Add_Edge_jmp[RTYPE])(v,e);
  }

  void MV_Add_Face(MVertex_ptr v, MFace_ptr f) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) v);
    (*MV_Add_Face_jmp[RTYPE])(v,f);
  }

  void MV_Add_Region(MVertex_ptr v, MRegion_ptr r) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) v);
    (*MV_Add_Region_jmp[RTYPE])(v,r);
  }

  void MV_Rem_Edge(MVertex_ptr v, MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) v);
    (*MV_Rem_Edge_jmp[RTYPE])(v,e);
  }

  void MV_Rem_Face(MVertex_ptr v, MFace_ptr f) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) v);
    (*MV_Rem_Face_jmp[RTYPE])(v,f);
  }

  void MV_Rem_Region(MVertex_ptr v, MRegion_ptr r) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) v);
    (*MV_Rem_Region_jmp[RTYPE])(v,r);
  }


#ifdef MSTK_HAVE_MPI

  PType MV_PType(MVertex_ptr v) {
    return MEnt_PType((MEntity_ptr) v);
  }

  void  MV_Set_PType(MVertex_ptr v, PType ptype) {
    MEnt_Set_PType((MEntity_ptr) v ,ptype);
  }

  int   MV_MasterParID(MVertex_ptr v) {
    return MEnt_MasterParID((MEntity_ptr) v ); 
  }

  void  MV_Set_MasterParID(MVertex_ptr v, int masterparid) {
    MEnt_Set_MasterParID((MEntity_ptr) v, masterparid) ;
  }

  int   MV_GlobalID(MVertex_ptr v) {
    return   MEnt_GlobalID((MEntity_ptr) v);
  }

  void  MV_Set_GlobalID(MVertex_ptr v, int globalid) {
    MEnt_Set_GlobalID((MEntity_ptr) v, globalid);
  }

  MVertex_ptr MV_GhostNew(Mesh_ptr mesh) {
    MVertex_ptr v;
    RepType RTYPE;
    v = (MVertex_ptr) MSTK_malloc(sizeof(MVertex));

    MEnt_Init_CmnData((MEntity_ptr) v);
    MEnt_Set_Mesh((MEntity_ptr) v,mesh);
    MEnt_Set_Dim((MEntity_ptr) v,0);
    MEnt_Set_GEntDim((MEntity_ptr) v,4);
    MEnt_Set_GEntID((MEntity_ptr) v,0);

    v->xyz[0] = v->xyz[1] = v->xyz[2] = 0.0;
    v->adj = NULL;

    RTYPE = mesh ? MESH_RepType(mesh) : F1;
    MV_Set_RepType(v,RTYPE);

    if (mesh) {
	MESH_Add_GhostVertex(mesh,v);
    }

    return v;
  } 
#endif  /* MSTK_HAVE_MPI */

#ifdef __cplusplus
}
#endif
