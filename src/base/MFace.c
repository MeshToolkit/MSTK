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
    RepType RTYPE;

    f = (MFace_ptr) MSTK_malloc(sizeof(MFace));

    MEnt_Init_CmnData(f);
    MEnt_Set_Mesh(f,mesh);
    MEnt_Set_Dim(f,2);
    MEnt_Set_GEntDim(f,4); /* nonsensical value as we don't know what it is */
    MEnt_Set_GEntID(f,0);

    f->upadj = (void *) NULL;
    f->sameadj = (void *) NULL;
    f->downadj = (void *) NULL;

    RTYPE = mesh ? MESH_RepType(mesh) : F1;
    MF_Set_RepType(f,RTYPE);
    
    if (mesh) MESH_Add_Face(mesh,f);

    return f;
  }

  void MF_Delete(MFace_ptr f, int keep) {
    RepType RTYPE = MEnt_RepType(f);
    Mesh_ptr mesh;

    (*MF_Delete_jmp[RTYPE])(f, keep);

    if (MEnt_Dim(f) != MDELETED) {
      mesh = MEnt_Mesh(f);
      MESH_Rem_Face(mesh,f);
      MEnt_Set_DelFlag(f);
    }

    if (!keep) {
      MEnt_Free_CmnData(f);
      MSTK_free(f);
    }
  }

  void MF_Restore(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType(f);
    Mesh_ptr mesh = MEnt_Mesh(f);

#ifdef DEBUG
    if (MEnt_Dim(f) != MDELETED) {
      MSTK_Report("MF_Restore",
		  "Trying to restore face that is not deleted",WARN);
      return;
    }
#endif

    MEnt_Rem_DelFlag(f);

    MESH_Add_Face(mesh,f);

    (*MF_Restore_jmp[RTYPE])(f);
  }

  void MF_Destroy_For_MESH_Delete(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType(f);

    (*MF_Destroy_For_MESH_Delete_jmp[RTYPE])(f);

    MEnt_Free_CmnData(f);

    MSTK_free(f);
  }

  void MF_Set_RepType(MFace_ptr f, RepType RTYPE) {
    MEnt_Set_RepType_Data(f,RTYPE);
    (*MF_Set_RepType_jmp[RTYPE])(f);
  }

  void MF_Set_GEntity(MFace_ptr f, GEntity_ptr gent) {
  }

  void MF_Set_GEntDim(MFace_ptr f, int gdim) {
    MEnt_Set_GEntDim(f,gdim);
  }

  void MF_Set_GEntID(MFace_ptr f, int gid) {
    MEnt_Set_GEntID(f,gid);
  }

  int MF_Set_GInfo_Auto(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType(f);
    return (*MF_Set_GInfo_Auto_jmp[RTYPE])(f);
  }


  void MF_Set_ID(MFace_ptr f, int id) {
    MEnt_Set_ID(f,id);
  }

  void MF_Set_Edges(MFace_ptr f, int n, MEdge_ptr *edges, int *dir) {
    RepType RTYPE = MEnt_RepType(f);
    (*MF_Set_Edges_jmp[RTYPE])(f,n,edges,dir);
  }

  void MF_Replace_Edges(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, 
			MEdge_ptr *nuedges) {
    RepType RTYPE = MEnt_RepType(f);
    (*MF_Replace_Edges_jmp[RTYPE])(f,nold,oldedges,nnu,nuedges);
  }

  void MF_Replace_Edges_i(MFace_ptr f, int nold, int i, int nnu, 
			  MEdge_ptr *nuedges) {
    RepType RTYPE = MEnt_RepType(f);
    (*MF_Replace_Edges_i_jmp[RTYPE])(f,nold,i,nnu,nuedges);
  }

  void MF_Set_Vertices(MFace_ptr f, int n, MVertex_ptr *verts) {
    RepType RTYPE = MEnt_RepType(f);
    (*MF_Set_Vertices_jmp[RTYPE])(f,n,verts);
  }

  void MF_Replace_Vertex(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv) {
    RepType RTYPE = MEnt_RepType(f);
    (*MF_Replace_Vertex_jmp[RTYPE])(f,v,nuv);
  }

  void MF_Replace_Vertex_i(MFace_ptr f, int i, MVertex_ptr nuv) {
    RepType RTYPE = MEnt_RepType(f);
    (*MF_Replace_Vertex_i_jmp[RTYPE])(f,i,nuv);
  }

  void MF_Insert_Vertex(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v) {
    RepType RTYPE = MEnt_RepType(f);
    (*MF_Insert_Vertex_jmp[RTYPE])(f,nuv,b4v);
  }

  void MF_Insert_Vertex_i(MFace_ptr f, MVertex_ptr nuv, int i) {
    RepType RTYPE = MEnt_RepType(f);
    (*MF_Insert_Vertex_i_jmp[RTYPE])(f,nuv,i);
  }

  Mesh_ptr MF_Mesh(MFace_ptr f) {
    return MEnt_Mesh(f);
  }

  int MF_ID(MFace_ptr f) {
    return MEnt_ID(f);
  }

  int MF_GEntDim(MFace_ptr f) {
    return MEnt_GEntDim(f);
  }

  int MF_GEntID(MFace_ptr f) {
    return MEnt_GEntID(f);
  }

  GEntity_ptr MF_GEntity(MFace_ptr f) {
    return MEnt_GEntity(f);
  }

  int MFs_AreSame(MFace_ptr f1, MFace_ptr f2) {
    RepType RTYPE = MEnt_RepType(f1);
    return (*MFs_AreSame_jmp[RTYPE])(f1,f2);
  }

  int MF_UsesEntity(MFace_ptr f, MEntity_ptr e, int etype) {
    RepType RTYPE = MEnt_RepType(f);

    switch (etype) {
    case 3:
      return 0;
    case 2:
      return (f == (MFace_ptr) e);
    case 1:
      return (*MF_UsesEdge_jmp[RTYPE])(f,(MEdge_ptr)e);
    case 0:
      return (*MF_UsesVertex_jmp[RTYPE])(f,(MVertex_ptr)e);
    default:
      MSTK_Report("MF_UsesEntity","Invalid entity type",ERROR);
    }
    return 0;
  }

  int MF_Num_Vertices(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType(f);
    return (*MF_Num_Vertices_jmp[RTYPE])(f);
  }

  int MF_Num_Edges(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType(f);
    return (*MF_Num_Edges_jmp[RTYPE])(f);
  }

  List_ptr MF_Vertices(MFace_ptr f, int dir, MVertex_ptr v) {
    RepType RTYPE = MEnt_RepType(f);
    return (*MF_Vertices_jmp[RTYPE])(f,dir, v);
  }

  List_ptr MF_Edges(MFace_ptr f, int dir, MVertex_ptr v) {
    RepType RTYPE = MEnt_RepType(f);
    return (*MF_Edges_jmp[RTYPE])(f,dir,v);
  }

  int MF_Num_AdjFaces(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType(f);
    return (*MF_Num_AdjFaces_jmp[RTYPE])(f);
  }

  int MF_EdgeDir(MFace_ptr f, MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType(f);
    return (*MF_EdgeDir_jmp[RTYPE])(f,e);
  }

  int MF_EdgeDir_i(MFace_ptr f, int i) {
    RepType RTYPE = MEnt_RepType(f);
    return (*MF_EdgeDir_i_jmp[RTYPE])(f,i);
  }

  List_ptr MF_AdjFaces(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType(f);
    return (*MF_AdjFaces_jmp[RTYPE])(f);
  }

  List_ptr MF_Regions(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType(f);
    return (*MF_Regions_jmp[RTYPE])(f);
  }

  MRegion_ptr MF_Region(MFace_ptr f, int side) {
    RepType RTYPE = MEnt_RepType(f);
    return (*MF_Region_jmp[RTYPE])(f,side);
  }

  void MF_Add_AdjFace(MFace_ptr f, int edgnum, MFace_ptr af) {
    RepType RTYPE = MEnt_RepType(f);
    (*MF_Add_AdjFace_jmp[RTYPE])(f,edgnum,af);
  }

  void MF_Rem_AdjFace(MFace_ptr f, int edgnum, MFace_ptr af) {
    RepType RTYPE = MEnt_RepType(f);
    (*MF_Rem_AdjFace_jmp[RTYPE])(f,edgnum,af);
  }

  void MF_Add_Region(MFace_ptr f, MRegion_ptr r, int side) {
    RepType RTYPE = MEnt_RepType(f);
    (*MF_Add_Region_jmp[RTYPE])(f,r,side);
  }

  void MF_Rem_Region(MFace_ptr f, MRegion_ptr r) {
    RepType RTYPE = MEnt_RepType(f);
    (*MF_Rem_Region_jmp[RTYPE])(f,r);
  }

  void MF_Dummy1(MFace_ptr f) {
    return;
  }

  void MF_Dummy2a(MFace_ptr f, int i) {
    return;
  }

  void MF_Dummy2b(MFace_ptr f, MRegion_ptr r, int side) {
    return;
  }

  MFace_ptr MVs_CommonFace(int nv, MVertex_ptr *fverts) {
    MFace_ptr common_face = NULL, vface;
    int i, j;
    int nvf0, contains_all;
    List_ptr vfaces0, fvtxlist;

    vfaces0 = MV_Faces(fverts[0]);
    if (!vfaces0)
      return NULL;

    nvf0 = List_Num_Entries(vfaces0);
    
    for (i = 0; i < nvf0; i++) {
      vface = List_Entry(vfaces0,i);
      fvtxlist = MF_Vertices(vface,1,0);

      contains_all = 1;
      for (j = 1; j < nv; j++) {
	if (!List_Contains(fvtxlist,fverts[j])) {
	  contains_all = 0;
	  break;
	}
      }      
      List_Delete(fvtxlist);

      if (contains_all) {
	common_face = vface;
	break;
      }
    }
    List_Delete(vfaces0);

    return common_face;
  }

  MFace_ptr MEs_CommonFace(int ne, MEdge_ptr *fedges) {
    MFace_ptr common_face = NULL, eface;
    int i, j;
    int nef0, contains_all;
    List_ptr efaces0, fedglist;

    efaces0 = ME_Faces(fedges[0]);
    if (!efaces0)
      return NULL;

    nef0 = List_Num_Entries(efaces0);
    
    for (i = 0; i < nef0; i++) {
      eface = List_Entry(efaces0,i);
      fedglist = MF_Edges(eface,1,0);

      contains_all = 1;
      for (j = 1; j < ne; j++) {
	if (!List_Contains(fedglist,fedges[j])) {
	  contains_all = 0;
	  break;
	}
      }      
      List_Delete(fedglist);

      if (contains_all) {
	common_face = eface;
	break;
      }
    }
    List_Delete(efaces0);

    return common_face;
  }

#ifdef __cplusplus
}
#endif

