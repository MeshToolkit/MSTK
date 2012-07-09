#define _H_MRegion_Private

#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

 int mrtype_nv[6]={0,4,5,6,8,-1};
 int mrtype_ne[6]={0,6,8,9,12,-1};
 int mrtype_nf[6]={0,4,5,5,6,-1};

  MRegion_ptr MR_New(Mesh_ptr mesh) {
    MRegion_ptr r;
    RepType RTYPE;

    r = (MRegion_ptr) MSTK_malloc(sizeof(MRegion));

    MEnt_Init_CmnData((MEntity_ptr) r);
    MEnt_Set_Mesh((MEntity_ptr) r,mesh);
    MEnt_Set_Dim((MEntity_ptr) r,MREGION);
    MEnt_Set_GEntDim((MEntity_ptr) r,3);
    MEnt_Set_GEntID((MEntity_ptr) r,0);

    r->adj = (void *) NULL;
    r->mrtype = RUNKNOWN;

    RTYPE = mesh ? MESH_RepType(mesh) : F1;
    MR_Set_RepType(r,RTYPE);

    if (mesh) MESH_Add_Region(mesh,r);

    return r;
  }

  void MR_Delete(MRegion_ptr r, int keep) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    Mesh_ptr mesh;

    (*MR_Delete_jmp[RTYPE])(r,keep);

    if (MEnt_Dim((MEntity_ptr) r) != MDELETED) {
      mesh = MEnt_Mesh((MEntity_ptr) r);

      if (MR_PType(r) == PGHOST)
	MESH_Rem_GhostRegion(mesh,r);
      else
	MESH_Rem_Region(mesh,r);

      MEnt_Set_DelFlag((MEntity_ptr) r);
    }

    if (!keep) {
      MEnt_Free_CmnData((MEntity_ptr) r);
      MSTK_free(r);
    }
  }

  void MR_Restore(MRegion_ptr r) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    Mesh_ptr mesh = MEnt_Mesh((MEntity_ptr) r);

#ifdef DEBUG
    if (MEnt_Dim((MEntity_ptr) r) != MDELETED) {
      MSTK_Report("MR_Restore",
		  "Trying to restore region that is not deleted",MSTK_WARN);
      return;
    }
#endif

    MEnt_Rem_DelFlag((MEntity_ptr) r);

    if (MR_PType(r) == PGHOST)
      MESH_Add_GhostRegion(mesh,r);
    else
      MESH_Add_Region(mesh,r);

    (*MR_Restore_jmp[RTYPE])(r);
  }


  /* Delete region without worrying about other entities that are
     pointing to it. TO BE USED ONLY BY MESH_DELETE FUNCTION WHERE ALL
     ELEMENTS WILL GET DELETED. BAD THINGS WILL HAPPEN IF USED
     ELSEWHERE !!!! */

  void MR_Destroy_For_MESH_Delete(MRegion_ptr r) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);

    (*MR_Destroy_For_MESH_Delete_jmp[RTYPE])(r);

    MEnt_Free_CmnData((MEntity_ptr) r);

    MSTK_free(r);
  }

  void MR_Set_RepType(MRegion_ptr r, RepType RTYPE) {
    MEnt_Set_RepType_Data((MEntity_ptr) r,RTYPE);
    (*MR_Set_RepType_jmp[RTYPE])(r);
  }

  void MR_Set_GEntity(MRegion_ptr r, GEntity_ptr gent) {
  }

  void MR_Set_GEntDim(MRegion_ptr r, int gdim) {
    MEnt_Set_GEntDim((MEntity_ptr) r,gdim);
  }

  void MR_Set_GEntID(MRegion_ptr r, int gid) {
    MEnt_Set_GEntID((MEntity_ptr) r,gid);
  }

  int MR_Set_GInfo_Auto(MRegion_ptr r) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    return (*MR_Set_GInfo_Auto_jmp[RTYPE])(r);
  }

  void MR_Set_ID(MRegion_ptr r, int id) {
    MEnt_Set_ID((MEntity_ptr) r,id);
  }

  void MR_Set_Faces(MRegion_ptr r, int nf, MFace_ptr *mfaces, int *dirs) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    (*MR_Set_Faces_jmp[RTYPE])(r,nf,mfaces,dirs);
  }

  void MR_Set_Vertices(MRegion_ptr r, int nv, MVertex_ptr *mvertices, int nf, int **template) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    (*MR_Set_Vertices_jmp[RTYPE])(r,nv,mvertices,nf,template);
  }

  Mesh_ptr MR_Mesh(MRegion_ptr r) {
    return MEnt_Mesh((MEntity_ptr) r);
  }

  int MR_ID(MRegion_ptr r) {
    return MEnt_ID((MEntity_ptr) r);
  }

  int MR_GEntDim(MRegion_ptr r) {
    return MEnt_GEntDim((MEntity_ptr) r);
  }

  int MR_GEntID(MRegion_ptr r) {
    return MEnt_GEntID((MEntity_ptr) r);
  }

  GEntity_ptr MR_GEntity(MRegion_ptr r) {
    return NULL;
  }


  MRType MR_ElementType(MRegion_ptr r) {
    return r->mrtype;
  }


  int MR_Num_Vertices(MRegion_ptr r) {
    List_ptr rverts;
    int nrv;

    rverts = MR_Vertices(r);
    nrv = List_Num_Entries(rverts);
    List_Delete(rverts);

#ifdef DEBUG
    {
      RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
      if (RTYPE != R1 || RTYPE != R2)
	MSTK_Report("MR_Num_Vertices",
		    "Inefficient to use this routine for this representation",
		    MSTK_WARN);
    }
#endif

    return nrv;
  }

  int MR_Num_Edges(MRegion_ptr r) {
    List_ptr redges;
    int nre;

    redges = MR_Edges(r);
    nre = List_Num_Entries(redges);
    List_Delete(redges);

#ifdef DEBUG
    {
      RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
      if (RTYPE != R1 || RTYPE != R2)
	MSTK_Report("MR_Num_Vertices",
		    "Inefficient to use this routine for this representation",
		    MSTK_WARN);
    }
#endif

    return nre;
  }

  int MR_Num_Faces(MRegion_ptr r) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    return (*MR_Num_Faces_jmp[RTYPE])(r);
  }

  /*
  int MR_Num_AdjRegions(MRegion_ptr r) {
    return (*MR_Num_AdjRegions_jmp[RTYPE])(r);
  }
  */

  List_ptr MR_Vertices(MRegion_ptr r) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    return (*MR_Vertices_jmp[RTYPE])(r);
  }

  List_ptr MR_Edges(MRegion_ptr r) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    return (*MR_Edges_jmp[RTYPE])(r);
  }

  List_ptr MR_Faces(MRegion_ptr r) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    return (*MR_Faces_jmp[RTYPE])(r);
  }

  List_ptr MR_AdjRegions(MRegion_ptr r) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    return (*MR_AdjRegions_jmp[RTYPE])(r);
  }

  int MR_FaceDir(MRegion_ptr r, MFace_ptr f) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    return (*MR_FaceDir_jmp[RTYPE])(r,f);
  }

  int MR_FaceDir_i(MRegion_ptr r, int i) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    return (*MR_FaceDir_i_jmp[RTYPE])(r,i);
  }

  int MR_Rev_FaceDir(MRegion_ptr r, MFace_ptr f) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    return (*MR_Rev_FaceDir_jmp[RTYPE])(r,f);
  }

  int MR_Rev_FaceDir_i(MRegion_ptr r, int i) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    return (*MR_Rev_FaceDir_i_jmp[RTYPE])(r,i);
  }

  void MR_Rem_Face(MRegion_ptr r, MFace_ptr f) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    (*MR_Rem_Face_jmp[RTYPE])(r,f);
  }

  void MR_Replace_Face(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf, int nudir) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    (*MR_Replace_Face_jmp[RTYPE])(r,f,nuf,nudir);
  }

  void MR_Replace_Vertex(MRegion_ptr r, MVertex_ptr v, MVertex_ptr nuv) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    (*MR_Replace_Vertex_jmp[RTYPE])(r,v,nuv);
  }

  void MR_Replace_Face_i(MRegion_ptr r, int i, MFace_ptr nuf, int nudir) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    (*MR_Replace_Face_i_jmp[RTYPE])(r,i,nuf,nudir);
  }

  void MR_Replace_Vertex_i(MRegion_ptr r, int i, MVertex_ptr nuv) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    (*MR_Replace_Vertex_i_jmp[RTYPE])(r,i,nuv);
  }

  void MR_Add_AdjRegion(MRegion_ptr r, int facenum, MRegion_ptr aregion) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    (*MR_Add_AdjRegion_jmp[RTYPE])(r,facenum,aregion);
  }

  void MR_Rem_AdjRegion(MRegion_ptr r, MRegion_ptr aregion) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    (*MR_Rem_AdjRegion_jmp[RTYPE])(r,aregion);
  }

  int MR_UsesEntity(MRegion_ptr r, MEntity_ptr e, int etype) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) r);
    
    switch (etype) {
    case 3:
      return (r == (MRegion_ptr) e);
    case 2:
      return (*MR_UsesFace_jmp[RTYPE])(r,(MFace_ptr)e);
    case 1:
      return (*MR_UsesEdge_jmp[RTYPE])(r,(MEdge_ptr)e);
    case 0:
      return (*MR_UsesVertex_jmp[RTYPE])(r,(MVertex_ptr)e);
    default:
      MSTK_Report("MR_UsesEntity","Invalid entity type",MSTK_ERROR);
    }
    return 0;
  }


  PType MR_PType(MRegion_ptr r) {
    return MEnt_PType((MEntity_ptr) r);
  }

  void  MR_Set_PType(MRegion_ptr r, PType ptype) {
    MEnt_Set_PType((MEntity_ptr) r, ptype);
  }

  int   MR_MasterParID(MRegion_ptr r) {
    return MEnt_MasterParID((MEntity_ptr) r ); 
  }

  void  MR_Set_MasterParID(MRegion_ptr r, int masterparid) {
    MEnt_Set_MasterParID((MEntity_ptr) r, masterparid) ;
  }

  int   MR_GlobalID(MRegion_ptr r) {
    return   MEnt_GlobalID((MEntity_ptr) r);
  }

  void  MR_Set_GlobalID(MRegion_ptr r, int globalid) {
    MEnt_Set_GlobalID((MEntity_ptr) r, globalid);
  }

  MRegion_ptr MR_GhostNew(Mesh_ptr mesh) {
    MRegion_ptr r;
    RepType RTYPE;
    r = (MRegion_ptr) MSTK_malloc(sizeof(MRegion));

    MEnt_Init_CmnData((MEntity_ptr) r);
    MEnt_Set_Mesh((MEntity_ptr) r,mesh);
    MEnt_Set_Dim((MEntity_ptr) r,3);
    MEnt_Set_GEntDim((MEntity_ptr) r,3);
    MEnt_Set_GEntID((MEntity_ptr) r,0);

    r->adj = (void *) NULL;

    RTYPE = mesh ? MESH_RepType(mesh) : F1;
    MR_Set_RepType(r,RTYPE);

    if (mesh) {
	MESH_Add_GhostRegion(mesh,r);
    }

    return r;
  }

  
#ifdef __cplusplus
}
#endif
