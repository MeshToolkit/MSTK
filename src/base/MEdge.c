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

    MEnt_Init_CmnData((MEntity_ptr) e);
    MEnt_Set_Mesh((MEntity_ptr) e,mesh);
    MEnt_Set_Dim((MEntity_ptr) e,1);
    MEnt_Set_GEntDim((MEntity_ptr) e,4); /* Nonsensical value as we don't 
                                            know what it is */
    MEnt_Set_GEntID((MEntity_ptr) e,0);

    e->adj = (void *) NULL;
    e->vertex[0] = e->vertex[1] = (MVertex_ptr) NULL;

    RTYPE = mesh ? MESH_RepType(mesh) : F1;
    ME_Set_RepType(e,RTYPE);

    if (mesh) MESH_Add_Edge(mesh,e);

    return e;
  }

  void ME_Delete(MEdge_ptr e, int keep) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) e);
    Mesh_ptr mesh;

    (*ME_Delete_jmp[RTYPE])(e,keep);

    if (MEnt_Dim((MEntity_ptr) e) != MDELETED) {
      mesh = MEnt_Mesh((MEntity_ptr) e);

#ifdef MSTK_HAVE_MPI
      if (ME_PType(e) == PGHOST)
	MESH_Rem_GhostEdge(mesh,e);
      else
#endif
	MESH_Rem_Edge(mesh,e);

      MEnt_Set_DelFlag((MEntity_ptr) e);
    }

    if (!keep) {
      MEnt_Free_CmnData((MEntity_ptr) e);
      MSTK_free((MEntity_ptr) e);
    }
  }

  void ME_Restore(MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) e);
    Mesh_ptr mesh = MEnt_Mesh((MEntity_ptr) e);

#ifdef DEBUG
    if (MEnt_Dim((MEntity_ptr) e) != MDELETED) {
      MSTK_Report("ME_Restore",
		  "Trying to restore edge that is not deleted",MSTK_WARN);
      return;
    }
#endif

    MEnt_Rem_DelFlag((MEntity_ptr) e);

#ifdef MSTK_HAVE_MPI
    if (ME_PType(e) == PGHOST)
      MESH_Add_GhostEdge(mesh,e);
    else
#endif
      MESH_Add_Edge(mesh,e);

    (*ME_Restore_jmp[RTYPE])(e);
  }

  void ME_Destroy_For_MESH_Delete(MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) e);

    (*ME_Destroy_For_MESH_Delete_jmp[RTYPE])(e);

    MEnt_Free_CmnData((MEntity_ptr) e);

    MSTK_free(e);
  }

  void ME_Set_RepType(MEdge_ptr e, RepType RTYPE) {
    MEnt_Set_RepType_Data((MEntity_ptr) e,RTYPE);
    (*ME_Set_RepType_jmp[RTYPE])(e);
  }

  void ME_Set_GEntity(MEdge_ptr e, GEntity_ptr gent) {
  }

  void ME_Set_GEntDim(MEdge_ptr e, int gdim) {
    MEnt_Set_GEntDim((MEntity_ptr) e,gdim);
  }

  void ME_Set_GEntID(MEdge_ptr e, int gid) {
    MEnt_Set_GEntID((MEntity_ptr) e,gid);
  }

  void ME_Set_ID(MEdge_ptr e, int id) {
    MEnt_Set_ID((MEntity_ptr) e,id);
  }

  void ME_Set_Vertex(MEdge_ptr e, int i, MVertex_ptr v) {
    e->vertex[i] = v;
    MV_Add_Edge(v,e);
  }

  void ME_Replace_Vertex(MEdge_ptr e, MVertex_ptr oldv, MVertex_ptr nuv) {
    if (e->vertex[0] == oldv) {
      e->vertex[0] = nuv;
      MV_Add_Edge(nuv,e);
      MV_Rem_Edge(oldv,e);
    }
    else if (e->vertex[1] == oldv) {
      e->vertex[1] = nuv;
      MV_Add_Edge(nuv,e);
      MV_Rem_Edge(oldv,e);
    }
    else
      MSTK_Report("ME_Replace_Vertex","Cannot find vertex in edge",MSTK_ERROR);
  }

  int ME_Set_GInfo_Auto(MEdge_ptr e) {
    MVertex_ptr v0, v1;
    int gdim0, gdim1, gid0, gid1, egdim, egid;

    v0 = e->vertex[0];
    gdim0 = MV_GEntDim(v0); gid0 = MV_GEntID(v0);
    v1 = e->vertex[1];
    gdim1 = MV_GEntDim(v1); gid1 = MV_GEntID(v1);

    egdim = 4;
    egid = -1;

    if (gdim0 == gdim1) {
      if (gid0 == gid1) {
	egdim = gdim0;
	egid = gid0;
      }
      else { /* Unknown classification */
	egdim = 4;
	egid = -1;
      }
    }
    else {
      if (gdim0 > gdim1) {
	egdim = gdim0;
	egid = gid0;
      }
      else {
	egdim = gdim1;
	egid = gid1;
      }
    }

    /* if the classification could not be determined, look at
       higher-dimensional entities connected to the edge to see if we
       can at least find the dimension of the model entity on which it
       is classified */

    if (egdim == 4) {
      List_ptr efaces;
      MFace_ptr eface;

      efaces = ME_Faces(e);
      if (!efaces) {
	egdim = 2;
	egid = -1;
      }
      else {
	int allinterior, allsame, idx, gfid0;

	/* First check if all faces connected to edge are interior faces */

	idx = 0;
	allinterior = 1;
	while ((eface = List_Next_Entry(efaces,&idx))) {
	  if (MEnt_GEntDim((MEntity_ptr) eface) != 3) {
	    allinterior = 0;
	    break;
	  }
	}

	if (allinterior) {
	  egdim = 3;
	  egid = MEnt_GEntID((MEntity_ptr) List_Entry(efaces,0));
	}
	else {

	  /* edge has some boundary faces connected to it. If all the
	     boundary faces are classified on the same model face,
	     then the edge is on that model face. If boundary faces
	     are on different model faces, the edge is on a model
	     edge, although we don't know which one */

	  gfid0 = -1;

	  allsame = 1;
	  idx = 0;
	  while ((eface = List_Next_Entry(efaces,&idx))) {
	    if (MEnt_GEntDim((MEntity_ptr) eface) != 2) continue;

	    if (gfid0 == -1) { /* first boundary face */
	      gfid0 = MEnt_GEntID((MEntity_ptr) eface);
	    }
	    else {
	      if (MEnt_GEntID((MEntity_ptr) eface) != gfid0) {
		allsame = 0;
		break;
	      }
	    }
	  }

	  if (allsame && gfid0 != -1) {
	    egdim = 2;
	    egid = gfid0;
	  }
	  else {
	    egdim = 1;
	    egid = -1;
	  }

	}

	List_Delete(efaces);
      }

    }

    MEnt_Set_GEntDim((MEntity_ptr) e,egdim);
    MEnt_Set_GEntID((MEntity_ptr) e,egid);

    if (egdim == 4)
      return 0;
    else
      return 1;

  }

  Mesh_ptr ME_Mesh(MEdge_ptr e) {
    return MEnt_Mesh((MEntity_ptr) e);
  }

  int ME_ID(MEdge_ptr e) {
    return MEnt_ID((MEntity_ptr) e);
  }

  int ME_GEntDim(MEdge_ptr e) {
    return MEnt_GEntDim((MEntity_ptr) e);
  }

  int ME_GEntID(MEdge_ptr e) {
    return MEnt_GEntID((MEntity_ptr) e);
  }

  GEntity_ptr ME_GEntity(MEdge_ptr e) {
    return MEnt_GEntity((MEntity_ptr) e);
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
    RepType RTYPE = MEnt_RepType((MEntity_ptr) e);
    return (*ME_Num_Faces_jmp[RTYPE])(e);
  }

  int ME_Num_Regions(MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) e);
    return (*ME_Num_Regions_jmp[RTYPE])(e);
  }

  List_ptr ME_Faces(MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) e);
    return (*ME_Faces_jmp[RTYPE])(e);
  }

  List_ptr ME_Regions(MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) e);
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
      MSTK_Report("ME_UsesEntity","Invalid entity type",MSTK_ERROR);
    }
    return 0;
  }

  void ME_Add_Face(MEdge_ptr e, MFace_ptr f) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) e);
    (*ME_Add_Face_jmp[RTYPE])(e,f);
  }

  void ME_Rem_Face(MEdge_ptr e, MFace_ptr f) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) e);
    (*ME_Rem_Face_jmp[RTYPE])(e,f);
  }

  void ME_Add_Region(MEdge_ptr e, MRegion_ptr r) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) e);
    (*ME_Add_Region_jmp[RTYPE])(e,r);
  }

  void ME_Rem_Region(MEdge_ptr e, MRegion_ptr r) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) e);
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

  /* Extra functionality for hash-tables */
  
  MEdge_ptr ME_NextInHash(MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) e);
    return ME_NextInHash_jmp[RTYPE](e);
  }
  
  void ME_Set_NextInHash(MEdge_ptr e, MEdge_ptr next) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) e);
    return ME_Set_NextInHash_jmp[RTYPE](e, next);
  }

  void ME_HashKey(MEdge_ptr e, unsigned int *pn, void* **pp) {
    *pn = 2;
    *pp = e->vertex;
  }
  
  void ME_Lock(MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) e);
    ME_Lock_jmp[RTYPE](e);
  }
  
  void ME_UnLock(MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) e);
    ME_UnLock_jmp[RTYPE](e);
  }
  
  int ME_IsLocked(MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) e);
    return ME_IsLocked_jmp[RTYPE](e);
  }

#ifdef MSTK_HAVE_MPI

  PType ME_PType(MEdge_ptr e) {
    return MEnt_PType((MEntity_ptr) e);
  }

  void  ME_Set_PType(MEdge_ptr e, PType ptype) {
    MEnt_Set_PType((MEntity_ptr) e, ptype);
  }

  int   ME_MasterParID(MEdge_ptr e) {
    return MEnt_MasterParID((MEntity_ptr) e ); 
  }

  /* Rename to ME_Set_MasterPartID? */

  void  ME_Set_MasterParID(MEdge_ptr e, int masterparid) {
    MEnt_Set_MasterParID((MEntity_ptr) e, masterparid) ;
  }

  int   ME_GlobalID(MEdge_ptr e) {
    return   MEnt_GlobalID((MEntity_ptr) e);
  }

  void  ME_Set_GlobalID(MEdge_ptr e, int globalid) {
    MEnt_Set_GlobalID((MEntity_ptr) e, globalid);
  }


  MEdge_ptr ME_GhostNew(Mesh_ptr mesh) {
    MEdge_ptr e;
    RepType RTYPE;
    e = (MEdge_ptr) MSTK_malloc(sizeof(MEdge));

    MEnt_Init_CmnData((MEntity_ptr) e);
    MEnt_Set_Mesh((MEntity_ptr) e,mesh);
    MEnt_Set_Dim((MEntity_ptr) e,1);
    MEnt_Set_GEntDim((MEntity_ptr) e,4); /* Nonsensical value as we don't 
                                            know what it is */
    MEnt_Set_GEntID((MEntity_ptr) e,0);

    e->adj = (void *) NULL;
    e->vertex[0] = e->vertex[1] = (MVertex_ptr) NULL;

    RTYPE = mesh ? MESH_RepType(mesh) : F1;
    ME_Set_RepType(e,RTYPE);

    if (mesh) {
	MESH_Add_GhostEdge(mesh,e);
    }
    return e;
  }

#endif

  
#ifdef __cplusplus
}
#endif

