#define _H_MFace_Private

#include <stdlib.h>
#include "MFace.h"
#include "MFace_jmp.h"
#include "MSTK_private.h"


#ifdef __cplusplus
extern "C" {
#endif

  MFace_ptr MF_New(Mesh_ptr mesh) {
    MFace_ptr f;
    RepType RTYPE;

    f = (MFace_ptr) malloc(sizeof(MFace));

    MEnt_Init_CmnData((MEntity_ptr) f);
    MEnt_Set_Mesh((MEntity_ptr) f,mesh);
    MEnt_Set_Dim((MEntity_ptr) f,MFACE);
    MEnt_Set_GEntDim((MEntity_ptr) f,4); /* nonsensical value as we don't know what it is */
    MEnt_Set_GEntID((MEntity_ptr) f,0);

    f->adj = (void *) NULL;

    RTYPE = mesh ? MESH_RepType(mesh) : F1;
    MF_Set_RepType(f,RTYPE);
    
    if (mesh) MESH_Add_Face(mesh,f);

    return f;
  }

  void MF_Delete(MFace_ptr f, int keep) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    Mesh_ptr mesh;

    (*MF_Delete_jmp[RTYPE])(f, keep);

    if (MEnt_Dim((MEntity_ptr) f) != MDELETED) {
      mesh = MEnt_Mesh((MEntity_ptr) f);

      if (mesh) {
#ifdef MSTK_HAVE_MPI
      if (MF_PType(f) == PGHOST) 
	MESH_Rem_GhostFace(mesh,f);
      else
#endif
	MESH_Rem_Face(mesh,f);
      }

      MEnt_Set_DelFlag((MEntity_ptr) f);
    }

    if (!keep) {
      MEnt_Free_CmnData((MEntity_ptr) f);
      free(f);
    }
  }

  void MF_Restore(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    Mesh_ptr mesh = MEnt_Mesh((MEntity_ptr) f);

#ifdef DEBUG
    if (MEnt_Dim((MEntity_ptr) f) != MDELETED) {
      MSTK_Report("MF_Restore",
		  "Trying to restore face that is not deleted",MSTK_WARN);
      return;
    }
#endif

    MEnt_Rem_DelFlag((MEntity_ptr) f);

#ifdef MSTK_HAVE_MPI
    if (MF_PType(f) == PGHOST)
      MESH_Add_GhostFace(mesh,f);
    else
#endif
      MESH_Add_Face(mesh,f);

    (*MF_Restore_jmp[RTYPE])(f);
  }

  void MF_Destroy_For_MESH_Delete(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);

    (*MF_Destroy_For_MESH_Delete_jmp[RTYPE])(f);

    MEnt_Free_CmnData((MEntity_ptr) f);

    free(f);
  }

  void MF_Set_RepType(MFace_ptr f, RepType RTYPE) {
    MEnt_Set_RepType_Data((MEntity_ptr) f,RTYPE);
    (*MF_Set_RepType_jmp[RTYPE])(f);
  }

  void MF_Set_GEntity(MFace_ptr f, GEntity_ptr gent) {
  }

  void MF_Set_GEntDim(MFace_ptr f, int gdim) {
    MEnt_Set_GEntDim((MEntity_ptr) f,gdim);
  }

  void MF_Set_GEntID(MFace_ptr f, int gid) {
    MEnt_Set_GEntID((MEntity_ptr) f,gid);
  }

  int MF_Set_GInfo_Auto(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    return (*MF_Set_GInfo_Auto_jmp[RTYPE])(f);
  }


  void MF_Set_ID(MFace_ptr f, int id) {
    MEnt_Set_ID((MEntity_ptr) f,id);
  }

  void MF_Set_Edges(MFace_ptr f, int n, MEdge_ptr *edges, int *dir) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    (*MF_Set_Edges_jmp[RTYPE])(f,n,edges,dir);
  }

  void MF_Rem_Edge(MFace_ptr f, MEdge_ptr edge) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    (*MF_Rem_Edge_jmp[RTYPE])(f,edge);
  }

  void MF_Replace_Edges(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, 
			MEdge_ptr *nuedges) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    (*MF_Replace_Edges_jmp[RTYPE])(f,nold,oldedges,nnu,nuedges);
  }

  void MF_Replace_Edges_i(MFace_ptr f, int nold, int i, int nnu, 
			  MEdge_ptr *nuedges) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    (*MF_Replace_Edges_i_jmp[RTYPE])(f,nold,i,nnu,nuedges);
  }

  void MF_Set_Vertices(MFace_ptr f, int n, MVertex_ptr *verts) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);

#ifdef DEBUG
    if (n < 3)
      MSTK_Report("MF_Set_Vertices","Trying to create face with only 2 vertices",MSTK_ERROR);
#endif

    (*MF_Set_Vertices_jmp[RTYPE])(f,n,verts);
  }

  void MF_Replace_Vertex(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    (*MF_Replace_Vertex_jmp[RTYPE])(f,v,nuv);
  }

  void MF_Replace_Vertex_i(MFace_ptr f, int i, MVertex_ptr nuv) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    (*MF_Replace_Vertex_i_jmp[RTYPE])(f,i,nuv);
  }

  void MF_Insert_Vertex(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    (*MF_Insert_Vertex_jmp[RTYPE])(f,nuv,b4v);
  }

  void MF_Insert_Vertex_i(MFace_ptr f, MVertex_ptr nuv, int i) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    (*MF_Insert_Vertex_i_jmp[RTYPE])(f,nuv,i);
  }


  int MF_Rev_EdgeDir(MFace_ptr f, MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    return (*MF_Rev_EdgeDir_jmp[RTYPE])(f,e);
  }

  int MF_Rev_EdgeDir_i(MFace_ptr f, int i) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    return (*MF_Rev_EdgeDir_i_jmp[RTYPE])(f,i);
  }


  Mesh_ptr MF_Mesh(MFace_ptr f) {
    return MEnt_Mesh((MEntity_ptr) f);
  }

  int MF_ID(MFace_ptr f) {
    return MEnt_ID((MEntity_ptr) f);
  }

  int MF_GEntDim(MFace_ptr f) {
    return MEnt_GEntDim((MEntity_ptr) f);
  }

  int MF_GEntID(MFace_ptr f) {
    return MEnt_GEntID((MEntity_ptr) f);
  }

  GEntity_ptr MF_GEntity(MFace_ptr f) {
    return MEnt_GEntity((MEntity_ptr) f);
  }

  int MFs_AreSame(MFace_ptr f1, MFace_ptr f2) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f1);
    return (*MFs_AreSame_jmp[RTYPE])(f1,f2);
  }

  int MF_UsesEntity(MFace_ptr f, MEntity_ptr e, int etype) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);

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
      MSTK_Report("MF_UsesEntity","Invalid entity type",MSTK_ERROR);
    }
    return 0;
  }


  MFType MF_ElementType(MFace_ptr f) {
    int nfv = MF_Num_Vertices(f);
    switch (nfv) {
    case 0: case 1: case 2:
      return FUNKNOWN;
      break;
    case 3:
      return TRI;
      break;
    case 4:
      return QUAD;
      break;
    default:
      return POLYGON;
    }

    return RUNKNOWN;
  }

  int MF_Num_Vertices(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    return (*MF_Num_Vertices_jmp[RTYPE])(f);
  }

  int MF_Num_Edges(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    return (*MF_Num_Edges_jmp[RTYPE])(f);
  }

  List_ptr MF_Vertices(MFace_ptr f, int dir, MVertex_ptr v) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    return (*MF_Vertices_jmp[RTYPE])(f,dir, v);
  }

  void MF_VertexIDs(MFace_ptr f, int dir, int v0, int *nfv,
                    int *fvertids) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    (*MF_VertexIDs_jmp[RTYPE])(f,dir,v0,nfv,fvertids);
  }

  List_ptr MF_Edges(MFace_ptr f, int dir, MVertex_ptr v) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    return (*MF_Edges_jmp[RTYPE])(f,dir,v);
  }

  void MF_EdgeIDs(MFace_ptr f, int dir, int v0, int *nfe,
                    int *fedgeids) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    (*MF_EdgeIDs_jmp[RTYPE])(f,dir,v0,nfe,fedgeids);
  }  

  int MF_Num_AdjFaces(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    return (*MF_Num_AdjFaces_jmp[RTYPE])(f);
  }

  int MF_EdgeDir(MFace_ptr f, MEdge_ptr e) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    return (*MF_EdgeDir_jmp[RTYPE])(f,e);
  }

  int MF_EdgeDir_i(MFace_ptr f, int i) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    return (*MF_EdgeDir_i_jmp[RTYPE])(f,i);
  }

  List_ptr MF_AdjFaces(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    return (*MF_AdjFaces_jmp[RTYPE])(f);
  }

  List_ptr MF_Regions(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    return (*MF_Regions_jmp[RTYPE])(f);
  }

  void MF_RegionIDs(MFace_ptr f, int *nfr, int *fregionids) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    (*MF_RegionIDs_jmp[RTYPE])(f,nfr,fregionids);
  }

  MRegion_ptr MF_Region(MFace_ptr f, int side) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    return (*MF_Region_jmp[RTYPE])(f,side);
  }

  int MF_RegionID(MFace_ptr f, int side) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    return (*MF_RegionID_jmp[RTYPE])(f,side);
  }

  void MF_Add_AdjFace(MFace_ptr f, int edgnum, MFace_ptr af) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    (*MF_Add_AdjFace_jmp[RTYPE])(f,edgnum,af);
  }

  void MF_Rem_AdjFace(MFace_ptr f, int edgnum, MFace_ptr af) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    (*MF_Rem_AdjFace_jmp[RTYPE])(f,edgnum,af);
  }

  void MF_Add_Region(MFace_ptr f, MRegion_ptr r, int side) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    (*MF_Add_Region_jmp[RTYPE])(f,r,side);
  }

  void MF_Rem_Region(MFace_ptr f, MRegion_ptr r) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
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

  
  /* NOTE: MVs_CommonFace can give a positive result if a face has all
     'nv' vertices but has other vertices not in the list */
  /* The implementation also assumes that fverts are given in order
     around the face regardless of direction - Otherwise this routine
     can be inefficient */

  MFace_ptr MVs_CommonFace(int nv, MVertex_ptr *fverts) {
    MFace_ptr common_face = NULL, eface;
    int i, j;
    int nef0, contains_all;
    List_ptr efaces0, fvtxlist;

    MEdge_ptr edge0 = MVs_CommonEdge(fverts[0],fverts[1]);
    if (!edge0)
      return NULL;

    efaces0 = ME_Faces(edge0);
    if (!efaces0)
      return NULL;

    nef0 = List_Num_Entries(efaces0);
    
    for (i = 0; i < nef0; i++) {
      eface = List_Entry(efaces0,i);
      int nfv = MF_Num_Vertices(eface);

      if (nfv != nv) continue;

      fvtxlist = MF_Vertices(eface,1,0);

      contains_all = 1;
      for (j = 0; j < nv; j++) {
        MVertex_ptr vtx = List_Entry(fvtxlist,j);

        if (vtx == fverts[0] || vtx == fverts[1]) continue;

        int k, found = 0;
        for (k = 2; k < nv; k++)
          if (vtx == fverts[k]) {
            found = 1;
            break;
          }

        if (!found) {
          contains_all = 0;
          break;
        }
      }  
    
      List_Delete(fvtxlist);

      if (contains_all) {
	common_face = eface;
	break;
      }
    }
    List_Delete(efaces0);

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

  /* Extra functionality for hash-tables */
  
  MFace_ptr MF_NextInHash(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    return MF_NextInHash_jmp[RTYPE](f);
  }
  
  void MF_Set_NextInHash(MFace_ptr f, MFace_ptr next) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    MF_Set_NextInHash_jmp[RTYPE](f, next);
  }

  void MF_HashKey(MFace_ptr f, unsigned int *pn, void* **pp) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    MF_HashKey_jmp[RTYPE](f, pn, pp);
  }

  void MF_Lock(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    MF_Lock_jmp[RTYPE](f);
  }
  
  void MF_UnLock(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    MF_UnLock_jmp[RTYPE](f);
  }
  
  int MF_IsLocked(MFace_ptr f) {
    RepType RTYPE = MEnt_RepType((MEntity_ptr) f);
    return MF_IsLocked_jmp[RTYPE](f);
  }


  int   MF_GlobalID(MFace_ptr f) {
#ifdef MSTK_HAVE_MPI
    return   MEnt_GlobalID((MEntity_ptr) f);
#else
    return MEnt_ID((MEntity_ptr) f);
#endif
  }

#ifdef MSTK_HAVE_MPI

  PType MF_PType(MFace_ptr f) {
    return MEnt_PType((MEntity_ptr) f);
  }

  void  MF_Set_PType(MFace_ptr f, PType ptype) {
    MEnt_Set_PType((MEntity_ptr) f, ptype);
  }

  int MF_OnParBoundary(MFace_ptr f) {
    return MEnt_OnParBoundary((MEntity_ptr) f); 
  }

  void MF_Flag_OnParBoundary(MFace_ptr f) {
    return MEnt_Flag_OnParBoundary((MEntity_ptr) f);
  }

  void MF_Unflag_OnParBoundary(MFace_ptr f) {
    return MEnt_Unflag_OnParBoundary((MEntity_ptr) f);
  }

  int   MF_MasterParID(MFace_ptr f) {
    return MEnt_MasterParID((MEntity_ptr) f ); 
  }

  void  MF_Set_MasterParID(MFace_ptr f, int masterparid) {
    MEnt_Set_MasterParID((MEntity_ptr) f, masterparid) ;
  }

  void  MF_Set_GlobalID(MFace_ptr f, int globalid) {
    MEnt_Set_GlobalID((MEntity_ptr) f, globalid);
  }

  MFace_ptr MF_GhostNew(Mesh_ptr mesh) {
    MFace_ptr f;
    RepType RTYPE;
    f = (MFace_ptr) malloc(sizeof(MFace));

    MEnt_Init_CmnData((MEntity_ptr) f);
    MEnt_Set_Mesh((MEntity_ptr) f,mesh);
    MEnt_Set_Dim((MEntity_ptr) f,2);
    MEnt_Set_GEntDim((MEntity_ptr) f,4); /* nonsensical value as we don't know what it is */
    MEnt_Set_GEntID((MEntity_ptr) f,0);

    f->adj = (void *) NULL;

    RTYPE = mesh ? MESH_RepType(mesh) : F1;
    MF_Set_RepType(f,RTYPE);
    
    if (mesh) {
	MESH_Add_GhostFace(mesh,f);
    }

    return f;
  }

#endif /* MSTK_HAVE_MPI */
  
#ifdef __cplusplus
}
#endif

