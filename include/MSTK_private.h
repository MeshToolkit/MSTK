#ifndef _H_MSTK_PRIVATE
#define _H_MSTK_PRIVATE

#include "MSTK_types.h"
#include "MSTK_externs.h"
#include "MSTK_util.h"
#include "MSTK_malloc.h"
#include "MSTK.h"

#ifdef   __cplusplus
extern "C" {
#endif

/* If MSTK_KEEP_DELETED is 1, then entities will only be marked as deleted */

extern int MSTK_KEEP_DELETED;

/* Don't want users to see this */

typedef enum MDelType {MDELREGION=-40, MDELFACE=-30, MDELEDGE=-20, MDELVERTEX=-10} MDelType;


/* THIS FILE HAS ADDITIONAL FUNCTIONS THAT THE NORMAL USER NEED NOT SEE */


  void       MESH_Add_Vertex(Mesh_ptr mesh, MVertex_ptr v);
  void       MESH_Add_Edge(Mesh_ptr mesh, MEdge_ptr e);
  void       MESH_Add_Face(Mesh_ptr mesh, MFace_ptr f);
  void       MESH_Add_Region(Mesh_ptr mesh, MRegion_ptr r);

  void       MESH_Rem_Vertex(Mesh_ptr mesh, MVertex_ptr v);
  void       MESH_Rem_Edge(Mesh_ptr mesh, MEdge_ptr e);
  void       MESH_Rem_Face(Mesh_ptr mesh, MFace_ptr f);
  void       MESH_Rem_Region(Mesh_ptr mesh, MRegion_ptr r);


  void       MESH_Add_Attrib(Mesh_ptr mesh, MAttrib_ptr attrib);
  void       MESH_Rem_Attrib(Mesh_ptr mesh, MAttrib_ptr attrib);


/*
  void MV_Set_RepType(MVertex_ptr v, RepType rtype);
*/
  void MV_Add_Edge(MVertex_ptr mvertex, MEdge_ptr medge);
  void MV_Add_Face(MVertex_ptr mvertex, MFace_ptr mface);
  void MV_Add_Region(MVertex_ptr mvertex, MRegion_ptr mregion);
  void MV_Rem_Edge(MVertex_ptr mvertex, MEdge_ptr medge);
  void MV_Rem_Face(MVertex_ptr mvertex, MFace_ptr mface);
  void MV_Rem_Region(MVertex_ptr mvertex, MRegion_ptr mregion);

/*
  void ME_Set_RepType(MEdge_ptr medge, RepType rtype*/

  void ME_Add_Face(MEdge_ptr medge, MFace_ptr mface);
  void ME_Add_Region(MEdge_ptr medge, MRegion_ptr mregion);
  void ME_Rem_Face(MEdge_ptr medge, MFace_ptr mface);
  void ME_Rem_Region(MEdge_ptr medge, MRegion_ptr mregion);



/*
  void MF_Set_RepType(MFace_ptr f, RepType rtype);
*/
  /* Add, Remove AdjFace can be automatically called when faces are
     being created or deleted. They do not need user invocation */

  void MF_Add_AdjFace(MFace_ptr f, int enbr, MFace_ptr af);
  void MF_Rem_AdjFace(MFace_ptr f, int enbr, MFace_ptr af);


  /* Add, Remove Region can be called internally when a region is
     being created from a list of faces */

  void MF_Add_Region(MFace_ptr f, MRegion_ptr r, int side);
  void MF_Rem_Region(MFace_ptr f, MRegion_ptr r);


  /*
  void MR_Set_RepType(MRegion_ptr r, RepType rtype);
  */

  /* Adjacent region info will be updated by private functions when
     regions are created or deleted. There is no need for invocation
     by users/applications */

  void MR_Add_AdjRegion(MRegion_ptr r, int facenum, MRegion_ptr ar);
  void MR_Rem_AdjRegion(MRegion_ptr r, MRegion_ptr ar);


  /* These are functions that MESH_Delete calls. These functions
     destroy entities without worrying about updating entities that
     point to them */

  void MR_Destroy_For_MESH_Delete(MRegion_ptr r);
  void MF_Destroy_For_MESH_Delete(MFace_ptr r);
  void ME_Destroy_For_MESH_Delete(MEdge_ptr r);
  void MV_Destroy_For_MESH_Delete(MVertex_ptr r);


  /* Init/Free data common to all entities */
  void MEnt_Init_CmnData(MEntity_ptr ent);
  void MEnt_Free_CmnData(MEntity_ptr ent);

  /* Set Entity dimension */
  void MEnt_Set_Dim(MEntity_ptr ent, MType dim);

  /* Set the mesh that entity belongs to */
  void MEnt_Set_Mesh(MEntity_ptr ent, Mesh_ptr mesh);

  /* Mark entity as volatile/temporary in reduced representations */
  void MEnt_Set_Volatile(MEntity_ptr ent);

  /* Mark an entity as deleted - 
     no update of topology or association with a mesh */

  void MEnt_Set_DelFlag(MEntity_ptr ent);
  void MEnt_Rem_DelFlag(MEntity_ptr ent);
 
  /* Setting RepType for entities */

  /* This just sets the data */
  void MEnt_Set_RepType_Data(MEntity_ptr mentity, RepType rtype);

  /* This calls M*_Set_RepType based on the dimension of the entity */
  void MEnt_Set_RepType(MEntity_ptr mentity, RepType rtype);
  

  /* Private functions for defining instance of attributes, i.e., objects 
     carry the actual attribute value for some mesh entity */

  MAttIns_ptr MAttIns_New(MAttrib_ptr attrib);
  MAttrib_ptr MAttIns_Attrib(MAttIns_ptr att);
  void        MAttIns_Set_Value(MAttIns_ptr att, int ival, double rval, void *pval);
  int         MAttIns_Get_Value(MAttIns_ptr att, int *ival, double *rval, void **pval);
  void        MAttIns_Delete(MAttIns_ptr att);

  /* Building mesh classification (only top level function in MSTK.h) */

  int MESH_BuildFaceClassfn(Mesh_ptr mesh);
  int MESH_BuildEdgeClassfn(Mesh_ptr mesh);
  int MESH_BuildVertexClassfn(Mesh_ptr mesh);


  /* Export */

  int MESH_Surf_ExportToFLAGX3D_Par(Mesh_ptr mesh, const char *filename, 
				      const int nparts, const int natt, 
				      const char **attnames, int *opts,
				      int *procids);
  int MESH_Vol_ExportToFLAGX3D_Par(Mesh_ptr mesh, const char *filename, 
				      const int nparts, const int natt, 
				      const char **attnames, int *opts,
				      int *procids);

  int MESH_ExportToDXBin(Mesh_ptr mesh, const char *filename);

  /* Adhoc routines for parallel output of MSTK files */

  int MESH_Init_ParAtts(Mesh_ptr mesh);

  int MEnt_NumProcs(MEntity_ptr ent);
  void MEnt_Set_ProcIDs(MEntity_ptr ent, int np, int *procids);
  int MEnt_ProcIDs(MEntity_ptr ent, int *np, int *procids);
  void MEnt_Set_LocalID(MEntity_ptr ent, int procid, int lnum);
  int MEnt_LocalID(MEntity_ptr ent, int procid);


  /* Extra functionality for List manipulation - risky for uninformed users */

  int List_Size_Raw(List_ptr l);
  void *List_Entry_Raw(List_ptr l, int i);
  int List_Remi_Raw(List_ptr l, int i);

  /* Extra functionality for hash-tables */


  
  Hash_ptr Hash_New(unsigned int inisize, int type);
  void Hash_Delete(Hash_ptr h);
  
  Hash_ptr Hash_Add(Hash_ptr h, void *entry, unsigned int np, void* *p);
  Hash_ptr Hash_ChknAdd(Hash_ptr h, void *entry, unsigned int np, void* *p);
  int      Hash_Rem(Hash_ptr h, void *entry, unsigned int np, void* *p);
  void    *Hash_Entry(Hash_ptr h, unsigned int np, void* *p);
  int      Hash_Num_Entries(Hash_ptr h);
  List_ptr Hash_Entries(Hash_ptr h);

  void     Hash_Print(Hash_ptr h);

  void Hash_Lock(int *plock);
  void Hash_UnLock(int *plock);
  int Hash_IsLocked(int lock);

  int  Hash_AutoRemove(Hash_ptr h);
  void Hash_Set_AutoRemove(Hash_ptr h, int t);
  unsigned int Hash_Remove_Unused(Hash_ptr h);

  MEdge_ptr ME_NextInHash(MEdge_ptr medge);
  void ME_Set_NextInHash(MEdge_ptr medge, MEdge_ptr next);
  void ME_HashKey(MEdge_ptr medge, unsigned int *pn, void* **pp);

  MFace_ptr MF_NextInHash(MFace_ptr mface);
  void MF_Set_NextInHash(MFace_ptr mface, MFace_ptr next);
  void MF_HashKey(MFace_ptr medge, unsigned int *pn, void* **pp);
  
  MEntity_ptr MEnt_NextInHash(MEntity_ptr ent);
  void MEnt_Set_NextInHash(MEntity_ptr ent, MEntity_ptr next);
  void MEnt_HashKey(MEntity_ptr ent, unsigned int *pn, void* **pp);
  
  int ME_IsLocked(MEdge_ptr e);
  int MF_IsLocked(MFace_ptr f);
  int MEnt_IsLocked(MEntity_ptr ent);



#ifdef MSTK_HAVE_MPI
  /* build ghost list */
  int        MESH_Build_GhostLists(Mesh_ptr mesh);

  /* create entity */
  /* no OvEntity create function since they are just pointers to regular Entity */
  MVertex_ptr MV_GhostNew(Mesh_ptr mesh);
  MEdge_ptr   ME_GhostNew(Mesh_ptr mesh);
  MFace_ptr   MF_GhostNew(Mesh_ptr mesh);
  MRegion_ptr MR_GhostNew(Mesh_ptr mesh);


  /* functions to update ghost info and attributes */
  /* Must be preceded by MESH_Update_ProcessorRel - Not recommended for users */
  int        MESH_UpdateAttr(Mesh_ptr mesh, const char *attr_name, int rank, 
			     int num,  MPI_Comm comm);
  int        MESH_Update_ProcessorRel(Mesh_ptr mesh, int rank, int num,  
				   MPI_Comm comm);
  int*       MESH_ProcessorRel(Mesh_ptr mesh);
  int*       MESH_LocalInfo(Mesh_ptr mesh);
  void       MESH_Set_ProcessorRel(Mesh_ptr mesh, int *proc_mesh_rel);
  void       MESH_Set_LocalInfo(Mesh_ptr mesh, int *local_info);


  /* Functions for entity sets */

  int        MESH_CopySet(Mesh_ptr mesh, Mesh_ptr submesh, MSet_ptr mset);
  int        MESH_SendMSet(Mesh_ptr mesh, const char *set_name, int rank, MPI_Comm comm);
  int        MESH_RecvMSet(Mesh_ptr mesh, const char *set_name, int rank, int recv_rank, MPI_Comm comm);


  /**** ????????? ******/
  /* the following 16 functions is not necessary, 
     we can operate on the whole list, and then call
     MESH_BuildGhostList()
  */
  void       MESH_Add_GhostVertex(Mesh_ptr mesh, MVertex_ptr v);
  void       MESH_Add_GhostEdge(Mesh_ptr mesh, MEdge_ptr e);
  void       MESH_Add_GhostFace(Mesh_ptr mesh, MFace_ptr f);
  void       MESH_Add_GhostRegion(Mesh_ptr mesh, MRegion_ptr r);

  void       MESH_Add_OverlapVertex(Mesh_ptr mesh, MVertex_ptr v);
  void       MESH_Add_OverlapEdge(Mesh_ptr mesh, MEdge_ptr e);
  void       MESH_Add_OverlapFace(Mesh_ptr mesh, MFace_ptr f);
  void       MESH_Add_OverlapRegion(Mesh_ptr mesh, MRegion_ptr r);

  void       MESH_Rem_GhostVertex(Mesh_ptr mesh, MVertex_ptr v);
  void       MESH_Rem_GhostEdge(Mesh_ptr mesh, MEdge_ptr e);
  void       MESH_Rem_GhostFace(Mesh_ptr mesh, MFace_ptr f);
  void       MESH_Rem_GhostRegion(Mesh_ptr mesh, MRegion_ptr r);

  
  void       MESH_Rem_OverlapVertex(Mesh_ptr mesh, MVertex_ptr v);
  void       MESH_Rem_OverlapEdge(Mesh_ptr mesh, MEdge_ptr e);
  void       MESH_Rem_OverlapFace(Mesh_ptr mesh, MFace_ptr f);
  void       MESH_Rem_OverlapRegion(Mesh_ptr mesh, MRegion_ptr r);

  int compareGlobalID(MEntity_ptr *a, MEntity_ptr *b);
#endif


#ifdef __cplusplus
	   }
#endif

#endif
