#ifndef _H_MSTK_PRIVATE
#define _H_MSTK_PRIVATE

#include "MSTK_types.h"
#include "MSTK_externs.h"
#include "MSTK_util.h"
#include "MSTK.h"

#ifdef MSTK_HAVE_MPI
#include "mpi.h"
#endif

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

  void       MESH_Clear_Attrib(Mesh_ptr mesh, MAttrib_ptr attrib);


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


  /* Update the element type of the region */

  void MR_Update_ElementType(MRegion_ptr r);

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
  void MF_Destroy_For_MESH_Delete(MFace_ptr f);
  void ME_Destroy_For_MESH_Delete(MEdge_ptr e);
  void MV_Destroy_For_MESH_Delete(MVertex_ptr v);

  void MAttrib_Destroy_For_MESH_Delete(MAttrib_ptr att);

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

  int MESH_BuildFaceClassfn(Mesh_ptr mesh, int use_geometry);
  int MESH_BuildEdgeClassfn(Mesh_ptr mesh, int use_geometry);
  int MESH_BuildVertexClassfn(Mesh_ptr mesh, int use_geometry);


  /* Export */
  /*
  int MESH_Surf_ExportToFLAGX3D_Par(Mesh_ptr mesh, const char *filename, 
				      const int nparts, const int natt, 
				      const char **attnames, int *opts,
				      int *procids);
  int MESH_Vol_ExportToFLAGX3D_Par(Mesh_ptr mesh, const char *filename, 
				      const int nparts, const int natt, 
				      const char **attnames, int *opts,
				      int *procids);
  */
  int MESH_ExportToDXBin(Mesh_ptr mesh, const char *filename);
  int MESH_ReadExodusII_Serial(Mesh_ptr mesh, const char *filename, const int rank);


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


  /* MSTK's internal quicksort functionality that can sort with respect
   to an auxilliary key */

  void mstk_quicksort(int array[], int key[], int i_beg, int i_end);


  /* MSTK's implementation of graph renumbering algorithm by Gibbs-Poole-Stockmeyer */

  int Graph_Renumber_GPS(int nnodes, int nstart, int nadj[], int adj[], 
                         int newmap[], int *depth, int *maxwidth);



  int compareINT(const void *a, const void *b);
  int compareGlobalID(const void *a, const void *b);
  int compareID(const void *a, const void *b);
  int compareVertexCoor(const void *a, const void *b);
  int compareCoorDouble(const void * a, const void * b);
  int compareEdgeINT(const void *a, const void *b);
  int compareEdgeID(const void *a, const void *b);
  int compareFaceINT(const void *a, const void *b);
  int compareFaceID(const void *a, const void *b);
  int compareValence(const void *a, const void *b);




#ifdef MSTK_HAVE_MPI

  /* If you call the routines to set master partition ID or Global ID
     without knowing what you are doing, you can shoot yourself in the
     foot. So if you are casual MSTK user, you are strongly advised
     against calling the Set_MasterParID and Set_GlobalID routines */  

  void  MV_Set_PType(MVertex_ptr v, PType ptype);
  void  MV_Flag_OnParBoundary(MVertex_ptr v);
  void  MV_Unflag_OnParBoundary(MVertex_ptr v);
  void  MV_Set_MasterParID(MVertex_ptr v, int masterpartid);
  void  MV_Set_GlobalID(MVertex_ptr v, int globalid);

  void  ME_Set_PType(MEdge_ptr e, PType ptype);
  void  ME_Flag_OnParBoundary(MEdge_ptr v);
  void  ME_Unflag_OnParBoundary(MEdge_ptr v);
  void  ME_Set_MasterParID(MEdge_ptr e, int masterparid);
  void  ME_Set_GlobalID(MEdge_ptr e, int globalid);

  void  MF_Set_PType(MFace_ptr f, PType ptype);
  void  MF_Flag_OnParBoundary(MFace_ptr v);
  void  MF_Unflag_OnParBoundary(MFace_ptr v);
  void  MF_Set_MasterParID(MFace_ptr f, int masterpartid);
  void  MF_Set_GlobalID(MFace_ptr f, int globalid);

  void  MR_Set_PType(MRegion_ptr r, PType ptype);
  void  MR_Set_MasterParID(MRegion_ptr r, int masterpartid);
  void  MR_Set_GlobalID(MRegion_ptr r, int globalid);

  void  MEnt_Set_PType(MEntity_ptr ent, PType ptype);
  void  MEnt_Set_MasterParID(MEntity_ptr ent, int masterpartid);
  void  MEnt_Set_GlobalID(MEntity_ptr ent, int globalid);


  void         MESH_Set_Prtn(Mesh_ptr mesh, unsigned int partition, 
                             unsigned int numpartitions);
  unsigned int MESH_Prtn(Mesh_ptr mesh);

  void         MESH_Flag_Has_Ghosts_From_Prtn(Mesh_ptr mesh, unsigned int prtn,
                                              MType mtype);  
  unsigned int MESH_Has_Ghosts_From_Prtn(Mesh_ptr mesh, unsigned int prtn, 
                                         MType mtype);
  void         MESH_Flag_Has_Overlaps_On_Prtn(Mesh_ptr mesh, unsigned int prtn, 
                                              MType mtype);
  unsigned int MESH_Has_Overlaps_On_Prtn(Mesh_ptr mesh, unsigned int prtn, 
                                         MType mtype);

  /* Given information about the partitions from with each partition
   * will receive info from, derive information about partitions to
   * which each partition has to send info to */

  int          MESH_Get_OverlapAdj_From_GhostAdj(Mesh_ptr mesh, MSTK_Comm comm);

  unsigned int MESH_Num_GhostPrtns(Mesh_ptr mesh);
  void         MESH_GhostPrtns(Mesh_ptr mesh, unsigned int *pnums);

  void         MESH_Init_Par_Recv_Info(Mesh_ptr mesh);
  void         MESH_Set_Num_Recv_From_Prtn(Mesh_ptr mesh, unsigned int prtn, 
                                           MType mtype, unsigned int numrecv);
  unsigned int MESH_Num_Recv_From_Prtn(Mesh_ptr mesh, unsigned int prtn,
                                       MType mtype);

  int MESH_Sort_GhostLists(Mesh_ptr mesh, 
                           int (*compfunc)(const void*, const void*));

  int MESH_BuildConnection(Mesh_ptr submesh, int topodim, MSTK_Comm comm);

  int MESH_Parallel_AddGhost(Mesh_ptr submesh, int topodim, MSTK_Comm comm);
  int MESH_AssignGlobalIDs_p2p(Mesh_ptr submesh, int topodim, MSTK_Comm comm);


  /* Mesh Partitioning Routines*/

  int        MESH_PartitionWithMetis(Mesh_ptr mesh, int nparts, int **part);
  int        MESH_PartitionWithZoltan(Mesh_ptr mesh, int nparts, int **part, 
                                      int noptions, char **options,
                                      MSTK_Comm comm);
  int        FixColumnPartitions(Mesh_ptr mesh, int *part, MSTK_Comm comm);
  int        FixColumnPartitions_IsSide(Mesh_ptr mesh, MRegion_ptr mr, MFace_ptr rf);
  int        FixColumnPartitions_UpDown(Mesh_ptr mesh, MRegion_ptr mr, MFace_ptr* up, MFace_ptr* dn);

  int        MESH_Partition(Mesh_ptr parentmesh, int num, int *part, Mesh_ptr *submeshes);
  int        MESH_Partition_Some(Mesh_ptr parentmesh, int num, int *part, int begin, int end, Mesh_ptr *submeshes);
  int        MESH_Partition_and_Send(Mesh_ptr parentmesh, int num, int *parts,
                                     int *toranks, int ring, int with_attr, 
                                     int del_inmesh, MSTK_Comm comm, 
                                     Mesh_ptr *mymesh);
  int        MESH_CopyAttr(Mesh_ptr mesh, int num, Mesh_ptr *submesh, const char *attr_name);

  /* build processor boundary */
  int        MESH_BuildPBoundary(Mesh_ptr mesh, Mesh_ptr submesh);
  int        MESH_AssignGlobalIDs(Mesh_ptr submesh, int topodim, int have_GIDs, MSTK_Comm comm);
  int        MESH_LabelPType(Mesh_ptr submesh, int topodim, MSTK_Comm comm);
  int        MESH_BuildSubMesh(Mesh_ptr mesh, int topodim, Mesh_ptr submesh);
  int        MESH_ConcatSubMesh(Mesh_ptr mesh, int topodim, int num, 
                                Mesh_ptr *submeshes);
  /* add ghost elements */
  int        MESH_AddGhost(Mesh_ptr mesh, Mesh_ptr submesh, int part_no, 
                           int ring);
  /* send and receive mesh */
  int        MESH_SendMesh(Mesh_ptr mesh, int torank, int with_attr,
                           MSTK_Comm comm,
                           int *numreq, int *maxreq, MPI_Request **requests,
                           int *numptrs2free, int *maxptrs2free,
                           void ***ptrs2free);
  int        MESH_Send_MetaData(Mesh_ptr mesh, int torank, MSTK_Comm comm,
                                int *numreq, int *maxreq, 
                                MPI_Request **requests,
                                int *numptrs2free, int *maxptrs2free,
                                void ***ptrs2free);
  int        MESH_Send_Vertices(Mesh_ptr mesh, int torank, MSTK_Comm comm, 
                                int *numreq, int *maxreq, 
                                MPI_Request **requests,
                                int *numptrs2free, int *maxptrs2free,
                                void ***ptrs2free);
  int        MESH_Send_VertexCoords(Mesh_ptr mesh, int torank, MSTK_Comm comm, 
                                    int *numreq, int *maxreq, 
                                    MPI_Request **requests,
                                    int *numptrs2free, int *maxptrs2free,
                                    void ***ptrs2free);
  int        MESH_Send_NonVertexEntities(Mesh_ptr mesh, int torank, 
                                         MSTK_Comm comm, 
                                         int *numreq, int *maxreq, 
                                         MPI_Request **requests,
                                         int *numptrs2free, int *maxptrs2free,
                                         void ***ptrs2free);
  int        MESH_SendAttributes(Mesh_ptr mesh, int torank, MSTK_Comm comm,
                                 int *numreq, int *maxreq, 
                                 MPI_Request **requests, int *numptrs2free, 
                                 int *maxptrs2free, void ***ptrs2free);
  int        MESH_Send_AttributeMetaData(Mesh_ptr mesh, int torank, 
                                         MSTK_Comm comm,
                                         int *numreq, int *maxreq, 
                                         MPI_Request **requests, 
                                         int *numptrs2free, 
                                         int *maxptrs2free, void ***ptrs2free);
  int        MESH_Send_Attribute(Mesh_ptr mesh, MAttrib_ptr attrib, int torank, 
                                 MSTK_Comm comm, int *numreq, int *maxreq, 
                                 MPI_Request **requests, int *numptrs2free, 
                                 int *maxptrs2free, void ***ptrs2free);
  int        MESH_Send_MSetMetaData(Mesh_ptr mesh, int torank, 
                                    MSTK_Comm comm,
                                    int *numreq, int *maxreq, 
                                    MPI_Request **requests, 
                                    int *numptrs2free, 
                                    int *maxptrs2free, void ***ptrs2free);
  int        MESH_SendMSets(Mesh_ptr mesh, int torank, MSTK_Comm comm,
                            int *numreq, int *maxreq, MPI_Request **requests, 
                            int *numptrs2free, int *maxptrs2free, 
                            void ***ptrs2free);
  int        MESH_Send_MSet(Mesh_ptr mesh, MSet_ptr attrib, int torank, 
                            MSTK_Comm comm, int *numreq, int *maxreq, 
                            MPI_Request **requests, int *numptrs2free, 
                            int *maxptrs2free, void ***ptrs2free);
  int        MESH_RecvMesh(Mesh_ptr mesh, int fromrank, int with_attr,
                           MSTK_Comm comm);
  int        MESH_Recv_MetaData(Mesh_ptr mesh, int fromrank, RepType *rtype,
                                int *nv, int *ne, int *nf, int *nr,
                                MSTK_Comm comm);
  int        MESH_Recv_Vertices(Mesh_ptr mesh, int fromrank, int nv,
                                MSTK_Comm comm);
  int        MESH_Recv_VertexCoords(Mesh_ptr mesh, int fromrank, int nv,
                                MSTK_Comm comm);
  int        MESH_Recv_NonVertexEntities(Mesh_ptr mesh, int fromrank, int ne,
                                         int nf, int nr, MSTK_Comm comm);
  int        MESH_RecvAttributes(Mesh_ptr mesh, int fromrank, MSTK_Comm comm);
  int        MESH_Recv_AttributeMetaData(Mesh_ptr mesh, int fromrank, 
                                         MSTK_Comm comm);
  int        MESH_Recv_Attribute(Mesh_ptr mesh, MAttrib_ptr attrib, 
                                 int fromrank, MSTK_Comm comm);
  int        MESH_RecvMSets(Mesh_ptr mesh, int fromrank, MSTK_Comm comm);
  int        MESH_Recv_MSetMetaData(Mesh_ptr mesh, int fromrank, 
                                    MSTK_Comm comm);
  int        MESH_Recv_MSet(Mesh_ptr mesh, MSet_ptr mset, int fromrank, 
                            MSTK_Comm comm);



  /* build ghost list */
  int        MESH_Build_GhostLists(Mesh_ptr mesh, int topodim);

  /* create entity */

  MVertex_ptr MV_GhostNew(Mesh_ptr mesh);
  MEdge_ptr   ME_GhostNew(Mesh_ptr mesh);
  MFace_ptr   MF_GhostNew(Mesh_ptr mesh);
  MRegion_ptr MR_GhostNew(Mesh_ptr mesh);


  /* functions to update ghost info and attributes */
  int        MESH_ParallelAdj_Current(Mesh_ptr mesh);
  void       MESH_Mark_ParallelAdj_Current(Mesh_ptr mesh);
  void       MESH_Mark_ParallelAdj_Stale(Mesh_ptr mesh);
  int        MESH_Update_ParallelAdj(Mesh_ptr mesh, MSTK_Comm comm);


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


  /* Functions for entity sets */

  int        MESH_CopySet(Mesh_ptr mesh, int num, Mesh_ptr *submeshes, 
                          MSet_ptr mset);

  /* Routine to exchange attributes between ghost and master entities 
     when a 1-to-1 mapping exists (edges in 2D and faces in 3D) */

  int        MESH_XchngEdgeFaceAttrib(Mesh_ptr mesh, MAttrib_ptr attrib,
                                      MSTK_Comm comm);

  /* Functions for improving searching for entities by global ID */

  void       MESH_Enable_GlobalIDSearch(Mesh_ptr mesh);
  void       MESH_Disable_GlobalIDSearch(Mesh_ptr mesh);
#endif /* MSTK_HAVE_MPI */  



#ifdef __cplusplus
	   }
#endif

#endif
