#ifndef _H_MSTK
#define _H_MSTK

/* Main Header file defining the MSTK API. 

   Author: Rao V Garimella
   Email: rao@lanl.gov
   Website: https://software.lanl.gov/MeshTools/hg/mstk

   The MSTK source is available for use under an LGPL license 
*/


#define MSTK_VERSION "2.00rc4"

#include <stdarg.h>

#include "MSTK_defines.h"
#include "MSTK_types.h"
#include "MSTK_externs.h"
#include "MSTK_util.h"
#include "MSTK_malloc.h"

#ifdef MSTK_HAVE_MPI
#include <mpi.h>
#endif


/* 
   IMPORTANT NOTE: 
   For serial builds, MSTK_types.h typedefs MSTK_Comm to be "void *" so
   that serial codes can pass in NULL for the comm argument. Also, 
   serial codes DO NOT have to include mpi.h or link with mpi.
*/


#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************/
void MSTK_Init(void);

/********************************************************************/
/*        MESH OBJECT OPERATORS                                     */
/********************************************************************/  

  Mesh_ptr    MESH_New(RepType type);
  void        MESH_Delete(Mesh_ptr mesh);

  /* Read from the native MSTK format. 'comm' can be NULL for serial codes */
  int         MESH_InitFromFile(Mesh_ptr mesh, const char *filename, MSTK_Comm comm);

  /* Initialize a mesh object on the fly from a general description of the mesh */
  int         MESH_InitFromGenDesc(Mesh_ptr mesh, int nv, double (*xyz)[3],
				  int nf, int *nfv, int **fvids, int nr, 
				  int *nrv, int **rvids, int *nrf, 
				  int ***rfvtemplate);

  /* Import mesh data into the mesh object from various
     formats. 'comm' can be NULL for serial codes */
  /* opts is an integer array with different meanings for different
     formats */

  int         MESH_ImportFromFile(Mesh_ptr mesh, const char *filename, 
                                  const char *format, int *opts, MSTK_Comm comm);
  int         MESH_ImportFromGMV(Mesh_ptr mesh, const char *filename, MSTK_Comm comm); 

  /*--------------------------------------------------------------------------*/
  /* Read an Exodus II file into MSTK */
  /* The input arguments parallel_opts indicate if we should build
     parallel adjacencies and some options for controlling that 

     if parallel_opts == NULL, then only processor 0 will read the mesh and
     do nothing else
     
     parallel_opts[0] = 1/0  ---- distribute/Don't distribute the mesh
     parallel_opts[1] = 0-3  ---- 0: read/partition on P0 and distribute
                                  1: read/partition on P0, reread local portion
                                     on each process and weave/connect
                                  2: read on multiple processors Pn where n is 
                                  an arithmetic progression (5,10,15 etc) and
                                  distribute to a subset of processors (like P5
     Opts 1,2,3 NOT ACTIVE        distributes to P0-P4, P10 to P6-P9 etc)
                                  3: read portions of the mesh on each processor
                                  and repartition
     parallel_opts[1] = N    ---- Number of ghost layers around mesh on proc
     parallel_opts[2] = 0/1  ---- partitioning method
                                  0: Metis
                                  1: Zoltan */                            
  /*--------------------------------------------------------------------------*/

  int         MESH_ImportFromExodusII(Mesh_ptr mesh, const char *filename, int *parallel_opts, MSTK_Comm comm);


  /*--------------------------------------------------------------------------*/
  /* Read an Nemesis I file into MSTK */
  /* The input arguments parallel_opts indicate if we should build
     parallel adjacencies and some options for controlling that 
     
     parallel_opts[0] = 1/0  ---- Weave/Don't weave the distributed meshes
                                  together to form parallel connections
     parallel_opts[1] = N    ---- Number of ghost layers around mesh  */
  /*--------------------------------------------------------------------------*/

  int         MESH_ImportFromNemesisI(Mesh_ptr mesh, const char *filename, int *parallel_opts, MSTK_Comm comm);



  int         MESH_ImportFromFLAGX3D(Mesh_ptr mesh, const char *filename, MSTK_Comm comm);


  /* Write mesh data into a file in the native MSTK format. 'comm' can
     be NULL for serial codes */
  int         MESH_WriteToFile(Mesh_ptr mesh, const char *filename, RepType rtype, MSTK_Comm comm);

  /* Export mesh data to various file formats. 'comm' can be NULL for
     serial codes */
  int         MESH_ExportToFile(Mesh_ptr mesh, const char *filename,
                                const char *format, const int natt, 
                                const char **attnames, const int *opts, MSTK_Comm comm);
  int         MESH_ExportToGMV(Mesh_ptr mesh, const char *filename, 
                               const int natt, const char **attnames,
                               const int *opts, MSTK_Comm comm);
  int         MESH_ExportToFLAGX3D(Mesh_ptr mesh, const char *filename, 
                                   const int natt, const char **attnames, 
                                   const int *opts, MSTK_Comm comm);
  int         MESH_ExportToExodusII(Mesh_ptr mesh, const char *filename, 
                                    const int natt, const char **attnames, 
                                    const int *opts, MSTK_Comm comm);

  int         MESH_ExportToSTL(Mesh_ptr mesh, const char *filename);
  int         MESH_ExportToDX(Mesh_ptr mesh, const char *filename, int binary);

  int         MESH_BuildClassfn(Mesh_ptr mesh);
  int         MESH_DelInterior(Mesh_ptr mesh);
  int         MESH_Tet2Hex(Mesh_ptr tetmesh, Mesh_ptr *hexmesh);

  /* Check if mesh topology is valid */
  int         MESH_CheckTopo(Mesh_ptr mesh);

  /* Check if mesh geometry is valid */
  /* int        MESH_CheckGeom(Mesh_ptr mesh); */

  /* Duplicate a mesh optionally with attributes and mesh sets */
  int         MESH_Copy(Mesh_ptr srcmesh, Mesh_ptr dstmesh, int with_attr, 
                        int with_sets);


  GModel_ptr  MESH_GModel(Mesh_ptr mesh);
  RepType     MESH_RepType(Mesh_ptr mesh);
  char       *MESH_RepType_Str(Mesh_ptr mesh);

  int         MESH_Num_Attribs(Mesh_ptr mesh);
  MAttrib_ptr MESH_Attrib(Mesh_ptr mesh, int i);
  MAttrib_ptr MESH_Next_Attrib(Mesh_ptr mesh, int *index);
  MAttrib_ptr MESH_AttribByName(Mesh_ptr mesh, const char *name);

  int         MESH_Num_MSets(Mesh_ptr mesh);
  MSet_ptr    MESH_MSet(Mesh_ptr mesh, int i);
  MSet_ptr    MESH_Next_MSet(Mesh_ptr mesh, int *index);
  MSet_ptr    MESH_MSetByName(Mesh_ptr mesh, const char *name);

  int         MESH_Num_Vertices(Mesh_ptr mesh);
  int         MESH_Num_Edges(Mesh_ptr mesh);
  int         MESH_Num_Faces(Mesh_ptr mesh);
  int         MESH_Num_Regions(Mesh_ptr mesh);

  MVertex_ptr MESH_Vertex(Mesh_ptr mesh, int i);
  MEdge_ptr   MESH_Edge(Mesh_ptr mesh, int i);
  MFace_ptr   MESH_Face(Mesh_ptr mesh, int i);
  MRegion_ptr MESH_Region(Mesh_ptr mesh, int i);

  MVertex_ptr MESH_Next_Vertex(Mesh_ptr mesh, int *index);
  MEdge_ptr   MESH_Next_Edge(Mesh_ptr mesh, int *index);
  MFace_ptr   MESH_Next_Face(Mesh_ptr mesh, int *index);
  MRegion_ptr MESH_Next_Region(Mesh_ptr mesh, int *index);


  MVertex_ptr MESH_VertexFromID(Mesh_ptr mesh, int i);
  MEdge_ptr   MESH_EdgeFromID(Mesh_ptr mesh, int i);
  MFace_ptr   MESH_FaceFromID(Mesh_ptr mesh, int i);
  MRegion_ptr MESH_RegionFromID(Mesh_ptr mesh, int i);
  MEntity_ptr MESH_EntityFromID(Mesh_ptr mesh, MType mtype, int i);

  
  int         MESH_SetRepType(Mesh_ptr mesh, RepType type);
  void        MESH_SetGModel(Mesh_ptr mesh, GModel_ptr geom);
  int         MESH_Change_RepType(Mesh_ptr mesh, int nurep);

  void        MESH_Renumber(Mesh_ptr mesh);


#ifdef MSTK_HAVE_MPI

  /* PRIMARY ROUTINES FOR PARALLEL APPLICATIONS */

  /* Read a mesh in, partition it and distribute it to 'num' processors */
  /* 'ring' indicates the number of ghost layers (can only do 1 for now)*/
  /* 'with_attr' distribute attributes to submeshes                     */
  /* 'method' indicates the partitioning method (0 - metis, 1 - zoltan) */  

  int         MSTK_Mesh_Read_Distribute(Mesh_ptr *recv_mesh, 
                                        const char* global_mesh_name, 
                                        int *topodim, int ring, int with_attr, 
                                        int method, MSTK_Comm comm);

  /* Partition an existing mesh (globalmesh) on processor 0 and
     distribute it to 'num' processors. The resulting mesh on my
     partition is in mymesh. If mymesh is already initialized, that
     mesh pointer is used as is and the information filled in. It is
     safest to set mymesh to NULL*/

  int         MSTK_Mesh_Distribute(Mesh_ptr globalmesh, Mesh_ptr *mymesh, 
                                   int *topodim, int ring, int with_attr, 
                                   int method, MSTK_Comm comm);



  /* 'Weave' a set of distributed meshes by building ghost layers */

  int         MSTK_Weave_DistributedMeshes(Mesh_ptr mesh, int topodim,
                                           int num_ghost_layers,
                                           int input_type, MSTK_Comm comm);

  /* Parallel update attribute values for ghost entities */

  int         MSTK_UpdateAttr(Mesh_ptr mesh, MSTK_Comm comm);


  /* Update vertex coordinates for ghost vertices */
  int         MESH_UpdateVertexCoords(Mesh_ptr mesh, MSTK_Comm comm);

  /* Check parallel consistency */

  int         MESH_Parallel_Check(Mesh_ptr mesh, MSTK_Comm comm);


  MVertex_ptr MESH_VertexFromGlobalID(Mesh_ptr mesh, int global_id);
  MEdge_ptr   MESH_EdgeFromGlobalID(Mesh_ptr mesh, int global_id);
  MFace_ptr   MESH_FaceFromGlobalID(Mesh_ptr mesh, int global_id);
  MRegion_ptr MESH_RegionFromGlobalID(Mesh_ptr mesh, int global_id);
  MEntity_ptr MESH_EntityFromGlobalID(Mesh_ptr mesh, MType mtype, int i);


  /* Get a partitioning for mesh using METIS (method=1) or ZOLTAN (method=2) */
  /* Doesn't actually partition the mesh or distribute it                    */

  int        MESH_Get_Partitioning(Mesh_ptr mesh, int method, int **part, 
                                   MPI_Comm);


#endif /* MSTK_HAVE_MPI */



  Hash_ptr    MESH_Hash_Edges(Mesh_ptr mesh);
  Hash_ptr    MESH_Hash_Faces(Mesh_ptr mesh);
  int         MESH_AutoLock(Mesh_ptr mesh);
  void        MESH_Set_AutoLock(Mesh_ptr mesh, int autolock);
  

/********************************************************************/
/*        MESH VERTEX OPERATORS                                     */
/********************************************************************/

  MVertex_ptr MV_New(Mesh_ptr mesh);
  void        MV_Delete(MVertex_ptr mvertex, int keep);
  void        MV_Restore(MVertex_ptr mvertex);
  void        MV_Set_RepType(MVertex_ptr mvertex, RepType reptype);
  void        MV_Set_Coords(MVertex_ptr mvertex, double *xyz);
  void        MV_Set_GEntity(MVertex_ptr mvertex, GEntity_ptr gent);
  void        MV_Set_GEntDim(MVertex_ptr mvertex, int gdim);
  void        MV_Set_GEntID(MVertex_ptr mvertex, int gid);
  void        MV_Add_AdjVertex(MVertex_ptr mvertex, MVertex_ptr adjvertex);
  void        MV_Rem_AdjVertex(MVertex_ptr mvertex, MVertex_ptr adjvertex);
  void        MV_Set_ID(MVertex_ptr mvertex, int id);

  Mesh_ptr    MV_Mesh(MVertex_ptr mv);
  int         MV_ID(MVertex_ptr mvertex);
  int         MV_GEntDim(MVertex_ptr mvertex);
  int         MV_GEntID(MVertex_ptr mvertex);
  GEntity_ptr MV_GEntity(MVertex_ptr mvertex);

  void        MV_Coords(MVertex_ptr mvertex, double *xyz);

  int         MV_Num_AdjVertices(MVertex_ptr mvertex);
  int         MV_Num_Edges(MVertex_ptr mvertex);
  int         MV_Num_Faces(MVertex_ptr mvertex);
  int         MV_Num_Regions(MVertex_ptr mvertex);
  List_ptr    MV_AdjVertices(MVertex_ptr mvertex);
  List_ptr    MV_Edges(MVertex_ptr mvertex);
  List_ptr    MV_Faces(MVertex_ptr mvertex);
  List_ptr    MV_Regions(MVertex_ptr mvertex);

#ifdef MSTK_HAVE_MPI

  PType       MV_PType(MVertex_ptr v);
  int         MV_MasterParID(MVertex_ptr v);
  int         MV_GlobalID(MVertex_ptr v);

#endif

/********************************************************************/
/*        MESH EDGE OPERATORS                                       */
/********************************************************************/

  MEdge_ptr   ME_New(Mesh_ptr mesh);
  void        ME_Delete(MEdge_ptr medge, int keep);
  void        ME_Restore(MEdge_ptr medge);
  void        ME_Set_RepType(MEdge_ptr medge, RepType reptype);
  void        ME_Set_GEntity(MEdge_ptr medge, GEntity_ptr gent);
  void        ME_Set_GEntDim(MEdge_ptr medge, int gdim);
  void        ME_Set_GEntID(MEdge_ptr medge, int gid);
  int         ME_Set_GInfo_Auto(MEdge_ptr medge);
  void        ME_Set_ID(MEdge_ptr medge, int id);

  void        ME_Set_Vertex(MEdge_ptr medge, int i, MVertex_ptr vertex);
  void        ME_Replace_Vertex(MEdge_ptr medge, MVertex_ptr vert, MVertex_ptr nuvert);

  Mesh_ptr    ME_Mesh(MEdge_ptr medge);
  int         ME_ID(MEdge_ptr medge);
  int         ME_GEntDim(MEdge_ptr medge);
  int         ME_GEntID(MEdge_ptr medge);
  GEntity_ptr ME_GEntity(MEdge_ptr medge);
  int         ME_Num_Faces(MEdge_ptr medge);
  int         ME_Num_Regions(MEdge_ptr medge);
  MVertex_ptr ME_Vertex(MEdge_ptr medge, int i);
  MVertex_ptr ME_OppVertex(MEdge_ptr medge, MVertex_ptr ov);
  int         ME_UsesEntity(MEdge_ptr medge, MEntity_ptr mentity, int etype);
  List_ptr    ME_Faces(MEdge_ptr medge);
  List_ptr    ME_Regions(MEdge_ptr medge);


  MEdge_ptr   MVs_CommonEdge(MVertex_ptr v1, MVertex_ptr v2);  

  double      ME_Len(MEdge_ptr e);
  double      ME_LenSqr(MEdge_ptr e);
  void        ME_Vec(MEdge_ptr e, double *evec);

  int         MEs_AreSame(MEdge_ptr e1, MEdge_ptr e2);

  void        ME_Lock(MEdge_ptr e);
  void        ME_UnLock(MEdge_ptr e);

#ifdef MSTK_HAVE_MPI

  PType       ME_PType(MEdge_ptr e);
  int         ME_MasterParID(MEdge_ptr e);
  int         ME_GlobalID(MEdge_ptr e);

#endif

/********************************************************************/
/*        MESH FACE OPERATORS                                       */
/********************************************************************/
  MFace_ptr   MF_New(Mesh_ptr mesh);
  void        MF_Delete(MFace_ptr mface, int keep);
  void        MF_Restore(MFace_ptr mface);
  void        MF_Set_RepType(MFace_ptr mface, RepType reptype);
  void        MF_Set_GEntity(MFace_ptr mface, GEntity_ptr gent);
  void        MF_Set_GEntDim(MFace_ptr mface, int gdim);
  void        MF_Set_GEntID(MFace_ptr mface, int gid);
  int         MF_Set_GInfo_Auto(MFace_ptr mface);
  void        MF_Set_ID(MFace_ptr mface, int id);

  /* These can be called by a user/application after calling MF_New() */
  void        MF_Set_Edges(MFace_ptr mface, int n, MEdge_ptr *edges, int *dirs);
  void        MF_Set_Vertices(MFace_ptr mface, int n, MVertex_ptr *verts);

  /* Can be called by higher level mesh modification operators */
  void        MF_Replace_Edges(MFace_ptr mface, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges);
  void        MF_Replace_Vertex(MFace_ptr mface, MVertex_ptr mvertex, MVertex_ptr nuvertex);
  void        MF_Replace_Edges_i(MFace_ptr mface, int nold, int i, int nnu, MEdge_ptr *nuedge);
  void MF_Rem_Edge(MFace_ptr mface, MEdge_ptr medge);
  void        MF_Replace_Vertex_i(MFace_ptr mface, int i, MVertex_ptr nuvertex);
  void        MF_Insert_Vertex(MFace_ptr mface, MVertex_ptr nuv, MVertex_ptr b4v);
  void        MF_Insert_Vertex_i(MFace_ptr mface, MVertex_ptr nuv, int i);
  int         MF_Rev_EdgeDir(MFace_ptr mface, MEdge_ptr medge);
  int         MF_Rev_EdgeDir_i(MFace_ptr mface, int i);

  Mesh_ptr    MF_Mesh(MFace_ptr mf);
  int         MF_ID(MFace_ptr mface);
  int         MF_GEntDim(MFace_ptr mface);
  int         MF_GEntID(MFace_ptr mface);
  GEntity_ptr MF_GEntity(MFace_ptr mface);

  MFType MF_ElementType(MFace_ptr mf);
  int         MF_Num_Vertices(MFace_ptr mface);
  int         MF_Num_Edges(MFace_ptr mface);
  int         MF_Num_AdjFaces(MFace_ptr mface);
  List_ptr    MF_Vertices(MFace_ptr mface, int dir, MVertex_ptr mvert);
  List_ptr    MF_Edges(MFace_ptr mface, int dir, MVertex_ptr mvert);
  List_ptr    MF_AdjFaces(MFace_ptr mface);

  /* Returns 1 or 0 */
  int         MF_EdgeDir(MFace_ptr mface, MEdge_ptr medge);
  int         MF_EdgeDir_i(MFace_ptr mface, int i);


  int         MF_UsesEntity(MFace_ptr mface, MEntity_ptr mentity, int type);
  List_ptr    MF_Regions(MFace_ptr mface);
  MRegion_ptr MF_Region(MFace_ptr mface, int side);

  MFace_ptr   MVs_CommonFace(int nv, MVertex_ptr *fverts);
  MFace_ptr   MEs_CommonFace(int ne, MEdge_ptr *fedges);

  int         MFs_AreSame(MFace_ptr f1, MFace_ptr f2);


  void        MF_Coords(MFace_ptr mface, int *n, double (*xyz)[3]);
       
  void        MF_Lock(MFace_ptr f);
  void        MF_UnLock(MFace_ptr f);

#ifdef MSTK_HAVE_MPI

  PType       MF_PType(MFace_ptr f);  
  int         MF_MasterParID(MFace_ptr f);
  int         MF_GlobalID(MFace_ptr f);

#endif

/********************************************************************/
/*        MESH REGN OPERATORS                                       */
/********************************************************************/

  MRegion_ptr MR_New(Mesh_ptr mesh);
  void        MR_Delete(MRegion_ptr mregion, int keep);
  void        MR_Restore(MRegion_ptr mregion);
  void        MR_Set_RepType(MRegion_ptr mregion, RepType reptype);
  void        MR_Set_GEntity(MRegion_ptr mregion, GEntity_ptr gent);
  void        MR_Set_GEntDim(MRegion_ptr mregion, int gdim);
  void        MR_Set_GEntID(MRegion_ptr mregion, int gid);
  int         MR_Set_GInfo_Auto(MRegion_ptr mregion);
  void        MR_Set_ID(MRegion_ptr mregion, int id);

  /* Can be called by user/application after creating region by MR_New(); */
  void        MR_Set_Faces(MRegion_ptr mregion, int nf, MFace_ptr *mfaces, int *dirs);
  void        MR_Set_Vertices(MRegion_ptr mregion, int nv, MVertex_ptr *mvertices, 
		       int nf, int **rfvtemplate);

  /* Can be called by mesh modification routines */
  void        MR_Replace_Face(MRegion_ptr mregion, MFace_ptr mface, MFace_ptr nuface, int dir);
  void        MR_Replace_Vertex(MRegion_ptr mregion, MVertex_ptr mvertex, MVertex_ptr nuvertex);
  void        MR_Replace_Face_i(MRegion_ptr mregion, int i, MFace_ptr mface, int dir);
  void        MR_Rem_Face(MRegion_ptr mregion, MFace_ptr mface);
  void        MR_Replace_Vertex_i(MRegion_ptr mregion, int i, MVertex_ptr mvertex);
  int         MR_Rev_FaceDir(MRegion_ptr mregion, MFace_ptr mface);
  int         MR_Rev_FaceDir_i(MRegion_ptr mregion, int i);
  


  Mesh_ptr    MR_Mesh(MRegion_ptr mregion);
  int         MR_ID(MRegion_ptr mregion);
  int         MR_GEntDim(MRegion_ptr mregion);
  int         MR_GEntID(MRegion_ptr mregion);
  GEntity_ptr MR_GEntity(MRegion_ptr mregion);

  MRType MR_ElementType(MRegion_ptr mregion);
  int         MR_Num_Vertices(MRegion_ptr mregion);
  int         MR_Num_Edges(MRegion_ptr mregion);
  int         MR_Num_Faces(MRegion_ptr mregion);
  int         MR_Num_AdjRegions(MRegion_ptr mregion);
  List_ptr    MR_Vertices(MRegion_ptr mregion);
  List_ptr    MR_Edges(MRegion_ptr mregion);
  List_ptr    MR_Faces(MRegion_ptr mregion);
  List_ptr    MR_AdjRegions(MRegion_ptr mregion);

  /* Returns 1 or 0 */
  int         MR_FaceDir(MRegion_ptr mregion, MFace_ptr mface);
  int         MR_FaceDir_i(MRegion_ptr mregion, int i);

  int         MR_UsesEntity(MRegion_ptr mregion, MEntity_ptr ment, int type);


  void        MR_Coords(MRegion_ptr mregion, int *n, double (*xyz)[3]);

#ifdef MSTK_HAVE_MPI

  PType       MR_PType(MRegion_ptr r);  
  int         MR_MasterParID(MRegion_ptr r);
  int         MR_GlobalID(MRegion_ptr r);

#endif

  /************************************************************************/
  /* GENERIC ENTITY OPERATORS                                             */
  /************************************************************************/

  void        MEnt_Set_GEntity(MEntity_ptr mentity, GEntity_ptr gent);
  void        MEnt_Set_GEntDim(MEntity_ptr mentity, int gdim);
  void        MEnt_Set_GEntID(MEntity_ptr mentity, int gid);
  void        MEnt_Set_ID(MEntity_ptr mentity, int id);
  void        MEnt_Set_RepType(MEntity_ptr mentity, RepType rtype);

  int         MEnt_ID(MEntity_ptr mentity);
  MType       MEnt_Dim(MEntity_ptr mentity);
  MType       MEnt_OrigDim(MEntity_ptr mentity);
  int         MEnt_IsVolatile(MEntity_ptr mentity);
  Mesh_ptr    MEnt_Mesh(MEntity_ptr mentity);
  int         MEnt_GEntDim(MEntity_ptr mentity);
  int         MEnt_GEntID(MEntity_ptr mentity);
  GEntity_ptr MEnt_GEntity(MEntity_ptr mentity);
  RepType     MEnt_RepType(MEntity_ptr mentity);
  
  void        MEnt_Delete(MEntity_ptr mentity, int keep);

  /* Attributes on entities */
  
  void        MEnt_Set_AttVal(MEntity_ptr ent, MAttrib_ptr attrib, int ival, 
			double lval, void *pval);
  void        MEnt_Rem_AttVal(MEntity_ptr ent, MAttrib_ptr attrib);
  int         MEnt_Get_AttVal(MEntity_ptr ent, MAttrib_ptr attrib, int *ival, 
			double *lval, void **pval);  
  void        MEnt_Clear_AttVal(MEntity_ptr ent, MAttrib_ptr attrib);
  void        MEnt_Rem_AllAttVals(MEntity_ptr);


#ifdef MSTK_HAVE_MPI

  PType       MEnt_PType(MEntity_ptr ent);
  int         MEnt_MasterParID(MEntity_ptr ent);
  int         MEnt_GlobalID(MEntity_ptr ent);

#endif


  /************************************************************************/
  /* ENTITY SET OPERATORS                                             */
  /************************************************************************/

  MSet_ptr    MSet_New(Mesh_ptr mesh, const char *set_name, MType entdim);
  char       *MSet_Name(MSet_ptr set, char *set_name);
  Mesh_ptr    MSet_Mesh(MSet_ptr set);
  MType       MSet_EntDim(MSet_ptr set);
  void        MSet_Delete(MSet_ptr set);

  MSet_ptr    MSet_Add(MSet_ptr set, void *entry);
  MSet_ptr    MSet_ChknAdd(MSet_ptr set, void *entry);
  MSet_ptr    MSet_Insert(MSet_ptr set, void *nuentry, void *b4entry);
  MSet_ptr    MSet_Inserti(MSet_ptr set, void *nuentry, int i);
  int         MSet_Rem(MSet_ptr set, void *entry);
  int         MSet_Remi(MSet_ptr set, int i);
  int         MSet_Replace(MSet_ptr set, void *entry, void *nuentry);
  int         MSet_Replacei(MSet_ptr set, int i, void *nuentry);
  int         MSet_Contains(MSet_ptr set, void *entry);
  int         MSet_Locate(MSet_ptr set, void *entry);
  void       *MSet_Entry(MSet_ptr set, int i);
  void       *MSet_Next_Entry(MSet_ptr set, int *i);
  int         MSet_Num_Entries(MSet_ptr set);
  MSet_ptr    MSet_Cat(MSet_ptr dest, MSet_ptr src);
  MSet_ptr    MSet_Copy(MSet_ptr oldset);

  MSet_ptr    MSets_Union(MSet_ptr s1, MSet_ptr s2);
  MSet_ptr    MSets_Intersect(MSet_ptr s1, MSet_ptr s2);
  MSet_ptr    MSets_Subtract(MSet_ptr s1, MSet_ptr s2);

#ifdef DEBUG
  void        MSet_Print(MSet_ptr set);
#endif



  /************************************************************************/
  /* ENTITY MARKING                                                       */
  /************************************************************************/

  int         MSTK_GetMarker(void);
  void        MSTK_FreeMarker(int mkr);
  void        MEnt_Mark(MEntity_ptr ent, int mkr);
  int         MEnt_IsMarked(MEntity_ptr ent, int mkr);
  void        MEnt_Unmark(MEntity_ptr ent, int mkr);
  void        List_Mark(List_ptr list, int mkr);
  void        List_Unmark(List_ptr list, int mkr);
  void        MSet_Mark(MSet_ptr set, int mkr);
  void        MSet_Unmark(MSet_ptr set, int mkr);


  /************************************************************************/
  /* ATTRIBUTE DEFINITION                                                 */
  /************************************************************************/

  /*
    MAttrib_ptr MAttrib_New(Mesh_ptr mesh, const char *att_name, 
    MAttType att_type, MType entdim, ...);

    When the attribute type is INT, DOUBLE or POINTER, the calling routine
    can end the argument list at entdim, like so:

    myatt = MAttrib_New(mesh,"anatt",INT,MVERTEX);

    When the attribute type is VECTOR OR DOUBLE, the calling routine has
    to specify the number of components as the last argument like so:	     

    myatt = MAttrib_New(mesh,"stress_tensor",TENSOR,MREGION,21);
  */

  MAttrib_ptr MAttrib_New(Mesh_ptr mesh, const char *att_name, 
			  MAttType att_type, MType entdim, ...);
  char       *MAttrib_Get_Name(MAttrib_ptr attrib, char *att_name);
  MAttType    MAttrib_Get_Type(MAttrib_ptr attrib);
  MType       MAttrib_Get_EntDim(MAttrib_ptr attrib);
  int         MAttrib_Get_NumComps(MAttrib_ptr attrib);
  void        MAttrib_Delete(MAttrib_ptr attrib);
  void        MAttrib_Clear(MAttrib_ptr attrib);


/*************************************************************************/
/*  SOME GENERIC OPERATORS FOR MESH REGIONS BASED ON TYPE AND CONVENTION */
/*************************************************************************/

  int         RType_NumVerts(MRType type);
  int         RType_NumFaces(MRType type);
  int         RType_NumFVerts(MRType type, int locfnum);
  int         RType_LocFVNums(MRType type, int locfnum, MVertex_ptr *lverts);


/***********************************************************************/
/* Mesh Modification Operators                                         */
/***********************************************************************/

  int         ME_Swap2D(MEdge_ptr e, MEdge_ptr *enew, MFace_ptr fnew[2]);
  MVertex_ptr MVs_Merge(MVertex_ptr v1, MVertex_ptr v2); /* v2 is deleted */
  MEdge_ptr   MEs_Merge(MEdge_ptr e1, MEdge_ptr e2); /* e2 is deleted */
  MFace_ptr   MFs_Merge(MFace_ptr f1, MFace_ptr f2); /* f2 is deleted */
  MFace_ptr   MFs_Join(MFace_ptr f1, MFace_ptr f2, MEdge_ptr e);
  MVertex_ptr ME_Collapse(MEdge_ptr e, MVertex_ptr ovkeep, int topoflag);
  

/**********************************************************************/
/* More parallel operators                                            */
/**********************************************************************/

#ifdef MSTK_HAVE_MPI

  /* ROUTINES FOR MORE FINE-GRAINED CONTROL OF PARALLEL APPLICATION  */
  /* IF YOU CALL THESE ROUTINES WITHOUT KNOWING YOUR WAY AROUND      */
  /* YOU WILL GET WHAT YOU DESERVE                                   */

  int         MESH_Num_GhostVertices(Mesh_ptr mesh);
  int         MESH_Num_GhostEdges(Mesh_ptr mesh);
  int         MESH_Num_GhostFaces(Mesh_ptr mesh);
  int         MESH_Num_GhostRegions(Mesh_ptr mesh);

  int         MESH_Num_OverlapVertices(Mesh_ptr mesh);
  int         MESH_Num_OverlapEdges(Mesh_ptr mesh);
  int         MESH_Num_OverlapFaces(Mesh_ptr mesh);
  int         MESH_Num_OverlapRegions(Mesh_ptr mesh);

  MVertex_ptr MESH_GhostVertex(Mesh_ptr mesh, int i);
  MEdge_ptr   MESH_GhostEdge(Mesh_ptr mesh, int i);
  MFace_ptr   MESH_GhostFace(Mesh_ptr mesh, int i);
  MRegion_ptr MESH_GhostRegion(Mesh_ptr mesh, int i);

  MVertex_ptr MESH_OverlapVertex(Mesh_ptr mesh, int i);
  MEdge_ptr   MESH_OverlapEdge(Mesh_ptr mesh, int i);
  MFace_ptr   MESH_OverlapFace(Mesh_ptr mesh, int i);
  MRegion_ptr MESH_OverlapRegion(Mesh_ptr mesh, int i);

  MVertex_ptr MESH_Next_GhostVertex(Mesh_ptr mesh, int *index);
  MEdge_ptr   MESH_Next_GhostEdge(Mesh_ptr mesh, int *index);
  MFace_ptr   MESH_Next_GhostFace(Mesh_ptr mesh, int *index);
  MRegion_ptr MESH_Next_GhostRegion(Mesh_ptr mesh, int *index);

  MVertex_ptr MESH_Next_OverlapVertex(Mesh_ptr mesh, int *index);
  MEdge_ptr   MESH_Next_OverlapEdge(Mesh_ptr mesh, int *index);
  MFace_ptr   MESH_Next_OverlapFace(Mesh_ptr mesh, int *index);
  MRegion_ptr MESH_Next_OverlapRegion(Mesh_ptr mesh, int *index);

  /*end for mpi */

#endif


/**********************************************************************/
/* Element quality evaluation                                         */
/**********************************************************************/

  void        MF_CondNums(MFace_ptr f, int *nfv, double *condnums);
  void        MR_CondNums(MRegion_ptr r, int *nrv, double *condnums);


/**********************************************************************/
/* MESH + GEOMETRY                                                    */
/**********************************************************************/

  double      MFs_DihedralAngle(MFace_ptr face1, MFace_ptr face2, MEdge_ptr edge);
  double      MEs_Angle(MEdge_ptr e1, MEdge_ptr e2);


/*******************************************/
/* UTILITIES                               */
/*******************************************/

  void        MSTK_Report(const char *module, const char *message, ErrType severity);


/*******************************************/
/* LISTS                                   */
/*******************************************/

  List_ptr    List_New(int inisize);
  void        List_Delete(List_ptr list);
  List_ptr    List_Compress(List_ptr list);
  List_ptr    List_Copy(List_ptr list);
  
  List_ptr    List_Add(List_ptr l, void *entry);
  List_ptr    List_ChknAdd(List_ptr l, void *entry);
  List_ptr    List_Insert(List_ptr l, void *nuentry, void *b4entry);
  List_ptr    List_Inserti(List_ptr l, void *nuentry, int i);
  int         List_Rem(List_ptr l, void *entry);
  int         List_Remi(List_ptr l, int i);
  int         List_RemSorted(List_ptr l, void *entry, int (*entry2int)(void *));
  int         List_Replace(List_ptr l, void *entry, void *nuentry);
  int         List_Replacei(List_ptr l, int i, void *nuentry);
  int         List_Contains(List_ptr l, void *entry);
  int         List_Locate(List_ptr l, void *entry);
  void       *List_Entry(List_ptr l, int i);
  void       *List_Next_Entry(List_ptr l, int *i);
  int         List_Num_Entries(List_ptr l);
  List_ptr    List_Cat(List_ptr dest, List_ptr src);
  /* Sort a list based on a user/application supplied comparison function */
  void        List_Sort(List_ptr l, size_t num, size_t size,
		     int(*comp)(const void *,const void *));
  /* Search a sorted list for a key based on a user/application supplied
     comparison function - unpredictable results for unsorted list */
  void       *List_Search(List_ptr l, const void *key, size_t num, size_t size,
                          int(*comp)(const void *,const void *));
#ifdef DEBUG
  void        List_Print(List_ptr l);
#endif

  /* Extra functionality for hash-tables */

  void*      *List_Entries(List_ptr l);


  /* Functions to print information about mesh entities. The argument
   level controls how detailed the information printed about the
   entity is. If level is 0, then minimal information is printed about
   the entity.  If it is 1, then classification information is also
   printed for the entity to the extent available. If it is 2 or
   higher, then adjacency information is printed for the entity */

  void        MV_Print(MVertex_ptr v, int level);
  void        ME_Print(MEdge_ptr e, int level);
  void        MF_Print(MFace_ptr f, int level);
  void        MR_Print(MRegion_ptr, int level);

  /* Mesh Statistics. If 'level' is 1, only entity counts are
     printed. If 'level' is 2 or more, then element quality numbers
     are also printed. If it is 3 or more, additional info about the
     worst quality element is printed. The routine returns silently,
     if 'level' is less than 1 */

  void        MESH_PrintStats(Mesh_ptr mesh, int level);
  

#ifdef __cplusplus
	   }
#endif

#endif
