/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#ifndef _H_MSTK
#define _H_MSTK

/* Main Header file defining the MSTK API. 

   Author: Rao V Garimella
   Email: rao@lanl.gov
   Website: https://software.lanl.gov/MeshTools/hg/mstk

   The MSTK source is available for use under an LGPL license 
*/



#include <stdarg.h>
#include <stddef.h>

#include "MSTK_defines.h"
#include "MSTK_types.h"
#include "MSTK_externs.h"
#include "MSTK_util.h"

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

  /* Retrieve version info                                           */
  void MSTK_Version(int *major_version, int *minor_version, int *patch_version,
                    char **version_string);

  
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

  /* Generate a structured mesh in 2D or in 3D */

  Mesh_ptr    MESH_Gen_Structured(double lx, double ly, double lz,
				  double ux, double uy, double uz,
				  int nx, int ny, int nz);
  
  
  /* Import mesh data into the mesh object from various
     formats. 'comm' can be NULL for serial codes. 'opts' is an
     integer array with different meanings for different
     formats. 'format' can be 'mstk' for the native format, 'gmv' for
     the General Mesh Viewer format
     (http://www.generalmeshviewer.com), 'exo' for the SEACAS/ExodusII
     format from Sandia National Laboratories
     (https://github.com/gsjaardema/seacas), 'par' for parallel Exodus
     files in the SEACAS/NemesisI format
     (https://github.com/gsjaardema/seacas) */

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
  /* Get element graph from Exodus II file                                    */
  /* adjbeg - pointer to array containing the offsets of the adjelems list at */
  /*          which the adjacent elements are enumerated for each element     */
  /*--------------------------------------------------------------------------*/
  
  int        ExodusII_GetElementGraph(const char *filename, int *ndim, int *nelems,
				      int **adjbeg, int **adjelems, int get_coords,
				      double (**elemcen)[3]);

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
     serial codes. 'format' is the mesh file format to be exported
     to. 'natt' indicates how many attributes to export (-1 means
     export all attributes present on the mesh), 'attnames' is the
     list of attributes to export (NULL if 'natt' is -1), 'opts' is a
     list of format specific output options (can be NULL). 'format'
     can be 'mstk' for the native format, 'gmv' for the General Mesh
     Viewer format (http://www.generalmeshviewer.com), 'exo' for the
     SEACAS/ExodusII format from Sandia National Laboratories
     (https://github.com/gsjaardema/seacas), 'x3d' for the FLAG X3D
     format, 'stl' for the STL format. If the number of ranks is
     greater than 1, selecting the 'exo' format will write out
     parallel Exodus files (or in other words Nemesis I files) with
     the extension .par.N.n */

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

  int         MESH_BuildClassfn(Mesh_ptr mesh, int use_geometry);
  int         MESH_DelInterior(Mesh_ptr mesh);
  int         MESH_Tet2Hex(Mesh_ptr tetmesh, Mesh_ptr *hexmesh);
  int         MESH_Tri2Quad(Mesh_ptr trimesh, Mesh_ptr *quadmesh);

  /*--------------------------------------------------------------------------*/
  /* Get element graph from mesh                                              */
  /* adjbeg - pointer to array containing the offsets of the adjelems list at */
  /*          which the adjacent elements are enumerated for each element     */
  /*--------------------------------------------------------------------------*/
  
  int        MESH_GetElementGraph(Mesh_ptr mesh, int *nelems,
				  int **adjbeg, int **adjelems);

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


  /* Search for entities by Local ID. If the mesh has been modified or
   * renumbered, it is recommended that MESH_Enable_LocalIDSearch be
   * called before a block of code that does such searches - the call
   * will sort lists for faster searching */

  void MESH_Enable_LocalIDSearch(Mesh_ptr mesh);

  MVertex_ptr MESH_VertexFromID(Mesh_ptr mesh, int i);
  MEdge_ptr   MESH_EdgeFromID(Mesh_ptr mesh, int i);
  MFace_ptr   MESH_FaceFromID(Mesh_ptr mesh, int i);
  MRegion_ptr MESH_RegionFromID(Mesh_ptr mesh, int i);
  MEntity_ptr MESH_EntityFromID(Mesh_ptr mesh, MType mtype, int i);

  
  int         MESH_SetRepType(Mesh_ptr mesh, RepType type);
  void        MESH_SetGModel(Mesh_ptr mesh, GModel_ptr geom);
  int         MESH_Change_RepType(Mesh_ptr mesh, int nurep);

  /* Renumber all mesh entities so that they have contiguous IDs */
  /* renum_type = 0 --- simple sequential numbering
                = 1 --- renumbering using Reverse Cuthill-McKee algorithm 
     mtype      = Type of entity to renumber (MVERTEX, MEDGE, MFACE, MREGION
                  or MALLTYPE)

     PER PROCESSOR RENUMBERING ONLY! See Mesh_RenumberGlobalIDs for
     consistent global ID renumbering of distributed meshes
  */

  void        MESH_Renumber(Mesh_ptr mesh, int renum_type, MType mtype);


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
     safest to set mymesh to NULL. If the del_inmesh flag is non-zero,
     the input mesh 'globalmesh' will be deleted soon after partitioning
     and creation of submeshes - this will save memory */

  int         MSTK_Mesh_Distribute(Mesh_ptr globalmesh, Mesh_ptr *mymesh, 
                                   int *topodim, int ring, int with_attr, 
                                   int method, int del_inmesh, MSTK_Comm comm);



  /* 'Weave' a set of distributed mesh partitions together to build
   parallel connections and ghost info
   
   input_type indicates what info is already present on the mesh

   0 -- we are give NO information about how these meshes are
   connected other than the knowledge that they come from partitioning
   of a single mesh

   1 -- we are given partitioned meshes with a unique global ID on
   each mesh vertex

   2 -- we are given parallel neighbor information, but no global ID
   on each mesh vertex
  */

  int         MSTK_Weave_DistributedMeshes(Mesh_ptr mesh, int topodim,
                                           int num_ghost_layers,
                                           int input_type, MSTK_Comm comm);

  /* Parallel update attribute values for ghost entities */


  int        MESH_UpdateAttributes(Mesh_ptr mesh, MSTK_Comm comm);
  int        MESH_Update1Attribute(Mesh_ptr mesh, MAttrib_ptr attrib,
                                   MSTK_Comm comm);

  /* Parallel gather-scatter of attributes */
  int        MESH_GathScat1Attribute(Mesh_ptr mesh, MAttrib_ptr attrib,
	                             MAttOpType optype, MSTK_Comm comm);

  /* Keep this version of the call around for backward compatibility */
  int        MSTK_UpdateAttr(Mesh_ptr mesh, MSTK_Comm comm); 


  /* Update vertex coordinates for ghost vertices */
  int         MESH_UpdateVertexCoords(Mesh_ptr mesh, MSTK_Comm comm);

  /* Check parallel consistency */

  int         MESH_Parallel_Check(Mesh_ptr mesh, MSTK_Comm comm);

  /* Query Global IDs - To search for entities by GlobalIDs, one has
   * to enable them first. Since it takes up additional storage, it
   * should be disabled as soon as the need is done. */

  void       MESH_Enable_GlobalIDSearch(Mesh_ptr mesh);
  void       MESH_Disable_GlobalIDSearch(Mesh_ptr mesh);

  MVertex_ptr MESH_VertexFromGlobalID(Mesh_ptr mesh, int global_id);
  MEdge_ptr   MESH_EdgeFromGlobalID(Mesh_ptr mesh, int global_id);
  MFace_ptr   MESH_FaceFromGlobalID(Mesh_ptr mesh, int global_id);
  MRegion_ptr MESH_RegionFromGlobalID(Mesh_ptr mesh, int global_id);
  MEntity_ptr MESH_EntityFromGlobalID(Mesh_ptr mesh, MType mtype, int i);

  /* Renumber global IDs in a parallel mesh */
  /* If mtype = MALLTYPE, all entity types are renumbered and made contiguous */
  /* Only method = 0 (sequential) is supported for now                        */
  /* Eventually, one could precompute some global IDs and just send it        */
  /* to this routine (NOT IMPLEMENTED)                                        */
  /* Cannot use mtype = MALLTYPE for preassigned GIDs                             */

  int         MESH_Renumber_GlobalIDs(Mesh_ptr mesh, MType mtype, int method,
                                      int *preassigned_gids, MSTK_Comm comm);

  /* Get a partitioning for mesh using METIS (method=1) or ZOLTAN (method=2) */
  /* Doesn't actually partition the mesh or distribute it                    */

  int        MESH_Get_Partitioning(Mesh_ptr mesh, int method, int **part, 
                                   MPI_Comm comm);
  int        MSTK_Mesh_Partition(Mesh_ptr mesh, int num, int *part,  int ring, 
	                         int with_attr, int del_inmesh,
                                 Mesh_ptr *submeshes);

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
  void        MV_AdjVertexIDs(MVertex_ptr mvertex, int *nvadj, int *adjvertids);
  List_ptr    MV_Edges(MVertex_ptr mvertex);
  void        MV_EdgeIDs(MVertex_ptr mvertex, int *nve, int *vedgeids);
  List_ptr    MV_Faces(MVertex_ptr mvertex);
  void        MV_FaceIDs(MVertex_ptr mvertex, int *nvf, int *vfaceids);
  List_ptr    MV_Regions(MVertex_ptr mvertex);
  void        MV_RegionIDs(MVertex_ptr mvertex, int *nvr, int *vregionids);

  int         MV_GlobalID(MVertex_ptr v);
  void        MV_Set_GlobalID(MVertex_ptr v, int gid);  /* ** CAUTION!! ** */

  PType       MV_PType(MVertex_ptr v);    /* PINTERIOR, PGHOST, POVERLAP */
  int         MV_OnParBoundary(MVertex_ptr v); 
  int         MV_MasterParID(MVertex_ptr v);

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
  int         ME_VertexID(MEdge_ptr medge, int i);
  MVertex_ptr ME_OppVertex(MEdge_ptr medge, MVertex_ptr v);
  int         ME_OppVertexID(MEdge_ptr medge, int vid);
  int         ME_UsesEntity(MEdge_ptr medge, MEntity_ptr mentity, int etype);
  List_ptr    ME_Faces(MEdge_ptr medge);
  void        ME_FaceIDs(MEdge_ptr medge, int *nef, int *efaceids);
  List_ptr    ME_Regions(MEdge_ptr medge);
  void        ME_RegionIDs(MEdge_ptr medge, int *ner, int *eregionids);


  MEdge_ptr   MVs_CommonEdge(MVertex_ptr v1, MVertex_ptr v2);  

  double      ME_Len(MEdge_ptr e);
  double      ME_LenSqr(MEdge_ptr e);
  void        ME_Vec(MEdge_ptr e, double *evec);

  int         MEs_AreSame(MEdge_ptr e1, MEdge_ptr e2);

  void        ME_Lock(MEdge_ptr e);
  void        ME_UnLock(MEdge_ptr e);

  int         ME_GlobalID(MEdge_ptr e);
  void        ME_Set_GlobalID(MEdge_ptr e, int gid);  /* ** CAUTION!! ** */

  PType       ME_PType(MEdge_ptr e);
  int         ME_OnParBoundary(MEdge_ptr e);
  int         ME_MasterParID(MEdge_ptr e);

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
  void        MF_Rem_Edge(MFace_ptr mface, MEdge_ptr medge);
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

  MFType      MF_ElementType(MFace_ptr mf);
  int         MF_Num_Vertices(MFace_ptr mface);
  int         MF_Num_Edges(MFace_ptr mface);
  int         MF_Num_AdjFaces(MFace_ptr mface);
  List_ptr    MF_Vertices(MFace_ptr mface, int dir, MVertex_ptr startvert);
  void        MF_VertexIDs(MFace_ptr mface, int dir, int startvertid, int *nfv, int *fvertids);
  List_ptr    MF_Edges(MFace_ptr mface, int dir, MVertex_ptr startvert);
  void        MF_EdgeIDs(MFace_ptr mface, int dir, int startvertid, int *nfe, int *fedgeids);
  List_ptr    MF_AdjFaces(MFace_ptr mface);
  /* void        MF_AdjFaceIDs(MFace_ptr mface, int *nadjf, int *adjfaceids); */

  /* Returns 1 or 0 for valid queries, -1 otherwise */
  int         MF_EdgeDir(MFace_ptr mface, MEdge_ptr medge);
  int         MF_EdgeDir_i(MFace_ptr mface, int i);


  int         MF_UsesEntity(MFace_ptr mface, MEntity_ptr mentity, int type);
  List_ptr    MF_Regions(MFace_ptr mface);
  void        MF_RegionIDs(MFace_ptr mface, int *nfr, int *fregids);
  MRegion_ptr MF_Region(MFace_ptr mface, int side);
  int         MF_RegionID(MFace_ptr mface, int side);

  MFace_ptr   MVs_CommonFace(int nv, MVertex_ptr *fverts);
  MFace_ptr   MEs_CommonFace(int ne, MEdge_ptr *fedges);

  int         MFs_AreSame(MFace_ptr f1, MFace_ptr f2);


  void        MF_Coords(MFace_ptr mface, int *n, double (*xyz)[3]);
       
  void        MF_Lock(MFace_ptr f);
  void        MF_UnLock(MFace_ptr f);

  int         MF_GlobalID(MFace_ptr f);
  void        MF_Set_GlobalID(MFace_ptr f, int gid);  /* ** CAUTION!! ** */

  PType       MF_PType(MFace_ptr f);  
  int         MF_OnParBoundary(MFace_ptr f);
  int         MF_MasterParID(MFace_ptr f);

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
  void        MR_Replace_Faces(MRegion_ptr mregion, int nold, MFace_ptr *oldf, 
                               int nnu, MFace_ptr *nuf, int *nudir);
  void        MR_Replace_Vertex(MRegion_ptr mregion, MVertex_ptr mvertex, MVertex_ptr nuvertex);
  void        MR_Replace_Face_i(MRegion_ptr mregion, int i, MFace_ptr mface, int dir);
  void        MR_Rem_Face(MRegion_ptr mregion, MFace_ptr mface);
  void        MR_Replace_Vertex_i(MRegion_ptr mregion, int i, MVertex_ptr mvertex);

  /* Reverse the direction in which the face is used by the region
   * (and the adjacent region, if it exists) */

  int         MR_Rev_FaceDir(MRegion_ptr mregion, MFace_ptr mface);
  int         MR_Rev_FaceDir_i(MRegion_ptr mregion, int i);


  /* Set the direction in which the face is used by the region. Also,
   * attach the region on the correct side of the face consistent with
   * the specified direction. If a region already exists on that side,
   * its overwritten */

  int         MR_Set_FaceDir(MRegion_ptr mregion, MFace_ptr mface, int dir);
  int         MR_Set_FaceDir_i(MRegion_ptr mregion, int i, int dir);


  Mesh_ptr    MR_Mesh(MRegion_ptr mregion);
  int         MR_ID(MRegion_ptr mregion);
  int         MR_GEntDim(MRegion_ptr mregion);
  int         MR_GEntID(MRegion_ptr mregion);
  GEntity_ptr MR_GEntity(MRegion_ptr mregion);

  MRType      MR_ElementType(MRegion_ptr mregion);
  int         MR_Num_Vertices(MRegion_ptr mregion);
  int         MR_Num_Edges(MRegion_ptr mregion);
  int         MR_Num_Faces(MRegion_ptr mregion);
  int         MR_Num_AdjRegions(MRegion_ptr mregion);
  List_ptr    MR_Vertices(MRegion_ptr mregion);
  void        MR_VertexIDs(MRegion_ptr mregion, int *nrv, int *rvertids);
  List_ptr    MR_Edges(MRegion_ptr mregion);
  void        MR_EdgeIDs(MRegion_ptr mregion, int *nre, int *redgeids);
  List_ptr    MR_Faces(MRegion_ptr mregion);
  void        MR_FaceIDs(MRegion_ptr mregion, int *nrf, int *rfaceids);
  List_ptr    MR_AdjRegions(MRegion_ptr mregion);
  void        MR_AdjRegionIDs(MRegion_ptr mregion, int *nradj, int *adjregids);


  /* Returns 1 or 0 for valid queries, -1 otherwise */
  int         MR_FaceDir(MRegion_ptr mregion, MFace_ptr mface);
  int         MR_FaceDir_i(MRegion_ptr mregion, int i);

  int         MR_UsesEntity(MRegion_ptr mregion, MEntity_ptr ment, int type);


  void        MR_Coords(MRegion_ptr mregion, int *n, double (*xyz)[3]);

  int         MR_GlobalID(MRegion_ptr r);
  void        MR_Set_GlobalID(MRegion_ptr r, int gid);  /* ** CAUTION!! ** */

  PType       MR_PType(MRegion_ptr r);  
  int         MR_MasterParID(MRegion_ptr r);

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
  void        MEnt_Set_IntAttVal(MEntity_ptr ent, MAttrib_ptr attrib, int ival);
  void        MEnt_Set_DblAttVal(MEntity_ptr ent, MAttrib_ptr attrib,
                                 double rval);
  void        MEnt_Set_PtrAttVal(MEntity_ptr ent, MAttrib_ptr attrib,
                                 void *pval);
  void        MEnt_Set_VecAttVal(MEntity_ptr ent, MAttrib_ptr attrib,
                                 double *vval);
  void        MEnt_Set_TnsrAttVal(MEntity_ptr ent, MAttrib_ptr attrib,
                                  double *tval);
  int         MEnt_Get_AttVal(MEntity_ptr ent, MAttrib_ptr attrib, int *ival,
                              double *lval, void **pval);
  int         MEnt_Get_IntAttVal(MEntity_ptr ent, MAttrib_ptr attrib);
  double      MEnt_Get_DblAttVal(MEntity_ptr ent, MAttrib_ptr attrib);
  void *      MEnt_Get_PtrAttVal(MEntity_ptr ent, MAttrib_ptr attrib);
  double *    MEnt_Get_VecAttVal(MEntity_ptr ent, MAttrib_ptr attrib);
  double *    MEnt_Get_TnsrAttVal(MEntity_ptr ent, MAttrib_ptr attrib);
  void        MEnt_Rem_AttVal(MEntity_ptr ent, MAttrib_ptr attrib);
  void        MEnt_Clear_AttVal(MEntity_ptr ent, MAttrib_ptr attrib);
  void        MEnt_Rem_AllAttVals(MEntity_ptr);


  int         MEnt_GlobalID(MEntity_ptr ent);
  void        MEnt_Set_GlobalID(MEntity_ptr ent, int gid);  /* caution */

  PType       MEnt_PType(MEntity_ptr ent);
  int         MEnt_OnParBoundary(MEntity_ptr ent);
  int         MEnt_MasterParID(MEntity_ptr ent);


  /************************************************************************/
  /* ENTITY SET OPERATORS                                             */
  /************************************************************************/

  MSet_ptr    MSet_New(Mesh_ptr mesh, const char *set_name, MType entdim);
  char       *MSet_Name(MSet_ptr set, char *set_name);
  Mesh_ptr    MSet_Mesh(MSet_ptr set);
  MType       MSet_EntDim(MSet_ptr set);
  void        MSet_Delete(MSet_ptr set);

  void        MSet_Rename(MSet_ptr set, char *newname);
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
  /*                                                                      */
  /* Entity marking can be used to avoid searching through lists when     */
  /* merging two lists, when walking through the mesh and other such      */
  /* situations. Entity marking is available only if the build is         */
  /* configured with MSTK_USE_MARKERS=On (Default: Off)                   */
  /*                                                                      */
  /* IN MULTI-THREADED SETTINGS THERE ARE 3 IMPORTANT POINTS TO NOTE:     */
  /* (1) ENTITY MARKING USES pthread_mutex_lock/unlock WHICH CAN KEEP THE */
  /* CODE FROM SCALING WELL (2) SOME MULTI-THREADING PARADIGMS MAY NOT    */
  /* USE PTHREADS (3) THERE IS A LIMIT ON THE NUMBER MARKERS THAT CAN BE  */
  /* IN USE AT A TIME AND IT IS EASY TO HIT THE LIMIT IN MULTI-THREADED   */
  /* RUNS. SO IT IS ADVISABLE TO USE MARKERS ONLY FOR SINGLE-THREADED     */
  /* APPLICATIONS                                                         */
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

  /* Merge vertices v1 and v2 - v2 is deleted */
  /* topoflag = 1 means respect model topology (as given by mesh entity
   classification) and do not allow dimensional reduction in the mesh.
   topoflag = 0 means the caller does not care about classification
   and would like the two faces to be merged as requested */

  MVertex_ptr MVs_Merge(MVertex_ptr v1, MVertex_ptr v2, int topoflag);

  /* Merge edges e1 and e2 - e2 is deleted */
  /* See MVs_Merge for meaning of topoflag */
  MEdge_ptr   MEs_Merge(MEdge_ptr e1, MEdge_ptr e2, int topoflag); 

  int         ME_Swap2D(MEdge_ptr e, MEdge_ptr *enew, MFace_ptr fnew[2]);

  /* Low level edge split - split edge only and incorporate into adjacent
     faces/regions - works for any type of mesh */
  MVertex_ptr ME_Split(MEdge_ptr esplit, double *xyz);

  /* Low level edge split with multiple points - split edge only and
     incorporate into adjacent faces/regions - works for any type of
     mesh. RETURN IS A LIST OF NEW EDGES. 'xyz' can be NULL in which
     case, split points are equispaced */

  List_ptr    ME_MultiSplit(MEdge_ptr esplit, int n, double (*xyz)[3]);

  /* High level edge split for tri/tet meshes - split all connected elements too */
  MVertex_ptr ME_Split_SimplexMesh(MEdge_ptr esplit, double *xyz);

  /* collapse out edge and delete any faces and regions connected to it */
  /* See MVs_Merge for meaning of topoflag                              */
  /* deleted_entities contains a list of all entities deleted due to    */
  /* this operation. merged_entity_pairs contains info about merged     */
  /* entities by storing the deleted and kept entity in each merged pair*/
  /* (in that order)                                                    */

  MVertex_ptr ME_Collapse(MEdge_ptr e, MVertex_ptr ovkeep, int topoflag,
                          List_ptr *deleted_entities,
                          List_ptr *merged_entity_pairs);

  /* Merge faces f1 and f2 - f2 is deleted                      */
  /* See MVs_Merge for meaning of topoflag                      */
  MFace_ptr   MFs_Merge(MFace_ptr f1, MFace_ptr f2, int topoflag); 


  /* join two faces along a common edge to form a larger face that is
     bounded by the edges of both faces except the common edge. The
     common edge and the two faces are deleted and a new face is
     returned */

  MFace_ptr   MFs_Join(MFace_ptr f1, MFace_ptr f2, MEdge_ptr e);

  /* join two regions along a common face to form a larger region that
     is bounded by the faces of both region except the common
     face. The common face and the second region are deleted and the
     first region is returned with the extra faces (NOTE DIFFERENCE
     FROM MFs_Join) */

  MRegion_ptr   MRs_Join(MRegion_ptr r1, MRegion_ptr r2, MFace_ptr f);

  /* Low level face split - split a face along edge connecting vertex vnew0
     and vnew1 and incorporate new faces into connected regions -
     works for any type of mesh */
  MEdge_ptr   MF_Split_with_Edge(MFace_ptr fsplit, MVertex_ptr vnew0, 
                                 MVertex_ptr vnew1);

  /* Low level face split - split a face at a given location, creating
     new triangular faces from the split vertex and each of the face's
     edges. Incorporate the new faces into the connected regions */

  MVertex_ptr MF_Split(MFace_ptr fsplit, double *xyz);

  /*  High level split for tri/tet meshes - split all connected elements too */
  MVertex_ptr MF_Split_SimplexMesh(MFace_ptr fsplit, double *splitxyz);

  /* Split a region into two given a loop of edges that would form an
     face intersecting the edge - works for any type of mesh */

  MFace_ptr   MR_Split_with_EdgeLoop(MRegion_ptr rsplit, int nfe, MEdge_ptr *fedges);



/**********************************************************************/
/* More parallel operators                                            */
/**********************************************************************/

#ifdef MSTK_HAVE_MPI

  /* ROUTINES FOR MORE FINE-GRAINED CONTROL OF PARALLEL APPLICATION  */
  /* IF YOU CALL THESE ROUTINES WITHOUT KNOWING YOUR WAY AROUND      */
  /* YOU WILL GET WHAT YOU DESERVE                                   */

  /* The MESH_Ghost_Vertex routine will return all vertices of ghost
     faces (2D) or ghost regions (3D) that are not owned by the
     current processors. This means that if the mesh has 1 layer of
     ghost elements, the ghost vertex list will include vertices that
     not on necessarily only on the partition boundaries. The same
     applies to edges (in 2D and 3D) and faces (in 3D) */

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
