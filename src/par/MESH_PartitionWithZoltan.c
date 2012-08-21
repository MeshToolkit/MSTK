#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "zoltan.h"

#include "MSTK.h"

/* Get partitioning of mesh using Zoltan. Doesn't actually do anything 
   the mesh */


#ifdef __cplusplus
extern "C" {
#endif

typedef struct{
  int numMyVertices; /* total vertices in in my partition */
  int numAllNbors;   /* total number of neighbors of my vertices */
  ZOLTAN_ID_TYPE *vertexGID;    /* global ID of each of my vertices */
  int *nborIndex;    /* nborIndex[i] is location of start of neighbors for vertex i */
  ZOLTAN_ID_TYPE *nborGID;      /* nborGIDs[nborIndex[i]] is first neighbor of vertex i */
  int *nborProc;     /* process owning each nbor in nborGID */
} GRAPH_DATA;
  
  /* Application defined query functions */
  
  static int get_number_of_vertices(void *data, int *ierr);
  static void get_vertex_list(void *data, int sizeGID, int sizeLID,
			      ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
			      int wgt_dim, float *obj_wgts, int *ierr);
  static void get_num_edges_list(void *data, int sizeGID, int sizeLID,
				 int num_obj,
				 ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
				 int *numEdges, int *ierr);                                                                                                           
  static void get_edge_list(void *data, int sizeGID, int sizeLID,                                                                                   
			    int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,                                                                               
			    int *num_edges,                                                                                                                           
			    ZOLTAN_ID_PTR nborGID, int *nborProc,                                                                                                     
			    int wgt_dim, float *ewgts, int *ierr);     

  
  int MESH_PartitionWithZoltan(Mesh_ptr mesh, int nparts, int **part) { 

  MEdge_ptr fedge;
  MFace_ptr mf, oppf, rface;
  MRegion_ptr mr, oppr;
  List_ptr fedges, efaces, rfaces, fregions;
  int  i, ncells;
  int  nv, ne, nf, nr, nfe, nef, nfr, nrf, idx, idx2;
  int  numflag, nedgecut, ipos;
  int  wtflag, metisopts[5] = {0,0,0,0,0};

  int rc;
  float ver;
  struct Zoltan_Struct *zz;
  GRAPH_DATA graph;
  int changes, numGidEntries, numLidEntries, numImport, numExport;
   ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
  int *importProcs, *importToPart, *exportProcs, *exportToPart;

  MPI_Comm comm = MSTK_Comm();
  int rank = MSTK_Comm_rank();
 
  printf("zoltan on processor %d\n",rank);
  rc = Zoltan_Initialize(0, NULL, &ver);

  if (rc != ZOLTAN_OK){
    printf("sorry...\n");
    MPI_Finalize();
    exit(0);
  }

  /******************************************************************
  ** Create a Zoltan library structure for this instance of partition 
  ********************************************************************/
  zz = Zoltan_Create(mpi_comm);

  if(rank == 0) {
    nv = MESH_Num_Vertices(mesh);
    ne = MESH_Num_Edges(mesh);
    nf = MESH_Num_Faces(mesh);
    nr = MESH_Num_Regions(mesh);

    ipos = 0;

    /* build vertices and neighbors list, similar as in partition with metis
       Assign processor 0 the whole mesh, assign other processors a NULL mesh */
  
    if (nr == 0) {
      if (nf == 0) {
	fprintf(stderr,"Cannot partition wire meshes with Zoltan\n");
	exit(-1);
      
      }

      graph.vertexGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * nf);
      graph.nborIndex = (int *)malloc(sizeof(int) * (nf + 1));
      graph.nborGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * 2*ne);
      graph.nborProc = (int *)malloc(sizeof(int) * 2*ne);
      
      graph.nborIndex[0] = 0;
      
      ncells = nf;
    
      /* Surface mesh */
      idx = 0; i = 0;
      while ((mf = MESH_Next_Face(mesh,&idx))) {
	graph.vertexGID[i] = MF_ID(mf);
	fedges = MF_Edges(mf,1,0);
	nfe = List_Num_Entries(fedges);
	
	idx2 = 0;
	while ((fedge = List_Next_Entry(fedges,&idx2))) {
	  
	  efaces = ME_Faces(fedge);
	  nef = List_Num_Entries(efaces);
	  
	  if (nef > 2) {
	    fprintf(stderr,"Non-manifold surface mesh. Exit!\n");
	    exit(-1);
	  }
	  else if (nef == 1) {
	    continue;          /* boundary edge; nothing to do */
	  }
	  else {
	    oppf = List_Entry(efaces,0);
	    if (oppf == mf)
	      oppf = List_Entry(efaces,1);
	    
	    graph.nborGID[ipos] = MF_ID(oppf);
	    /* initially set all vertices on processor 0 */
	    graph.nborProc[ipos] = 0;
	    ipos++;
	  }
	  
	  List_Delete(efaces);
	  
	}
	
	List_Delete(fedges);
	i++;
	graph.nborIndex[i] = ipos;
      }
      graph.numMyVertices = i;
      graph.numAllNbors = ipos;
    }
    else {
      graph.vertexGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * nr);
      graph.nborIndex = (int *)malloc(sizeof(int) * (nr + 1));
      graph.nborGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * 2*nf);
      graph.nborProc = (int *)malloc(sizeof(int) * 2*nf);
      
      graph.nborIndex[0] = 0;
      
      ncells = nr;
      
      /* Volume mesh */
      
      idx = 0; i = 0;
      while ((mr = MESH_Next_Region(mesh,&idx))) {
	graph.vertexGID[i] = MR_ID(mr);
	rfaces = MR_Faces(mr);
	nrf = List_Num_Entries(rfaces);
      
	idx2 = 0;
	while ((rface = List_Next_Entry(rfaces,&idx2))) {
	  
	  fregions = MF_Regions(rface);
	  nfr = List_Num_Entries(fregions);
	  
	  if (nfr == 1) {
	    continue;          /* boundary face; nothing to do */
	  }
	  else {
	    oppr = List_Entry(fregions,0);
	    if (oppr == mr)
	      oppr = List_Entry(fregions,1);
	    
	    graph.nborGID[ipos] = MR_ID(oppr);
	    /* initially set all vertices on processor 0 */
	    graph.nborProc[ipos] = 0;
	    ipos++;
	  }
	  
	  List_Delete(fregions);
	  
	}
	
	List_Delete(rfaces);
	
	i++;
	graph.nborIndex[i] = ipos;
      }
      graph.numMyVertices = i;
      graph.numAllNbors = ipos;
    }
    /* other processors have a NULL graph */
  }
  else
    {
      graph.numMyVertices = 0;
      graph.numAllNbors = 0;
    }
  

  /* General parameters for Zoltan */

  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
  Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");

  /* Graph parameters */

  Zoltan_Set_Param(zz, "CHECK_GRAPH", "2");
  Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", ".35");  /* 0-remove all, 1-remove none */

  /* Query functions - defined in simpleQueries.h */

  Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, &graph);
  Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, &graph);
  Zoltan_Set_Num_Edges_Multi_Fn(zz, get_num_edges_list, &graph);
  Zoltan_Set_Edge_List_Multi_Fn(zz, get_edge_list, &graph);

  /* Partition the graph */
  /******************************************************************                                                                             
   ** Zoltan can now partition the graph.                                                                                                   
   ** We assume the number of partitions is                                                                                
   ** equal to the number of processes.  Process rank 0 will own                                                                                   
   ** partition 0, process rank 1 will own partition 1, and so on.                                                                                 
   ******************************************************************/
  rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
			   &changes,        /* 1 if partitioning was changed, 0 otherwise */
			   &numGidEntries,  /* Number of integers used for a global ID */
			   &numLidEntries,  /* Number of integers used for a local ID */
			   &numImport,      /* Number of vertices to be sent to me */
			   &importGlobalGids,  /* Global IDs of vertices to be sent to me */
			   &importLocalGids,   /* Local IDs of vertices to be sent to me */
			   &importProcs,    /* Process rank for source of each incoming vertex */
			   &importToPart,   /* New partition for each incoming vertex */
			   &numExport,      /* Number of vertices I must send to other processes*/
			   &exportGlobalGids,  /* Global IDs of the vertices I must send */
			   &exportLocalGids,   /* Local IDs of the vertices I must send */
			   &exportProcs,    /* Process to which I send each of the vertices */
			   &exportToPart);  /* Partition to which each vertex will belong */

  if (rc != ZOLTAN_OK){
    printf("sorry...\n");
    MPI_Finalize();
    Zoltan_Destroy(&zz);
    exit(0);
  }

  if(rank == 0) {
    *part = (int *) malloc(ncells*sizeof(int));
    for ( i = 0; i < ncells; i++ ) 
      (*part)[i] = 0;
    for ( i = 0; i < numExport; i++ ) {
      (*part)[exportGlobalGids[i]-1] = exportToPart[i];
    }
    free(graph.vertexGID);
    free(graph.nborIndex);
    free(graph.nborGID);
    free(graph.nborProc);
  }
  else { 
    *part = NULL;
  }


  Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids,                                                                                       		      &exportProcs, &exportToPart);                                                                                               
  Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids,                                                                                                              &importProcs, &importToPart);                                                                                               Zoltan_Destroy(&zz);                                                                                                                           
  return 1;
}

  /* Application defined query functions */

  static int get_number_of_vertices(void *data, int *ierr)
  {
    GRAPH_DATA *graph = (GRAPH_DATA *)data;
    *ierr = ZOLTAN_OK;
    return graph->numMyVertices;
  }

  static void get_vertex_list(void *data, int sizeGID, int sizeLID,
			      ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
			      int wgt_dim, float *obj_wgts, int *ierr)
  {
    int i;

    GRAPH_DATA *graph = (GRAPH_DATA *)data;
    *ierr = ZOLTAN_OK;

    /* In this example, return the IDs of our vertices, but no weights.                                                                             
     * Zoltan will assume equally weighted vertices.                                                                                                
     */

    for (i=0; i<graph->numMyVertices; i++){
      globalID[i] = graph->vertexGID[i];
      localID[i] = i;
    }
  }

  static void get_num_edges_list(void *data, int sizeGID, int sizeLID,
				 int num_obj,
				 ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
				 int *numEdges, int *ierr)
  {
    int i, idx;

    GRAPH_DATA *graph = (GRAPH_DATA *)data;

    if ( (sizeGID != 1) || (sizeLID != 1) || (num_obj != graph->numMyVertices)){
      *ierr = ZOLTAN_FATAL;
      return;
    }

    for (i=0;  i < num_obj ; i++){
      idx = localID[i];
      numEdges[i] = graph->nborIndex[idx+1] - graph->nborIndex[idx];
    }

    *ierr = ZOLTAN_OK;
    return;
  }

  static void get_edge_list(void *data, int sizeGID, int sizeLID,
			    int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
			    int *num_edges,
			    ZOLTAN_ID_PTR nborGID, int *nborProc,
			    int wgt_dim, float *ewgts, int *ierr)
  {
    int i, j, from, to;
    int *nextProc;
    ZOLTAN_ID_TYPE *nextNbor;

    GRAPH_DATA *graph = (GRAPH_DATA *)data;
    *ierr = ZOLTAN_OK;

    if ( (sizeGID != 1) || (sizeLID != 1) ||
	 (num_obj != graph->numMyVertices)||
	 (wgt_dim != 0)){
      *ierr = ZOLTAN_FATAL;
      return;
    }

    nextNbor = nborGID;
    nextProc = nborProc;

    for (i=0; i < num_obj; i++){
      /*                                                                                                                                            
       * In this example, we are not setting edge weights.  Zoltan will                                                                             
       * set each edge to weight 1.0.                                                                                                               
       */

      to = graph->nborIndex[localID[i]+1];
      from = graph->nborIndex[localID[i]];
      if ((to - from) != num_edges[i]){
	*ierr = ZOLTAN_FATAL;
	return;
      }

      for (j=from; j < to; j++){

	*nextNbor++ = graph->nborGID[j];
	*nextProc++ = graph->nborProc[j];
      }
    }
    return;

  }


#ifdef __cplusplus
  }
#endif

