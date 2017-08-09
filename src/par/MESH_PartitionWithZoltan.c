#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "zoltan.h"

#include "MSTK.h"
#include "MSTK_private.h"

/* Get partitioning of mesh using Zoltan. Doesn't actually do anything 
   the mesh */

/* When options[0] is LB_PARTITION=RCB in a N-dimensional mesh, then
   it is assumed that we want to ignore the last or N'th dimension of
   the 3D mesh and partition only in first N-1 dimensions */


#ifdef __cplusplus
extern "C" {
#endif

typedef struct{
  int numMyNodes;          /* total nodes in in my partition of the graph */
  int numAllNbors;            /* total number of neighbors of my graph nodes */  
  ZOLTAN_ID_TYPE *nodeGID;  /* global ID of each of my graph nodes */
  double *nodeCoords;      /* coordinates of my graph nodes = 3*numMyNodes */
  int *nborIndex;    /* nborIndex[i] is location of start of neighbors for graph node i */
  ZOLTAN_ID_TYPE *nborGID;   /* nborGIDs[nborIndex[i]] is first neighbor of node i */
  int *nborProc;             /* process owning each nbor in nborGID */
} GRAPH_DATA;
  
int NDIM_4_ZOLTAN = 3;

  /* Application defined query callback functions for ZOLTAN */
  
  static int get_number_of_nodes(void *data, int *ierr);
  static void get_node_list(void *data, int sizeGID, int sizeLID,
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

  static int get_num_dimensions_reduced(void *data, int *ierr);
  static void get_element_centers_reduced(void *data, int sizeGID, int sizeLID, int num_obj,
                                          ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                                          int num_dim, double *coord_vec, int *ierr);
  
  int MESH_PartitionWithZoltan(Mesh_ptr mesh, int nparts, int **part, int noptions, 
                               char **options, MSTK_Comm comm) { 

  MEdge_ptr fedge;
  MFace_ptr mf, oppf, rface;
  MRegion_ptr mr, oppr;
  List_ptr fedges, efaces, rfaces, fregions;
  int  i, j, k, id;
  int  nv, ne, nf, nr=0, nfe, nef, nfr, nrf, idx, idx2;
  int  numflag, nedgecut, ipos;
  int  wtflag;

  int rc;
  float ver;
  struct Zoltan_Struct *zz;
  GRAPH_DATA graph;
  int changes, numGidEntries, numLidEntries, numImport, numExport;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
  int *importProcs, *importToPart, *exportProcs, *exportToPart;

  int rank;
  MPI_Comm_rank(comm,&rank);
 
  rc = Zoltan_Initialize(0, NULL, &ver);

  if (rc != ZOLTAN_OK){
    MSTK_Report("MESH_PartitionWithZoltan","Could not initialize Zoltan",MSTK_FATAL);
    MPI_Finalize();
    exit(0);
  }

  /******************************************************************
  ** Create a Zoltan library structure for this instance of partition 
  ********************************************************************/
  zz = Zoltan_Create(comm);

  /*****************************************************************
   ** Figure out partitioning method
   *****************************************************************/
  
  char partition_method_str[32];
  strcpy(partition_method_str,"RCB");  /* Default - Recursive Coordinate Bisection */
  if (noptions) {
    for (i = 0; i < noptions; i++) {
      if (strncmp(options[i],"LB_PARTITION",12) == 0) {
        char *result = NULL, instring[256];
        strcpy(instring,options[i]);
        result = strtok(instring,"=");
        result = strtok(NULL," ");
        strcpy(partition_method_str,result);
      }
    }
  }
  
  if (rank == 0) {
    char mesg[256];
    sprintf(mesg,"Using partitioning method %s for ZOLTAN\n",partition_method_str);
    MSTK_Report("MESH_PartitionWithZoltan",mesg,MSTK_MESG);
  }

  /* General parameters for Zoltan */
  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param(zz, "LB_METHOD", partition_method_str);
  Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");


  graph.numMyNodes = 0;
  graph.numAllNbors = 0;
  graph.nodeGID = NULL;
  graph.nodeCoords = NULL;
  graph.nborIndex = NULL;
  graph.nborGID = NULL;
  graph.nborProc = NULL;

  if (strcmp(partition_method_str,"RCB") == 0) {
    if (rank == 0) {
      nr = MESH_Num_Regions(mesh);
      nf = MESH_Num_Faces(mesh);

      if (!nf && !nr)
        MSTK_Report("MESH_PartitionWithZoltan","Cannot partition wire meshes",
                    MSTK_FATAL);

      if (nr == 0) { /* Surface or planar mesh */

        int ndim = 2;       /* assume mesh is planar */
        idx = 0;
        MVertex_ptr mv;
        while ((mv = MESH_Next_Vertex(mesh,&idx))) {
          double vxyz[3];
          MV_Coords(mv,vxyz);
          if (vxyz[2] != 0.0) {
            ndim = 3;  /* non-planar or planar with non-zero z */
            break;
          }
        }
        NDIM_4_ZOLTAN = ndim-1;  /* ignore last dimension to avoid partitioning in that dimension */

        graph.numMyNodes = nf;

        graph.nodeGID = (ZOLTAN_ID_TYPE *) malloc(sizeof(ZOLTAN_ID_TYPE) * nf);
        graph.nodeCoords = (double *) malloc(sizeof(double) * NDIM_4_ZOLTAN * nf);

        idx = 0;
        while ((mf = MESH_Next_Face(mesh,&idx))) {
          double fxyz[MAXPV2][3], cen[3];
          int nfv;

          MF_Coords(mf,&nfv,fxyz);
          cen[0] = cen[1] = cen[2] = 0.0;
          for (j = 0; j < nfv; j++)
            for (k = 0; k < NDIM_4_ZOLTAN; k++) 
              cen[k] += fxyz[j][k];              
          for (k = 0; k < NDIM_4_ZOLTAN; k++) cen[k] /= nfv;

          id = MF_ID(mf);
          graph.nodeGID[id-1] = id;
          memcpy(&(graph.nodeCoords[NDIM_4_ZOLTAN*(id-1)]),cen,NDIM_4_ZOLTAN*sizeof(double));
        }

      }
      else { /* Volume mesh */

        int ndim = 3;
        NDIM_4_ZOLTAN = ndim-1;  /* ignore last dimension  to avoid partitioning in that dimension */
        graph.numMyNodes = nr;

        graph.nodeGID = (ZOLTAN_ID_TYPE *) malloc(sizeof(ZOLTAN_ID_TYPE) * nr);
        graph.nodeCoords = (double *) malloc(sizeof(double) * NDIM_4_ZOLTAN * nr);

        idx = 0;
        while ((mr = MESH_Next_Region(mesh,&idx))) {
          double rxyz[MAXPV3][3], cen[3];
          int nrv;
          
          MR_Coords(mr,&nrv,rxyz);
          cen[0] = cen[1] = cen[2] = 0.0;
          for (j = 0; j < nrv; j++)
            for (k = 0; k < NDIM_4_ZOLTAN; k++)
              cen[k] += rxyz[j][k];
          for (k = 0; k < NDIM_4_ZOLTAN; k++) cen[k] /= nrv;
          for (k = 0; k < NDIM_4_ZOLTAN; k++) 
            if (fabs(cen[k]) < 1.0e-10) cen[k] = 0.0; 

          id = MR_ID(mr);
          graph.nodeGID[id-1] = id;
          memcpy(&(graph.nodeCoords[NDIM_4_ZOLTAN*(id-1)]),cen,NDIM_4_ZOLTAN*sizeof(double));
        }

      }
    }

    MPI_Bcast(&NDIM_4_ZOLTAN,1,MPI_INT,0,comm);

    /* Set some default values */
    Zoltan_Set_Param(zz, "RCB_RECTILINEAR_BLOCKS","1");
    //    Zoltan_Set_Param(zz, "AVERAGE_CUTS", "1");

    if (noptions > 1) {
      for (i = 1; i < noptions; i++) {
        char *paramstr = NULL, *valuestr = NULL, instring[256];
        strcpy(instring,options[i]);
        paramstr = strtok(instring,"=");
        valuestr = strtok(NULL," ");
        Zoltan_Set_Param(zz,paramstr,valuestr);
      }
    }

    /* Query functions - defined in simpleQueries.h */

    Zoltan_Set_Num_Obj_Fn(zz, get_number_of_nodes, &graph);
    Zoltan_Set_Obj_List_Fn(zz, get_node_list, &graph);
    Zoltan_Set_Num_Geom_Fn(zz, get_num_dimensions_reduced, &graph);    /* reduced dimensions */
    Zoltan_Set_Geom_Multi_Fn(zz, get_element_centers_reduced, &graph); /* reduced dimension centers */

  }
  else if (strcmp(partition_method_str,"GRAPH") == 0) {

    if(rank == 0) {
      nv = MESH_Num_Vertices(mesh);
      ne = MESH_Num_Edges(mesh);
      nf = MESH_Num_Faces(mesh);
      nr = MESH_Num_Regions(mesh);
      
      ipos = 0;
      
      /* build nodes and neighbors list, similar as in partition with metis
         Assign processor 0 the whole mesh, assign other processors a NULL mesh */
  
      if (nr == 0) {
        if (nf == 0) {
          MSTK_Report("MESH_PartitionWithZoltan",
                      "Cannot partition wire meshes with Zoltan",MSTK_FATAL);
          exit(-1);
      
        }

        graph.nodeGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * nf);
        graph.nborIndex = (int *)malloc(sizeof(int) * (nf + 1));
        graph.nborGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * 2*ne);
        graph.nborProc = (int *)malloc(sizeof(int) * 2*ne);
      
        graph.nborIndex[0] = 0;
      
        /* Surface mesh */
        idx = 0; i = 0;
        while ((mf = MESH_Next_Face(mesh,&idx))) {
          graph.nodeGID[i] = MF_ID(mf);
          fedges = MF_Edges(mf,1,0);
          nfe = List_Num_Entries(fedges);
	
          idx2 = 0;
          while ((fedge = List_Next_Entry(fedges,&idx2))) {
	  
            efaces = ME_Faces(fedge);
            nef = List_Num_Entries(efaces);
	  
            if (nef == 1) {
              continue;          /* boundary edge; nothing to do */
            } else {
              int j;
              for (j = 0; j < nef; j++) {
                oppf = List_Entry(efaces,j);
                if (oppf == mf) {
                  graph.nborGID[ipos] = MF_ID(oppf);
                  /* initially set all nodes on processor 0 */
                  graph.nborProc[ipos] = 0;
                  ipos++;
                }
              }
            }
	  
            List_Delete(efaces);
	  
          }
	
          List_Delete(fedges);
          i++;
          graph.nborIndex[i] = ipos;
        }
        graph.numMyNodes = i;
        graph.numAllNbors = ipos;
      }
      else {
        graph.nodeGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * nr);
        graph.nborIndex = (int *)malloc(sizeof(int) * (nr + 1));
        graph.nborGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * 2*nf);
        graph.nborProc = (int *)malloc(sizeof(int) * 2*nf);
      
        graph.nborIndex[0] = 0;
      
        /* Volume mesh */
      
        idx = 0; i = 0;
        while ((mr = MESH_Next_Region(mesh,&idx))) {
          graph.nodeGID[i] = MR_ID(mr);
          rfaces = MR_Faces(mr);
          nrf = List_Num_Entries(rfaces);
      
          idx2 = 0;
          while ((rface = List_Next_Entry(rfaces,&idx2))) {
	  
            fregions = MF_Regions(rface);
            nfr = List_Num_Entries(fregions);
	  
            if (nfr > 1) {
              oppr = List_Entry(fregions,0);
              if (oppr == mr)
                oppr = List_Entry(fregions,1);
	    
              graph.nborGID[ipos] = MR_ID(oppr);
              /* initially set all nodes on processor 0 */
              graph.nborProc[ipos] = 0;
              ipos++;
            }
	  
            List_Delete(fregions);
	  
          }
	
          List_Delete(rfaces);
	
          i++;
          graph.nborIndex[i] = ipos;
        }
        graph.numMyNodes = i;
        graph.numAllNbors = ipos;
      }
    }

    /* Graph parameters */

    /* Zoltan_Set_Param(zz, "CHECK_GRAPH", "2"); */
    Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", ".35");  /* 0-remove all, 1-remove none */

    /* Query functions - defined in simpleQueries.h */

    Zoltan_Set_Num_Obj_Fn(zz, get_number_of_nodes, &graph);
    Zoltan_Set_Obj_List_Fn(zz, get_node_list, &graph);
    Zoltan_Set_Num_Edges_Multi_Fn(zz, get_num_edges_list, &graph);
    Zoltan_Set_Edge_List_Multi_Fn(zz, get_edge_list, &graph);    
  }

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
			   &numImport,      /* Number of nodes to be sent to me */
			   &importGlobalGids,  /* Global IDs of nodes to be sent to me */
			   &importLocalGids,   /* Local IDs of nodes to be sent to me */
			   &importProcs,    /* Process rank for source of each incoming node */
			   &importToPart,   /* New partition for each incoming node */
			   &numExport,      /* Number of nodes I must send to other processes*/
			   &exportGlobalGids,  /* Global IDs of the nodes I must send */
			   &exportLocalGids,   /* Local IDs of the nodes I must send */
			   &exportProcs,    /* Process to which I send each of the nodes */
			   &exportToPart);  /* Partition to which each node will belong */

  if (rc != ZOLTAN_OK){
    if (rank == 0)
      MSTK_Report("MESH_PartitionWithZoltan","Could not partition mesh with ZOLTAN",
                  MSTK_ERROR);
    Zoltan_Destroy(&zz);
    MPI_Finalize();
    return 0;
  }

  if(rank == 0) {
    *part = (int *) calloc(graph.numMyNodes,sizeof(int));
    for ( i = 0; i < numExport; i++ ) {
      (*part)[exportGlobalGids[i]-1] = exportToPart[i];
    }
    if (graph.nodeGID) free(graph.nodeGID);
    if (graph.nodeCoords) free(graph.nodeCoords);
    if (graph.nborIndex) free(graph.nborIndex);
    if (graph.nborGID) free(graph.nborGID);
    if (graph.nborProc) free(graph.nborProc);
  }
  else { 
    *part = NULL;
  }


  Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);
  Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, &importProcs, &importToPart);
  Zoltan_Destroy(&zz);                

  return 1;
}

  /* Application defined query functions */

  static int get_number_of_nodes(void *data, int *ierr)
  {
    GRAPH_DATA *graph = (GRAPH_DATA *)data;
    *ierr = ZOLTAN_OK;
    return graph->numMyNodes;
  }

  static void get_node_list(void *data, int sizeGID, int sizeLID,
			      ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
			      int wgt_dim, float *obj_wgts, int *ierr)
  {
    int i;

    GRAPH_DATA *graph = (GRAPH_DATA *)data;
    *ierr = ZOLTAN_OK;

    /* In this example, return the IDs of our nodes, but no weights.                                                                             
     * Zoltan will assume equally weighted nodes.                                                                                                
     */

    for (i=0; i<graph->numMyNodes; i++){
      globalID[i] = graph->nodeGID[i];
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

    if ( (sizeGID != 1) || (sizeLID != 1) || (num_obj != graph->numMyNodes)){
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
	 (num_obj != graph->numMyNodes)||
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


  static int get_num_dimensions_reduced(void *data, int *ierr) {
    return NDIM_4_ZOLTAN;
  }

  static void get_element_centers_reduced(void *data, int sizeGID, int sizeLID, int num_obj,
                                  ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                                  int num_dim, double *coord_vec, int *ierr) {
    
    GRAPH_DATA *graph = (GRAPH_DATA *) data;
    *ierr = ZOLTAN_OK;

    if (graph->nodeCoords == NULL) return;
    memcpy(coord_vec,graph->nodeCoords,num_dim*num_obj*sizeof(double));
  }

#ifdef __cplusplus
  }
#endif

