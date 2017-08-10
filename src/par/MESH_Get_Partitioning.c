#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"
#include "MSTK_private.h"


#ifdef __cplusplus
extern "C" {
#endif

  /* 
     This function derives a partitioning of the mesh on processor 0

     Author(s): Duo Wang, Rao Garimella

  */


int MESH_Get_Partitioning(Mesh_ptr mesh, PartitionMethod method,
                          int **part, MSTK_Comm comm) {
  int i, ok = 1, nf, nr, ncells=0;
  int noptions;
  char **options = (char **) malloc(10*sizeof(char *));
  for (i = 0; i < 10; i++) options[i] = (char *) malloc(64*sizeof(char));

  int rank, num;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&num);

  *part = NULL;  

  /* basic mesh information */
  
  if ( rank == 0 ) {
    nf = MESH_Num_Faces(mesh);
    nr = MESH_Num_Regions(mesh);
    ncells = (nr) ? nr : nf; 

#ifdef DEBUG
    if (ncells == 0) {
      MSTK_Report("MESH_Get_Partition",
                  "Mesh has no cells to partition",MSTK_ERROR);
      exit(-1);
    }
#endif
  }
  
  if (num == 1) {
    *part = (int *) calloc(ncells,sizeof(int)); /* will initialize to 0 */
    for (i = 0; i < 10; ++i) free(options[i]);
    free(options);
    return 1;
  }


  switch (method) {
  case METIS: {

    /* Call the partitioner only on processor zero */

    if (rank == 0) {
#ifdef _MSTK_HAVE_METIS
      ok = MESH_PartitionWithMetis(mesh, num, part);
#else
      MSTK_Report("MESH_Partition","Metis not enabled",MSTK_FATAL);
#endif
    }

    break;
  }
  case ZOLTAN_GRAPH: {
    /* This invokes the graph partitioner in Zoltan */
    
    /* Even though we assign all entities to processor zero and
       ask Zoltan to partition the mesh, we have to invoke the 
       Zoltan partitioner on all processors */
    
#ifdef _MSTK_HAVE_ZOLTAN
    noptions = 1;
    strcpy(options[0],"LB_PARTITION=GRAPH");
    ok = MESH_PartitionWithZoltan(mesh, num, part, noptions, options, comm);
#else
    MSTK_Report("MESH_Partition","Zoltan not enabled",MSTK_FATAL);
#endif
    
    break;
  }
  case ZOLTAN_RCB: {
    /* This invokes the Recursive Coordinate Bisection partitioner in Zoltan */

    /* Even though we assign all entities to processor zero and
       ask Zoltan to partition the mesh, we have to invoke the 
       Zoltan partitioner on all processors */

#ifdef _MSTK_HAVE_ZOLTAN
    noptions = 1;
    strcpy(options[0],"LB_PARTITION=RCB");
    ok = MESH_PartitionWithZoltan(mesh, num, part, noptions, options, comm);
#else
    MSTK_Report("MESH_Partition","Zoltan not enabled",MSTK_FATAL);
#endif

    break;
  }
  case METIS_COLUMNAR: {
      
    /* Call the partitioner only on processor zero */
    /* For volume meshes, this also assumes that the mesh has a columnar structure
     and does a redistribution to ensure that all regions in the column are in
     the same partition */
    if (rank == 0) {
#ifdef _MSTK_HAVE_METIS
      ok = MESH_PartitionWithMetis(mesh, num, part);
      if (ok)
        ok = FixColumnPartitions(mesh, *part, comm);
#else
      MSTK_Report("MESH_Partition","Metis not enabled",MSTK_FATAL);
#endif
    }
      
    break;
  }
  case ZOLTAN_GRAPH_COLUMNAR: {
    /* This invokes the graph partitioner in Zoltan */
    /* For volume meshes, this also assumes that the mesh has a columnar structure
     and does a redistribution to ensure that all regions in the column are in
     the same partition */
    /* Even though we assign all entities to processor zero and
     ask Zoltan to partition the mesh, we have to invoke the
     Zoltan partitioner on all processors */
      
#ifdef _MSTK_HAVE_ZOLTAN
    noptions = 1;
    strcpy(options[0],"LB_PARTITION=GRAPH");
    ok = MESH_PartitionWithZoltan(mesh, num, part, noptions, options, comm);
    if (ok && (rank == 0))
      ok = FixColumnPartitions(mesh, *part, comm);
#else
    MSTK_Report("MESH_Partition","Zoltan not enabled",MSTK_FATAL);
#endif
      
    break;
  }
  case ZOLTAN_RCB_COLUMNAR: {
    /* This invokes the Recursive Coordinate Bisection partitioner in Zoltan */
    /* For volume meshes, this also assumes that the mesh has a columnar structure
       and does a redistribution to ensure that all regions in the column are in 
       the same partition */
    /* Even though we assign all entities to processor zero and
     ask Zoltan to partition the mesh, we have to invoke the
     Zoltan partitioner on all processors */
      
#ifdef _MSTK_HAVE_ZOLTAN
    noptions = 1;
    strcpy(options[0],"LB_PARTITION=RCB");
    ok = MESH_PartitionWithZoltan(mesh, num, part, noptions, options, comm);
    if (ok && (rank == 0))
      ok = FixColumnPartitions(mesh, *part, comm);
#else
    MSTK_Report("MESH_Partition","Zoltan not enabled",MSTK_FATAL);
#endif
      
    break;
  }
  default:
    MSTK_Report("MESH_Get_Partition","Unknown partitioning method",MSTK_FATAL);
  }


  for (i = 0; i < 10; i++)
    free(options[i]);
  free(options);

  return ok;
}
    

  
#ifdef __cplusplus
}
#endif

