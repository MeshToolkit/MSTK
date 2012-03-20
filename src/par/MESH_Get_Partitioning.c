#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif

  /* 
     This function derives a partitioning of the mesh on processor 0

     Author(s): Duo Wang, Rao Garimella

  */


int MESH_Get_Partitioning(Mesh_ptr mesh, int num, int method, int rank, MPI_Comm comm, int **part) {
  int ok = 1, nf, nr, ncells;
  
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

    *part = (int *) MSTK_malloc(ncells*sizeof(int));
  }
  else
    *part = NULL;

  switch (method) {
  case 0:

    /* Call the partitioner only on processor zero */

    if (rank == 0) {
#ifdef _MSTK_WITH_METIS
      ok = MESH_PartitionWithMetis(mesh, num, part);
#else
      MSTK_Report("MESH_Partition","Metis not enabled",MSTK_FATAL);
#endif
    }

    break;
  case 1:
    /* Even though we assign all entities to processor zero and
       ask Zoltan to partition the mesh, we have to invoke the 
       Zoltan partitioner on all processors */

#ifdef _MSTK_WITH_ZOLTAN
    ok = MESH_PartitionWithZoltan(mesh, num, part,rank,comm);
#else
    MSTK_Report("MESH_Partition","Zoltan not enabled",MSTK_FATAL);
#endif

    break;
  default:
    MSTK_Report("MESH_Get_Partition","Unknown partitioning method",MSTK_FATAL);
  }

  return ok;
}
    

  
#ifdef __cplusplus
}
#endif

