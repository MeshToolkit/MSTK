#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif


  /* Partition a given mesh and distribute it to 'num' processors 

     Authors: Rao Garimella
              Duo Wang
  */


  int MSTK_Mesh_Distribute(Mesh_ptr parentmesh, Mesh_ptr *mysubmesh, int *dim, 
			   int ring, int with_attr, int method, 
			   int del_inmesh, MSTK_Comm comm) {
    int i, a, m, n, recv_dim;
    int *send_dim, *part=NULL;
    int rank, numprocs, *toranks;
    int DebugWait=0;
    MAttrib_ptr attrib;
    MSet_ptr mset;

    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&numprocs);
    recv_dim = rank+5;          /* ??? */

    while (DebugWait) ;

#ifdef DEBUG
    double elapsed_time;
    double t0 = MPI_Wtime();
#endif

    send_dim = (int *) malloc(numprocs*sizeof(int));
    for (i = 0; i < numprocs; i++) send_dim[i] = *dim;

    MPI_Scatter(send_dim, 1, MPI_INT, &recv_dim, 1, MPI_INT, 0, comm);
    free(send_dim);
    if (rank != 0)
      *dim = recv_dim;


    /* PARTITIONING OF THE MESH GRAPH AND GETTING THE PARTITION
       NUMBERS FOR EACH ELEMENT */

    MESH_Get_Partitioning(parentmesh, method, &part, comm);


    if (rank == 0) {
      toranks = (int *) malloc(numprocs*sizeof(int));
      for (i = 0; i < numprocs; i++) toranks[i] = i;

      *mysubmesh = MESH_New(MESH_RepType(parentmesh));
      MESH_Partition_and_Send(parentmesh, numprocs, part, toranks, ring, 
                              with_attr, del_inmesh, comm, mysubmesh);

      free(toranks);

      fprintf(stderr,"Finished partioning and sending on rank 0\n");
    }

    if (part) free(part);

    
    /* RECEIVING THE INFORMATION ON OTHER PROCESSORS - THIS IS
       STRAIGHTFORWARD AND USES BLOCKING MPI RECEIVES BECAUSE EACH
       STEP NEEDS TO BE COMPLETED BEFORE THE NEXT */


    if (rank > 0) { /* Receive the mesh from processor 0 */

      int fromrank = 0;
      int nv, ne, nf, nr;
      RepType rtype;

      *mysubmesh = MESH_New(UNKNOWN_REP);
      MESH_RecvMesh(*mysubmesh, fromrank, with_attr, comm);

    }



    /* FINISH UP */
    /* Build the sorted lists of ghost and overlap entities */

    MESH_Build_GhostLists(*mysubmesh,*dim);


    /* Some additional parallel information update to indicate which
       processors communicate with which others */
    
    MESH_Set_Prtn(*mysubmesh,rank,numprocs);

    MESH_Update_ParallelAdj(*mysubmesh,comm);


    MESH_Disable_GlobalIDSearch(*mysubmesh);


#ifdef DEBUG
    elapsed_time = MPI_Wtime() - t0;
    fprintf(stderr,"Elapsed time after mesh distribution on processor %-d is %lf s\n",rank,elapsed_time);
#endif


    /* Put a barrier so that distribution of meshes takes place one at a time 
       in a simulation that may have multiple mesh objects on each processor */

    MPI_Barrier(comm);

    return 1;
  }



#ifdef __cplusplus
}
#endif

