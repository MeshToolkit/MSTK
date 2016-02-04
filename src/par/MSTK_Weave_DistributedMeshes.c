#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif



  /* Weave a set of distributed mesh partitions together to build the
     parallel connections and ghost info.

     input_type indicates what info is already present on the mesh
     
     0 -- we are given NO information about how these meshes are connected
          other than the knowledge that they come from the partitioning of
          a single mesh

     1 -- we are given partitioned meshes with a unique global ID on 
          each mesh vertex

     2 -- we are given parallel neighbor information, but no global ID on 
          each mesh vertex


  */
     


  int MSTK_Weave_DistributedMeshes(Mesh_ptr mesh, int topodim,
                                   int num_ghost_layers, int input_type,
                                   MSTK_Comm comm) {

    int have_GIDs = 0;
    int rank, num;

    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&num);

    if (num_ghost_layers > 1)
      MSTK_Report("MSTK_Weave_DistributedMeshes", "Only 1 ghost layer supported currently", MSTK_FATAL);

    if (input_type > 2) 
      MSTK_Report("MSTK_Weave_DistributedMeshes","Unrecognized input type for meshes", MSTK_WARN);

    // This partition does not have a mesh or has an empty mesh which is ok

    if (mesh == NULL)
      MSTK_Report("MSTK_Weave_DistributedMeshes","MESH is null on this processor",MSTK_FATAL);


    MESH_Set_Prtn(mesh, rank, num);
    
    if (input_type == 0)
      have_GIDs = 0;
    else if (input_type == 1)
      have_GIDs = 1;

    if (input_type == 0 || input_type == 1) {
      /* MESH_MatchEnts_ParBdry(mesh, have_GIDs, rank, num, comm); */
      MESH_AssignGlobalIDs(mesh, topodim, have_GIDs, comm);
      MESH_BuildConnection(mesh, topodim, comm);
    }
    else if (input_type == 2) 
      MESH_AssignGlobalIDs_p2p(mesh, topodim, comm);

    MESH_LabelPType(mesh, topodim, comm);

    MESH_Parallel_AddGhost(mesh, topodim, comm);

    MESH_Build_GhostLists(mesh, topodim);

    MESH_Update_ParallelAdj(mesh, comm);
    return 1;
  }


#ifdef __cplusplus
}
#endif

