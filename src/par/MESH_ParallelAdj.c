#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* 
     this function updates (or intializes) inter-processor
     relationships of the mesh


     for applications that do not involve topology change,
     this routine can be called only once during initialization

     mesh_par_adj_flags is an integer array of size num
     mesh_par_adj_flags[i] on a processor j indicates the relationship of the mesh 
     on processor j with the mesh on processor i.

     1 and 2 bit(from right) 
     indicate relation on vertex
     0(00) no relation,
     1(01) has ghost entities related to processor j,
     2(10) has overlap entities related to processor j,
     3(11) both
     3 and 4 bit indicate relation on edge
     5 and 6 bit indicate relation on face
     7 and 8 bit indicate relation on region

     for example: mesh_par_adj_flags[3] = 51 (00110011) means mesh on this
     processor has ghost vertices and ghost faces whose master
     processor is 3 and has overlap vertices and overlap faces which
     are ghost entities on processor 3.



     mesh_par_adj_info stores the number of overlap entities for allocating
     recv buffer 

     mesh_par_adj_info[0]: number of processors that have overlap
     with mesh on a particular processor

     mesh_par_adj_info[1]-mesh_par_adj_info[num_recv_rank+1]: processor ids that
     have overlap with mesh on particular processor (What is num_recv_rank?)

     Then num_recv_rank*4 array stores the number of overlap
     entities. (DOES THAT MEAN THAT THE SUBSEQUENT ENTRIES INDICATE
     THE NUMBER OF OVERLAP ENTITIES OF EACH TYPE FOR EACH OF THE
     RELEVANT PROCESSORS?)

     Author(s): Duo Wang, Rao Garimella

  */


int MESH_Update_ParallelAdj(Mesh_ptr mesh, int rank, int num,  MPI_Comm comm) {
  int i,idx,nv,ne,nf,nr,local_ov_num[4];
  int ebit, ov_index, ov_index2;
  MVertex_ptr mv;
  MEdge_ptr me;
  MFace_ptr mf;
  MRegion_ptr mr;
  int *mesh_par_adj_info, *mesh_par_adj_flags;


  int *global_ov_num = (int *) MSTK_malloc(4*num*sizeof(int));
  int *global_ranks = (int *) MSTK_malloc(num*num*sizeof(int));

  
  int *local_ranks = (int *)MSTK_malloc(num*sizeof(int));
  for(i = 0; i < num; i++)
    local_ranks[i] = 0;

  /* set ghost adjacencies */

  idx = 0;
  while(mv = MESH_Next_GhostVertex(mesh,&idx))
    MESH_Flag_Has_GhostEnts_From_Proc(mesh,MVERTEX,MV_MasterParID(mv));
  idx = 0;
  while(me = MESH_Next_GhostEdge(mesh,&idx))
    MESH_Flag_Has_GhostEnts_From_Proc(mesh,MVERTEX,ME_MasterParID(me));
  idx = 0;
  while(mf = MESH_Next_GhostFace(mesh,&idx))
    MESH_Flag_Has_GhostEnts_From_Proc(mesh,MFACE,MF_MasterParID(mf));
  idx = 0;
  while(mr = MESH_Next_GhostRegion(mesh,&idx))
    MESH_Flag_Has_GhostEnts_From_Proc(mesh,MREGION,MR_MasterParID(mr));

  /* derive which processors this processor has overlaps with */

  MESH_Update_OverlapAdj(mesh, rank, num, comm);




  /* local ov num */
  local_ov_num[0] = MESH_Num_OverlapVertices(mesh);
  local_ov_num[1] = MESH_Num_OverlapEdges(mesh);
  local_ov_num[2] = MESH_Num_OverlapFaces(mesh);
  local_ov_num[3] = MESH_Num_OverlapRegions(mesh);

  /* allgather ov num info */
  MPI_Allgather(local_ov_num,4,MPI_INT,global_ov_num,4,MPI_INT,comm);

  for(i = 0; i < num; i++) {
    if (MESH_Has_OverlapEnts_On_Proc(mesh,MANYTYPE,i)) {
      for (mtype = MVERTEX; mtype < MREGION; mtype++) 
        MESH_Set_Num_OverlapEnts_On_Proc(mesh,mtype,i,global_ov_num[4*i+mtype]);
    }
  }

  MSTK_free(global_ov_num);


 return 1;
}



void MESH_Update_OverlapAdj(Mesh_ptr mesh, unsigned int myrank, unsigned int numprocs, MPI_Comm comm) {

  int *local_par_adj = (int *) MSTK_malloc(numprocs*sizeof(int));
  int *global_par_adj = (int *) MSTK_malloc(numprocs*numprocs*sizeof(int));

  /* NOTE: While we can build a direct accessor to mesh_par_adj_flags
     variable of the mesh object, we avoid doing that in order to
     preserve the encapsulation of the mesh object data (More object
     oriented). For now we will duplicate the packing strategy of the
     mesh object to store/exchange the parallel ghost adjacency info,
     but the data separation will mean that we can change the internal
     representation of the mesh but still keep this routine as
     is. Since this operation will not be done too often, we do not
     expect a penalty */

  for (i = 0; i < numprocs; i++) {
    local_par_adj[i] = 0;

    for (mtype = MVERTEX; mtype < MREGION; mtype++) {

      j = MESH_Has_GhostEnts_From_Proc(mesh,mtype,i);

      local_par_adj[i] |= j<<(2*mtype);

    }
  }
     
  /* At this point, it is assumed that this processor ('rank') has
     knowledge of all the processors that it has ghost entities from
     and what type of entities they are. We do an MPI_Allgather so
     that the processor can find out the reverse info, i.e., which
     processors are expecting ghost entities from this processor and
     what type of entities. This info then goes in as the overlap
     entity info for this processor */

  MPI_Allgather(local_par_adj,MPI_INT,global_par_adj,numprocs,MPI_INT,comm);

  /* Now set overlap adjacency flags */

  ovnum = 0;
  procnums = (unsigned int *) malloc(numprocs*sizeof(unsigned int));
  for (i = 0; i < numprocs; i++) {
    for (mtype = MVERTEX; mtype < MREGION; mtype++) {

      j = global_par_adj[i*numprocs + myrank] & 1<<(2*mtype);

      if (j) {
        MESH_Flag_Has_OverlapEnts_On_Proc(mesh,mtype,i);

        if (ovnum == 0 || 
            (ovnum > 0 && procnums[ovnum-1] != i)) {
          procnums[ovnum] = i;
          ovnum++;
        }
      } /* if my proc (rank) has ghosts from proc i */
    }
  }

  /* Set the number of processors with which this processor has overlap */

  MESH_Set_OverlapProcs(mesh,ovnum,procnums);

}

  
#ifdef __cplusplus
}
#endif

