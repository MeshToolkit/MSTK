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
     this function updates (or intialize) global info associated with mesh
     this is a collective operation called by every process

     for applications that do not involve topology change,
     this routine can be called only once during initialization

     global_info is an integer array of size num
     global_info[i] on a processor j indicates the relationship of the mesh 
     on processor j with the mesh on processor i.

     1 and 2 bit(from left - ERRRR, SHOULD THAT BE FROM THE RIGHT?? - RVG) 
     indicate relation on vertex
     0(00) no relation,
     1(01) has ghost entities related to processor j,
     2(10) has overlap entities related to processor j,
     3(11) both
     3 and 4 bit indicate relation on edge
     5 and 6 bit indicate relation on face
     7 and 8 bit indicate relation on region

     for example: global_info[3] = 51 (00110011) means mesh on this
     processor has ghost vertices and ghost faces whose master
     processor is 3 and has overlap vertices and overlap faces which
     are ghost entities on processor 3.

     PERHAPS GLOBAL INFO SHOULD BE CALLED MY_MESH_CON OR MESH_PROC_CON
     BECAUSE IT IS DESCRIBING HOW THE MESH ON THIS PROCESSOR IS
     CONNECTED TO THE MESH ON OTHER PROCESSORS. USING THE PHRASE
     GLOBAL_INFO LEADS US TO BELIEVE THAT THIS INFO IS UNIFORM ACROSS
     ALL PROCESSORS BUT THAT IS UNTRUE.



     local_info stores the number of overlap entities for allocating
     recv buffer 

     local_info[0]: number of processors that have overlap
     with mesh on a particular processor

     local_info[1]-local_info[num_recv_rank+1]: processor ids that
     have overlap with mesh on particular processor (What is num_recv_rank?)

     Then num_recv_rank*4 array stores the number of overlap
     entities. (DOES THAT MEAN THAT THE SUBSEQUENT ENTRIES INDICATE
     THE NUMBER OF OVERLAP ENTITIES OF EACH TYPE FOR EACH OF THE
     RELEVANT PROCESSORS?)

     THIS IS ALSO SOME SORT OF CONNECTIVITY INFO - WHAT SHOULD WE RENAME IT?

     Author(s): Duo Wang, Rao Garimella

  */


int MESH_UpdateGlobalInfo(Mesh_ptr mesh, int rank, int num,  MPI_Comm comm) {
  int i,idx,nv,ne,nf,nr,local_ov_num[4];
  int ebit, ov_index, ov_index2;
  MVertex_ptr mv;
  MEdge_ptr me;
  MFace_ptr mf;
  MRegion_ptr mr;
  int *local_info, *global_info;


  int *global_ov_num = (int *) MSTK_malloc(4*num*sizeof(int));
  int *global_ranks = (int *) MSTK_malloc(num*num*sizeof(int));

  
  global_info = (int *) MSTK_malloc(num*sizeof(int));
  
  /* local ov num */
  local_ov_num[0] = MESH_Num_OverlapVertices(mesh);
  local_ov_num[1] = MESH_Num_OverlapEdges(mesh);
  local_ov_num[2] = MESH_Num_OverlapFaces(mesh);
  local_ov_num[3] = MESH_Num_OverlapRegions(mesh);

  /* allgather ov num info */
  MPI_Allgather(local_ov_num,4,MPI_INT,global_ov_num,4,MPI_INT,comm);

  nv = MESH_Num_GhostVertices(mesh);
  ne = MESH_Num_GhostEdges(mesh);
  nf = MESH_Num_GhostFaces(mesh);
  nr = MESH_Num_GhostRegions(mesh);
  int *local_ranks = (int *)MSTK_malloc(num*sizeof(int));
  for(i = 0; i < num; i++)
    local_ranks[i] = 0;

  /* set bit to 1 */
  if(nv) {
    idx = 0;
    while(mv = MESH_Next_GhostVertex(mesh,&idx))
      local_ranks[MV_MasterParID(mv)] |= 1;
  }
  if(ne) {
    idx = 0;
    while(me = MESH_Next_GhostEdge(mesh,&idx))
      local_ranks[ME_MasterParID(me)] |= (1<<2);
  }
  if(nf) {
    idx = 0;
    while(mf = MESH_Next_GhostFace(mesh,&idx))
      local_ranks[MF_MasterParID(mf)] |= (1<<4);
  }
  if(nr) {
    idx = 0;
    while(mr = MESH_Next_GhostRegion(mesh,&idx))
      local_ranks[MR_MasterParID(mr)] |= (1<<6);
  }

  /* allgather info */
  MPI_Allgather(local_ranks,num,MPI_INT,global_ranks,num,MPI_INT,comm);

  /* check which processor has ghost entity on rank processor */
  for(i = 0; i < num; i++) 
    global_info[i] = local_ranks[i] | (global_ranks[i*num+rank] << 1);
  


  /* get ov info */
  ov_index = 0;
  for(i = 0; i < num; i++)
    if(local_ranks[i])
      ov_index++;
  /*
    local_info:
    local_info[0]: number of master partitions of ghost
    then is the partition id;
    then is the number of ov vertices, edges, faces and regions
  */
  local_info = (int *) MSTK_malloc((5*(ov_index)+1)*sizeof(int));      
  for(i = 0; i < 5*(ov_index)+1; i++)
    local_info[i] = 0;
  local_info[0] = ov_index;
  ov_index2 = 0;
  for(i = 0; i < num; i++) {
    ebit = local_ranks[i];
    if(ebit) {
      local_info[ov_index2+1] = i;
      if(ebit & 1)
	local_info[ov_index+ov_index2*4+1] = global_ov_num[4*i];
      if( (ebit>>2) & 1)
	local_info[ov_index+ov_index2*4+2] = global_ov_num[4*i+1];
      if( (ebit>>4) & 1)
	local_info[ov_index+ov_index2*4+3] = global_ov_num[4*i+2];
      if( (ebit>>6) & 1)
	local_info[ov_index+ov_index2*4+4] = global_ov_num[4*i+3];
      ov_index2++;
    }
  }


  MESH_Set_GlobalInfo(mesh,global_info);
  MESH_Set_LocalInfo(mesh,local_info);

  MSTK_free(local_ranks);
  MSTK_free(global_ranks);
  MSTK_free(global_ov_num);


 return 1;
}
  
#ifdef __cplusplus
}
#endif

