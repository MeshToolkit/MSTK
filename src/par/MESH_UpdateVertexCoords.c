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
     Update vertex coordinates across processors
     
     Author(s): Duo Wang, Rao Garimella
  */

  int MESH_UpdateVertexCoords(Mesh_ptr mesh) {
    int i, j, k, *loc, count;
    int num_ghost=0, num_ov=0, send_size,recv_size;
    int *gid_list_send, *gid_list_recv, recv_index, num_recv_procs;
    int global_id, par_id;
    double xyz[3];
    MType mtype;
    MAttType att_type;
    MVertex_ptr mv;
    MPI_Status status;

    MPI_Comm comm = MSTK_Comm();
    int myrank = MSTK_Comm_rank();
    int numprocs = MSTK_Comm_size();

    MPI_Barrier(comm);

    num_ghost = MESH_Num_GhostVertices(mesh);
    num_ov = MESH_Num_OverlapVertices(mesh);


    /* start collecting attribute info to send */

    gid_list_send = (int *)MSTK_malloc(num_ov*sizeof(int));
  
    double *coords_send = (double *)MSTK_malloc(num_ov*3*sizeof(double));
  
    /* collect vertex global ids and coordinates from overlap vertices */
    for (j = 0; j < num_ov; j++) {
      mv = MESH_OverlapVertex(mesh,j);
      MV_Coords(mv,xyz);
      gid_list_send[j] = MEnt_GlobalID(mv);
      for(k = 0; k < 3; k++)
	coords_send[3*j+k] = xyz[k];
    }


    /* Info about how many processors we will receive info from, their
       IDs, number of entities they will send etc */

    num_recv_procs = MESH_Num_GhostPrtns(mesh);

    /* this stores the global to local processor id map */
    int *rank_g2l = (int*) MSTK_malloc((numprocs+1)*sizeof(int));

    unsigned int *recv_procs = (unsigned int *) MSTK_malloc(num_recv_procs*sizeof(unsigned int));
    int *recv_pos = (int*) MSTK_malloc((num_recv_procs+1)*sizeof(int));


    MESH_GhostPrtns(mesh,recv_procs);
    
    recv_pos[0] = 0;
    for (i = 0; i < num_recv_procs; i++) {
      rank_g2l[recv_procs[i]] = i;
      recv_pos[i+1] = recv_pos[i] + MESH_Num_Recv_From_Prtn(mesh,recv_procs[i],MVERTEX);
    }

    /* Allocate storage for receiving attribute info */

    gid_list_recv = (int *)MSTK_malloc(recv_pos[num_recv_procs]*sizeof(int));
    double *coords_recv = (double *)MSTK_malloc(recv_pos[num_recv_procs]*3*sizeof(double));
  


    /* check whether need to send or recv info */
    for (i = 0; i < numprocs; i++) {
      int found;
      int has_ghosts, has_overlaps;

      has_ghosts = MESH_Has_Ghosts_From_Prtn(mesh,i,MVERTEX);
      has_overlaps = MESH_Has_Overlaps_On_Prtn(mesh,i,MVERTEX);

      /* How many entries will we send to processor 'i' ? */
      send_size = num_ov;

      /* search for the index of the processor in recv_procs list */
      found = 0;
      for (k = 0; k < num_recv_procs; k++) {
	if (recv_procs[k] == i) { 	  /* found the processor */
	  recv_index = k;
	  found = 1;
	  break;
	}
      }

      /* How many entries will we receiver from processor 'i' */
      recv_size = found ? MESH_Num_Recv_From_Prtn(mesh,recv_procs[recv_index],MVERTEX) : 0;

      /* if the current rank is lower, send first and receive */
      if (myrank < i) {

	if ( has_overlaps && send_size ) {
	  MPI_Send(gid_list_send,send_size,MPI_INT,i,i,comm);
#ifdef DEBUG_MAX
	  fprintf(stderr,"send %d attr to processor %d on rank %d\n",send_size,i,myrank);
#endif

	  MPI_Send(coords_send,send_size*3,MPI_DOUBLE,i,i,comm);
	}


	if( (has_ghosts) ) {

	/* Check for recv_size which can be zero if the two processors
	   do not need to exchange entities of mtype but only lower
	   order entities */

	  if (recv_size) {
	    MPI_Recv(&gid_list_recv[recv_pos[recv_index]],recv_size,MPI_INT,i,
		     myrank,comm, &status);
	    MPI_Get_count(&status,MPI_INT,&count);
#ifdef DEBUG_MAX
	    fprintf(stderr,"receive %d attr from processor %d on rank %d\n",count,i,myrank);
#endif

	    MPI_Recv(&coords_recv[recv_pos[recv_index]*3],
		     count*3,MPI_DOUBLE,i,myrank,comm,&status);
	  }
	}
      }

      /* if the current rank is higher, recv first and send */
      if ( myrank > i ) {
	if ( (has_ghosts) ) {
	  if (recv_size) {
	    MPI_Recv(&gid_list_recv[recv_pos[recv_index]],recv_size,MPI_INT,i,
		     myrank,comm, &status);
	    MPI_Get_count(&status,MPI_INT,&count);
#ifdef DEBUG_MAX
	    fprintf(stderr,"receive %d attr from processor %d on rank %d\n",count,i,myrank);
#endif
	    MPI_Recv(&coords_recv[recv_pos[recv_index]*3],
		     count*3,MPI_DOUBLE,i,myrank,comm,&status);
	  }
	}
	if ( has_overlaps && send_size ) {
	  MPI_Send(gid_list_send,send_size,MPI_INT,i,i,comm);
#ifdef DEBUG_MAX
	  fprintf(stderr,"send %d attr to processor %d on rank %d\n",send_size,i,myrank);
#endif
	  MPI_Send(coords_send,send_size*3,MPI_DOUBLE,i,i,comm);
	}
      }

    }

    /* Data has been received. Assign it to ghosts */
    for(j = 0; j < num_ghost; j++) {
      int nent;

      mv = MESH_GhostVertex(mesh,j);
    
      global_id = MEnt_GlobalID(mv);
    
      par_id = MEnt_MasterParID(mv);

      nent = MESH_Num_Recv_From_Prtn(mesh, par_id, MVERTEX);

      /* since the ov list is already sorted, use binary search */
      loc = (int *)bsearch(&global_id,
			   &gid_list_recv[recv_pos[rank_g2l[par_id]]],
                           nent,
			   sizeof(int),
			   compareINT);
      /* get the index */
      i = (int)(loc - &gid_list_recv[0]);
      for(k = 0; k < 3; k++)
	xyz[k] = coords_recv[i*3+k];
      MV_Set_Coords(mv,xyz);
    }

    /* release the send buffer */
    MSTK_free(gid_list_send);
    MSTK_free(coords_send);
    MSTK_free(gid_list_recv);
    MSTK_free(coords_recv);
    MSTK_free(recv_pos);    
    MSTK_free(rank_g2l);    
    return 1;
}
  
#ifdef __cplusplus
}
#endif
  
  
