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
     this function sends ov attributes and update ghost attributes
     
     must call MESH_Update_ParallelAdj() first
     called by every process
     
     attr_name: the attribute name

     Author(s): Duo Wang, Rao Garimella
  */

  int MESH_UpdateAttr(Mesh_ptr mesh, const char *attr_name, MSTK_Comm comm) {
    int i, j, k, *loc, ebit, count;
    int num_ghost=0, num_ov=0, ncomp, ival,send_size,recv_size;
    double rval;
    void *pval;
    double *rval_arr=NULL;
    int *list_info_send, *list_info_recv, recv_index=0, num_recv_procs;
    int global_id, par_id;
    MType mtype;
    MAttType att_type;
    MEntity_ptr ment=NULL;
    MAttrib_ptr attrib;
    MPI_Status status;

    int myrank, numprocs;
    MPI_Comm_rank(comm,&myrank);
    MPI_Comm_size(comm,&numprocs);


    attrib = MESH_AttribByName(mesh,attr_name);
    /* if there is no such attribute */
    if (!attrib) {
      MSTK_Report("MESH_UpdateAttr()","There is no such attribute",MSTK_ERROR);
      return 0;
    }

    /* get attribute properties */
    att_type = MAttrib_Get_Type(attrib);
    ncomp = MAttrib_Get_NumComps(attrib);
    mtype = MAttrib_Get_EntDim(attrib);

    if (att_type == POINTER) {
      char mesg[1024], attname[256];
      MAttrib_Get_Name(attrib,attname);
#ifdef DEBUG
      sprintf(mesg,"Meaningless to update pointer attributes across processors (attribute: %s", attname);
      MSTK_Report("MESH_UpdateAttr()",mesg,MSTK_WARN);
#endif
      return 0;
    }

    /*  printf("attrib type is %d, ncomp is %d, mtype is %d\n",att_type,ncomp,mtype); */
    /* attribute entity type */
    switch (mtype) {
    case MVERTEX:
      num_ghost = MESH_Num_GhostVertices(mesh);
      num_ov = MESH_Num_OverlapVertices(mesh);
      break;
    case MEDGE:
      num_ghost = MESH_Num_GhostEdges(mesh);
      num_ov = MESH_Num_OverlapEdges(mesh);
      break;
    case MFACE:
      num_ghost = MESH_Num_GhostFaces(mesh);
      num_ov = MESH_Num_OverlapFaces(mesh);
      break;
    case MREGION:
      num_ghost = MESH_Num_GhostRegions(mesh);
      num_ov = MESH_Num_OverlapRegions(mesh);
      break;
    case MALLTYPE:
      MSTK_Report("MESH_UpdateAttr","Not implemented for MALLTYPE",MSTK_WARN);
      return 0;
      break;
    default:
      MSTK_Report("MESH_UpdateAttr()","Invalid entity type",MSTK_FATAL);
    }
  



    /* start collecting attribute info to send */

    list_info_send = (int *)MSTK_malloc(num_ov*sizeof(int));
  
      int *list_value_int_send = (int *)MSTK_malloc(num_ov*ncomp*sizeof(int));
    double *list_value_double_send = (double *)MSTK_malloc(num_ov*ncomp*sizeof(double));
  
    /* collect attribute info from overlap entities */
    for (j = 0; j < num_ov; j++) {
      switch (mtype) {
      case MVERTEX:
	ment = MESH_OverlapVertex(mesh,j);
	break;
      case MEDGE:
	ment = MESH_OverlapEdge(mesh,j);
	break;
      case MFACE:
	ment = MESH_OverlapFace(mesh,j);
	break;
      case MREGION:
	ment = MESH_OverlapRegion(mesh,j);
	break;
      default:
	MSTK_Report("MESH_SendAttr()","Invalid entity type",MSTK_FATAL);
      }
      MEnt_Get_AttVal(ment,attrib,&ival,&rval,&pval);
      list_info_send[j] = MEnt_GlobalID(ment);
      if (att_type == INT) {
	list_value_int_send[j] = ival;
      }
      else {
	if (ncomp == 1)
	  list_value_double_send[j] = rval;
	if (ncomp > 1) {
	  rval_arr = (double *)pval;
	  for(k = 0; k < ncomp; k++)
	    list_value_double_send[ncomp*j+k] = rval_arr[k];
	}
      }
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
      recv_pos[i+1] = recv_pos[i] + MESH_Num_Recv_From_Prtn(mesh,recv_procs[i],mtype);
    }

    /* Allocate storage for receiving attribute info */

    list_info_recv = (int *)MSTK_malloc(recv_pos[num_recv_procs]*sizeof(int));
    int *list_value_int_recv = (int *)MSTK_malloc(recv_pos[num_recv_procs]*ncomp*sizeof(int));
    double *list_value_double_recv = (double *)MSTK_malloc(recv_pos[num_recv_procs]*ncomp*sizeof(double));
  


    /* check whether need to send or recv info */
    for (i = 0; i < numprocs; i++) {
      int found;
      int has_ghosts, has_overlaps;

      /* shift to proper bit based on the type of entity attribute lives on */
      /* ebit = mesh_par_adj_flags[i] >> 2*mtype; */

      has_ghosts = MESH_Has_Ghosts_From_Prtn(mesh,i,mtype);
      has_overlaps = MESH_Has_Overlaps_On_Prtn(mesh,i,mtype);

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

      recv_size = found ? MESH_Num_Recv_From_Prtn(mesh,recv_procs[recv_index],mtype) : 0;

      /* if the current rank is lower, send first and receive */
      if (myrank < i) {

	if ( has_overlaps && send_size ) {
	  MPI_Send(list_info_send,send_size,MPI_INT,i,i,comm);
#ifdef DEBUG_MAX
	  fprintf(stderr,"send %d attr to processor %d on rank %d\n",send_size,i,myrank);
#endif

	  if (att_type == INT)
	    MPI_Send(list_value_int_send,send_size*ncomp,MPI_INT,i,i,comm);
	  else
	    MPI_Send(list_value_double_send,send_size*ncomp,MPI_DOUBLE,i,i,comm);
	}


	if( (has_ghosts) ) {

	/* Check for recv_size which can be zero if the two processors
	   do not need to exchange entities of mtype but only lower
	   order entities */

	  if (recv_size) {
	    MPI_Recv(&list_info_recv[recv_pos[recv_index]],recv_size,MPI_INT,i,
		     myrank,comm, &status);
	    MPI_Get_count(&status,MPI_INT,&count);
#ifdef DEBUG_MAX
	    fprintf(stderr,"receive %d attr from processor %d on rank %d\n",count,i,myrank);
#endif
	    if (att_type == INT) 
	      MPI_Recv(&list_value_int_recv[recv_pos[recv_index]*ncomp],
		       count*ncomp,MPI_INT,i,myrank,comm,&status);
	    else
	      MPI_Recv(&list_value_double_recv[recv_pos[recv_index]*ncomp],
		       count*ncomp,MPI_DOUBLE,i,myrank,comm,&status);
	  }
	}
      }

      /* if the current rank is higher, recv first and send */
      if ( myrank > i ) {
	if ( (has_ghosts) ) {
	  if (recv_size) {
	    MPI_Recv(&list_info_recv[recv_pos[recv_index]],recv_size,MPI_INT,i,
		     myrank,comm, &status);
	    MPI_Get_count(&status,MPI_INT,&count);
#ifdef DEBUG_MAX
	    fprintf(stderr,"receive %d attr from processor %d on rank %d\n",count,i,myrank);
#endif
	    if (att_type == INT) 
	      MPI_Recv(&list_value_int_recv[recv_pos[recv_index]*ncomp],
		       count*ncomp,MPI_INT,i,myrank,comm,&status);
	    else
	      MPI_Recv(&list_value_double_recv[recv_pos[recv_index]*ncomp],
		       count*ncomp,MPI_DOUBLE,i,myrank,comm,&status);
	  }
	}
	if ( has_overlaps && send_size ) {
	  MPI_Send(list_info_send,send_size,MPI_INT,i,i,comm);
#ifdef DEBUG_MAX
	  fprintf(stderr,"send %d attr to processor %d on rank %d\n",send_size,i,myrank);
#endif
	  if (att_type == INT)
	    MPI_Send(list_value_int_send,send_size*ncomp,MPI_INT,i,i,comm);
	  else
	    MPI_Send(list_value_double_send,send_size*ncomp,MPI_DOUBLE,i,i,comm);
	}
      }

    }

    /* Data has been received. Assign it to ghosts */
    for(j = 0; j < num_ghost; j++) {
      int nent;

      switch (mtype) {
      case MVERTEX:
	ment = MESH_GhostVertex(mesh,j);
	break;
      case MEDGE:
	ment = MESH_GhostEdge(mesh,j);
	break;
      case MFACE:
	ment = MESH_GhostFace(mesh,j);
	break;
      case MREGION:
	ment = MESH_GhostRegion(mesh,j);
	break;
      default:
	MSTK_Report("MESH_UpdateAttr()","Invalid entity type",MSTK_FATAL);
      }
    
      global_id = MEnt_GlobalID(ment);
    
      par_id = MEnt_MasterParID(ment);

      nent = MESH_Num_Recv_From_Prtn(mesh, par_id, mtype);

      /* since the ov list is already sorted, use binary search */
      loc = (int *)bsearch(&global_id,
			   &list_info_recv[recv_pos[rank_g2l[par_id]]],
                           nent,
			   sizeof(int),
			   compareINT);

      if (loc == NULL)
        MSTK_Report("MESH_UpdateAttr","Cannot find global ID in list",MSTK_ERROR);

      /* get the index */
      i = (int)(loc - &list_info_recv[0]);
      if (att_type == INT)
	ival = list_value_int_recv[i];
      else {
	if(ncomp == 1) 
	  rval = list_value_double_recv[i];
	if(ncomp > 1)
	  {
            int allzero = 1;
            for (k = 0; k < ncomp; k++)
              if (list_value_double_recv[i*ncomp+k] != 0.0) {
                allzero = 0;
                break;
              }
            if (!allzero) {
              rval_arr = (double *)MSTK_malloc(ncomp*sizeof(double));
              for(k = 0; k < ncomp; k++)
                rval_arr[k] = list_value_double_recv[i*ncomp+k];
            }	
            else
              rval_arr = NULL;
          }
      }
      if ((att_type == INT && ival) ||
          (att_type == DOUBLE && rval) ||
          ((att_type == VECTOR || att_type == TENSOR) && rval_arr))
        MEnt_Set_AttVal(ment,attrib,ival,rval,(void*) rval_arr);
    }

    /* release the send buffer */
    MSTK_free(list_info_send);
    MSTK_free(list_value_int_send);
    MSTK_free(list_value_double_send);
    MSTK_free(list_info_recv);
    MSTK_free(list_value_int_recv);
    MSTK_free(list_value_double_recv);
    MSTK_free(recv_pos);    
    MSTK_free(rank_g2l);    
    return 1;
}
  
#ifdef __cplusplus
}
#endif
  
  
