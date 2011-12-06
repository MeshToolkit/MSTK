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
     
     must call MESH_UpdateGlobalInfo() first
     called by every process
     
     attr_name: the attribute name

     Author(s): Duo Wang, Rao Garimella
  */

  int ov_compare(const void * a, const void * b) {
    if ( *(int*)a > *(int*)b ) 
      return 1;
    else if ( *(int*)a < *(int*)b ) 
      return -1;
    else 
      return 0;
  }


  int MESH_UpdateAttr(Mesh_ptr mesh, const char *attr_name, int rank, int num,  MPI_Comm comm) {
    int i, j, k, *loc, ebit, count;
    int num_ghost=0, num_ov=0, ncomp, ival,send_size,recv_size;
    double rval;
    void *pval;
    double *rval_arr=NULL;
    int *list_info_send, *list_info_recv, recv_index, num_recv_rank;
    int global_id, par_id;
    MType mtype;
    MAttType att_type;
    MEntity_ptr ment=NULL;
    MAttrib_ptr attrib;
    MPI_Status status;

    int *local_info = MESH_LocalInfo(mesh);
    int *global_info = MESH_GlobalInfo(mesh);
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
  
    /* this stores the global to local processor id map */
    int *rank_g2l = (int*) MSTK_malloc((num+1)*sizeof(int));
    num_recv_rank = local_info[0];
    int *recv_pos = (int*) MSTK_malloc((num_recv_rank+1)*sizeof(int));
    recv_pos[0] = 0;
    for (i = 0; i < num_recv_rank; i++) {
      rank_g2l[local_info[i+1]] = i;
      recv_pos[i+1] = recv_pos[i]+local_info[4*i+num_recv_rank+mtype+1];
    }


    /* attribute index and global id */ 
    list_info_send = (int *)MSTK_malloc(num_ov*sizeof(int));
    list_info_recv = (int *)MSTK_malloc(recv_pos[num_recv_rank]*sizeof(int));
  
    /* attribute values */
    int *list_value_int_send = (int *)MSTK_malloc(num_ov*ncomp*sizeof(int));
    double *list_value_double_send = (double *)MSTK_malloc(num_ov*ncomp*sizeof(double));
  
    int *list_value_int_recv = (int *)MSTK_malloc(recv_pos[num_recv_rank]*ncomp*sizeof(int));
    double *list_value_double_recv = (double *)MSTK_malloc(recv_pos[num_recv_rank]*ncomp*sizeof(double));
  
    /* collect ov data */
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

    /* check whether need to send or recv info */
    recv_index = 0;
    for (i = 0; i < num; i++) {

      /* shift to proper bit based on the type of entity attribute lives on */
      ebit = global_info[i] >> 2*mtype;

      /* send and recv size */
      send_size = num_ov;
      recv_size = local_info[4*recv_index+num_recv_rank+mtype+1];

      /* if the current rank is lower, send first and receive */
      if (rank < i) {
	if ( (ebit & 2) && send_size ) {
	  MPI_Send(list_info_send,send_size,MPI_INT,i,i,comm);
#ifdef DEBUG_MAX
	  fprintf(stderr,"send %d attr to processor %d on rank %d\n",send_size,i,rank);
#endif

	  if (att_type == INT)
	    MPI_Send(list_value_int_send,send_size*ncomp,MPI_INT,i,i,comm);
	  else
	    MPI_Send(list_value_double_send,send_size*ncomp,MPI_DOUBLE,i,i,comm);
	}


	if( (ebit & 1) ) {

	/* Check for recv_size which can be zero if the two processors
	   do not need to exchange entities of mtype but only lower
	   order entities */

	  if (recv_size) {
	    MPI_Recv(&list_info_recv[recv_pos[recv_index]],recv_size,MPI_INT,i,
		     rank,comm, &status);
	    MPI_Get_count(&status,MPI_INT,&count);
#ifdef DEBUG_MAX
	    fprintf(stderr,"receive %d attr from processor %d on rank %d\n",count,i,rank);
#endif
	    if (att_type == INT) 
	      MPI_Recv(&list_value_int_recv[recv_pos[recv_index]*ncomp],
		       count*ncomp,MPI_INT,i,rank,comm,&status);
	    else
	      MPI_Recv(&list_value_double_recv[recv_pos[recv_index]*ncomp],
		       count*ncomp,MPI_DOUBLE,i,rank,comm,&status);
	  }
	}
      }

      /* if the current rank is higher, recv first and send */
      if ( rank > i ) {
	if ( (ebit & 1) ) {
	  if (recv_size) {
	    MPI_Recv(&list_info_recv[recv_pos[recv_index]],recv_size,MPI_INT,i,
		     rank,comm, &status);
	    MPI_Get_count(&status,MPI_INT,&count);
#ifdef DEBUG_MAX
	    fprintf(stderr,"receive %d attr from processor %d on rank %d\n",count,i,rank);
#endif
	    if (att_type == INT) 
	      MPI_Recv(&list_value_int_recv[recv_pos[recv_index]*ncomp],
		       count*ncomp,MPI_INT,i,rank,comm,&status);
	    else
	      MPI_Recv(&list_value_double_recv[recv_pos[recv_index]*ncomp],
		       count*ncomp,MPI_DOUBLE,i,rank,comm,&status);
	  }
	}
	if ( (ebit & 2) && send_size ) {
	  MPI_Send(list_info_send,send_size,MPI_INT,i,i,comm);
#ifdef DEBUG_MAX
	  fprintf(stderr,"send %d attr to processor %d on rank %d\n",send_size,i,rank);
#endif
	  if (att_type == INT)
	    MPI_Send(list_value_int_send,send_size*ncomp,MPI_INT,i,i,comm);
	  else
	    MPI_Send(list_value_double_send,send_size*ncomp,MPI_DOUBLE,i,i,comm);
	}
      }

      if (i != rank) recv_index++;
    }

    /* assign ghost data */
    for(j = 0; j < num_ghost; j++) {
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
      /* since the ov list is already sorted, use binary search */
      loc = (int *)bsearch(&global_id,
			   &list_info_recv[recv_pos[rank_g2l[par_id]]],
			   local_info[4*rank_g2l[par_id]+num_recv_rank+mtype+1],
			   sizeof(int),
			   ov_compare);
      /* get the index */
      i = (int)(loc - &list_info_recv[0]);
      if (att_type == INT)
	ival = list_value_int_recv[i];
      else {
	if(ncomp == 1) 
	  rval = list_value_double_recv[i];
	if(ncomp > 1)
	  {
	    rval_arr = (double *)MSTK_malloc(ncomp*sizeof(double));
	    for(k = 0; k < ncomp; k++)
	      rval_arr[k] = list_value_double_recv[i*ncomp+k];
	  }	
      }
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
  
  
