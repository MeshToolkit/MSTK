#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* 
     this function sends ov attributes and update ghost attributes
     
     attr_name: the attribute name

     Author(s): Rao Garimella
  */

  int MESH_GathScat1Attribute(Mesh_ptr mesh, MAttrib_ptr attrib,
                              MAttOpType optype, MSTK_Comm comm)  {
    int i, j, k, *loc, count;
    int num_ghost=0, num_overlap=0, ival;
    double rval;
    void *pval;
    double *rval_arr=NULL;
    int *list_info_send, recv_index=0, num_recv_procs;
    int global_id, par_id;
    MEntity_ptr ment=NULL;
    MPI_Status status;    


    int myrank, numprocs;
    MPI_Comm_rank(comm,&myrank);
    MPI_Comm_size(comm,&numprocs);

    if (numprocs == 1) return 0;  /* Nothing to do */

    /* get attribute properties */
    MAttType att_type = MAttrib_Get_Type(attrib);
    int ncomp = MAttrib_Get_NumComps(attrib);
    MType mtype = MAttrib_Get_EntDim(attrib);

    if (att_type == POINTER) {
      char mesg[1024], attname[256];
      MAttrib_Get_Name(attrib,attname);
#ifdef DEBUG2
      sprintf(mesg,"Meaningless to update pointer attributes across processors (attribute: %s)", attname);
      MSTK_Report("MESH_GathScat1Attribute",mesg,MSTK_WARN);
#endif
      return 0;
    }

    if (att_type == INT && optype == ATTOP_AVG) {
      char mesg[1024], attname[256];
      MAttrib_Get_Name(attrib,attname);
      sprintf(mesg,"Averaging integer attributes will result in real valued attribute (attribute: %s)",attname);
      MSTK_Report("MESH_GathScat1Attribute",mesg,MSTK_ERROR);
      return 0;
    }


    int glob_paradj_status=0;
    int loc_paradj_status = MESH_ParallelAdj_Current(mesh);
    MPI_Allreduce(&loc_paradj_status,&glob_paradj_status,1,MPI_INT,MPI_MIN,comm);
    if (glob_paradj_status == 0)
      MESH_Update_ParallelAdj(mesh,comm);


    
    /* Get the ghosts to transmit their values back to their masters */

    /* attribute entity type */
    switch (mtype) {
      case MVERTEX:
        num_ghost = MESH_Num_GhostVertices(mesh);
        num_overlap = MESH_Num_OverlapVertices(mesh);
        break;
      case MEDGE:
        num_ghost = MESH_Num_GhostEdges(mesh);
        num_overlap = MESH_Num_OverlapEdges(mesh);
        break;
      case MFACE:
        num_ghost = MESH_Num_GhostFaces(mesh);
        num_overlap = MESH_Num_OverlapFaces(mesh);
        break;
      case MREGION:
        num_ghost = MESH_Num_GhostRegions(mesh);
        num_overlap = MESH_Num_OverlapRegions(mesh);
        break;
      default:
        MSTK_Report("MESH_GathScat1Attribute","Invalid entity type for this operation",MSTK_FATAL);
    }


     /* First we have to communicate how many entities we are sending
       to each processors - and correspondingly, record how many
       entities we are receiving from each processor */
    
    int *send_size = (int *) calloc(numprocs,sizeof(int));
    int *recv_size = (int *) calloc(numprocs,sizeof(int));

    for (i = 0; i < numprocs; i++) {
      if (i == myrank) continue;

      int found;
      int has_ghosts, has_overlaps;

      has_ghosts = MESH_Has_Ghosts_From_Prtn(mesh,i,mtype);
      has_overlaps = MESH_Has_Overlaps_On_Prtn(mesh,i,mtype);

      /* How many entries will we send to processor 'i' ? */
      for (j = 0; j < num_ghost; j++) {
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
          default: /* Case filtered earlier. Won't come here */ 
            MSTK_Report("MESH_GathScat1Attribute","Invalid entity type for this operation",MSTK_FATAL);
	}
	
	if (MEnt_MasterParID(ment) == i) 
	  send_size[i]++;
      }

      /* if the current rank is lower, send first and receive */
      if (myrank < i) {
	if (has_ghosts)
	  MPI_Send(&(send_size[i]),1,MPI_INT,i,i,comm);
	if (has_overlaps)
	  MPI_Recv(&(recv_size[i]),1,MPI_INT,i,myrank,comm,&status);
      }

      /* if the current rank is higher, recv first and send */
      if ( myrank > i ) {
	if (has_overlaps)
	  MPI_Recv(&(recv_size[i]),1,MPI_INT,i,myrank,comm,&status);
	if (has_ghosts)
	  MPI_Send(&(send_size[i]),1,MPI_INT,i,i,comm);
      }

    } 



    /* Info about how many processors we will receive info from, their
       IDs, number of entities they will send etc */

    unsigned int *recv_procs = (unsigned int *) malloc(numprocs*sizeof(unsigned int));

    num_recv_procs = 0;
    for (i = 0; i < numprocs; ++i) {
      if (i == myrank) continue;
      if (MESH_Has_Overlaps_On_Prtn(mesh,i,mtype))
	recv_procs[num_recv_procs++] = i;
    }

    /* this stores the global to local processor id map */
    int *localrank = (int*) malloc((numprocs+1)*sizeof(int));
    for (i = 0; i < num_recv_procs; i++)
      localrank[recv_procs[i]] = i;

    /* data from different processors will be stored in a single
       array end-to-end; so calculate offsets where the data for a 
       processor starts */

    int tot_recv_size = 0;
    int *recv_pos = (int*) malloc((num_recv_procs+1)*sizeof(int));
    recv_pos[0] = 0;
    for (i = 0; i < num_recv_procs; i++) {
      recv_pos[i+1] = recv_pos[i] + recv_size[recv_procs[i]];
      tot_recv_size += recv_size[recv_procs[i]];
    }


    /* Allocate storage for receiving attribute info */

    int *list_info_recv = NULL;
    int *list_value_int_recv = NULL;
    double *list_value_double_recv = NULL;

    if (tot_recv_size) { /* can be zero if it has no master entities */
      list_info_recv = (int *) malloc(tot_recv_size*sizeof(int));
      list_value_int_recv = (int *) malloc(tot_recv_size*sizeof(int));
      list_value_double_recv = (double *) malloc(tot_recv_size*ncomp*sizeof(double));
    }



    /* check whether need to send or recv info */
    for (i = 0; i < numprocs; i++) {
      if (i == myrank) continue;
      int num_ghost = 0;
      switch (mtype) {
        case MVERTEX:
          num_ghost = MESH_Num_GhostVertices(mesh);
          break;
        case MEDGE:
          num_ghost = MESH_Num_GhostEdges(mesh);
          break;
        case MFACE:
          num_ghost = MESH_Num_GhostFaces(mesh);
          break;
        case MREGION:
          num_ghost = MESH_Num_GhostRegions(mesh);
          break;
        default: /* Case filtered earlier. Won't come here */ 
          MSTK_Report("MESH_GathScat1Attribute",
                      "Invalid entity type for this operation",MSTK_FATAL);
      }

      int found;
      int has_ghosts, has_overlaps;
      int *list_info_send;
      int *list_value_int_send;
      double *list_value_double_send;

      has_ghosts = MESH_Has_Ghosts_From_Prtn(mesh,i,mtype);
      has_overlaps = MESH_Has_Overlaps_On_Prtn(mesh,i,mtype);

      if (has_ghosts && send_size[i]) {

	/* collect attribute info to send */
	
	list_info_send = (int *) malloc(send_size[i]*sizeof(int));
	
	list_value_int_send = (int *) malloc(send_size[i]*sizeof(int));
	list_value_double_send = (double *) malloc(send_size[i]*ncomp*sizeof(double));
	
	/* collect attribute info from overlap entities */
        int n = 0;
	for (j = 0; j < num_ghost; j++) {
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
              MSTK_Report("MESH_GathScat1Attribute",
                          "Invalid entity type for this operation",MSTK_FATAL);
	  }
	  if (MEnt_MasterParID(ment) == i) {
	    MEnt_Get_AttVal(ment,attrib,&ival,&rval,&pval);
	    list_info_send[n] = MEnt_GlobalID(ment);
	    if (att_type == INT) {
	      list_value_int_send[n] = ival;
	    }
	    else {
	      if (ncomp == 1)
		list_value_double_send[n] = rval;
	      if (ncomp > 1) {
		rval_arr = (double *)pval;
		for(k = 0; k < ncomp; k++)
		  list_value_double_send[ncomp*n+k] = rval_arr[k];
	      }
	    }
            n++;
            if (n > send_size[i]) 
              MSTK_Report("MESH_GathScat1Attribute",
                          "More values to send then accounted for",MSTK_FATAL);
	  }
	} /* for (j = 0; j < num_ghost; ++j) */
      }
      else {
        list_info_send = NULL;
        list_value_int_send = NULL;
        list_value_double_send = NULL;
      }

      /* if the current rank is lower, send first and receive */
      if (myrank < i) {

	if (has_ghosts && send_size[i]) {

          /* Send the list of global IDs of entities whose values we
           * are sending out */
	  MPI_Send(list_info_send,send_size[i],MPI_INT,i,i,comm);

          /* Send the values themselves */
	  if (att_type == INT)
	    MPI_Send(list_value_int_send,send_size[i],MPI_INT,i,i,comm);
	  else
	    MPI_Send(list_value_double_send,send_size[i]*ncomp,MPI_DOUBLE,i,i,comm);
	}


        /* Note: check recv_size > 0 as two processors may only be
         * sharing entites of lower order than mtype. But doesn't
         * 'has_overlaps' or MESH_Has_Overlaps_On_Prtn check for that
         * above? Anyway, it can't hurt */

	if (has_overlaps && recv_size[i]) {
	  int offset = recv_pos[localrank[i]];

          /* Receive the list of global IDs of entities whose values
           * we are receiving */
          MPI_Recv(&list_info_recv[offset],recv_size[i],MPI_INT,
                   i,myrank,comm,&status);
          MPI_Get_count(&status,MPI_INT,&count);
          if (count > recv_size[i])
            MSTK_Report("MESH_GathScat1Attribute",
                        "Received more entities than expected. Memory corruption will result",MSTK_FATAL);
          
          /* Receive the values themselves - why are we using count
           * and not recv_size[i]? */
          if (att_type == INT) 
            MPI_Recv(&list_value_int_recv[offset],
                     count,MPI_INT,i,myrank,comm,&status);
          else
            MPI_Recv(&list_value_double_recv[offset*ncomp],
                     count*ncomp,MPI_DOUBLE,i,myrank,comm,&status);
        }
      } /* if (myrank < i) */


      /* if the current rank is higher, recv first and send */

      if (myrank > i) {

	if (has_overlaps && recv_size[i]) {
          int offset = recv_pos[localrank[i]];

          /* Receive the list of global IDs of entities whose values
           * we are receiving */
          MPI_Recv(&list_info_recv[offset],recv_size[i],MPI_INT,
                   i,myrank,comm, &status);
          MPI_Get_count(&status,MPI_INT,&count);
          if (count > recv_size[i])
            MSTK_Report("MESH_GathScat1Attribute",
                        "Received more entities than expected. Memory corruption will result",MSTK_FATAL);

          /* Receive the values themselves - why are we using count
           * and not recv_size[i]? */
          if (att_type == INT) 
            MPI_Recv(&list_value_int_recv[offset],
                     count,MPI_INT,i,myrank,comm,&status);
          else
            MPI_Recv(&list_value_double_recv[offset*ncomp],
                     count*ncomp,MPI_DOUBLE,i,myrank,comm,&status);
	}

	if (has_ghosts && send_size[i]) {

          /* Send the list of global IDs of entities whose values we
           * are sending out */
	  MPI_Send(list_info_send,send_size[i],MPI_INT,i,i,comm);

          /* Send the values themselves */
	  if (att_type == INT)
	    MPI_Send(list_value_int_send,send_size[i],MPI_INT,i,i,comm);
	  else
	    MPI_Send(list_value_double_send,send_size[i]*ncomp,MPI_DOUBLE,i,i,comm);
	}
      } /* if (myrank > i) */

      if (list_info_send) free(list_info_send);
      if (list_value_int_send) free(list_value_int_send);
      if (list_value_double_send) free(list_value_double_send);
      
    } /* for (i = 0; i < numprocs; i++) */


    /* Data has been received from ghosts on other processors. Assign
     * it to the master/overlap entities */
    for (j = 0; j < num_overlap; j++) {
      int nent;

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
        default: /* Case filtered earlier. Won't come here */ 
          MSTK_Report("MESH_GathScat1Attribute",
                      "Invalid entity type for this operation",MSTK_FATAL);
      }
    
      global_id = MEnt_GlobalID(ment);
   
      /* Search for the global ID in the list_info_recv array. NOTE
       * that there may be multiple instances and we have to find all
       * of them. The 'optype' variable will inform us how to
       * reconcile them */

      int nfound = 0;
      int final_ival;
      double final_rval;
      double *final_rval_arr = (double *) malloc(ncomp*sizeof(double));

      MEnt_Get_AttVal(ment,attrib,&final_ival,&final_rval,
                      (void *) final_rval_arr);

      for (i = 0; i < num_recv_procs; i++) {
        int offset = recv_pos[i];   /* recv_pos array is num_recv_procs long */
        nent = recv_size[recv_procs[i]]; /* recv_size array is numprocs long */


        /* Search for the global ID in the list_info_recv array */
        /* The search is conducted in chunks corresponding to the
           globalIDs received from each processor because they are
           guaranteed to be sorted, making a binary search possible */
        
        loc = (int *)bsearch(&global_id,
                             &list_info_recv[offset],
                             nent,
                             sizeof(int),
                             compareINT);
        
        if (loc == NULL) continue;

        /* get the index */
        int ind = (int)(loc - &list_info_recv[0]);
        if (att_type == INT) {
          ival = list_value_int_recv[ind];
          switch (optype) {
            case ATTOP_MAX:
              if (ival > final_ival) final_ival = ival;
              break;
            case ATTOP_MIN:
              if (ival < final_ival) final_ival = ival;
              break;
            case ATTOP_SUM: case ATTOP_AVG:
              final_ival += ival;
              break;
            case ATTOP_UNDEF:
              final_ival = ival;
              break;
          }
        }
        else {
          if(ncomp == 1) {
            rval = list_value_double_recv[ind];
            switch (optype) {
              case ATTOP_MAX:
                if (rval > final_rval) final_rval = rval;
                break;
              case ATTOP_MIN:
                if (rval < final_rval) final_rval = rval;
                break;
              case ATTOP_SUM: case ATTOP_AVG:
                final_rval += rval;
                break;
              case ATTOP_UNDEF:
                final_rval = rval;
                break;
            }
          }
          if(ncomp > 1) {
            int allzero = 1;
            for (k = 0; k < ncomp; k++)
              if (list_value_double_recv[ind*ncomp+k] != 0.0) {
                allzero = 0;
                break;
              }
            if (!allzero) {
              rval_arr = (double *)malloc(ncomp*sizeof(double));
              for(k = 0; k < ncomp; k++)
                rval_arr[k] = list_value_double_recv[ind*ncomp+k];
              double final_norm = 0.0, norm = 0.0;
              switch (optype) {
                case ATTOP_MAX:
                  for (k = 0; k < ncomp; k++) 
                    final_norm += final_rval_arr[k]*final_rval_arr[k];
                  for (k = 0; k < ncomp; k++)
                    norm += rval_arr[k]*rval_arr[k];
                  if (norm > final_norm)
                    memcpy(rval_arr,final_rval_arr,ncomp*sizeof(double));
                  break;
                case ATTOP_MIN:
                  for (k = 0; k < ncomp; k++) 
                    final_norm += final_rval_arr[k]*final_rval_arr[k];
                  for (k = 0; k < ncomp; k++)
                    norm += rval_arr[k]*rval_arr[k];
                  if (norm < final_norm)
                    memcpy(rval_arr,final_rval_arr,ncomp*sizeof(double));
                case ATTOP_SUM: case ATTOP_AVG:
                  for (k = 0; k < ncomp; k++)
                    final_rval_arr[k] += rval_arr[k];
                  break;
                default:
                  for (k = 0; k < ncomp; k++)
                    final_rval_arr[k] = rval_arr[k];
                  break;
              }
            }	
            else
              rval_arr = NULL;
          }
        }

        nfound++;
      } /* for (i = 0; i < num_recv_procs; i++) */

      if (!nfound)
        MSTK_Report("MESH_GathScat1Attribute",
                    "Cannot find global ID in list",MSTK_ERROR);
      
      /* if we are asked to compute an average divide by number of
       * found items + 1 (for the value on this processor) */

      if (optype == ATTOP_AVG) {
        if (att_type == DOUBLE) final_rval /= nfound+1;
        if ((att_type == VECTOR || att_type == TENSOR) && rval_arr)
          for (k = 0; k < ncomp; k++)
            final_rval_arr[k] /= nfound+1;
      }
      
      if ((att_type == INT && final_ival) ||
          (att_type == DOUBLE && final_rval) ||
          ((att_type == VECTOR || att_type == TENSOR) && final_rval_arr))
        MEnt_Set_AttVal(ment,attrib,final_ival,final_rval,(void*) final_rval_arr);

    }
    
        
    /* release the send buffer */
    if (list_info_recv) free(list_info_recv);
    if (list_value_int_recv) free(list_value_int_recv);
    if (list_value_double_recv) free(list_value_double_recv);
    free(recv_procs);
    free(recv_pos);    
    free(localrank);  
    free(send_size);
    free(recv_size);


    /* Now call the routine to get masters to update the ghosts */

    MESH_Update1Attribute(mesh,attrib,comm);

    
    return 1;
}
  
#ifdef __cplusplus
}
#endif
  
  
