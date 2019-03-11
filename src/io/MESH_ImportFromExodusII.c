/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MSTK.h"
#include "MSTK_private.h"

#include "exodusII.h"
#ifdef EXODUSII_4
#include "exodusII_ext.h"
#endif

#ifdef _MSTK_HAVE_METIS
#include "metis.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif


  /* Read an Exodus II file into MSTK */
  /* 
     The input arguments parallel_opts indicate if we should build
     parallel adjacencies and some options for controlling that 

     if parallel_opts == NULL, then only processor 0 will read the mesh and
     do nothing else
     
     parallel_opts[0] = 1/0  ---- distribute/Don't distribute the mesh
     parallel_opts[1] = 0-3  ---- 0: read/partition on P0 and distribute
                                  1: read/partition on P0, reread local portion
                                     on each process and weave/connect
                                  2: read on multiple processors Pn where n is 
                                  an arithmetic progression (5,10,15 etc) and
                                  distribute to a subset of processors (like P5
     Opts 2,3 NOT ACTIVE          distributes to P0-P4, P10 to P6-P9 etc)
                                  3: read portions of the mesh on each processor
                                  and repartition
     parallel_opts[2] = N    ---- Number of ghost layers around mesh on proc
     parallel_opts[3] = 0/1  ---- partitioning method
                                  0: Metis
                                  1: Zoltan
  */                            

  
  /* Right now we are creating an attribute for material sets, side
     sets and node sets. We are ALSO creating meshsets for each of
     these entity sets. We are also setting mesh geometric entity IDs
     - Should we pick one of the first two? */

  int MESH_ImportFromExodusII(Mesh_ptr mesh, const char *filename,
			      int *parallel_opts, 
                              MSTK_Comm comm) {

  char mesg[256], funcname[32]="MESH_ImportFromExodusII";
  int distributed=0;
  
  ex_init_params exopar;

  FILE *fp;
  char basename[256], modfilename[256], *ext;

  
#ifdef MSTK_HAVE_MPI

  int rank=0, numprocs=1;

  if (comm) {
    MPI_Comm_size(comm,&numprocs);
    MPI_Comm_rank(comm,&rank);
  }

  if (numprocs > 1 && parallel_opts) {

    if (parallel_opts[0] == 1) { /* distribute the mesh */

      int num_ghost_layers = parallel_opts[2];
      
      int have_zoltan = 0, have_metis = 0;
#ifdef _MSTK_HAVE_ZOLTAN
      have_zoltan = 1;
#endif
#ifdef _MSTK_HAVE_METIS
      have_metis = 1;
#endif


      int part_method = parallel_opts[3];
      
      if (part_method == 0) {
        if (!have_metis) {
          MSTK_Report(funcname,"Metis not available. Trying Zoltan",MSTK_WARN);
          part_method = 1;
          if (!have_zoltan) 
            MSTK_Report(funcname,"No partitioner defined",MSTK_FATAL);
        }
      }
      else if (part_method == 1 || part_method == 2) {
        if (!have_zoltan) {
          MSTK_Report(funcname,"Zoltan not available. Trying Metis",MSTK_WARN);
          part_method = 0;
          if (!have_metis) 
            MSTK_Report(funcname,"No partitioner defined",MSTK_FATAL);
        }
      }
      
      if (parallel_opts[1] == 0) {
        
      
        /* Read on processor 0 and distribute to all other processors */
        
        Mesh_ptr globalmesh;
        int topodim;
        
        if (rank == 0) { /* Read only on processor 0 */
          
          globalmesh = MESH_New(MESH_RepType(mesh));
          int read_status = MESH_ReadExodusII_Serial(globalmesh,filename,rank);
          
          if (!read_status) {
            sprintf(mesg,"Could not read Exodus II file %s successfully\n",
                    filename);
            MSTK_Report(funcname,mesg,MSTK_FATAL);
          }
          
          topodim = (MESH_Num_Regions(globalmesh) == 0) ? 2 : 3;
        }
        else
          globalmesh = NULL;
        
        int with_attr = 1;      
        int del_inmesh = 1;
        int dist_status = MSTK_Mesh_Distribute(globalmesh, &mesh, &topodim, 
                                               num_ghost_layers,
                                               with_attr, part_method, 
                                               del_inmesh, comm);
        if (!dist_status)
          MSTK_Report(funcname,
                      "Could not distribute meshes to other processors",
                      MSTK_FATAL);
        
      }
      else if (parallel_opts[1] == 1) {
	int dim;
	int nelems, num_myelems;
	int *myelems=NULL, *elemgraphoff=NULL, *elemgraphadj=NULL;
	int get_coords = 0;
	double (*elemcen)[3];
	
	if (part_method == 0) {
	  if (rank == 0) {
	    
	    ExodusII_GetElementGraph(filename, &dim, &nelems, &elemgraphoff,
				     &elemgraphadj, get_coords, NULL);
	    
	    /* Metis needs us to send it the graph using idx_t arrays - so
	       copy the arrays verbatim */
	  
	    idx_t *xadj = (idx_t *) malloc((nelems+1)*sizeof(idx_t));
	    idx_t *adjncy = (idx_t *) malloc(elemgraphoff[nelems]*sizeof(idx_t));
	
	    memcpy(xadj, elemgraphoff, (nelems+1)*sizeof(idx_t));
	    memcpy(adjncy, elemgraphadj, elemgraphoff[nelems]*sizeof(idx_t));

	    idx_t wtflag = 0;        /* No weights are specified */
	    idx_t *vwgt = NULL;
	    idx_t *adjwgt = NULL;
	  
	    idx_t numflag = 0;    /* C style numbering of elements (nodes of
				     the dual graph) */
	    idx_t ngraphvtx = nelems; 
	    idx_t numparts = numprocs;  

	    idx_t *idxpart = (idx_t *) malloc(nelems*sizeof(idx_t));

	    idx_t ncons = 1;  /* Number of constraints */
	    idx_t nedgecut;
	    idx_t *vsize = NULL;  
	    real_t *tpwgts = NULL;
	    real_t *ubvec = NULL;
	    idx_t metisopts[METIS_NOPTIONS];

	    METIS_SetDefaultOptions(metisopts);
	    metisopts[METIS_OPTION_NUMBERING] = 0;

	    if (numprocs <= 8)
	      METIS_PartGraphRecursive(&ngraphvtx, &ncons, xadj, adjncy, vwgt,
				       vsize, adjwgt, &numparts, tpwgts, ubvec,
				       metisopts, &nedgecut, idxpart);
	    else
	      METIS_PartGraphKway(&ngraphvtx, &ncons, xadj, adjncy, vwgt, vsize,
				  adjwgt, &numparts, tpwgts, ubvec, metisopts,
				  &nedgecut, idxpart);

	    free(xadj);
	    free(adjncy);


	    // Collect the elements on each partition

	    int *part_num_elems = (int *) calloc(numprocs, sizeof(int));
	    int **part_elem_gids = (int **) malloc(numprocs*sizeof(int *));
	    for (int i = 0; i < nelems; i++)
	      part_num_elems[idxpart[i]]++;
	    for (int i = 0; i < numprocs; i++)
	      part_elem_gids[i] = (int *) malloc(part_num_elems[i]*sizeof(int));

	    for (int i = 0; i < numprocs; i++)
	      part_num_elems[i] = 0;  // reset so we can use as counters
	    for (int i = 0; i < nelems; i++) {
	      int p = idxpart[i];
	      int n = part_num_elems[p];
	      part_elem_gids[p][n] = i+1;
	      part_num_elems[p]++;
	    }

	    
	    // straight copy into this partition's arrays
	    num_myelems = part_num_elems[0]; 
	    myelems = (int *) malloc(part_num_elems[0]*sizeof(int));	  
	    memcpy(myelems, part_elem_gids[0], part_num_elems[0]*sizeof(int));

	    // Now send them to all other partitions
	    MPI_Request *requests =
	      (MPI_Request *) malloc(3*(numparts-1)*sizeof(MPI_Request));
	      
	    for (int i = 1; i < numprocs; i++) {
	      MPI_Isend(&dim, 1, MPI_INT, i, 3*i, comm, &(requests[3*(i-1)]));
	      MPI_Isend(&part_num_elems[i], 1, MPI_INT, i, 3*i+1, comm,
			&(requests[3*(i-1)+1]));
	      MPI_Isend(part_elem_gids[i], part_num_elems[i], MPI_INT, i, 3*i+2,
			comm, &(requests[3*(i-1)+2]));
	    }
	    if (MPI_Waitall(3*(numparts-1), requests, MPI_STATUSES_IGNORE) !=
		MPI_SUCCESS)
	      MSTK_Report("MSTK_ImportFromExodusII",
			  "Cound not distribute partition info", MSTK_FATAL);
	    
	    for (int i = 0; i < numparts; i++)
	      free(part_elem_gids[i]);
	    free(part_elem_gids);
	    free(part_num_elems);
	    free(requests);
	  
	  } else {
	    // Get partiton info from rank 0
	    MPI_Request request;

	    MPI_Irecv(&dim, 1, MPI_INT, 0, 3*rank, comm, &request);
	    if (MPI_Wait(&request, MPI_STATUS_IGNORE) != MPI_SUCCESS)
	      MSTK_Report("MSTK_ImportFromExodusII",
			  "Could not receive mesh dimension on partiion",
			  MSTK_FATAL);

	    
	    MPI_Irecv(&num_myelems, 1, MPI_INT, 0, 3*rank+1, comm, &request);
	    if (MPI_Wait(&request, MPI_STATUS_IGNORE) != MPI_SUCCESS)
	      MSTK_Report("MSTK_ImportFromExodusII",
			  "Could not receive num elems on partiion",
			  MSTK_FATAL);

	    myelems = (int *) malloc(num_myelems*sizeof(int));
	    MPI_Irecv(myelems, num_myelems, MPI_INT, 0, 3*rank+2, comm, &request);
	    if (MPI_Wait(&request, MPI_STATUS_IGNORE) != MPI_SUCCESS)
	      MSTK_Report("MSTK_ImportFromExodusII",
			  "Could not receive elem GIDs on partiton",
			  MSTK_FATAL);
	  }
	} else {
	  MSTK_Report("MESH_ImportFromExodusII",
		      "Zoltan partitioning for ExodusII graph not implemented",
		      MSTK_FATAL);
	}

	MESH_ReadExodusII_Partial(mesh, filename, rank, num_myelems, myelems);


	/* Weave the meshes together to establish interprocessor connectivity */
	int num_ghost_layers = 1;
	int input_type = 1;
	MSTK_Weave_DistributedMeshes(mesh, dim, num_ghost_layers, input_type,
				     comm);
	
      }
      else if (parallel_opts[1] == 2) {
        
      }
      else if (parallel_opts[1] == 3) {
        
      }

      int parallel_check = MESH_Parallel_Check(mesh,comm);        
      
      if (!parallel_check)
        MSTK_Report(funcname, "Parallel mesh checks failed", MSTK_FATAL);
    }
    else { /* Read the mesh on rank 0 but don't distribute */
      if (rank == 0)
        return (MESH_ReadExodusII_Serial(mesh,filename,0));
    }
  
  } /* if numprocs > 1 */
  else {
    if (rank == 0) 
      return (MESH_ReadExodusII_Serial(mesh,filename,0));
    else
      return 1;
  }

#else /* Serial process */
  return (MESH_ReadExodusII_Serial(mesh,filename,0));
#endif 
  
  return 1;

}


#ifdef __cplusplus
}
#endif
