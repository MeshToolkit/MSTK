#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* 
     MSTK high level parallel functions for user/application

     Author(s): Duo Wang, Rao Garimella

  */


  /* Partition a given mesh into 'num' submeshes, adding a 'ring'
     layers of ghost elements around each partition. If 'with_attr' is
     1, attributes from the mesh are copied onto the submeshes. This
     routine does not send the meshes to other partitions */

  int MSTK_Mesh_Partition(Mesh_ptr mesh, int num, int *part,  int ring, 
			  int with_attr, int del_inmesh, Mesh_ptr *submeshes) {
    int i, j, natt, nset, idx, ival;
    double rval;
    char attr_name[256];
    MAttrib_ptr attrib, g2latt, l2gatt;
    MSet_ptr mset;
    MVertex_ptr mv;
    MEdge_ptr me;
    MFace_ptr mf;
    MRegion_ptr mr;
    List_ptr g2llist;


    /* Create an attribute to keep track of the connections from
       entities of the global mesh to entities of submeshes */
    /* This attribute value will be populated with lists in
       MESH_Partition along with another attribute called Local2Global
       and used in MESH_BuildPBoundary and MESH_AddGhost. They are not
       needed subsequently */

    g2latt = MAttrib_New(mesh,"Global2Local",POINTER,MALLTYPE);


    /* Split the mesh into 'num' submeshes */

    MESH_Partition(mesh, num, part, submeshes);

    for(i = 0; i < num; i++) {
      /* Tag entities as being in the partition interior or on the
         partition boundary */

      MESH_BuildPBoundary(mesh,submeshes[i]);

      /* Add ghost layers */

      MESH_AddGhost(mesh,submeshes[i],i,ring);
    }

    if (with_attr) {

      /* Also transmit attributes to the partitions */
      
      natt = MESH_Num_Attribs(mesh);
      for(j = 0; j < natt; j++) {
	attrib = MESH_Attrib(mesh,j);
	if (attrib == g2latt) continue;
	
	MAttrib_Get_Name(attrib,attr_name);
        MESH_CopyAttr(mesh,num,submeshes,attr_name);
      }
      
      /* Split the sets into subsets and send them to the partitions */
      nset = MESH_Num_MSets(mesh);
      for (j = 0; j < nset; j++) {
	mset = MESH_MSet(mesh,j);
	MESH_CopySet(mesh,num,submeshes,mset);
      }
      
    }

    /* Delete the lists associated with g2latt attribute but don't
       remove attribute itself from each of these entities - it will
       get deleted when the mesh gets deleted */

    idx = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx))) {
      MEnt_Get_AttVal(mv,g2latt,&ival,&rval,&g2llist);
      if (g2llist) List_Delete(g2llist);
    }
	
    idx = 0;
    while ((me = MESH_Next_Edge(mesh,&idx))) {
      MEnt_Get_AttVal(me,g2latt,&ival,&rval,&g2llist);
      if (g2llist) List_Delete(g2llist);
    }
	
    idx = 0;
    while ((mf = MESH_Next_Face(mesh,&idx))) {
      MEnt_Get_AttVal(mf,g2latt,&ival,&rval,&g2llist);
      if (g2llist) List_Delete(g2llist);
    }
	
    idx = 0;
    while ((mr = MESH_Next_Region(mesh,&idx))) {
      MEnt_Get_AttVal(mr,g2latt,&ival,&rval,&g2llist);
      if (g2llist) List_Delete(g2llist);
    }	
	
    if (del_inmesh) MESH_Delete(mesh);
    /*
    printf("global mesh has been partitioned into %d parts\n",num);
    */
    return 1;
  }


  /* Send a mesh to a processor 'torank' with or without attributes */

  int MSTK_SendMesh(Mesh_ptr mesh, int torank, int with_attr, MSTK_Comm comm,
                    int *numreq, int *maxreq, MPI_Request **requests,
                    int *numptrs2free, int *maxptrs2free,
                    void ***ptrs2free) {
    int i, natt, nset;
    char attr_name[256], mset_name[256];
    MAttrib_ptr attrib;
    MSet_ptr mset;
    
    if (requests == NULL)
      MSTK_Report("MSTK_SendMesh","Invalid MPI request buffer",MSTK_FATAL);
    
    if (*maxreq == 0) {
      *maxreq = 25;
      *requests = (MPI_Request *) malloc(*maxreq*sizeof(MPI_Request));
      *numreq = 0;
    }
    else if (*maxreq < (*numreq) + 2) {
      *maxreq *= 2;
      *requests = (MPI_Request *) realloc(*requests,*maxreq*sizeof(MPI_Request));
    }
    
    if (ptrs2free == NULL)
      MSTK_Report("MSTK_SendMesh","Invalid ptrs2free buffer",MSTK_FATAL);
    
    /* Send mesh only */
    
    MESH_SendMesh(mesh,torank,comm,numreq,maxreq,requests,numptrs2free,
                  maxptrs2free,ptrs2free);
    
    
    if (!with_attr) 
      return 1;
    
    /* Send attributes as well */
    
    natt = MESH_Num_Attribs(mesh);
    for(i = 0; i < natt; i++) {
      attrib = MESH_Attrib(mesh,i);
      MAttrib_Get_Name(attrib,attr_name);
      /*
        fprintf(stderr,"Sending attribute %s of type %d of entity dim %d \n",
        attr_name, MAttrib_Get_Type(attrib),
        MAttrib_Get_EntDim(attrib));
      */
      MESH_SendAttr(mesh,attr_name,torank,comm,numreq,maxreq,requests,
                    numptrs2free,maxptrs2free,ptrs2free);
    }
    
    /* Distribute entity sets */
    
    nset = MESH_Num_MSets(mesh);
    for(i = 0; i < nset; i++) {
      mset = MESH_MSet(mesh,i);
      MSet_Name(mset,mset_name);
      MESH_SendMSet(mesh,mset_name,torank,comm,numreq,maxreq,requests,
                    numptrs2free,maxptrs2free,ptrs2free);
    }
    
   return 1;
  }

  /* Receive a mesh of dimension 'dim' from processor 'send_rank' with
     or without attributes onto processor 'rank' */

  int MSTK_RecvMesh(Mesh_ptr mesh, int dim, int fromrank, int with_attr, MSTK_Comm comm) {
    int i, natt, nset;
    char attr_name[256], mset_name[256];
    MAttrib_ptr attrib;
    MSet_ptr mset;

    /* Receive the mesh only */

    MESH_RecvMesh(mesh,dim,fromrank,comm);

    /* Build the sorted lists of ghost and overlap entities */

    MESH_Build_GhostLists(mesh,dim);

    if (!with_attr)
      return 1;

    /* Receive attributes as well */

    natt = MESH_Num_Attribs(mesh);
    for(i = 0; i < natt; i++) {
      attrib = MESH_Attrib(mesh,i);
      MAttrib_Get_Name(attrib,attr_name);
      MESH_RecvAttr(mesh,attr_name,fromrank,comm);
    }


    /* Receive mesh sets as well */

    nset = MESH_Num_MSets(mesh);
    for(i = 0; i < nset; i++) {
      mset = MESH_MSet(mesh,i);
      MSet_Name(mset,mset_name);
      MESH_RecvMSet(mesh,mset_name,fromrank,comm);
    }

    return 1;
  }
    

  /* Read a mesh in, partition it and distribute it to 'num' processors */
  /* The rank, num and comm arguments will go away soon                 */

  int MSTK_Mesh_Read_Distribute(Mesh_ptr *recv_mesh, 
				const char* global_mesh_name, 
				int *dim, int ring, int with_attr,
                                int method, MSTK_Comm comm) {

    int i, ok;

    int rank;
    MPI_Comm_rank(comm,&rank);

    Mesh_ptr mesh=NULL;
    if(rank == 0) {
      mesh = MESH_New(UNKNOWN_REP);
      ok = MESH_InitFromFile(mesh,global_mesh_name,comm);
      if (!ok) {
	fprintf(stderr,"Cannot open input file %s\n\n\n",global_mesh_name);
	exit(-1);
      }
      fprintf(stdout,"mstk mesh %s read in successfully\n",global_mesh_name);

      *dim = MESH_Num_Regions(mesh) ? 3 : 2;
    }

    int del_inmesh = 1;
    MSTK_Mesh_Distribute(mesh, recv_mesh, dim, ring, with_attr, method, 
			 del_inmesh, comm);

    return 1;
  }


  /* Partition a given mesh and distribute it to 'num' processors */
  /* The rank, num and comm arguments will go away soon           */

  int MSTK_Mesh_Distribute(Mesh_ptr globalmesh, Mesh_ptr *mymesh, int *dim, 
			   int ring, int with_attr, int method, 
			   int del_inmesh, MSTK_Comm comm) {
    int i, recv_dim;
    int *send_dim, *part;
    int rank, num;
    int DebugWait=0;

    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&num);
    recv_dim = rank+5;

    while (DebugWait) ;

#ifdef DEBUG
    double t0 = MPI_Wtime();
#endif

    send_dim = (int *) malloc(num*sizeof(int));
    for (i = 0; i < num; i++) send_dim[i] = *dim;

    MPI_Scatter(send_dim, 1, MPI_INT, &recv_dim, 1, MPI_INT, 0, comm);
    free(send_dim);

    if (rank != 0)
      *dim = recv_dim;
    MESH_Get_Partitioning(globalmesh, method, &part, comm);

#ifdef DEBUG
    double elapsed_time = MPI_Wtime() - t0;

    fprintf(stderr,"Elapsed time after calculating partitioning on processor %-d is %lf s\n",rank,elapsed_time);
#endif

    if(rank == 0) {    

      /* Partition the mesh */
      Mesh_ptr *submeshes = (Mesh_ptr *) MSTK_malloc((num)*sizeof(Mesh_ptr));
      MSTK_Mesh_Partition(globalmesh, num, part, ring, with_attr, 
			  del_inmesh, submeshes);

#ifdef DEBUG
      fprintf(stderr,"Finished partitioning\n");
#endif
#ifdef DEBUG
      elapsed_time = MPI_Wtime() - t0;

      fprintf(stderr,"Elapsed time after creating partitions on processor %d is %lf s\n",rank,elapsed_time);
#endif
    
      int numreq = 0;
      int maxreq = 25; /* should be 17*(num-1) to avoid realloc */
      MPI_Request *requests = (MPI_Request *) malloc(maxreq*sizeof(MPI_Request));
      int numptrs2free = 0;
      int maxptrs2free = 25; /* should be about 12*(num-1) to avoid realloc */
      void ** ptrs2free = (void **) malloc(maxptrs2free*sizeof(void *));

      for(i = 1; i < num; i++) {

        int torank = i;
	MSTK_SendMesh(submeshes[i],torank,with_attr,comm,&numreq,&maxreq,
                      &requests,&numptrs2free,&maxptrs2free,&ptrs2free);

	MESH_Delete(submeshes[i]);  

#ifdef DEBUG
	elapsed_time = MPI_Wtime() - t0;

	fprintf(stderr,"Elapsed time after sending out mesh #%-d on processor %-d is %lf s\n",i, rank, elapsed_time);
#endif

	
	int flag;

	MPI_Testall(numreq,requests,&flag,MPI_STATUSES_IGNORE);
	if (flag) {
	  /* all queued requests are completed - reset the requests array 
	     and free up memory allocated in the MESH_Send* routines */

	  numreq = 0;
	  
	  /* free all the memory allocated in the send subroutines */

	  int p;
	  for (p = 0; p < numptrs2free; ++p) free(ptrs2free[p]); 
	  numptrs2free = 0;
	}
	else {
	  /* check if we buffered too many requests or this is the
	   last loop iteration - if so, we wait until all the data is
	   sent out; if not, we continue. One can control how
	   frequently we do a blocking wait for the send requests by
	   adjusting maxpendreq */

	  int maxpendreq = 200;
	  if (numreq > maxpendreq || i == num-1) {
	    if (MPI_Waitall(numreq,requests,MPI_STATUSES_IGNORE) != MPI_SUCCESS)
	      MSTK_Report("MSTK_Mesh_Distribute","Could not send mesh",
			  MSTK_FATAL);
	    numreq = 0;

	    int p;
	    for (p = 0; p < numptrs2free; ++p) free(ptrs2free[p]);
	    numptrs2free = 0;
	  }
	    
	}
      }
      if (maxptrs2free) free(ptrs2free);
      if (maxreq) free(requests);

      if (*mymesh == NULL)
        *mymesh = submeshes[0];
      else {
        MESH_Copy(submeshes[0],*mymesh,1,1);
        MESH_Delete(submeshes[0]);
      }


#ifdef DEBUG
      fprintf(stderr,"Sent meshes to all partitions\n");
#endif

      MESH_Build_GhostLists(*mymesh,*dim);

      free(submeshes);

    }

    if( rank > 0) {

      if (!(*mymesh)) *mymesh = MESH_New(UNKNOWN_REP);
      int fromrank = 0;
      MSTK_RecvMesh(*mymesh,*dim,fromrank,with_attr,comm);

#ifdef DEBUG
      fprintf(stderr,"Received mesh on partition %d\n",rank);
#endif

#ifdef DEBUG
      elapsed_time = MPI_Wtime() - t0;
      fprintf(stderr,"Elapsed time after receiving mesh on partition %d is %lf s\n",rank,elapsed_time);
#endif
      
    }

    MESH_Set_Prtn(*mymesh,rank,num);

    MESH_Update_ParallelAdj(*mymesh,comm);

#ifdef DEBUG
    elapsed_time = MPI_Wtime() - t0;
    fprintf(stdout,"Elapsed time after MESH_Update_ParallelAdj on processor %-d is %lf s\n",rank,elapsed_time);
#endif

    MESH_Disable_GlobalIDSearch(*mymesh);


#ifdef DEBUG
    fprintf(stderr,"Updated parallel adjacencies. Exiting mesh distribution\n");
#endif
#ifdef DEBUG
    elapsed_time = MPI_Wtime() - t0;
    fprintf(stderr,"Elapsed time after mesh distribution on processor %-d is %lf s\n",rank,elapsed_time);
#endif


    if (part)
      free(part);

    /* Put a barrier so that distribution of meshes takes place one at a time 
       in a simulation that may have multiple mesh objects on each processor */

    MPI_Barrier(comm);

    return 1;
  }


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

    if (input_type > 1) 
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


  /* Update attributes on partitions */

  int MSTK_UpdateAttr(Mesh_ptr mesh, MSTK_Comm comm) {
    int i, natt;
    char attr_name[256];
    MAttrib_ptr attrib;

    natt = MESH_Num_Attribs(mesh);
    for(i = 0; i < natt; i++) {
      attrib = MESH_Attrib(mesh,i);
      MAttrib_Get_Name(attrib,attr_name);
      MPI_Barrier(comm);
      MESH_UpdateAttr(mesh,attr_name,comm);
    }

    return 1;
  }

#ifdef __cplusplus
}
#endif

