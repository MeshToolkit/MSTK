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
			  int with_attr, Mesh_ptr *submeshes) {
    int i, j, natt, nset, idx;
    char attr_name[256];
    MAttrib_ptr attrib, g2latt, l2gatt;
    MSet_ptr mset;
    MVertex_ptr mv;
    MEdge_ptr me;
    MFace_ptr mf;
    MRegion_ptr mr;


    /* Create an attribute to keep track of the connections from
       entities of the global mesh to entities of submeshes */

    g2latt = MAttrib_New(mesh,"Global2Local",POINTER,MALLTYPE);


    /* Split the mesh into 'num' submeshes */

    MESH_Partition(mesh, num, part, submeshes);

    for(i = 0; i < num; i++) {
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
	
	for (i = 0; i < num; i++)
	  MESH_CopySet(mesh,submeshes[i],mset);
      }
      
    }


    /* Don't waste time removing these attributes */
    /* */
    /* /\* Remove these attributes from submeshes *\/ */

    /* for (i = 0; i < num; i++) { */
      
    /*   l2gatt = MESH_AttribByName(submeshes[i],"Local2Global"); */
      
    /*   if (l2gatt) { */
	
    /*     idx = 0; */
    /*     while ((mv = MESH_Next_Vertex(submeshes[i],&idx))) */
    /*       MEnt_Rem_AttVal(mv,l2gatt); */
	
    /*     idx = 0; */
    /*     while ((me = MESH_Next_Edge(submeshes[i],&idx))) */
    /*       MEnt_Rem_AttVal(me,l2gatt); */
	
    /*     idx = 0; */
    /*     while ((mf = MESH_Next_Face(submeshes[i],&idx))) */
    /*       MEnt_Rem_AttVal(mf,l2gatt); */
	
    /*     idx = 0; */
    /*     while ((mr = MESH_Next_Region(submeshes[i],&idx))) */
    /*       MEnt_Rem_AttVal(mr,l2gatt); */
	
    /*     MAttrib_Delete(l2gatt); */
	
    /*   } /\* if (l2gatt) *\/ */


    /*   g2latt = MESH_AttribByName(submeshes[i],"Global2Local"); */
      
    /*   if (g2latt) { */
	
    /*     idx = 0; */
    /*     while ((mv = MESH_Next_Vertex(submeshes[i],&idx))) */
    /*       MEnt_Rem_AttVal(mv,g2latt); */
	
    /*     idx = 0; */
    /*     while ((me = MESH_Next_Edge(submeshes[i],&idx))) */
    /*       MEnt_Rem_AttVal(me,g2latt); */
	
    /*     idx = 0; */
    /*     while ((mf = MESH_Next_Face(submeshes[i],&idx))) */
    /*       MEnt_Rem_AttVal(mf,g2latt); */
	
    /*     idx = 0; */
    /*     while ((mr = MESH_Next_Region(submeshes[i],&idx))) */
    /*       MEnt_Rem_AttVal(mr,g2latt); */
	
    /*     MAttrib_Delete(g2latt); */
	
    /*   } /\* if (g2latt) *\/ */

    /* } /\* for (i = 0; i < num; i++) *\/ */


    /*
    printf("global mesh has been partitioned into %d parts\n",num);
    */
    return 1;
  }


  /* Send a mesh to a processor 'torank' with or without attributes */

  int MSTK_SendMesh(Mesh_ptr mesh, int torank, int with_attr, MSTK_Comm comm) {
    int i, natt, nset;
    char attr_name[256], mset_name[256];
    MAttrib_ptr attrib;
    MSet_ptr mset;

    /* Send mesh only */

    MESH_SendMesh(mesh,torank,comm);

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
      MESH_SendAttr(mesh,attr_name,torank,comm);
    }

    /* Distribute entity sets */

    nset = MESH_Num_MSets(mesh);
    for(i = 0; i < nset; i++) {
      mset = MESH_MSet(mesh,i);
      MSet_Name(mset,mset_name);
      MESH_SendMSet(mesh,mset_name,torank,comm);
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

    if(!with_attr)
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

    MSTK_Mesh_Distribute(mesh, recv_mesh, dim, ring, with_attr, method, comm);

    if (rank == 0)
      MESH_Delete(mesh);

    return 1;
  }


  /* Partition a given mesh and distribute it to 'num' processors */
  /* The rank, num and comm arguments will go away soon           */

  int MSTK_Mesh_Distribute(Mesh_ptr globalmesh, Mesh_ptr *mymesh, int *dim, 
			   int ring, int with_attr, int method, MSTK_Comm comm) {
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

      /* Partition the mesh*/
      Mesh_ptr *submeshes = (Mesh_ptr *) MSTK_malloc((num)*sizeof(Mesh_ptr));
      MSTK_Mesh_Partition(globalmesh, num, part, ring, with_attr, submeshes);

#ifdef DEBUG
      fprintf(stderr,"Finished partitioning\n");
#endif
#ifdef DEBUG
      elapsed_time = MPI_Wtime() - t0;

      fprintf(stderr,"Elapsed time after creating partitions on processor %d is %lf s\n",rank,elapsed_time);
#endif
      
      for(i = 1; i < num; i++) {

        int torank = i;
	MSTK_SendMesh(submeshes[i],torank,with_attr,comm);

	MESH_Delete(submeshes[i]);  
#ifdef DEBUG
      elapsed_time = MPI_Wtime() - t0;

      fprintf(stderr,"Elapsed time after sending out mesh #%-d on processor %-d is %lf s\n",i, rank, elapsed_time);
#endif
      }

      if (*mymesh == NULL)
        *mymesh = submeshes[0];
      else
        MESH_Copy(submeshes[0],*mymesh,1,1);


#ifdef DEBUG
      fprintf(stderr,"Sent meshes to all partitions\n");
#endif

      MESH_Build_GhostLists(*mymesh,*dim);
     
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

