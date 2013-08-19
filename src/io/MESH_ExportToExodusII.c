#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"
#include "exodusII.h"


#ifdef __cplusplus
extern "C" {
#endif

  /* this function collects element block inforamtion based on element type       */
  void MESH_Get_Element_Block_Info(Mesh_ptr mesh, int *num_element_blocks,\
				   List_ptr **element_blocks, 
				   int **element_block_ids, 
				   char ***element_block_types,
                                   MAttrib_ptr block_type_att);

  /* this function collects side set information                                  */
  int MESH_Num_Face_Set(List_ptr esides, int *side_set_id, int enable_set);
  int MESH_Num_Edge_Set(List_ptr esides, int *side_set_id, int enable_set);
  void MESH_Get_Side_Set_Info(Mesh_ptr mesh, int *num_side_set, 
			      List_ptr **side_sets, int **side_set_ids);
  void MESH_Get_Node_Set_Info(Mesh_ptr mesh, int *num_node_set, 
			      List_ptr **node_sets, int **node_set_ids);


  /* this function get the local side number */
  int MF_LocalID_in_Region(MFace_ptr mf, MRegion_ptr mr);
  int ME_LocalID_in_Face(MEdge_ptr me, MFace_ptr mf);
  void itoa(int n, char s[]);
  void reverse(char s[]);


  /*this defines the maximum number of side set per element block*/

  /* Function to export MSTK mesh to EXODUSII format               

     Export only owned not ghost elements in distributed meshes

     Global Parameter information:                                     
     num_node, num_elem, num_elem_blk, num_node_set, num_side_set  
     Quality assurance information(optional)                       
     The nodes(vertices) list and elements list  
     For nodes, there is only one attribute associate with them    
     opts[0] = 1 --- verbose (print stat information)
     opts[1] = 1 --- enable_set (sidesets, nodesets based on GEntID, GEntDim)
  */

#ifdef MSTK_HAVE_MPI
  int ownedmk;
#endif

  int MESH_ExportToExodusII(Mesh_ptr mesh, const char *filename, 
			    const int natt, const char **attnames, 
			    const int *opts, MSTK_Comm comm) {

    int enable_set, verbose;
    int i, j, k, idx, idx2;
    int nv, ne, nf, nr, nvall, neall, nfall, nrall;
    int num_element_block, num_side_set, num_node_set;
    int nvfblock, block_id, *connect, *nnpe;
    MVertex_ptr mv;
    MEdge_ptr me;
    MFace_ptr mf;
    MRegion_ptr mr;
    RepType reptype;
    char MESH_rtype_str[5][3] = {"F1\0","F4\0","R1\0","R2\0","R4\0"};
    int jv, vid;
    double xyz_coord[3];
    int boundary_dim;  /* For pure 2D mesh, boundary_dim is 1, for surface mesh, it is 2 */
    int *element_block_ids, *side_set_ids, *node_set_ids;
    int num_face_block;
    char **element_block_types, block_name[256];
    List_ptr *element_blocks, *side_sets, *node_sets, face_block;
    List_ptr fverts, rverts;
    MAttrib_ptr vidatt=NULL, fidatt=NULL;
    char modfilename[256], funcname[64]="MESH_ExportToExodusII";
    int status;
    double rval;
    void *pval;

    int mstk2exo_facemap[5][6]={{0,0,0,0,0,0},
				{4,1,2,3,0,0}, /* TET */
				{0,0,0,0,0,0}, /* PYRAMID, no support in exo*/
			        {4,5,1,2,3,0}, /* PRISM, must verify nums */
			        {5,6,1,2,3,4}};/* HEX */


    reptype = MESH_RepType(mesh);


#ifdef MSTK_HAVE_MPI
    int rank, numprocs;
    MPI_Comm_size(comm,&numprocs);
    MPI_Comm_rank(comm,&rank);

    if (numprocs > 1) {
      char basename[256];

      int ndigits = 0;
      int div = 1;
      while (numprocs/div) {div *= 10; ndigits++;}
 
      strcpy(basename,filename);
      char *ext = strstr(basename,".par");
      
      if (!ext) {
        ext = strstr(basename,".exo");
        if (!ext)
          MSTK_Report(funcname,
                      "Extension for parallel Exodus II output should be .exo or .par",
                      MSTK_FATAL);
        else {          
          ext[0] = '\0';
          sprintf(modfilename,"%s.exo.%d.%0*d",basename,numprocs,ndigits,rank);
        }
      }
      else {
        ext[0] = '\0';
        sprintf(modfilename,"%s.par.%d.%0*d",basename,numprocs,ndigits,rank);
      }

    }
    else
      strcpy(modfilename,filename);
#else
    strcpy(modfilename,filename);
#endif    



    verbose = opts ? opts[0] : 0;
    enable_set = opts ? opts[1] : 1;

    if (verbose)
      fprintf(stdout,"\nMesh representation type is %s\n", 
	      MESH_rtype_str[reptype]);
    nv = nvall = MESH_Num_Vertices(mesh);
    ne = neall = MESH_Num_Edges(mesh);
    nf = nfall = MESH_Num_Faces(mesh);
    nr = nrall = MESH_Num_Regions(mesh);
    
    if (nv == 0) {
      fprintf(stdout,"No vertices information \n");
      exit(-1);
    }


    int element_dim = -1, side_dim = -1;
    int num_element = 0, num_side = 0;

    vidatt = MAttrib_New(mesh,"vidatt",INT,MVERTEX);

    if (nrall) {
      element_dim = 3; side_dim = 2; boundary_dim = 2;

#ifdef MSTK_HAVE_MPI

      /* Mark only the owned regions */

      ownedmk = MSTK_GetMarker();

      nr = 0; nf = 0; ne = 0; nv = 0;
      idx = 0;
      while ((mr = MESH_Next_Region(mesh,&idx))) {
        if (MR_PType(mr) != PGHOST) {
          MEnt_Mark(mr,ownedmk);
          nr++;

          List_ptr rfaces = MR_Faces(mr);
          idx2 = 0;
          while ((mf = (MFace_ptr) List_Next_Entry(rfaces,&idx2))) {
            if (!MEnt_IsMarked(mf,ownedmk)) {
              MEnt_Mark(mf,ownedmk);
              nf++;

              List_ptr fedges = MF_Edges(mf,1,0);
              int idx3 = 0;
              while ((me = (MEdge_ptr) List_Next_Entry(fedges,&idx3))) {
                if (!MEnt_IsMarked(me,ownedmk)) {
                  MEnt_Mark(me,ownedmk);
                  ne++;

                  int kk;
                  for (kk = 0; kk < 2; kk++) {
                    mv = ME_Vertex(me,kk);
                    if (!MEnt_IsMarked(mv,ownedmk)) {
                      MEnt_Mark(mv,ownedmk);
                      nv++;
                      MEnt_Set_AttVal(mv,vidatt,nv,0.0,NULL);
                    }
                  }
                }
              } /* while (me...) */
              List_Delete(fedges);
            }
          } /* while (mf...) */
          List_Delete(rfaces);
        }
      }

#endif
        
      
      num_element = nr; num_side = nf;

      if(verbose)
	fprintf(stdout,"\nThis is a 3D volume mesh with %d Vertices, %d Edges, %d Faces and %d Regions\n",\
		nv,ne,nf,nr);

    }
    else if (nfall) {

      element_dim = 2; side_dim = 1; boundary_dim = 1;

#ifdef MSTK_HAVE_MPI

      /* Mark only the owned faces */

      ownedmk = MSTK_GetMarker();

      nr = 0; nf = 0; ne = 0; nv = 0;
      idx2 = 0;
      while ((mf = (MFace_ptr) MESH_Next_Face(mesh,&idx2))) {
        if (MF_PType(mf) != PGHOST) {
          MEnt_Mark(mf,ownedmk);
          nf++;
          
          List_ptr fedges = MF_Edges(mf,1,0);
          int idx3 = 0;
          while ((me = (MEdge_ptr) List_Next_Entry(fedges,&idx3))) {
            if (!MEnt_IsMarked(me,ownedmk)) {
              MEnt_Mark(me,ownedmk);
              ne++;
              
              int kk;
              for (kk = 0; kk < 2; kk++) {
                mv = ME_Vertex(me,kk);
                if (!MEnt_IsMarked(mv,ownedmk)) {
                  MEnt_Mark(mv,ownedmk);
                  nv++;
                  MEnt_Set_AttVal(mv,vidatt,nv,0.0,NULL);
                }
              }
            }
          } /* while (me...) */
          List_Delete(fedges);
        }
      } /* while (mf...) */

#endif

      num_element = nf; num_side = ne;

      if(verbose)
      fprintf(stdout,"\nThis is a 3D surface mesh with %d Vertices, %d Edges and %d Faces\n", nv,ne,nf);

    }
    else {
      MSTK_Report(funcname,"Writing of wire meshes or point clouds not supported",MSTK_ERROR);
      exit(-1);
    }
    
    /* decide if it is 2d planar mesh or 3D surface mesh */

    if (element_dim == 2) {
      double xyz_coord0[3];
      mv = MESH_Vertex(mesh,0);
      MV_Coords(mv,xyz_coord0);

      idx = 0;
      while ((mv = MESH_Next_Vertex(mesh,&idx))) {
	MV_Coords(mv,xyz_coord);
	/* if z coordinate is not all the same, it is a surface mesh */
	if ((xyz_coord[2]!=xyz_coord0[2])) {
	  boundary_dim = 2;
	  break;
	}
      }
      if (verbose) {
	if (boundary_dim == 1)
	  fprintf(stdout,"\nThis is a 2D planar mesh\n");
      }
    }



    /* create ExodusII file */

    int CPU_word_size=sizeof(double);
    int IO_word_size=sizeof(double);
    int exoid, err;
    exoid = ex_create(modfilename,EX_CLOBBER,&CPU_word_size,&IO_word_size);
    if (exoid < 0) {
      fprintf(stderr, "after ex_create, error = %d\n", exoid);
      exit(-1);
    }


    
    if(verbose) {
      if (enable_set)
	fprintf(stdout,"\nElement block and side set based on geometric id is enabled.\n");
      else
	fprintf(stdout,"\nElement block and side set based on geometric id is disabled.\n");
    }


    /* COLLECT ELEMENT BLOCK INFO */

    /* Get element block information(based on element type, and
       geometric entity id) 
       In 3D,
       element_block_id = geometric_entity_id<<16 | MR_Num_Vertices<<8 | MR_Num_Faces  
       In 2D,
       element_block_id = geometric_entity_id<<16 | MF_Num_Vertices

    */

    if(verbose) {
      if (enable_set)
	fprintf(stdout,"\nElement block and side set based on geometric id is enabled.\n");
      else
	fprintf(stdout,"\nElement block and side set based on geometric id is disabled.\n");
    }

    /* Even though there is a native element type for every mesh
       region or face we have to consider what the element type is
       considered as for Exodus II output. For example, a HEX element
       in a mixed mesh will be written as a NFACED (polyhedral)
       element */

    MAttrib_ptr block_type_att = 
      MAttrib_New(mesh,"block_type",INT,MALLTYPE); 

    MESH_Get_Element_Block_Info(mesh, &num_element_block, &element_blocks,
                                &element_block_ids, &element_block_types, 
                                block_type_att);



    int num_element_block_glob;
    int *element_block_ids_glob;
    char **element_block_types_glob;
    List_ptr *element_blocks_glob;

#ifdef MSTK_HAVE_MPI
    
    int maxnum=0, maxnum1=0;
    MPI_Allreduce(&num_element_block,&maxnum,1,MPI_INT,MPI_MAX,comm);

    int *ebids_array_loc = (int *) MSTK_calloc(maxnum,sizeof(int));

    for (i = 0; i < num_element_block; i++)
      ebids_array_loc[i] = element_block_ids[i];

    int *ebids_array_glob;
    ebids_array_glob = (int *) MSTK_calloc(maxnum*numprocs,sizeof(int));
    
    MPI_Gather(ebids_array_loc,maxnum,MPI_INT,ebids_array_glob,maxnum,MPI_INT,0,comm);

    if (rank == 0) {

      /* Make a unique list of block IDs */

      maxnum1 = num_element_block;

      for (i = 1; i < numprocs; i++) {

        int offset = maxnum*i;
        for (j = 0; j < maxnum; j++) {
          if (ebids_array_glob[offset+j] == 0) continue;

          /* Is this element block ID in the list of global element
             block IDs being collected at the head of the
             ebids_array_glob? */

          int found = 0;
          for (k = 0; k < maxnum1; k++) {
            if (ebids_array_glob[offset+j] == ebids_array_glob[k]) {
              found = 1;
              break;
            }
          }

          if (!found) {
            
            /* Found element block ID on processor i that is not in
               the global list. Put it in. */

            ebids_array_glob[maxnum1] = ebids_array_glob[offset+j];
            maxnum1++;

          }
        }
      } /* for (i = 1; ....) */

      num_element_block_glob = maxnum1;

    } /* if (rank == 0) */

    /* Tell everyone how many element blocks there really are */

    MPI_Bcast(&num_element_block_glob,1,MPI_INT,0,comm);

    /* Send everyone the IDs of these element blocks */

    element_block_ids_glob = (int *) MSTK_malloc(num_element_block_glob*sizeof(int));
    if (rank == 0)
      for (i = 0; i < num_element_block_glob; i++)
        element_block_ids_glob[i] = ebids_array_glob[i];

    MPI_Bcast(element_block_ids_glob,num_element_block_glob,MPI_INT,
              0,comm);

    /* Populate the global element block data on each processor */

    element_blocks_glob = (List_ptr *) MSTK_calloc(num_element_block_glob,
                                                   sizeof(List_ptr));
    element_block_types_glob = (char **) MSTK_calloc(num_element_block_glob,
                                                   sizeof(char *));
    for (i = 0; i < num_element_block_glob; i++)
      element_block_types_glob[i] = (char *) MSTK_malloc(16*sizeof(char));

    for (i = 0; i < num_element_block_glob; i++) {
      int found = 0;
      for (j = 0; j < num_element_block; j++) {
        if (element_block_ids[j] == element_block_ids_glob[i]) {
          element_blocks_glob[i] = element_blocks[j];
          strcpy(element_block_types_glob[i],element_block_types[j]);
          found = 1;
          break;
        }        
      }      
      if (!found) {
        element_blocks_glob[i] = List_New(1); /* dummy list */
        if (element_dim == 3)
          strcpy(element_block_types_glob[i],"HEX");
        else if (element_dim == 2)
          strcpy(element_block_types_glob[i],"QUAD");
        else
          strcpy(element_block_types_glob[i],"UNKNOWN");
      }
    }

    if (verbose) {
      fprintf(stdout,"\nElement block information:\n");
      fprintf(stdout,"total %d elements\n", num_element); 
      fprintf(stdout,"total %d element blocks:\n\n", num_element_block);
    }

    MSTK_free(ebids_array_loc);
    MSTK_free(ebids_array_glob);

#else

    num_element_block_glob = num_element_block;
    element_block_ids_glob = element_block_ids;
    element_blocks_glob = element_blocks;
    element_block_types_glob = element_block_types;

#endif

    /* Compute the ID of each element according to the Exodus file */

    int maxid=1, id;
    if (MESH_Num_Regions(mesh)) {
      idx = 0;
      while ((mr = MESH_Next_Region(mesh,&idx))) {
        id = MR_ID(mr);
        if (maxid < id) maxid = id;
      }
    }
    else if (MESH_Num_Faces(mesh)) {
      idx = 0;
      while ((mf = MESH_Next_Face(mesh,&idx))) {
        id = MF_ID(mf);
        if (maxid < id) maxid = id;
      }
    }
    int *elem_id = (int *) calloc(maxid,sizeof(int));

    int offset = 0;
    for (i = 0; i < num_element_block; i++) {
      int nelem = List_Num_Entries(element_blocks[i]);
      for (k = 0; k < nelem; k++) {
        MEntity_ptr ment = List_Entry(element_blocks[i],k);
        int entid = MEnt_ID(ment);
        elem_id[entid-1] = offset+k+1;
      }

      offset += nelem;
    }
    

    /* COLLECT FACE BLOCK FOR POLYHEDRAL MESHES */

    /* Put all faces of polyhedral elements in one block */

    face_block = List_New(10);
    num_face_block = 1;
    nvfblock = 0;
    int mkid1 = MSTK_GetMarker();
    fidatt = MAttrib_New(mesh,"local_face_ids",INT,MFACE);
    

    j = 0;
    for (i = 0; i < num_element_block; i++) {
      if (strncasecmp(element_block_types[i],"NFACED",6) == 0) {

	idx = 0;
	while ((mr = List_Next_Entry(element_blocks[i],&idx))) {
	  List_ptr rfaces = MR_Faces(mr);
	  MFace_ptr rf;

	  idx2 = 0;
	  while ((rf = List_Next_Entry(rfaces,&idx2))) {
	    if (!MEnt_IsMarked(rf,mkid1)) {
	      List_Add(face_block,rf);
	      MEnt_Mark(rf,mkid1);
	      MEnt_Set_AttVal(rf,fidatt,++j,0.0,NULL);
	      nvfblock += MF_Num_Vertices(rf);
	    }
	  }
	  List_Delete(rfaces);
	}

      }
    }

    List_Unmark(face_block,mkid1);
    MSTK_FreeMarker(mkid1);
   
    if (List_Num_Entries(face_block) == 0) {
      List_Delete(face_block);
      face_block = NULL;
      num_face_block = 0;
    }
    


    /* COLLECT SIDE SET INFO */

    MESH_Get_Side_Set_Info(mesh, &num_side_set, &side_sets, &side_set_ids);


    int num_side_set_glob;
    int *side_set_ids_glob;
    List_ptr *side_sets_glob;

#ifdef MSTK_HAVE_MPI
    
    maxnum=0, maxnum1=0;
    MPI_Allreduce(&num_side_set,&maxnum,1,MPI_INT,MPI_MAX,comm);

    int *ssids_array_loc = (int *) MSTK_calloc(maxnum,sizeof(int));

    for (i = 0; i < num_side_set; i++)
      ssids_array_loc[i] = side_set_ids[i];

    int *ssids_array_glob = NULL;
    ssids_array_glob = (int *) MSTK_calloc(maxnum*numprocs,sizeof(int));
    
    MPI_Gather(ssids_array_loc,maxnum,MPI_INT,ssids_array_glob,maxnum,MPI_INT,0,comm);

    if (rank == 0) {

      /* Make a unique list of block IDs */

      maxnum1 = num_side_set;

      for (i = 1; i < numprocs; i++) {

        int offset = maxnum*i;
        for (j = 0; j < maxnum; j++) {
          if (ssids_array_glob[offset+j] == 0) continue;

          /* Is this element block ID in the list of global sideset
             IDs being collected at the head of the
             ssids_array_glob? */

          int found = 0;
          for (k = 0; k < maxnum1; k++) {
            if (ssids_array_glob[offset+j] == ssids_array_glob[k]) {
              found = 1;
              break;
            }
          }

          if (!found) {
            
            /* Found sideset ID on processor i that is not in
               the global list. Put it in. */

            ssids_array_glob[maxnum1] = ssids_array_glob[offset+j];
            maxnum1++;

          }
        }
      } /* for (i = 1; ....) */

      num_side_set_glob = maxnum1;
      
    } /* if (rank == 0) */

    /* Tell everyone how many sidesets there really are */

    MPI_Bcast(&num_side_set_glob,1,MPI_INT,0,comm);

    /* Send everyone the IDs of these sidesets */

    side_set_ids_glob = (int *) MSTK_malloc(num_side_set_glob*sizeof(int));
    if (rank == 0)
      for (i = 0; i < num_side_set_glob; i++)
        side_set_ids_glob[i] = ssids_array_glob[i];

    MPI_Bcast(side_set_ids_glob,num_side_set_glob,MPI_INT,0,comm);

    /* Populate the global sideset data on each processor */

    side_sets_glob = (List_ptr *) MSTK_calloc(num_side_set_glob,
                                                   sizeof(List_ptr));

    for (i = 0; i < num_side_set_glob; i++) {
      int found = 0;
      for (j = 0; j < num_side_set; j++) {
        if (side_set_ids[j] == side_set_ids_glob[i]) {
          side_sets_glob[i] = side_sets[j];
          found = 1;
          break;
        }
      }      
      if (!found)
        side_sets_glob[i] = List_New(1); /* dummy list */
    }

    if (verbose) {
      fprintf(stdout,"\nSide set information:\n");
      fprintf(stdout,"total %d side sets:\n\n", num_side_set);
    }

    MSTK_free(ssids_array_loc);
    MSTK_free(ssids_array_glob);

#else

    num_side_set_glob = num_side_set;
    side_set_ids_glob = side_set_ids;
    side_sets_glob = side_sets;

#endif





    /* COLLECT NODE SET INFO */

    MESH_Get_Node_Set_Info(mesh, &num_node_set, &node_sets, &node_set_ids);


    int num_node_set_glob;
    int *node_set_ids_glob;
    char **node_set_types_glob;
    List_ptr *node_sets_glob;

#ifdef MSTK_HAVE_MPI
    
    maxnum=0, maxnum1=0;
    MPI_Allreduce(&num_node_set,&maxnum,1,MPI_INT,MPI_MAX,comm);

    int *nsids_array_loc = (int *) MSTK_calloc(maxnum,sizeof(int));

    for (i = 0; i < num_node_set; i++)
      nsids_array_loc[i] = node_set_ids[i];

    int *nsids_array_glob = NULL;
    nsids_array_glob = (int *) MSTK_calloc(maxnum*numprocs,sizeof(int));
    
    MPI_Gather(nsids_array_loc,maxnum,MPI_INT,nsids_array_glob,maxnum,MPI_INT,0,comm);

    if (rank == 0) {

      /* Make a unique list of block IDs */

      maxnum1 = num_node_set;

      for (i = 1; i < numprocs; i++) {

        int offset = maxnum*i;
        for (j = 0; j < maxnum; j++) {
          if (nsids_array_glob[offset+j] == 0) continue;

          /* Is this element block ID in the list of global nodeset
             IDs being collected at the head of the
             nsids_array_glob? */

          int found = 0;
          for (k = 0; k < maxnum1; k++) {
            if (nsids_array_glob[offset+j] == nsids_array_glob[k]) {
              found = 1;
              break;
            }
          }

          if (!found) {
            
            /* Found nodeset ID on processor i that is not in
               the global list. Put it in. */

            nsids_array_glob[maxnum1] = nsids_array_glob[offset+j];
            maxnum1++;

          }
        }
      } /* for (i = 1; ....) */

      num_node_set_glob = maxnum1;

    } /* if (rank == 0) */

    /* Tell everyone how many nodesets there really are */

    MPI_Bcast(&num_node_set_glob,1,MPI_INT,0,comm);

    /* Send everyone the IDs of these nodesets */

    node_set_ids_glob = (int *) MSTK_malloc(num_node_set_glob*sizeof(int));
    if (rank == 0)
      for (i = 0; i < num_node_set_glob; i++)
        node_set_ids_glob[i] = nsids_array_glob[i];
  

    MPI_Bcast(node_set_ids_glob,num_node_set_glob,MPI_INT,
              0,comm);

    /* Populate the global nodeset data on each processor */

    node_sets_glob = (List_ptr *) MSTK_calloc(num_node_set_glob,
                                                   sizeof(List_ptr));

    for (i = 0; i < num_node_set_glob; i++) {
      int found = 0;
      for (j = 0; j < num_node_set; j++) {
        if (node_set_ids[j] == node_set_ids_glob[i]) {
          node_sets_glob[i] = node_sets[j];
          found = 1;
          break;
        }
      }      
      if (!found)
        node_sets_glob[i] = List_New(1); /* dummy list */
    }

    if (verbose) {
      fprintf(stdout,"\nNode set information:\n");
      fprintf(stdout,"total %d node sets:\n\n", num_node_set);
    }

    MSTK_free(nsids_array_loc);
    MSTK_free(nsids_array_glob);

#else

    num_node_set_glob = num_node_set;
    node_set_ids_glob = node_set_ids;
    node_sets_glob = node_sets;

#endif






    /* write global parameters */
    
    ex_init_params par;
    strcpy(par.title,filename);
    par.num_nodes = nv;

    if (element_dim == 3 || boundary_dim == 2)
      par.num_dim = 3;
    else 
      par.num_dim = 2;

    if (element_dim == 3) {
      par.num_edge = 0;
      par.num_face = face_block ? List_Num_Entries(face_block) : 0;
      par.num_face_blk = num_face_block;
      par.num_edge_blk = 0;
      par.num_elem = nr;
      par.num_elem_blk = num_element_block_glob;  
      par.num_node_sets = num_node_set_glob;
      par.num_edge_sets = 0;
      par.num_face_sets = 0;
      par.num_side_sets = num_side_set_glob;
      par.num_elem_sets = 0;

#ifdef MSTK_HAVE_MPI
      if (numprocs > 1) {
        par.num_node_maps = 1;
        par.num_edge_maps = 0;
        par.num_face_maps = 0;
        par.num_elem_maps = 1;
      }
      else {
        par.num_node_maps = 0;
        par.num_edge_maps = 0;
        par.num_face_maps = 0;
        par.num_elem_maps = 0;
      }
#else
      par.num_node_maps = 0;
      par.num_edge_maps = 0;
      par.num_face_maps = 0;
      par.num_elem_maps = 0;
#endif
    }
    else {
      par.num_edge = 0;
      par.num_face = 0;
      par.num_face_blk = 0;
      par.num_edge_blk = 0;
      par.num_elem = nf;
      par.num_elem_blk = num_element_block_glob;  
      par.num_node_sets = num_node_set_glob;
      par.num_edge_sets = 0;
      par.num_face_sets = 0;
      par.num_side_sets = num_side_set_glob;
      par.num_elem_sets = 0;

#ifdef MSTK_HAVE_MPI
      if (numprocs > 1) {
        par.num_node_maps = 1;
        par.num_edge_maps = 0;
        par.num_face_maps = 0;
        par.num_elem_maps = 1;
      }
      else {
        par.num_node_maps = 0;
        par.num_edge_maps = 0;
        par.num_face_maps = 0;
        par.num_elem_maps = 0;
      }
#else
      par.num_node_maps = 0;
      par.num_edge_maps = 0;
      par.num_face_maps = 0;
      par.num_elem_maps = 0;
#endif
    }

    ex_put_init_ext(exoid,&par);
  
    if (verbose)
      fprintf(stdout,"Global parameters wrote into %s\n", filename);


    /* write coordinate values */

    double *xcoord, *ycoord, *zcoord, vxyz[3];
    xcoord = (double *) MSTK_malloc(nv*sizeof(double));
    ycoord = (double *) MSTK_malloc(nv*sizeof(double));
    zcoord = (double *) MSTK_malloc(nv*sizeof(double));

    idx = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx))) {

#ifdef MSTK_HAVE_MPI 
      if (!MEnt_IsMarked(mv,ownedmk)) continue;
      MEnt_Get_AttVal(mv,vidatt,&vid,&rval,&pval);
#else
      vid = MV_ID(mv);
#endif
      MV_Coords(mv,vxyz);
      xcoord[vid-1] = vxyz[0];
      ycoord[vid-1] = vxyz[1];
      zcoord[vid-1] = vxyz[2];
    }
    err = ex_put_coord (exoid, xcoord, ycoord, zcoord);

    MSTK_free(xcoord);
    MSTK_free(ycoord);
    MSTK_free(zcoord);
    

    /* add coordinate names(optional) */

    char *coord_names[3] = {"xcoor","ycoor","zcoor"};
    err = ex_put_coord_names (exoid, coord_names);
    if (verbose)
      fprintf(stdout,"Coordinate values wrote into %s\n", filename);
    


    /* Write out face block if any */

    if (face_block) {
      strcpy(block_name,"face_block");
      block_id = 9999999;

      ex_put_block(exoid, EX_FACE_BLOCK, block_id, "nsided",
		   List_Num_Entries(face_block),nvfblock,
		   0,0,0);
      ex_put_name(exoid, EX_FACE_BLOCK, block_id, block_name);


      connect = (int *) calloc(nvfblock,sizeof(int));
      nnpe = (int *) calloc(List_Num_Entries(face_block),sizeof(int));

      idx = 0; i = 0; k = 0;
      while ((mf = List_Next_Entry(face_block,&idx))) {
	List_ptr fverts = MF_Vertices(mf,1,0);
	int idx2 = 0;
	MVertex_ptr fv;

	while ((fv = List_Next_Entry(fverts,&idx2)))
	  connect[k++] = MV_ID(fv);

	nnpe[i++] = List_Num_Entries(fverts);

	List_Delete(fverts);
      }

      /* Output number of nodes per face in face block */
      
      ex_put_entity_count_per_polyhedra(exoid, EX_FACE_BLOCK, block_id, nnpe);

      /* Output connectivity of faces in face block */

      ex_put_conn(exoid, EX_FACE_BLOCK, block_id, connect, NULL, NULL);

      free(connect);
      free(nnpe);
    }

    

    /* write element blocks, side sets and node sets information */

    if (verbose) {
      fprintf(stdout,"\nSide sets and node sets information:\n");
      fprintf(stdout,"total %d side sets:\n",num_side_set);
    }
    
    for (i = 0; i < num_element_block_glob; i++) { 

      char block_name[256];

      block_id = element_block_ids_glob[i];

      if (strncasecmp(element_block_types_glob[i],"NFACED",6) == 0) {
	
	int nfrblock = 0;

	nnpe = (int *) MSTK_calloc(List_Num_Entries(element_blocks_glob[i]),
				   sizeof(int));	
	idx = 0; j = 0;
	while ((mr = List_Next_Entry(element_blocks_glob[i],&idx))) {	  
	  nnpe[j] = MR_Num_Faces(mr);
	  nfrblock += nnpe[j];
	  j++;
	}

	connect = (int *) MSTK_calloc(nfrblock,sizeof(int));
	idx = 0; j = 0;
	while ((mr = List_Next_Entry(element_blocks_glob[i],&idx))) {
	  List_ptr rfaces = MR_Faces(mr);
	  MFace_ptr rf;
	  idx2 = 0;
	  while ((rf = List_Next_Entry(rfaces,&idx2))) {
	    int fid;
	    MEnt_Get_AttVal(rf,fidatt,&fid,NULL,NULL);
	    connect[j++] = fid;
	  }
	  List_Delete(rfaces);
	}


	sprintf(block_name,"POLYHEDRA_BLOCK_%-d",block_id);
	
	ex_put_block(exoid, EX_ELEM_BLOCK, block_id, "NFACED",
		     List_Num_Entries(element_blocks_glob[i]),
		     0, /* nodes */
		     0, /* edges */
		     nfrblock, /* faces */
		     0); /* attribute count */

	ex_put_name(exoid, EX_ELEM_BLOCK, block_id, block_name);

	ex_put_entity_count_per_polyhedra(exoid, EX_ELEM_BLOCK, block_id, nnpe);

	ex_put_conn(exoid, EX_ELEM_BLOCK, block_id, NULL, NULL, connect);

	MSTK_free(nnpe);

	MSTK_free(connect);

      }
      else if (strncasecmp(element_block_types_glob[i],"NSIDED",6) == 0) {

	int nvfblock = 0;

	nnpe = (int *) MSTK_calloc(List_Num_Entries(element_blocks_glob[i]),
				   sizeof(int));	
	idx = 0; j = 0;
	while ((mf = List_Next_Entry(element_blocks_glob[i],&idx))) {	  
	  nnpe[j] = MF_Num_Vertices(mf);
	  nvfblock += nnpe[j];
	  j++;
	}

	connect = (int *) MSTK_calloc(nvfblock,sizeof(int));
	idx = 0; j = 0;
	while ((mf = List_Next_Entry(element_blocks_glob[i],&idx))) {
	  List_ptr fverts = MF_Vertices(mf,1,0);
	  MVertex_ptr fv;
	  idx2 = 0;
	  while ((fv = List_Next_Entry(fverts,&idx2))) {
#ifdef MSTK_HAVE_MPI      
            MEnt_Get_AttVal(fv,vidatt,&vid,&rval,&pval);
#else
            vid = MV_ID(fv);
#endif
	    connect[j++] = vid;
          }
	  List_Delete(fverts);
	}


	sprintf(block_name,"POLYGON_BLOCK_%-d",block_id);
	
	ex_put_block(exoid, EX_ELEM_BLOCK, block_id, "NSIDED",
		     List_Num_Entries(element_blocks_glob[i]),
		     nvfblock, /* nodes */
		     0, /* edges */
		     0, /* faces */
		     0); /* attribute count */

	ex_put_name(exoid, EX_ELEM_BLOCK, block_id, block_name);

	ex_put_entity_count_per_polyhedra(exoid, EX_ELEM_BLOCK, block_id, nnpe);

	ex_put_conn(exoid, EX_ELEM_BLOCK, block_id, connect, NULL, NULL);

	MSTK_free(nnpe);

	MSTK_free(connect);


      }
      else {

	int nelem, nelnodes=0;

	sprintf(block_name,"BLOCK_%-d",block_id);
	nelem = List_Num_Entries(element_blocks_glob[i]);
	
	if (MESH_Num_Regions(mesh)) {
	  if (strncasecmp(element_block_types_glob[i],"TETRA",5) == 0) {
	    nelnodes = 4;
	    ex_put_elem_block(exoid, block_id, "TETRA", nelem, nelnodes, 1);
	  }
	  else if (strncasecmp(element_block_types_glob[i],"WEDGE",5) == 0) {
	    nelnodes = 6;
	    ex_put_elem_block(exoid, block_id, "WEDGE", nelem, nelnodes, 1);
	  }
	  else if (strncasecmp(element_block_types_glob[i],"HEX",3) == 0) {
	    nelnodes = 8;
	    ex_put_elem_block(exoid, block_id, "HEX", nelem, nelnodes, 1);
	  }
	  else
	    MSTK_Report("MESH_ExportToEXODUSII",
			"Element type unsupported by EXODUS II format",MSTK_FATAL);
	  
	  
	  connect = (int *) MSTK_malloc(nelnodes*nelem*sizeof(int));
	  
	  idx = 0; k = 0;
	  while ((mr = List_Next_Entry(element_blocks_glob[i],&idx))) {
	    rverts = MR_Vertices(mr);
	    for (j = 0; j < nelnodes; j++) {
#ifdef MSTK_HAVE_MPI      
              MEnt_Get_AttVal(List_Entry(rverts,j),vidatt,&vid,&rval,&pval);
#else
              vid = MV_ID(List_Entry(rverts,j));
#endif
	      connect[nelnodes*k+j] = vid;
            }
	    List_Delete(rverts);	    
	    k++;
	  }
	  
	  ex_put_elem_conn(exoid, block_id, connect);

	  MSTK_free(connect);
	  
	}
	else if (MESH_Num_Faces(mesh)) {
	  if (strncasecmp(element_block_types_glob[i],"TRIANGLE",8) == 0) {
	    nelnodes = 3;
	    ex_put_elem_block(exoid, block_id, "TRIANGLE", nelem, nelnodes, 1);
	  }
	  else if (strncasecmp(element_block_types_glob[i],"QUAD",4) == 0) {
	    nelnodes = 4;
	    ex_put_elem_block(exoid, block_id, "QUAD", nelem, nelnodes, 1);
	  }
	  else
	    MSTK_Report("MESH_ExportToEXODUSII",
			"Element type unsupported by EXODUS II format",MSTK_FATAL);

	  connect = (int *) MSTK_malloc(nelnodes*nelem*sizeof(int));

	  idx = 0; k = 0;
	  while ((mf = List_Next_Entry(element_blocks_glob[i],&idx))) {
	    fverts = MF_Vertices(mf,1,0);
	    for (j = 0; j < nelnodes; j++) {
#ifdef MSTK_HAVE_MPI      
              MEnt_Get_AttVal(List_Entry(fverts,j),vidatt,&vid,&rval,&pval);
#else
              vid = MV_ID(List_Entry(fverts,j);
#endif
	      connect[nelnodes*k+j] = vid;
            }
	    List_Delete(fverts);
	    k++;
	  }

	  ex_put_elem_conn(exoid, block_id, connect);

	  MSTK_free(connect);
	}

      }

    }

    if (verbose)
      fprintf(stdout,"Element block information written into %s\n", filename);



    /* Write out the node set information */

   
    for (i = 0; i < num_node_set_glob; i++) {
      int nnodes = List_Num_Entries(node_sets_glob[i]);
      ex_put_node_set_param(exoid, node_set_ids_glob[i], nnodes, 0);

      int *node_list = (int *) MSTK_malloc(nnodes*sizeof(int));

      idx = 0; j = 0;
      while ((mv = List_Next_Entry(node_sets_glob[i],&idx))) {
#ifdef MSTK_HAVE_MPI      
        MEnt_Get_AttVal(mv,vidatt,&vid,&rval,&pval);
#else
        vid = MV_ID(mv);
#endif
	node_list[j++] = vid;
      }

      ex_put_node_set(exoid, node_set_ids_glob[i], node_list);

      MSTK_free(node_list);
    }

    if (verbose)
      fprintf(stdout,"Node sets information written into %s\n", filename);



    
    /* Write out the side set information */


    for (i = 0; i < num_side_set_glob; i++) {

      int nsides = List_Num_Entries(side_sets_glob[i]);
      ex_put_side_set_param(exoid, side_set_ids_glob[i], nsides, 0);

      int *elem_list = (int *) MSTK_malloc(nsides*sizeof(int));
      int *side_list = (int *) MSTK_malloc(nsides*sizeof(int));

      if (nr) {
	idx = 0; j = 0;
	while ((mf = List_Next_Entry(side_sets_glob[i],&idx))) {

	  List_ptr fregs = MF_Regions(mf);

	  if (!fregs || List_Num_Entries(fregs) == 0)
	    MSTK_Report("MESH_ExportToEXODUSII",
			"Standalone face with no regions in side set",
			MSTK_FATAL);
	  mr = List_Entry(fregs,0);

#ifdef MSTK_HAVE_MPI
          if (MR_PType(mr) == PGHOST)
            mr = List_Entry(fregs,1);  /* at least 1 region should be owned */
#endif

          elem_list[j] = elem_id[MR_ID(mr)-1];
          if (elem_list[j] == 0)
            MSTK_Report("MESH_ExportToExodusII","Invalid Exodus II element ID",MSTK_ERROR);

	  
          /* Rather than look at the real element type of the region,
             we need to see what it will be classified as in the
             Exodus II file */

          MRType mrtype;

          MEnt_Get_AttVal(mr,block_type_att,&mrtype,NULL,NULL);

          int lid;
          if (mrtype == TET || mrtype == PRISM || mrtype == HEX) {
            lid = MF_LocalID_in_Region(mf,mr);
            side_list[j] = mstk2exo_facemap[mrtype][lid];
          }
          else { 
            List_ptr rfaces = MR_Faces(mr);
            int k, nrf = List_Num_Entries(rfaces);
            for (k = 0, lid = -1; k < nrf && lid < 0; k++)
              if (List_Entry(rfaces,k) == mf) lid = k;
            List_Delete(rfaces);
            if (lid == -1)
              MSTK_Report("MESH_ExportToExodusII",
              "Messed up mesh. Face not found in region connected to face",
                          MSTK_FATAL);
              
            side_list[j] = lid+1;
          }

	  List_Delete(fregs);
	  j++;
	}

	ex_put_side_set(exoid, side_set_ids_glob[i], elem_list, side_list);

      }
      else {

	idx = 0; j = 0;
	while ((me = List_Next_Entry(side_sets_glob[i],&idx))) {
	  List_ptr efaces = ME_Faces(me);

	  if (!efaces || List_Num_Entries(efaces) == 0)
	    MSTK_Report("MESH_ExportToEXODUSII",
			"Standalone edge with no faces in side set",
			MSTK_FATAL);
	  mf = List_Entry(efaces,0);

#ifdef MSTK_HAVE_MPI
          if (MF_PType(mf) == PGHOST)
            mf = List_Entry(efaces,1);  /* at least 1 face should be owned */
#endif

	  /* Since elements will be read back from the Exodus II file
	     element block by element block, we have to jump through
	     hoops to determine what the element number will be when
	     it is read */

          elem_list[j] = elem_id[MF_ID(mf)-1];
          if (elem_list[j] == 0)
            MSTK_Report("MESH_ExportToExodusII","Invalid Exodus II element ID",MSTK_ERROR);
	  
          int lid = ME_LocalID_in_Face(me,mf);
          side_list[j] = lid+1;

	  List_Delete(efaces);
	  j++;
	}

	ex_put_side_set(exoid, side_set_ids_glob[i], elem_list, side_list);

      }

      MSTK_free(elem_list);
      MSTK_free(side_list);

    }

    if (verbose)
      fprintf(stdout,"Side set information written into %s\n", filename);
    
    if (verbose)
      fprintf(stdout,"Node sets information written into %s\n", filename);


#ifdef MSTK_HAVE_MPI

    if (numprocs > 1) {

      /* Write out node map (global IDs of nodes) */

      int *node_map = (int *) malloc(nv*sizeof(int));
      
      idx = 0; i = 0;
      while ((mv = MESH_Next_Vertex(mesh,&idx)))
        if (MEnt_IsMarked(mv,ownedmk))
          node_map[i++] = MV_GlobalID(mv);

      status = ex_put_node_num_map(exoid, node_map);
      if (status < 0)
        MSTK_Report(funcname,"Error while writing node map in Exodus II file",
                    MSTK_FATAL);

      free(node_map);


      /* Write out element map (global IDs of elements) */

      int *elem_map;

      if (nr) {
        elem_map = (int *) malloc(nr*sizeof(int));
        idx = 0; i = 0;
        while ((mr = MESH_Next_Region(mesh,&idx)))
          if (MEnt_IsMarked(mr,ownedmk))
            elem_map[i++] = MR_GlobalID(mr);
      }
      else { // assume surface mesh
        elem_map = (int *) malloc(nf*sizeof(int));
        idx = 0; i = 0;
        while ((mf = MESH_Next_Face(mesh,&idx)))
          if (MEnt_IsMarked(mf,ownedmk))
            elem_map[i++] = MF_GlobalID(mf);
      }

      status = ex_put_elem_num_map(exoid, elem_map);
      if (status < 0)
        MSTK_Report(funcname,"Error while writing element map in Exodus II file",
                    MSTK_FATAL);

      free(elem_map);
    }

#endif /* MSTK_HAVE_MPI */

     
    /* write quality assurance information optional*/

    char *qa_record[1][4], date_str[256], time_str[256];
    qa_record[0][0] = "MSTK to ExodusII converter";
    qa_record[0][1] = "D Wang, R Garimella";
    time_t ctime;
    ctime = time(&ctime);
    strftime(date_str,sizeof(date_str),"%m/%d/%Y",localtime(&ctime));
    strftime(time_str,sizeof(time_str),"%H:%M:%S",localtime(&ctime));
    qa_record[0][2]= date_str;
    qa_record[0][3]= time_str;
    err = ex_put_qa(exoid, 1, qa_record);
    if (verbose)
      fprintf(stdout,"Quality assurance information wrote into %s\n", filename);

    ex_close(exoid);


    MSTK_free(element_block_ids_glob);
    for (i = 0; i < num_element_block_glob; i++)
      List_Delete(element_blocks_glob[i]);
    MSTK_free(element_blocks_glob);
    for (i = 0; i < num_element_block_glob; i++)
      MSTK_free(element_block_types_glob[i]);
    MSTK_free(element_block_types_glob);
    MSTK_free(elem_id);


    if (num_side_set_glob) {
      MSTK_free(side_set_ids_glob);
      for (i = 0; i < num_side_set_glob; i++)
        List_Delete(side_sets_glob[i]);
      MSTK_free(side_sets_glob);
    }
    
    if (num_node_set_glob) {
      MSTK_free(node_set_ids_glob);
      for (i = 0; i < num_node_set_glob; i++)
        List_Delete(node_sets_glob[i]);
      MSTK_free(node_sets_glob);
    }

#ifdef MSTK_HAVE_MPI
    idx = 0;
    while ((mr = MESH_Next_Region(mesh,&idx)))
      MEnt_Unmark(mr,ownedmk);
    idx = 0;
    while ((mf = MESH_Next_Face(mesh,&idx)))
      MEnt_Unmark(mf,ownedmk);
    idx = 0;
    while ((me = MESH_Next_Edge(mesh,&idx)))
      MEnt_Unmark(me,ownedmk);
    idx = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx)))
      MEnt_Unmark(mv,ownedmk);

    MSTK_FreeMarker(ownedmk);


    /* In parallel case, free the local lists also */
    /* No need to delete the lists because they are not duplicated */

    MSTK_free(element_block_ids);
    MSTK_free(element_blocks);
    for (i = 0; i < num_element_block; i++)
      MSTK_free(element_block_types[i]);
    MSTK_free(element_block_types);


    MSTK_free(side_set_ids);
    MSTK_free(side_sets);
      
    
    MSTK_free(node_set_ids);
    MSTK_free(node_sets);    

#endif


    idx = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx)))
      MEnt_Rem_AttVal(mv,vidatt);
    MAttrib_Delete(vidatt);
    if (fidatt) {
      idx = 0;
      while ((mf = MESH_Next_Face(mesh,&idx)))
        MEnt_Rem_AttVal(mf,fidatt);
      MAttrib_Delete(fidatt);
    }
      
    return 1;

  }



  /* Categorize the elements into blocks based on their geometric
     entity ID and the element type */

  void MESH_Get_Element_Block_Info(Mesh_ptr mesh, int *num_element_block, 
                                   List_ptr **element_blocks, 
                                   int **element_block_ids, 
                                   char ***element_block_types,
                                   MAttrib_ptr block_type_att) {
    MEdge_ptr me;
    MFace_ptr mf;
    MRegion_ptr mr;
    int i, nb, bid, found, nballoc;

    nb = 0; nballoc=10;
    *element_block_ids = (int *) calloc(nballoc,sizeof(int));
    *element_block_types = (char **) malloc(nballoc*sizeof(char *));
    *element_blocks = (List_ptr *) calloc(nballoc,sizeof(List_ptr));

    int idx = 0;
    
    if (MESH_Num_Regions(mesh)) {

      while ((mr = MESH_Next_Region(mesh, &idx))) {

#ifdef MSTK_HAVE_MPI
        if (!MEnt_IsMarked(mr,ownedmk)) continue;  /* Ghost element */
#endif

	List_ptr rverts = MR_Vertices(mr);
	int nrv = List_Num_Entries(rverts);
        List_Delete(rverts);
	int nrf = MR_Num_Faces(mr);
        MRType mrtype = MR_ElementType(mr);

        /*
	if (mrtype == TET || mrtype == PRISM || mrtype == HEX)
	  bid = (MEnt_GEntID(mr)<<16) | nrv<<8 | nrf;
	else
	  bid = MEnt_GEntID(mr)<<16;
        */

        bid = MEnt_GEntID(mr);
	
	found = 0; 
	i = 0;
	while (!found && i < nb) {
	  if (bid == (*element_block_ids)[i]) {           
	    found = 1;
	    List_Add((*element_blocks)[i],mr);    
            MEnt_Set_AttVal(mr,block_type_att,mrtype,0.0,NULL);
	    break;
	  }
	  i++;
	}

	if (!found) {
	  /* New element block */

	  if (nballoc == nb) {
	    nballoc *= 2;
	    *element_block_ids = (int *) realloc(*element_block_ids,nballoc*sizeof(int));
	    *element_blocks = (List_ptr *) realloc(*element_blocks,nballoc*sizeof(List_ptr));
	    *element_block_types = (char **) realloc(*element_block_types,nballoc*sizeof(char *));
	  }

	  (*element_block_ids)[nb] = bid;
	  (*element_blocks)[nb] = List_New(10);
	  (*element_block_types)[nb] = (char *) malloc(16*sizeof(char));

          if (mrtype == TET)
            strcpy((*element_block_types)[nb],"TETRA");
          else if (mrtype == PRISM)
            strcpy((*element_block_types)[nb],"WEDGE");
          else if (mrtype == HEX)
            strcpy((*element_block_types)[nb],"HEX");
          else
            strcpy((*element_block_types)[nb],"NFACED");

          MEnt_Set_AttVal(mr,block_type_att,mrtype,0.0,NULL);

	  List_Add((*element_blocks)[nb],mr);
	  nb++;
	}
      }
    }
    else if (MESH_Num_Faces(mesh)) {

      while ((mf = MESH_Next_Face(mesh,&idx))) {

#ifdef MSTK_HAVE_MPI
        if (!MEnt_IsMarked(mf,ownedmk)) continue; /* Ghost element */
#endif

	int nfv = MF_Num_Vertices(mf);

        /*
	if (nfv == 3 || nfv == 4)
	  bid = (MEnt_GEntID(mf)<<16) | nfv;
	else
	  bid = MEnt_GEntID(mf)<<16;
        */
        
        bid = MEnt_GEntID(mf);

	found = 0; i = 0;
	while (!found && i < nb) {
	  if (bid == (*element_block_ids)[i]) {
	    found = 1;
	    List_Add((*element_blocks)[i],mf);    
	    break;
	  }
	  i++;
	}

	if (!found) {
	  /* New element block */

	  if (nballoc == nb) {
	    nballoc *= 2;
	    *element_block_ids = (int *) realloc(*element_block_ids,nballoc*sizeof(int));
	    *element_blocks = (List_ptr *) realloc(*element_blocks,nballoc*sizeof(List_ptr));
	    *element_block_types = (char **) realloc(*element_block_types,nballoc*sizeof(char *));
	  }

	  (*element_block_ids)[nb] = bid;
	  (*element_blocks)[nb] = List_New(10);
	  (*element_block_types)[nb] = (char *) malloc(16*sizeof(char));
	  if (nfv > 4) {
	    strcpy((*element_block_types)[nb],"NSIDED");
            MEnt_Set_AttVal(mf,block_type_att,POLYGON,0.0,NULL);
          }
	  else {
	    if (nfv == 4) {
	      strcpy((*element_block_types)[nb],"QUAD");
              MEnt_Set_AttVal(mf,block_type_att,QUAD,0.0,NULL);
            }
	    else if (nfv == 3) {
	      strcpy((*element_block_types)[nb],"TRIANGLE");
              MEnt_Set_AttVal(mf,block_type_att,TRI,0.0,NULL);
            }
	  }

	  List_Add((*element_blocks)[nb],mf);
	  nb++;
	}
      }
    }

    *num_element_block = nb;

    /* Post process so that if the same element block has different
       types of elements, the block type will be labeled as NSIDED or
       NFACED */

    if (MESH_Num_Regions(mesh)) {
      for (i = 0; i < nb; i++) {
        MRType mrtype0, mrtype;
        MRegion_ptr mr;
        
        mr = List_Entry((*element_blocks)[i],0);
        mrtype0 = MR_ElementType(mr);
        
        int found = 0, idx = 0;
        while (!found && ((mr = List_Next_Entry((*element_blocks)[i],&idx)))) {
          if (MR_ElementType(mr) != mrtype0)
            found = 1;
        }
        
        if (found) { 
          /* found different element types in element block */
          /* label it as a block containing polyhedra       */
          
          strcpy((*element_block_types)[i],"NFACED");

          /* Attach attribute to each of its elements indicating that
             the element type according to Exodus II will be NSIDED */

          idx = 0;
          while ((mr = List_Next_Entry((*element_blocks)[i],&idx))) 
            MEnt_Set_AttVal(mr,block_type_att,POLYHED,0.0,NULL);
        }
      }
    }
    else if (MESH_Num_Faces(mesh)) {
      for (i = 0; i < nb; i++) {
        int nfv0;
        MFace_ptr mf;

        mf = List_Entry((*element_blocks)[i],0);
        nfv0 = MF_Num_Vertices(mf);

        int found = 0, idx = 0;
        while (!found && ((mf = List_Next_Entry((*element_blocks)[i],&idx)))) {
          if (MF_Num_Vertices(mf) != nfv0)
            found = 1;
        }

        if (found) {
          /* found different element types in element block */
          /* label it as a block containing polyhedra       */
          
          strcpy((*element_block_types)[i],"NSIDED");
          
          /* Attach attribute to each of its elements indicating that
             the element type according to Exodus II will be NSIDED */

          idx = 0;
          while ((mf = List_Next_Entry((*element_blocks)[i],&idx))) 
            MEnt_Set_AttVal(mf,block_type_att,POLYGON,0.0,NULL);
        }
      }
    }

  }


 
  void MESH_Get_Side_Set_Info(Mesh_ptr mesh, int *num_side_set, 
			      List_ptr **side_sets, int **side_set_ids) {
    MFace_ptr mf;
    MEdge_ptr me;
    int idx = 0, i = 0, found;
    int ns, nr, nf, ne, sid, nsalloc;


    ns = 0;
    nsalloc = 10;
    *side_sets = (List_ptr *) MSTK_malloc(nsalloc*sizeof(List_ptr));
    *side_set_ids = (int *) MSTK_malloc(nsalloc*sizeof(int));

    nr = MESH_Num_Regions(mesh);
    nf = MESH_Num_Faces(mesh);
    ne = MESH_Num_Edges(mesh);

    if (nr) {
      while ((mf = MESH_Next_Face(mesh,&idx))) {
	if (MF_GEntDim(mf) == 3) continue; /* Internal face */

#ifdef MSTK_HAVE_MPI
        if (!MEnt_IsMarked(mf,ownedmk)) continue; /* Face cnctd to only ghost elements */
#endif

	sid = MF_GEntID(mf);

	found = 0;
	i = 0;
	while (!found && i < ns) {
	  if ((*side_set_ids)[i] == sid) {
	    found = 1;
	    List_Add((*side_sets)[i],mf);
	    break;
	  }
	  i++;
	}

	if (!found) {
	  if (ns == nsalloc) {
	    nsalloc *= 2;
	    *side_sets = (List_ptr *) MSTK_realloc(*side_sets,nsalloc*sizeof(List_ptr));
	    *side_set_ids = (int *) MSTK_realloc(*side_set_ids,nsalloc*sizeof(int));
	  }
	  (*side_sets)[ns] = List_New(10);
	  List_Add((*side_sets)[ns],mf);
	  (*side_set_ids)[ns] = sid;
	  ns++;
	}
      }
    }
    else if (nf) {
      while ((me = MESH_Next_Edge(mesh,&idx))) {
	if (ME_GEntDim(me) != 1) continue;

#ifdef MSTK_HAVE_MPI
        if (!MEnt_IsMarked(me,ownedmk)) continue; /* Edge cnctd to only ghost elements */
#endif

	sid = ME_GEntID(me);

	found = 0;
	i = 0;
	while (!found && i < ns) {
	  if ((*side_set_ids)[i] == sid) {
	    found = 1;
	    List_Add((*side_sets)[i],me);
	    break;
	  }
	  i++;
	}

	if (!found) {
	  if (ns == nsalloc) {
	    nsalloc *= 2;
	    *side_sets = (List_ptr *) MSTK_realloc(*side_sets,nsalloc*sizeof(List_ptr));
	    *side_set_ids = (int *) MSTK_realloc(*side_set_ids,nsalloc*sizeof(int));
	  }
	  (*side_sets)[ns] = List_New(10);
	  List_Add((*side_sets)[ns],me);
	  (*side_set_ids)[ns] = sid;
	  ns++;
	}
      }
    }

    *num_side_set = ns;
  }
  

 
  void MESH_Get_Node_Set_Info(Mesh_ptr mesh, int *num_node_set, 
			      List_ptr **node_sets, int **node_set_ids) {
    MVertex_ptr mv;
    MEdge_ptr me, me2, ve;
    MFace_ptr mf, mf2, ef;
    int idx, idx2, idx3, idx4, i, found, nnalloc;
    int nn, nid, nid2, mkid1, mkid2;
    List_ptr fverts, fedges, vedges, efaces, elist, flist, vlist;


    nn = 0;
    nnalloc = 10;
    *node_sets = (List_ptr *) MSTK_malloc(nnalloc*sizeof(List_ptr));
    *node_set_ids = (int *) MSTK_malloc(nnalloc*sizeof(int));



    /* Node sets formed by vertices on geometric model vertices */

    idx = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx))) {
      if (MV_GEntDim(mv) != 0) continue;
      

#ifdef MSTK_HAVE_MPI
      if (!MEnt_IsMarked(mv,ownedmk)) continue; /* Vtx cnctd to only ghost elements */
#endif

      nid = MV_GEntID(mv);
      
      found = 0;
      i = 0;
      while (!found && i < nn) {
	if ((*node_set_ids)[i] == nid) {
	  found = 1;
	  List_Add((*node_sets)[i],mv);
	  break;
	}
	i++;
      }
      
      if (!found) {
	if (nn == nnalloc) {
	  nnalloc *= 2;
	  *node_sets = (List_ptr *) MSTK_realloc(*node_sets,nnalloc*sizeof(List_ptr));
	  *node_set_ids = (int *) MSTK_realloc(*node_set_ids,nnalloc*sizeof(int));
	}
	(*node_sets)[nn] = List_New(10);
	List_Add((*node_sets)[nn],mv);
	(*node_set_ids)[nn] = nid;
	nn++;
      }
    }




    /* Node sets formed by vertices on the closure of geometric model edges */

    mkid1 = MSTK_GetMarker();
    mkid2 = MSTK_GetMarker();

    idx = 0;
    while ((me = MESH_Next_Edge(mesh,&idx))) {
      if (ME_GEntDim(me) != 1) continue;

      if (MEnt_IsMarked(me,mkid1)) continue;

      /* Found a mesh edge on an unprocessed geometric model edge */	

      elist = List_New(10);

      MEnt_Mark(me,mkid1);
      MEnt_Mark(me,mkid2);
      List_Add(elist,me);

      int egid = ME_GEntID(me);
      nid = 10000+egid;

      /* Make a list of the mesh edges on this geometric model edge */
      
      idx2 = 0;
      while ((me2 = List_Next_Entry(elist,&idx2))) {
	for (i = 0; i < 2; i++) {
	  mv = ME_Vertex(me2,i);
	  vedges = MV_Edges(mv);
	  
	  idx3 = 0;
	  while ((ve = List_Next_Entry(vedges,&idx3))) {

	    if (ME_GEntDim(ve) != 1 || ME_GEntID(ve) != egid) continue; /* Not on the same model edge */

	    if (!MEnt_IsMarked(ve,mkid2)) {
	      List_Add(elist,ve);
	      MEnt_Mark(ve,mkid1);
	      MEnt_Mark(ve,mkid2);
	    }
	  }
	  List_Delete(vedges);
	}
      }


      /* Get the unique vertices of edges on this model edge */

      vlist = List_New(10);
      idx2 = 0;
      while ((me2 = List_Next_Entry(elist,&idx2))) {
	for (i = 0; i < 2; i++) {
	  mv = ME_Vertex(me2,i);

#ifdef MSTK_HAVE_MPI
          if (!MEnt_IsMarked(mv,ownedmk)) continue; /* Vtx cnctd to only ghost elements */
#endif
	  if (!MEnt_IsMarked(mv,mkid2)) {
	    MEnt_Mark(mv,mkid2);
	    List_Add(vlist,mv);
	  }
	}
      }
      List_Unmark(vlist,mkid2);

      List_Unmark(elist,mkid2);
      List_Delete(elist);

      if (List_Num_Entries(vlist) == 0) {
        List_Delete(vlist);
        continue;
      }

      /* Make a nodeset from the vertices on the closure of this model edge */

      if (nn == nnalloc) {
	nnalloc *= 2;
	*node_sets = (List_ptr *) MSTK_realloc(*node_sets,nnalloc*sizeof(List_ptr));
	*node_set_ids = (int *) MSTK_realloc(*node_set_ids,nnalloc*sizeof(int));
      }

      (*node_sets)[nn] = vlist; /* don't delete vlist */
      (*node_set_ids)[nn] = nid;
      nn++;
    }

    idx = 0;
    while ((me = MESH_Next_Edge(mesh,&idx)))
      MEnt_Unmark(me,mkid1);


    MSTK_FreeMarker(mkid1);
    MSTK_FreeMarker(mkid2);

    *num_node_set = nn;
  }


  /* Local ID of face in region according to MSTK conventions - Zero based */

  int MF_LocalID_in_Region(MFace_ptr mf, MRegion_ptr mr) {
    MRType mrtype = MR_ElementType(mr);

    switch (mrtype) {
    case TET: case PYRAMID: case PRISM: case HEX: {
      int i, j;
      int allfound = 0;

      List_ptr rverts = MR_Vertices(mr);
      List_ptr fverts = MF_Vertices(mf,1,0);

      int nrf = MSTK_nrf_template[mrtype];
      for (i = 0; i < nrf; i++) {
        int nrfv = MSTK_rfv_template[mrtype][i][0];

        allfound = 1;
        for (j = 0; j < nrfv; j++) {
          MVertex_ptr rv = List_Entry(rverts,MSTK_rfv_template[mrtype][i][j+1]);
          if (!List_Contains(fverts,rv)) {
            allfound = 0;
            break;
          }
        }
        if (allfound)
          break;
      }
      
      List_Delete(rverts);
      List_Delete(fverts);

      if (allfound)
        return i;
      else {
        MSTK_Report("MF_LocalID_in_Region","Face not found in region",MSTK_ERROR);
        return -1;
      }
      break;
    }
    case POLYHED: case RUNKNOWN: {
      List_ptr rfaces = MR_Faces(mr);
      int i, found = 0;
      for (i = 0; i < List_Num_Entries(rfaces); i++) {
        if (mf == List_Entry(rfaces,i)) { /* will this work for reduced representations? */
          found = 1;
          break;
        }
      }
      List_Delete(rfaces);
      if (found)
        return i;
      else {
        MSTK_Report("MF_LocalID_in_Region","Face not found in region",MSTK_ERROR);
        return -1;
      }
      break;
    }
    case RDELETED: default:
      return -1;
    }
  }
	    

  /* Local ID of edge in face according to MSTK conventions - Zero based */

  int ME_LocalID_in_Face(MEdge_ptr me, MFace_ptr mf) {
    MFType mftype = MF_ElementType(mf);
    
    switch (mftype) {
    case TRI: case QUAD: case POLYGON: case FUNKNOWN: {
      List_ptr fedges = MF_Edges(mf,1,0);
      int i, found = 0;
      for (i = 0; i < List_Num_Entries(fedges); i++) {
        if (me == List_Entry(fedges,i)) { /* will this work for reduced representations? */
          found = 1;
          break;
        }
      }
      List_Delete(fedges);
      if (found)
        return i;
      else {
        MSTK_Report("ME_LocalID_in_Face","Edge not found in face",MSTK_ERROR);
        return -1;
      }
      break;
    }
    case FDELETED: default:
      return -1;
    }
  }


#ifdef __cplusplus
}
#endif

