#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"
#include "exodusII.h"


#ifdef __cplusplus
extern "C" {
#endif

  /* this function collects element block inforamtion based on element type */
  void MESH_Get_Element_Block_Info(Mesh_ptr mesh, 
                                   int *num_element_blocks_glob,
				   MSet_ptr **element_blocks_glob, 
				   int **element_block_ids_glob, 
				   char ***element_block_types_glob,
                                   MAttrib_ptr block_type_att,
                                   MSTK_Comm comm);

  /* this function collects side set information */
  void MESH_Get_Side_Set_Info(Mesh_ptr mesh, 
                              int enable_geometry_sets,
                              int *num_side_set_glob, 
			      MSet_ptr **side_sets_glob, 
                              int **side_set_ids_glob,
                              MSTK_Comm comm);
  void MESH_Get_Node_Set_Info(Mesh_ptr mesh, 
                              int enable_geometry_sets,
                              int *num_node_set_glob, 
			      MSet_ptr **node_sets_glob, 
                              int **node_set_ids_glob,
                              MSTK_Comm comm);
  void MESH_Get_Element_Set_Info(Mesh_ptr mesh, 
                              int *num_node_set_glob, 
			      MSet_ptr **node_sets_glob, 
                              int **node_set_ids_glob,
                              MSTK_Comm comm);
  void MESH_Get_Attribute_Info(Mesh_ptr mesh,
                               int *num_element_atts,
                               char ***element_att_names,
                               int *num_node_atts,
                               char ***node_att_names,
                               int *num_sideset_atts,
                               char ***sideset_att_names,
                               MSTK_Comm comm);

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
     opts[0] = 0/1 --- verbose (print stat information - default 0)
     opts[1] = 0/1 --- enable_geometry_sets (sidesets, nodesets based on 
                       GEntID, GEntDim - default 0)
     opts[2] = 0/1 --- Don't Enforce/Enforce contiguous entity numbering
                       (Default = 1 ** ENFORCE CONTIGUOUS ORDERING **)
                       (NOTE: This means that if you read in an Exodus II
                        file and write it out you will get differences 
                        in node numbering)
  */

  int ownedmk = 0;
  MAttrib_ptr ownedatt = NULL;

  int MESH_ExportToExodusII(Mesh_ptr mesh, const char *filename, 
			    const int natt, const char **attnames, 
			    const int *opts, MSTK_Comm comm) {

    int enable_geometry_sets, verbose, enforce_contiguous_ids;
    int i, j, k, idx, idx2, ival;
    int nvowned, neowned, nfowned, nrowned, nvall, neall, nfall, nrall;
    int num_element_block_glob, num_side_set_glob, num_node_set_glob,
      num_element_set_glob;
    int num_element_atts_glob, num_node_atts_glob, num_sideset_atts_glob;    
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
    int *element_block_ids_glob, *side_set_ids_glob, *node_set_ids_glob,
      *element_set_ids_glob;
    int num_face_block;
    char **element_block_types_glob, block_name[256];
    char **element_att_names_glob, **node_att_names_glob, 
      **sideset_att_names_glob, **elementset_att_names_glob;
    MSet_ptr *element_blocks_glob, *side_sets_glob, *node_sets_glob, 
      *element_sets_glob;
    List_ptr face_block;
    List_ptr fverts, rverts;
    MAttrib_ptr vidatt=NULL, fidatt=NULL, ridatt=NULL;
    char modfilename[256], funcname[64]="MESH_ExportToExodusII";
    int status;
    double rval;
    void *pval;
    MAttrib_ptr elatt, ndatt; /* Just for testing only */

    int mstk2exo_facemap[5][6]={{0,0,0,0,0,0},
				{4,1,2,3,0,0}, /* TET */
				{0,0,0,0,0,0}, /* PYRAMID, no support in exo*/
			        {4,5,1,2,3,0}, /* PRISM, must verify nums */
			        {5,6,1,2,3,4}};/* HEX */
    int rank, numprocs;


    int meshdim;
    if (MESH_Num_Regions(mesh))
      meshdim = 3;
    else if (MESH_Num_Faces(mesh))
      meshdim = 2;
    else if (MESH_Num_Edges(mesh))
      meshdim = 1;
    else 
      meshdim = 0;


    reptype = MESH_RepType(mesh);


#ifdef MSTK_HAVE_MPI
    if (comm != 0 && comm != MPI_COMM_NULL) {
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
    }
    else {
      rank = 0;
      numprocs = 1;
      strcpy(modfilename,filename);      
    }
#else
    strcpy(modfilename,filename);
#endif    


    verbose = opts ? opts[0] : 0;  // default is 0
    enable_geometry_sets = opts ? opts[1] : 0;  // default is 0
    enforce_contiguous_ids = opts ? opts[1] : 1; // default is 1

    if (verbose)
      fprintf(stdout,"\nMesh representation type is %s\n", 
	      MESH_rtype_str[reptype]);
    nvall = MESH_Num_Vertices(mesh);
    neall = MESH_Num_Edges(mesh);
    nfall = MESH_Num_Faces(mesh);
    nrall = MESH_Num_Regions(mesh);
    
    if (nvall == 0) {
      MSTK_Report("MESH_ExportToExodusII","No vertices in mesh?\n",MSTK_ERROR);
      return 0;
    }


    int element_dim = -1, side_dim = -1;
    int num_element = 0, num_side = 0;

    vidatt = MAttrib_New(mesh,"vidatt",INT,MVERTEX);

    if (nrall) {
      element_dim = 3; side_dim = 2; boundary_dim = 2;


      /* Mark the regions and their lower dimensional entities. In the
         case of, parallel meshes, mark only owned regions - This has
         the dual purpose of marking owned entities as well as giving
         the entities a temporary renumbering to ensure contiguous
         IDs */
        
#ifdef MSTK_USE_MARKERS
      ownedmk = MSTK_GetMarker();
#else
      ownedatt = MAttrib_New(mesh, "ownedatt", INT, MALLTYPE);
#endif
        
      nrowned = 0; nfowned = 0; neowned = 0; nvowned = 0;
      idx = 0;
      while ((mr = MESH_Next_Region(mesh,&idx))) {

#ifdef MSTK_HAVE_MPI
        if (comm == 0 || comm == MPI_COMM_NULL || MR_PType(mr) != PGHOST) {
#endif

          // If we created some degenerate elements (like boundary ghosts)
          // ignore them

          if (MR_ElementType(mr) == RUNKNOWN) continue;

#ifdef MSTK_USE_MARKERS
          MEnt_Mark(mr,ownedmk);
#else
          MEnt_Set_AttVal(mr, ownedatt, 1, 0.0, NULL);
#endif
          nrowned++;
            
          List_ptr rfaces = MR_Faces(mr);
          idx2 = 0;
          while ((mf = (MFace_ptr) List_Next_Entry(rfaces,&idx2))) {
            int fowned;
#ifdef MSTK_USE_MARKERS
            fowned = MEnt_IsMarked(mf,ownedmk);
#else
            MEnt_Get_AttVal(mf, ownedatt, &fowned, &rval, &pval);
#endif
            if (!fowned) {
#ifdef MSTK_USE_MARKERS
              MEnt_Mark(mf,ownedmk);
#else
              MEnt_Set_AttVal(mf, ownedatt, 1, 0.0, NULL);
#endif
              nfowned++;
              
              List_ptr fedges = MF_Edges(mf,1,0);
              int idx3 = 0;
              while ((me = (MEdge_ptr) List_Next_Entry(fedges,&idx3))) {
                int eowned;
#ifdef MSTK_USE_MARKERS
                eowned = MEnt_IsMarked(me,ownedmk);
#else
                MEnt_Get_AttVal(me, ownedatt, &eowned, &rval, &pval);
#endif
                if (!eowned) {
#ifdef MSTK_USE_MARKERS
                  MEnt_Mark(me,ownedmk);
#else
                  MEnt_Set_AttVal(me, ownedatt, 1, 0.0, NULL);
#endif
                  neowned++;
                  
                  int kk;
                  for (kk = 0; kk < 2; kk++) {
                    mv = ME_Vertex(me,kk);
                    int vowned;
#ifdef MSTK_USE_MARKERS
                    vowned = MEnt_IsMarked(mv,ownedmk);
#else
                    MEnt_Get_AttVal(mv, ownedatt, &vowned, &rval, &pval);
#endif
                    if (!vowned) {
#ifdef MSTK_USE_MARKERS
                      MEnt_Mark(mv,ownedmk);
#else
                      MEnt_Set_AttVal(mv, ownedatt, 1, 0.0, NULL);
#endif
                      nvowned++;
                      if (enforce_contiguous_ids)
                        MEnt_Set_AttVal(mv,vidatt,nvowned,0.0,NULL);
                      else
                        MEnt_Set_AttVal(mv,vidatt,MV_ID(mv),0.0,NULL);
                    }
                  }
                }
              } /* while (me...) */
              List_Delete(fedges);
            }
          } /* while (mf...) */
          List_Delete(rfaces);
          
#ifdef MSTK_HAVE_MPI
        }
#endif
      }      

      
      num_element = nrowned; num_side = nfowned;

      if(verbose)
	fprintf(stdout,"\nThis is a 3D volume mesh with %d Vertices, %d Edges, %d Faces and %d Regions\n",\
		nvowned,neowned,nfowned,nrowned);

    }
    else if (nfall) {
      
      element_dim = 2; side_dim = 1; boundary_dim = 1;

      /* Mark the faces and their lower dimensional entities. In the
         case of, parallel meshes, mark only owned faces - This has
         the dual purpose of marking owned entities as well as giving
         the entities a temporary renumbering to ensure contiguous
         IDs */
        
#ifdef MSTK_USE_MARKERS
      ownedmk = MSTK_GetMarker();
#else
      ownedatt = MAttrib_New(mesh, "ownedatt", INT, MALLTYPE);
#endif
      
      nrowned = 0; nfowned = 0; neowned = 0; nvowned = 0;
      idx2 = 0;
      while ((mf = (MFace_ptr) MESH_Next_Face(mesh,&idx2))) {

#ifdef MSTK_HAVE_MPI
        if (comm == 0 || comm == MPI_COMM_NULL || MF_PType(mf) != PGHOST) {
#endif

          // If we created some degenerate elements (like boundary ghosts)
          // ignore them

          if (MF_ElementType(mf) == FUNKNOWN) continue; 

#ifdef MSTK_USE_MARKERS
          MEnt_Mark(mf,ownedmk);
#else
          MEnt_Set_AttVal(mf, ownedatt, 1, 0, NULL);
#endif
          nfowned++;
          
          List_ptr fedges = MF_Edges(mf,1,0);
          int idx3 = 0;
          while ((me = (MEdge_ptr) List_Next_Entry(fedges,&idx3))) {
            int eowned;
#ifdef MSTK_USE_MARKERS
            eowned = MEnt_IsMarked(me, ownedmk);
#else
            MEnt_Get_AttVal(me, ownedatt, &eowned, &rval, &pval);
#endif
            if (!eowned) {
#ifdef MSTK_USE_MARKERS
              MEnt_Mark(me, ownedmk);
#else
              MEnt_Set_AttVal(me, ownedatt, 1, 0.0, NULL);
#endif
              neowned++;
              
              int kk;
              for (kk = 0; kk < 2; kk++) {
                mv = ME_Vertex(me,kk);
                int vowned;
#ifdef MSTK_USE_MARKERS
                vowned = MEnt_IsMarked(mv, ownedmk);
#else
                MEnt_Get_AttVal(mv, ownedatt, &vowned, &rval, &pval);
#endif
                if (!vowned) {
#ifdef MSTK_USE_MARKERS
                  MEnt_Mark(mv,ownedmk);
#else
                  MEnt_Set_AttVal(mv, ownedatt, 1, 0, NULL);
#endif
                  nvowned++;
                  if (enforce_contiguous_ids)
                    MEnt_Set_AttVal(mv,vidatt,nvowned,0.0,NULL);
                  else
                    MEnt_Set_AttVal(mv,vidatt,MV_ID(mv),0.0,NULL);
                }
              }
            }
          } /* while (me...) */
          List_Delete(fedges);

#ifdef MSTK_HAVE_MPI
        }
#endif

      } /* while (mf...) */
        

      num_element = nfowned; num_side = neowned;

      if(verbose)
      fprintf(stdout,"\nThis is a 3D surface mesh with %d Vertices, %d Edges and %d Faces\n", nvowned,neowned,nfowned);

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
      MSTK_Report("MESH_ExportToExodusII", "Error in creation of Exodus file",
                  MSTK_ERROR);
      return 0;
    }


    
    if(verbose) {
      if (enable_geometry_sets)
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
      if (enable_geometry_sets)
	fprintf(stdout,"\nElement block and side set based on geometric id is enabled.\n");
      else
	fprintf(stdout,"\nElement block and side set based on geometric id is disabled.\n");
    }

    /* Even though there is a native element type for every mesh
       region or face we have to consider what the element type is
       considered as for Exodus II output. For example, a HEX element
       in a mixed mesh will be written as a NFACED (polyhedral)
       element */

    MAttrib_ptr block_type_att = MAttrib_New(mesh,"block_type",INT,MALLTYPE); 

    MESH_Get_Element_Block_Info(mesh, &num_element_block_glob, 
                                &element_blocks_glob,
                                &element_block_ids_glob, 
                                &element_block_types_glob, 
                                block_type_att,
                                comm);


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
    for (i = 0; i < num_element_block_glob; i++) {
      int nelem = MSet_Num_Entries(element_blocks_glob[i]);
      for (k = 0; k < nelem; k++) {
        MEntity_ptr ment = MSet_Entry(element_blocks_glob[i],k);
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
#ifdef MSTK_USE_MARKERS
    int mkid1 = MSTK_GetMarker();
#endif
    fidatt = MAttrib_New(mesh,"local_face_ids",INT,MFACE);
    

    j = 0;
    for (i = 0; i < num_element_block_glob; i++) {
      if (strncasecmp(element_block_types_glob[i],"NFACED",6) == 0) {

	idx = 0;
	while ((mr = MSet_Next_Entry(element_blocks_glob[i],&idx))) {
	  List_ptr rfaces = MR_Faces(mr);
	  MFace_ptr rf;

	  idx2 = 0;
	  while ((rf = List_Next_Entry(rfaces,&idx2))) {
            int inlist;
#ifdef MSTK_USE_MARKERS
            inlist = MEnt_IsMarked(rf, mkid1);
#else
            MEnt_Get_AttVal(rf, fidatt, &ival, &rval, &pval);
            inlist = ival > 0 ? 1 : 0;
#endif
	    if (!inlist) {
	      List_Add(face_block,rf);
	      MEnt_Set_AttVal(rf,fidatt,++j,0.0,NULL);
	      nvfblock += MF_Num_Vertices(rf);
#ifdef MSTK_USE_MARKERS
	      MEnt_Mark(rf,mkid1);
#endif
	    }
	  }
	  List_Delete(rfaces);
	}

      }
    }

#ifdef MSTK_USE_MARKERS
    List_Unmark(face_block,mkid1);
    MSTK_FreeMarker(mkid1);
#endif
   
    if (List_Num_Entries(face_block) == 0) {
      List_Delete(face_block);
      face_block = NULL;
      num_face_block = 0;
    }


    /* COLLECT ELEMENT SET INFO - Different from element blocks - no
       requirement that their union cover the entire domain and their
       intersection be null - An element can occur in multiple sets or
       none at all */

    MESH_Get_Element_Set_Info(mesh, &num_element_set_glob, &element_sets_glob,
                              &element_set_ids_glob, comm);

    /* COLLECT SIDE SET INFO */

    MESH_Get_Side_Set_Info(mesh, enable_geometry_sets,
                           &num_side_set_glob, &side_sets_glob, 
                           &side_set_ids_glob, comm);



    /* COLLECT NODE SET INFO */

    MESH_Get_Node_Set_Info(mesh, enable_geometry_sets,
                           &num_node_set_glob, &node_sets_glob, 
                           &node_set_ids_glob, comm);

    

    /* INFORMATION GATHERING IS LARGELY DONE - WRITE THE FILE */

    /* write global parameters */
    
    ex_init_params par;
    strcpy(par.title,filename);
    par.num_nodes = nvowned;

    if (element_dim == 3 || boundary_dim == 2)
      par.num_dim = 3;
    else 
      par.num_dim = 2;

    if (element_dim == 3) {
      par.num_edge = 0;
      par.num_face = face_block ? List_Num_Entries(face_block) : 0;
      par.num_face_blk = num_face_block;
      par.num_edge_blk = 0;
      par.num_elem = nrowned;
      par.num_elem_blk = num_element_block_glob;  
      par.num_node_sets = num_node_set_glob;
      par.num_edge_sets = 0;
      par.num_face_sets = 0;
      par.num_side_sets = num_side_set_glob;
      par.num_elem_sets = num_element_set_glob;

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
      par.num_elem = nfowned;
      par.num_elem_blk = num_element_block_glob;  
      par.num_node_sets = num_node_set_glob;
      par.num_edge_sets = 0;
      par.num_face_sets = 0;
      par.num_side_sets = num_side_set_glob;
      par.num_elem_sets = num_element_set_glob;

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
    xcoord = (double *) malloc(nvowned*sizeof(double));
    ycoord = (double *) malloc(nvowned*sizeof(double));
    zcoord = (double *) malloc(nvowned*sizeof(double));

    idx = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx))) {
      int vowned;
#ifdef MSTK_HAVE_MPI
#ifdef MSTK_USE_MARKERS 
      vowned = MEnt_IsMarked(mv,ownedmk);
#else
      MEnt_Get_AttVal(mv, ownedatt, &vowned, &rval, &pval);
#endif
      if (!vowned) continue;
#endif

      MV_Coords(mv,vxyz);
      MEnt_Get_AttVal(mv,vidatt,&vid,&rval,&pval);
      xcoord[vid-1] = vxyz[0];
      ycoord[vid-1] = vxyz[1];
      zcoord[vid-1] = vxyz[2];
    }
    err = ex_put_coord (exoid, xcoord, ycoord, zcoord);

    free(xcoord);
    free(ycoord);
    free(zcoord);
    

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

	while ((fv = List_Next_Entry(fverts,&idx2))) {
          MEnt_Get_AttVal(fv,vidatt,&vid,&rval,&pval);
	  connect[k++] = vid;
        }

	nnpe[i++] = List_Num_Entries(fverts);

	List_Delete(fverts);
      }

      /* Output number of nodes per face in face block */
      
      ex_put_entity_count_per_polyhedra(exoid, EX_FACE_BLOCK, block_id, nnpe);

      /* Output connectivity of faces in face block */

      ex_put_conn(exoid, EX_FACE_BLOCK, block_id, connect, NULL, NULL);

      free(connect);
      free(nnpe);
      
      List_Delete(face_block);
    }

    

    /* write element blocks info */

    for (i = 0; i < num_element_block_glob; i++) { 

      char block_name[256];

      block_id = element_block_ids_glob[i];

      if (strncasecmp(element_block_types_glob[i],"NFACED",6) == 0) {
	
	int nfrblock = 0;

	nnpe = (int *) calloc(List_Num_Entries(element_blocks_glob[i]),
				   sizeof(int));	
	idx = 0; j = 0;
	while ((mr = MSet_Next_Entry(element_blocks_glob[i],&idx))) {	  
	  nnpe[j] = MR_Num_Faces(mr);
	  nfrblock += nnpe[j];
	  j++;
	}

	connect = (int *) calloc(nfrblock,sizeof(int));
	idx = 0; j = 0;
	while ((mr = MSet_Next_Entry(element_blocks_glob[i],&idx))) {
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
		     MSet_Num_Entries(element_blocks_glob[i]),
		     0, /* nodes */
		     0, /* edges */
		     nfrblock, /* faces */
		     0); /* attribute count */

	ex_put_name(exoid, EX_ELEM_BLOCK, block_id, block_name);

	ex_put_entity_count_per_polyhedra(exoid, EX_ELEM_BLOCK, block_id, nnpe);

	ex_put_conn(exoid, EX_ELEM_BLOCK, block_id, NULL, NULL, connect);

	free(nnpe);

	free(connect);

      }
      else if (strncasecmp(element_block_types_glob[i],"NSIDED",6) == 0) {

	int nvfblock = 0;

	nnpe = (int *) calloc(MSet_Num_Entries(element_blocks_glob[i]),
				   sizeof(int));	
	idx = 0; j = 0;
	while ((mf = MSet_Next_Entry(element_blocks_glob[i],&idx))) {	  
	  nnpe[j] = MF_Num_Vertices(mf);
	  nvfblock += nnpe[j];
	  j++;
	}

	connect = (int *) calloc(nvfblock,sizeof(int));
	idx = 0; j = 0;
	while ((mf = MSet_Next_Entry(element_blocks_glob[i],&idx))) {
	  List_ptr fverts = MF_Vertices(mf,1,0);
	  MVertex_ptr fv;
	  idx2 = 0;
	  while ((fv = List_Next_Entry(fverts,&idx2))) {
            MEnt_Get_AttVal(fv,vidatt,&vid,&rval,&pval);
	    connect[j++] = vid;
          }
	  List_Delete(fverts);
	}


	sprintf(block_name,"POLYGON_BLOCK_%-d",block_id);
	
	ex_put_block(exoid, EX_ELEM_BLOCK, block_id, "NSIDED",
		     MSet_Num_Entries(element_blocks_glob[i]),
		     nvfblock, /* nodes */
		     0, /* edges */
		     0, /* faces */
		     0); /* attribute count */

	ex_put_name(exoid, EX_ELEM_BLOCK, block_id, block_name);

	ex_put_entity_count_per_polyhedra(exoid, EX_ELEM_BLOCK, block_id, nnpe);

	ex_put_conn(exoid, EX_ELEM_BLOCK, block_id, connect, NULL, NULL);

	free(nnpe);

	free(connect);


      }
      else {

	int nelem, nelnodes=0;

	sprintf(block_name,"BLOCK_%-d",block_id);
	nelem = MSet_Num_Entries(element_blocks_glob[i]);
	
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
	  
	  
	  connect = (int *) malloc(nelnodes*nelem*sizeof(int));
	  
	  idx = 0; k = 0;
	  while ((mr = MSet_Next_Entry(element_blocks_glob[i],&idx))) {
	    rverts = MR_Vertices(mr);
	    for (j = 0; j < nelnodes; j++) {
              MEnt_Get_AttVal(List_Entry(rverts,j),vidatt,&vid,&rval,&pval);
	      connect[nelnodes*k+j] = vid;
            }
	    List_Delete(rverts);	    
	    k++;
	  }
	  
	  ex_put_elem_conn(exoid, block_id, connect);

	  free(connect);
	  
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

	  connect = (int *) malloc(nelnodes*nelem*sizeof(int));

	  idx = 0; k = 0;
	  while ((mf = MSet_Next_Entry(element_blocks_glob[i],&idx))) {
	    fverts = MF_Vertices(mf,1,0);
	    for (j = 0; j < nelnodes; j++) {
              MEnt_Get_AttVal(List_Entry(fverts,j),vidatt,&vid,&rval,&pval);
	      connect[nelnodes*k+j] = vid;
            }
	    List_Delete(fverts);
	    k++;
	  }

	  ex_put_elem_conn(exoid, block_id, connect);

	  free(connect);
	}

      }

    }

    if (verbose)
      fprintf(stdout,"Element block information written into %s\n", filename);



    /* Write out the node set information */

   
    for (i = 0; i < num_node_set_glob; i++) {
      int nnodes = MSet_Num_Entries(node_sets_glob[i]);
      ex_put_node_set_param(exoid, node_set_ids_glob[i], nnodes, 0);

      int *node_list = (int *) malloc(nnodes*sizeof(int));

      idx = 0; j = 0;
      while ((mv = MSet_Next_Entry(node_sets_glob[i],&idx))) {
        MEnt_Get_AttVal(mv,vidatt,&vid,&rval,&pval);
	node_list[j++] = vid;
      }

      ex_put_node_set(exoid, node_set_ids_glob[i], node_list);

      free(node_list);
    }

    if (verbose)
      fprintf(stdout,"Node sets information written into %s\n", filename);



    
    /* Write out the side set information */


    for (i = 0; i < num_side_set_glob; i++) {

      int nsides = MSet_Num_Entries(side_sets_glob[i]);
      ex_put_side_set_param(exoid, side_set_ids_glob[i], nsides, 0);

      int *elem_list = (int *) malloc(nsides*sizeof(int));
      int *side_list = (int *) malloc(nsides*sizeof(int));

      if (nrowned) {
	idx = 0; j = 0;
	while ((mf = MSet_Next_Entry(side_sets_glob[i],&idx))) {
          if (MF_PType(mf) == PGHOST &&  !MF_OnParBoundary(mf)) continue;

	  List_ptr fregs = MF_Regions(mf);

	  if (!fregs || List_Num_Entries(fregs) == 0)
	    MSTK_Report("MESH_ExportToEXODUSII",
			"Standalone face with no regions in side set",
			MSTK_FATAL);
	  mr = List_Entry(fregs,0);

#ifdef MSTK_HAVE_MPI
          if (comm && MR_PType(mr) == PGHOST) {
            if (List_Num_Entries(fregs) > 1)
              mr = List_Entry(fregs,1);  /* at least 1 region should be owned */
            else  // First check in loop should prevent it from coming here
              MSTK_Report("MESH_ExportToExodusII", "Trying to write out ghost side set face thats not on the parallel boundary", MSTK_FATAL);
          }
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
	while ((me = MSet_Next_Entry(side_sets_glob[i],&idx))) {
	  List_ptr efaces = ME_Faces(me);

	  if (!efaces || List_Num_Entries(efaces) == 0)
	    MSTK_Report("MESH_ExportToEXODUSII",
			"Standalone edge with no faces in side set",
			MSTK_FATAL);
	  mf = List_Entry(efaces,0);

#ifdef MSTK_HAVE_MPI
          if (comm && MF_PType(mf) == PGHOST)
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

      free(elem_list);
      free(side_list);

    }

    if (verbose)
      fprintf(stdout,"Side set information written into %s\n", filename);


    /* Write out element set information */
  
    for (i = 0; i < num_element_set_glob; i++) {
      int nelements = MSet_Num_Entries(element_sets_glob[i]);
      ex_put_set_param(exoid, EX_ELEM_SET, element_set_ids_glob[i], nelements, 0);

      int *element_list = (int *) malloc(nelements*sizeof(int));

      MEntity_ptr ment;
      idx = 0; j = 0;
      while ((ment = MSet_Next_Entry(element_sets_glob[i],&idx)))
        element_list[j++] = elem_id[MEnt_ID(ment)-1];

      ex_put_set(exoid, EX_ELEM_SET, element_set_ids_glob[i], element_list, NULL);

      free(element_list);
    }
    
    if (verbose)
      fprintf(stdout,"Node sets information written into %s\n", filename);


#ifdef MSTK_HAVE_MPI

    if (numprocs > 1) {

      /* Write out node map (global IDs of nodes) */

      int *node_map = (int *) malloc(nvowned*sizeof(int));
      
      idx = 0; i = 0;
      while ((mv = MESH_Next_Vertex(mesh,&idx))) {
        int vowned;
#ifdef MSTK_USE_MARKERS
        vowned = MEnt_IsMarked(mv,ownedmk);
#else
        MEnt_Get_AttVal(mv, ownedatt, &vowned, &rval, &pval);
#endif
        if (vowned)
          node_map[i++] = MV_GlobalID(mv);
      }

      status = ex_put_node_num_map(exoid, node_map);
      if (status < 0)
        MSTK_Report(funcname,"Error while writing node map in Exodus II file",
                    MSTK_FATAL);

      free(node_map);

      
      /* Write out element map (global IDs of elements) */

      int *elem_map;

      if (nrowned) {
        elem_map = (int *) malloc(nrowned*sizeof(int));
        idx = 0; i = 0;
        while ((mr = MESH_Next_Region(mesh,&idx))) {
          int rowned;
#ifdef MSTK_USE_MARKERS
          rowned = MEnt_IsMarked(mr, ownedmk);
#else
          MEnt_Get_AttVal(mr, ownedatt, &rowned, &rval, &pval);
#endif
          if (rowned)
            elem_map[i++] = MR_GlobalID(mr);
        }
      }
      else { // assume surface mesh
        elem_map = (int *) malloc(nfowned*sizeof(int));
        idx = 0; i = 0;
        while ((mf = MESH_Next_Face(mesh,&idx))) {
          int fowned;
#ifdef MSTK_USE_MARKERS
          fowned = MEnt_IsMarked(mf, ownedmk);
#else
          MEnt_Get_AttVal(mf, ownedatt, &fowned, &rval, &pval);
#endif
          if (fowned)
            elem_map[i++] = MF_GlobalID(mf);
        }
        
        status = ex_put_elem_num_map(exoid, elem_map);
        if (status < 0)
          MSTK_Report(funcname,"Error while writing element map in Exodus II file",
                      MSTK_FATAL);
        
        free(elem_map);
      }
    }

#endif /* MSTK_HAVE_MPI */

     
    /* Write out "results" variables tied to elements, nodes and sidesets */
    /* Even though they are called "results" they can be any quantity     */
    /* PS: Do we HAVE to output the time value? */


    /* COLLECT ATTRIBUTE/FIELD INFO */

    MESH_Get_Attribute_Info(mesh, &num_element_atts_glob, 
                            &element_att_names_glob,
                            &num_node_atts_glob, &node_att_names_glob,
                            &num_sideset_atts_glob, &sideset_att_names_glob,
                            comm);


    char veckey[16] = "_veccomp";
    int keylen = strlen(veckey);
    
    if (num_node_atts_glob) {

      /* If there are some vector/tensor attributes, each of their
         components has to be written out separately by suffixing the
         attribute name with a keyword and component index */

      int num_node_atts_out = 0;      
      char **node_att_names_out =  malloc(num_node_atts_glob*6*sizeof(char *));

      for (j = 0; j < num_node_atts_glob; j++) {
        MAttrib_ptr att = MESH_AttribByName(mesh,node_att_names_glob[j]);
        if (MAttrib_Get_Type(att) == VECTOR || MAttrib_Get_Type(att) == TENSOR) {
          int ncomp = MAttrib_Get_NumComps(att);
          for (k = 0; k < ncomp; k++) {
            node_att_names_out[num_node_atts_out] = malloc(256*sizeof(char));
            sprintf(node_att_names_out[num_node_atts_out],
                    "%s%s%1d",node_att_names_glob[j],veckey,k);
            num_node_atts_out++;
          }
        }
        else {
          node_att_names_out[num_node_atts_out] = malloc(256*sizeof(char));
          strcpy(node_att_names_out[num_node_atts_out],node_att_names_glob[j]);
          num_node_atts_out++;
        }
      }

      status = ex_put_var_param(exoid, "n", num_node_atts_out);
      if (status < 0)
        MSTK_Report(funcname,"Error while writing node variable parameters",
                    MSTK_FATAL);

      status = ex_put_var_names(exoid, "n", num_node_atts_out,
                                node_att_names_out);
      if (status < 0)
        MSTK_Report(funcname,"Error while writing node variable parameters",
                    MSTK_FATAL);

      /* Now write out each variable */

      double *node_vars = (double *) malloc(nvowned*sizeof(double));
      int attid = 1;
      for (j = 0; j < num_node_atts_glob; ++j) {
        MAttrib_ptr att = MESH_AttribByName(mesh,node_att_names_glob[j]);
        
        if (MAttrib_Get_Type(att) == VECTOR || MAttrib_Get_Type(att) == TENSOR) {
          /* Write out each component separately */

          int ncomp = MAttrib_Get_NumComps(att);
          int n;
          for (n = 0; n < ncomp; n++) {
            idx = 0;
            while ((mv = MESH_Next_Vertex(mesh,&idx))) {
              int vowned, vid;
#ifdef MSTK_HAVE_MPI
#ifdef MSTK_USE_MARKERS 
              vowned = MEnt_IsMarked(mv,ownedmk);
#else
              MEnt_Get_AttVal(mv, ownedatt, &vowned, &rval, &pval);
#endif
              if (!vowned) continue;
#endif

              MEnt_Get_AttVal(mv,vidatt,&vid,&rval,&pval);
              MEnt_Get_AttVal(mv,att,&ival,&rval,&pval);
              node_vars[vid-1] = ((double *) pval)[n];
            }
             
            status = ex_put_nodal_var(exoid, 1, attid, nvowned, node_vars);
            if (status < 0)
              MSTK_Report(funcname,"Error while writing node variable",
                          MSTK_FATAL);

            attid++;
          }
        }
        else {        
          idx = 0; k = 0;
          while ((mv = MESH_Next_Vertex(mesh,&idx))) {
            int vowned, vid;
#ifdef MSTK_HAVE_MPI
#ifdef MSTK_USE_MARKERS 
            vowned = MEnt_IsMarked(mv,ownedmk);
#else
            MEnt_Get_AttVal(mv, ownedatt, &vowned, &rval, &pval);
#endif
            if (!vowned) continue;
#endif

            MEnt_Get_AttVal(mv,vidatt,&vid,&rval,&pval);
            MEnt_Get_AttVal(mv,att,&ival,&rval,&pval);
            node_vars[vid-1] = rval;
          }
          
          status = ex_put_nodal_var(exoid, 1, attid, nvowned, node_vars);
          if (status < 0)
            MSTK_Report(funcname,"Error while writing node variable",
                        MSTK_FATAL);

          attid++;
        }
      }
      free(node_vars);
     
      for (i = 0; i < num_node_atts_out; i++)
        free(node_att_names_out[i]);
      free(node_att_names_out);
    }



    if (num_element_atts_glob) {
      /* If there are some vector/tensor attributes, each of their
         components has to be written out separately by suffixing the
         attribute name with a keyword and component index */

      int num_element_atts_out = 0;      
      char **element_att_names_out = 
          malloc(num_element_atts_glob*6*sizeof(char *));

      for (j = 0; j < num_element_atts_glob; j++) {
        MAttrib_ptr att = MESH_AttribByName(mesh,element_att_names_glob[j]);
        if (MAttrib_Get_Type(att) == VECTOR || MAttrib_Get_Type(att) == TENSOR) {
          int ncomp = MAttrib_Get_NumComps(att);
          for (k = 0; k < ncomp; k++) {
            element_att_names_out[num_element_atts_out] = malloc(256*sizeof(char));
            sprintf(element_att_names_out[num_element_atts_out],
                    "%s%s%1d",element_att_names_glob[j],veckey,k);
            num_element_atts_out++;
          }
        }
        else {
          element_att_names_out[num_element_atts_out] = malloc(256*sizeof(char));
          strcpy(element_att_names_out[num_element_atts_out],
                 element_att_names_glob[j]);
          num_element_atts_out++;
        }
      }

      status = ex_put_var_param(exoid, "e", num_element_atts_out);
      if (status < 0)
        MSTK_Report(funcname,"Error while writing variable parameters",
                    MSTK_FATAL);

      status = ex_put_var_names(exoid, "e", num_element_atts_out,
                                element_att_names_out);
      if (status < 0)
        MSTK_Report(funcname,"Error while writing element variable parameters",
                    MSTK_FATAL);

      /* for now we write out a value of a variable for every element
         of every block, so the truth table has all 1s */

      int *truth_table = 
        (int *) malloc(num_element_block_glob*num_element_atts_out*sizeof(int));
      for (i = 0, k = 0; i < num_element_block_glob; ++i)
        for (j = 0; j < num_element_atts_out; ++j)
          truth_table[k++] = 1;

      status = ex_put_elem_var_tab(exoid, num_element_block_glob, 
                                   num_element_atts_out, truth_table);
      if (status < 0)
        MSTK_Report(funcname,"Error while writing element variable truth table",
                    MSTK_FATAL);

      
      /* Now write out each variable */

      int attid = 1;
      for (j = 0; j < num_element_atts_glob; ++j) {
        MAttrib_ptr att = MESH_AttribByName(mesh,element_att_names_glob[j]);

        if (MAttrib_Get_Type(att) == VECTOR || MAttrib_Get_Type(att) == TENSOR) {
          int ncomp = MAttrib_Get_NumComps(att);
          int n;
          for (n = 0; n < ncomp; n++) {            
            for (i = 0; i < num_element_block_glob; ++i) {
              int nelem = MSet_Num_Entries(element_blocks_glob[i]);
              double *elem_vars = (double *) malloc(nelem*sizeof(double));
              
              MEntity_ptr ment;
              idx = 0; k = 0;
              while ((ment = MSet_Next_Entry(element_blocks_glob[i],&idx))) {
                MEnt_Get_AttVal(ment,att,&ival,&rval,&pval);
                elem_vars[k++] = ((double *) pval)[n];
              }
              
              status = ex_put_elem_var(exoid, 1, attid, 
                                       element_block_ids_glob[i], nelem, 
                                       elem_vars);
              if (status < 0)
                MSTK_Report(funcname,"Error while writing element variable",
                            MSTK_FATAL);
              
              free(elem_vars);
            }
            attid++;
          }
        }
        else {
          for (i = 0; i < num_element_block_glob; ++i) {
            int nelem = MSet_Num_Entries(element_blocks_glob[i]);
            double *elem_vars = (double *) malloc(nelem*sizeof(double));
            
            MEntity_ptr ment;
            idx = 0; k = 0;
            while ((ment = MSet_Next_Entry(element_blocks_glob[i],&idx))) {
              int entid = MEnt_ID(ment);
              MEnt_Get_AttVal(ment,att,&ival,&rval,&pval);
              elem_vars[k++] = rval;
            }
                        
            status = ex_put_elem_var(exoid, 1, attid, element_block_ids_glob[i],
                                     nelem, elem_vars);
            if (status < 0)
              MSTK_Report(funcname,"Error while writing element variable",
                          MSTK_FATAL);
            
            free(elem_vars);
          }
          attid++;
        }
      }

      for (i = 0; i < num_element_atts_out; i++)
        free(element_att_names_out[i]);
      free(element_att_names_out);
    }

    if (num_sideset_atts_glob) {
      /* MSTK_Report(funcname,"Variable writing not implement for sidesets",
         MSTK_WARN); */

      /*
      status = ex_put_var_param(exoid, "s", num_sideset_atts_glob);
      if (status < 0)
        MSTK_Report(funcname,"Error while writing variable parameters",
                    MSTK_FATAL);
      */
    }


    



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


    char msetname[256];

    free(element_block_ids_glob);
    for (i = 0; i < num_element_block_glob; i++) {
      MSet_Name(element_blocks_glob[i],msetname);
      if (strncmp(msetname,"TEMPORARY_",10) == 0)
        MSet_Delete(element_blocks_glob[i]);
    }
    free(element_blocks_glob);
    for (i = 0; i < num_element_block_glob; i++)
      free(element_block_types_glob[i]);
    free(element_block_types_glob);
    free(elem_id);


    if (num_element_set_glob) {
      free(element_set_ids_glob);
      for (i = 0; i < num_element_set_glob; i++) {
        MSet_Name(element_sets_glob[i],msetname);
        if (strncmp(msetname,"TEMPORARY_",10) == 0) 
          MSet_Delete(element_sets_glob[i]);
      }
      free(element_sets_glob);
    }
    
    if (num_side_set_glob) {
      free(side_set_ids_glob);
      for (i = 0; i < num_side_set_glob; i++) {
        MSet_Name(side_sets_glob[i],msetname);
        if (strncmp(msetname,"TEMPORARY_",10) == 0)
          MSet_Delete(side_sets_glob[i]);
      }
      free(side_sets_glob);
    }
    
    if (num_node_set_glob) {
      free(node_set_ids_glob);
      for (i = 0; i < num_node_set_glob; i++) {
        MSet_Name(node_sets_glob[i],msetname);
        if (strncmp(msetname,"TEMPORARY_",10) == 0)
          MSet_Delete(node_sets_glob[i]);
      }
      free(node_sets_glob);
    }

    if (num_element_atts_glob) {
      for (i = 0; i < num_element_atts_glob; ++i)
        free(element_att_names_glob[i]);
      free(element_att_names_glob);
    }
      
    if (num_node_atts_glob) {
      for (i = 0; i < num_node_atts_glob; ++i)
        free(node_att_names_glob[i]);
      free(node_att_names_glob);
    }
      
    if (num_sideset_atts_glob) {
      for (i = 0; i < num_sideset_atts_glob; ++i)
        free(sideset_att_names_glob[i]);
      free(sideset_att_names_glob);
    }
      

#ifdef MSTK_HAVE_MPI
    if (comm) {
#ifdef MSTK_USE_MARKERS
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
#else
      MAttrib_Delete(ownedatt);
#endif
    }
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

  void MESH_Get_Element_Block_Info(Mesh_ptr mesh, 
                                   int *num_element_block_glob, 
                                   MSet_ptr **element_blocks_glob, 
                                   int **element_block_ids_glob, 
                                   char ***element_block_types_glob,
                                   MAttrib_ptr block_type_att,
                                   MSTK_Comm comm) {
    MEdge_ptr me;
    MFace_ptr mf;
    MRegion_ptr mr;
    int i, j, k, nb, bid, found, nballoc;
    int rank=0, numprocs=1;
    double rval;
    void *pval;

#ifdef MSTK_HAVE_MPI
    if (comm) {
      MPI_Comm_size(comm,&numprocs);
      MPI_Comm_rank(comm,&rank);
    }
#endif

    nb = 0; nballoc=10;
    int *element_block_ids = (int *) calloc(nballoc,sizeof(int));
    char **element_block_types = (char **) malloc(nballoc*sizeof(char *));
    MSet_ptr *element_blocks = (MSet_ptr *) calloc(nballoc,sizeof(MSet_ptr));

    int idx = 0;
    
    if (MESH_Num_Regions(mesh)) {

      while ((mr = MESH_Next_Region(mesh, &idx))) {
        int rowned;
#ifdef MSTK_HAVE_MPI
#ifdef MSTK_USE_MARKERS
        rowned = MEnt_IsMarked(mr,ownedmk);  
#else
        MEnt_Get_AttVal(mr, ownedatt, &rowned, &rval, &pval);
#endif
        if (!rowned) continue;  /* Ghost element */
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
	  if (bid == element_block_ids[i]) {           
	    found = 1;
	    MSet_Add(element_blocks[i],mr);    
            MEnt_Set_AttVal(mr,block_type_att,mrtype,0.0,NULL);
	    break;
	  }
	  i++;
	}

	if (!found) {
	  /* New element block */

	  if (nballoc == nb) {
	    nballoc *= 2;
	    element_block_ids = (int *) realloc(element_block_ids,nballoc*sizeof(int));
	    element_blocks = (MSet_ptr *) realloc(element_blocks,nballoc*sizeof(MSet_ptr));
	    element_block_types = (char **) realloc(element_block_types,nballoc*sizeof(char *));
	  }

          /* create a set whose name starts with the string
             TEMPORARY_element_block so that we know we should delete
             it when we are done with it */

          char element_block_name[256];
          sprintf(element_block_name,"TEMPORARY_element_block_%-d",bid);
	  element_blocks[nb] = MSet_New(mesh,element_block_name,MREGION);
	  element_block_ids[nb] = bid;
	  element_block_types[nb] = (char *) malloc(16*sizeof(char));

          if (mrtype == TET)
            strcpy(element_block_types[nb],"TETRA");
          else if (mrtype == PRISM)
            strcpy(element_block_types[nb],"WEDGE");
          else if (mrtype == HEX)
            strcpy(element_block_types[nb],"HEX");
          else
            strcpy(element_block_types[nb],"NFACED");

          MEnt_Set_AttVal(mr,block_type_att,mrtype,0.0,NULL);

	  MSet_Add(element_blocks[nb],mr);
	  nb++;
	}
      }
    }
    else if (MESH_Num_Faces(mesh)) {

      while ((mf = MESH_Next_Face(mesh,&idx))) {
        int fowned;
#ifdef MSTK_HAVE_MPI
#ifdef MSTK_USE_MARKERS
        fowned = MEnt_IsMarked(mf,ownedmk);
#else
        MEnt_Get_AttVal(mf, ownedatt, &fowned, &rval, &pval);
#endif        
        if (!fowned) continue; /* Ghost element */
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
	  if (bid == element_block_ids[i]) {
	    found = 1;
	    MSet_Add(element_blocks[i],mf);    
	    break;
	  }
	  i++;
	}

	if (!found) {
	  /* New element block */

	  if (nballoc == nb) {
	    nballoc *= 2;
	    element_block_ids = (int *) realloc(element_block_ids,
                                                nballoc*sizeof(int));
	    element_blocks = (MSet_ptr *) realloc(element_blocks,
                                                  nballoc*sizeof(MSet_ptr));
	    element_block_types = (char **) realloc(element_block_types,
                                                    nballoc*sizeof(char *));
	  }

          /* Create a set whose name begins with the string
             TEMPORARY_element_block so that we know it should be
             deleted when we are finished with it */

          char element_block_name[256];
          sprintf(element_block_name,"TEMPORARY_element_block_%-d",bid);
	  element_blocks[nb] = MSet_New(mesh,element_block_name,MFACE);
	  element_block_ids[nb] = bid;
	  element_block_types[nb] = (char *) malloc(16*sizeof(char));
	  if (nfv > 4) {
	    strcpy(element_block_types[nb],"NSIDED");
            MEnt_Set_AttVal(mf,block_type_att,POLYGON,0.0,NULL);
          }
	  else {
	    if (nfv == 4) {
	      strcpy(element_block_types[nb],"QUAD");
              MEnt_Set_AttVal(mf,block_type_att,QUAD,0.0,NULL);
            }
	    else if (nfv == 3) {
	      strcpy(element_block_types[nb],"TRIANGLE");
              MEnt_Set_AttVal(mf,block_type_att,TRI,0.0,NULL);
            }
	  }

          MSet_Add(element_blocks[nb],mf);
	  nb++;
	}
      }
    }

    /* Post process so that if the same element block has different
       types of elements, the block type will be labeled as NSIDED or
       NFACED */

    if (MESH_Num_Regions(mesh)) {
      for (i = 0; i < nb; i++) {
        MRType mrtype0, mrtype;
        MRegion_ptr mr;
        
        mr = MSet_Entry(element_blocks[i],0);
        mrtype0 = MR_ElementType(mr);
        
        int found = 0, idx = 0;
        while (!found && ((mr = MSet_Next_Entry(element_blocks[i],&idx)))) {
          if (MR_ElementType(mr) != mrtype0)
            found = 1;
        }
        
        if (found) { 
          /* found different element types in element block */
          /* label it as a block containing polyhedra       */
          
          strcpy(element_block_types[i],"NFACED");

          /* Attach attribute to each of its elements indicating that
             the element type according to Exodus II will be NSIDED */

          idx = 0;
          while ((mr = MSet_Next_Entry(element_blocks[i],&idx))) 
            MEnt_Set_AttVal(mr,block_type_att,POLYHED,0.0,NULL);
        }
      }
    }
    else if (MESH_Num_Faces(mesh)) {
      for (i = 0; i < nb; i++) {
        int nfv0;
        MFace_ptr mf;

        mf = MSet_Entry(element_blocks[i],0);
        nfv0 = MF_Num_Vertices(mf);

        int found = 0, idx = 0;
        while (!found && ((mf = MSet_Next_Entry(element_blocks[i],&idx)))) {
          if (MF_Num_Vertices(mf) != nfv0)
            found = 1;
        }

        if (found) {
          /* found different element types in element block */
          /* label it as a block containing polyhedra       */
          
          strcpy(element_block_types[i],"NSIDED");
          
          /* Attach attribute to each of its elements indicating that
             the element type according to Exodus II will be NSIDED */

          idx = 0;
          while ((mf = MSet_Next_Entry(element_blocks[i],&idx))) 
            MEnt_Set_AttVal(mf,block_type_att,POLYGON,0.0,NULL);
        }
      }
    }


#ifdef MSTK_HAVE_MPI


    // Gather over all processors


    int maxnum=0, maxnum1=0;
    if (comm) {
      MPI_Allreduce(&nb,&maxnum,1,MPI_INT,MPI_MAX,comm);
      
      int *ebids_array_loc = (int *) calloc(maxnum,sizeof(int));
      
      for (i = 0; i < nb; i++)
        ebids_array_loc[i] = element_block_ids[i];
      
      int *ebids_array_glob;
      ebids_array_glob = (int *) calloc(maxnum*numprocs,sizeof(int));
      
      MPI_Gather(ebids_array_loc,maxnum,MPI_INT,ebids_array_glob,maxnum,MPI_INT,0,comm);
      
      if (rank == 0) {
        
        /* Make a unique list of block IDs */
        
        maxnum1 = nb;
        
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
        
        *num_element_block_glob = maxnum1;
        
      } /* if (rank == 0) */

      /* Tell everyone how many element blocks there really are */
      
      MPI_Bcast(num_element_block_glob,1,MPI_INT,0,comm);
      
      /* Send everyone the IDs of these element blocks */
      
      *element_block_ids_glob = 
        (int *) malloc((*num_element_block_glob)*sizeof(int));
      if (rank == 0)
        for (i = 0; i < (*num_element_block_glob); i++)
          (*element_block_ids_glob)[i] = ebids_array_glob[i];
      
      MPI_Bcast(*element_block_ids_glob,*num_element_block_glob,MPI_INT,
                0,comm);
      
      /* Populate the global element block data on each processor */
      
      *element_blocks_glob = (MSet_ptr *) calloc((*num_element_block_glob),
                                                 sizeof(MSet_ptr));
      *element_block_types_glob = (char **) calloc((*num_element_block_glob),
                                                   sizeof(char *));
      for (i = 0; i < (*num_element_block_glob); i++)
        (*element_block_types_glob)[i] = (char *) malloc(16*sizeof(char));
      
      for (i = 0; i < (*num_element_block_glob); i++) {
        int found = 0;
        for (j = 0; j < nb; j++) {
          if (element_block_ids[j] == (*element_block_ids_glob)[i]) {

            /* This global element block ID is present on this partition.
               So, make sure the global element block on this partition
               refers to the correct list of elements in this block. */

            (*element_blocks_glob)[i] = element_blocks[j];
            strcpy((*element_block_types_glob)[i],element_block_types[j]);
            found = 1;
            break;
          }        
        }      
        if (!found) { 
          
          /* This global element block ID is NOT present on this
             partition.  Put in a dummy list - create a set whose name
             begins with the string TEMPORARY_element_block so that we
             know we should delete it when we are finished with it */

          char element_block_name[256];
          sprintf(element_block_name,"TEMPORARY_element_block_%-d",
                  (*element_block_ids_glob)[i]);
          if (MESH_Num_Regions(mesh)) {
            (*element_blocks_glob)[i] = MSet_New(mesh,element_block_name,MREGION);
            strcpy((*element_block_types_glob)[i],"HEX");
          }
          else if (MESH_Num_Faces(mesh)) {
            (*element_blocks_glob)[i] = MSet_New(mesh,element_block_name,MFACE);            
            strcpy((*element_block_types_glob)[i],"QUAD");
          }
          else {
            (*element_blocks_glob)[i] = MSet_New(mesh,element_block_name,MUNKNOWNTYPE);
            strcpy((*element_block_types_glob)[i],"UNKNOWN");
          }
        }
      }
      
      free(ebids_array_loc);
      free(ebids_array_glob);
    }
    else { /* No MPI communicator given - treat it as a serial run */

      *num_element_block_glob = nb;

      *element_block_ids_glob = (int *) malloc(nb*sizeof(int));
      memcpy(*element_block_ids_glob,element_block_ids,nb*sizeof(int));

      *element_blocks_glob = (MSet_ptr *) calloc(nb,sizeof(MSet_ptr));      
      memcpy(*element_blocks_glob,element_blocks,nb*sizeof(MSet_ptr));

      *element_block_types_glob = (char **) calloc(nb,sizeof(char *));
      for (i = 0; i < nb; i++) {
        (*element_block_types_glob)[i] = (char *) malloc(16*sizeof(char));
        strncpy((*element_block_types_glob)[i],element_block_types[i],16);
      }
    }

#else

    *num_element_block_glob = nb;
    
    *element_block_ids_glob = (int *) malloc(nb*sizeof(int));
    memcpy(*element_block_ids_glob,element_block_ids,nb*sizeof(int));
    
    *element_blocks_glob = (MSet_ptr *) calloc(nb,sizeof(MSet_ptr));      
    memcpy(*element_blocks_glob,element_blocks,nb*sizeof(MSet_ptr));
    
    *element_block_types_glob = (char **) malloc(nb*sizeof(char *));
    for (i = 0; i < nb; ++i) {
      (*element_block_types_glob)[i] = (char *) malloc(16*sizeof(char));
      strncpy((*element_block_types_glob)[i],element_block_types[i],16);
    }
#endif


  }



  /* Collect side sets (sets of boundary entities that are one
     dimension less than the dimension of elements in the mesh, e.g.,
     faces in solid mesh)
  */

 
  void MESH_Get_Side_Set_Info(Mesh_ptr mesh, 
                              int enable_geometry_sets,
                              int *num_side_set_glob, 
			      MSet_ptr **side_sets_glob, 
                              int **side_set_ids_glob,
                              MSTK_Comm comm) {
    MSet_ptr mset;
    MFace_ptr mf;
    MEdge_ptr me;
    int idx, idx2, i, j, k, found, dim;
    int nsideset, nr, nf, sid, nsalloc, maxnum, maxnum1;
    char mset_name[256];
    int meshdim=3;
    int rank=0, numprocs=1;
    double rval;
    void *pval;

#ifdef MSTK_HAVE_MPI
    if (comm) {
      MPI_Comm_size(comm,&numprocs);
      MPI_Comm_rank(comm,&rank);
    }
#endif


    nsideset = 0;
    nsalloc = 10;
    MSet_ptr *side_sets = (MSet_ptr *) malloc(nsalloc*sizeof(MSet_ptr));
    int *side_set_ids = (int *) malloc(nsalloc*sizeof(int));

    nr = MESH_Num_Regions(mesh);
    nf = MESH_Num_Faces(mesh);

    if (nr) {
      meshdim = 3;

      /* first check if there are any mesh sets with entity type
         MFACE.  If so, write out these mesh sets as
         sidesets. Otherwise, use geometric classification of boundary
         faces to form and write sidesets */

      idx = 0;
      while ((mset = MESH_Next_MSet(mesh,&idx))) {
        MSet_Name(mset,mset_name);
        dim = MSet_EntDim(mset);

        if (dim != MFACE) continue; 
        if (strncmp(mset_name,"sideset_",8) != 0) continue;

        sscanf(mset_name+8,"%d",&sid);

        if (nsideset == nsalloc) {
          nsalloc *= 2;
          side_sets = (MSet_ptr *) realloc(side_sets,nsalloc*sizeof(MSet_ptr));
          side_set_ids = (int *) realloc(side_set_ids,nsalloc*sizeof(int));
        }

        side_set_ids[nsideset] = sid;
        side_sets[nsideset] = mset;
        nsideset++;
      }
    
      if (enable_geometry_sets) { 

        /* create sidesets from groups of faces classified on the same
           geometric model face */

        idx = 0;
        while ((mf = MESH_Next_Face(mesh,&idx))) {
          if (MF_GEntDim(mf) == 3) continue; /* Internal face */
          
          int fowned;
#ifdef MSTK_HAVE_MPI
#ifdef MSTK_USE_MARKERS
          fowned = MEnt_IsMarked(mf,ownedmk);
#else
          MEnt_Get_AttVal(mf, ownedatt, &fowned, &rval, &pval);
#endif
          if (!fowned) continue; /* Face cnctd to only ghost elements */
#endif
          
          sid = MF_GEntID(mf);          
          found = 0;
          i = 0;
          while (!found && i < nsideset) {
            if (side_set_ids[i] == sid) {
              found = 1;
              MSet_Add(side_sets[i],mf);
              break;
            }
            i++;
          }
          
          if (!found) {
            if (nsideset == nsalloc) {
              nsalloc *= 2;
              side_sets = (MSet_ptr *) realloc(side_sets,nsalloc*sizeof(MSet_ptr));
              side_set_ids = (int *) realloc(side_set_ids,nsalloc*sizeof(int));
            }
            
            /* create a sideset whose name starts with the string
               TEMPORARY so that we know to delete it when we are
               finished with it */

            char sidesetname[256];
            sprintf(sidesetname,"TEMPORARY_sideset_%-d",sid);
            side_sets[nsideset] = MSet_New(mesh,sidesetname,MFACE);
            MSet_Add(side_sets[nsideset],mf);
            side_set_ids[nsideset] = sid;
            nsideset++;
          }
        }

      }

    }
    else if (nf) {

      meshdim = 2;

      /* first check if there are any mesh sets with entity type
         MEDGE. If so, write out these mesh sets as
         sidesets. Otherwise, use geometric classification of boundary
         faces to form and write sidesets */

      idx = 0;
      while ((mset = MESH_Next_MSet(mesh,&idx))) {
        MSet_Name(mset,mset_name);
        dim = MSet_EntDim(mset);

        if (dim != MEDGE) continue; 
        if (strncmp(mset_name,"sideset_",8) != 0) continue;

        sscanf(mset_name+8,"%d",&sid);

        if (nsideset == nsalloc) {
          nsalloc *= 2;
          side_sets = (MSet_ptr *) realloc(side_sets,nsalloc*sizeof(MSet_ptr));
          side_set_ids = (int *) realloc(side_set_ids,nsalloc*sizeof(int));
        }

        side_set_ids[nsideset] = sid;
        side_sets[nsideset] = mset;
        nsideset++;
      }

      if (enable_geometry_sets) {

        /* create sidesets from sets of edges classified on the same
           geometric model edge */

        while ((me = MESH_Next_Edge(mesh,&idx))) {
          if (ME_GEntDim(me) != 1) continue;
          
          int eowned;

#ifdef MSTK_HAVE_MPI
#ifdef MSTK_USE_MARKERS
          eowned = MEnt_IsMarked(me,ownedmk);
#else
          MEnt_Get_AttVal(me, ownedatt, &eowned, &rval, &pval);
#endif
          if (!eowned) continue; /* Edge cnctd to only ghost elements */
#endif
          
          sid = ME_GEntID(me);
          
          found = 0;
          i = 0;
          while (!found && i < nsideset) {
            if (side_set_ids[i] == sid) {
              found = 1;
              MSet_Add(side_sets[i],me);
              break;
            }
            i++;
          }
          
          if (!found) {
            if (nsideset == nsalloc) {
              nsalloc *= 2;
              side_sets = (MSet_ptr *) realloc(side_sets,nsalloc*sizeof(MSet_ptr));
              side_set_ids = (int *) realloc(side_set_ids,nsalloc*sizeof(int));
            }

            /* create a sideset whose name starts with the string
               TEMPORARY so that we know to delete it when we are
               finished with it */

            char sidesetname[256];
            sprintf(sidesetname,"TEMPORARY_sideset_%-d",sid);
            side_sets[nsideset] = MSet_New(mesh,sidesetname,MEDGE);
            MSet_Add(side_sets[nsideset],me);
            side_set_ids[nsideset] = sid;
            nsideset++;
          }
        }
      }
    }


    /* Reconcile and aggregate over all processors */
      
#ifdef MSTK_HAVE_MPI
    
    if (comm) {

      maxnum=0, maxnum1=0;
      MPI_Allreduce(&nsideset,&maxnum,1,MPI_INT,MPI_MAX,comm);
      
      int *ssids_array_loc = (int *) calloc(maxnum,sizeof(int));
      
      for (i = 0; i < nsideset; i++)
        ssids_array_loc[i] = side_set_ids[i];
      
      int *ssids_array_glob = (int *) calloc(maxnum*numprocs,sizeof(int));
      
      MPI_Gather(ssids_array_loc,maxnum,MPI_INT,ssids_array_glob,maxnum,
                 MPI_INT,0,comm);
      
      if (rank == 0) {
        
        /* Make a unique list of block IDs */
        
        maxnum1 = nsideset;
        
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
        
        *num_side_set_glob = maxnum1;
        
      } /* if (rank == 0) */
      
      /* Tell everyone how many sidesets there really are globally */
      
      MPI_Bcast(num_side_set_glob,1,MPI_INT,0,comm);
      
      /* Send everyone the IDs of these sidesets */
      
      *side_set_ids_glob = (int *) malloc((*num_side_set_glob)*sizeof(int));
      if (rank == 0)
        for (i = 0; i < (*num_side_set_glob); i++)
          (*side_set_ids_glob)[i] = ssids_array_glob[i];
      
      MPI_Bcast(*side_set_ids_glob,*num_side_set_glob,MPI_INT,0,comm);
      
      /* Populate the global sideset data on each processor */
      
      *side_sets_glob = (MSet_ptr *) calloc((*num_side_set_glob),
                                            sizeof(MSet_ptr));
      
      for (i = 0; i < (*num_side_set_glob); i++) {
        int found = 0;
        for (j = 0; j < nsideset; j++) {
          if (side_set_ids[j] == (*side_set_ids_glob)[i]) {
            /* Found a global side set ID on local partition */

            (*side_sets_glob)[i] = side_sets[j];
            found = 1;
            break;
          }
        }      
        if (!found) {
          /* Sideset in global partition is not present on local partition */
          /* Put in a dummy set - create a sideset whose whose name
             starts with the string TEMPORARY so that we know to
             delete it when we are finished with it */

          char sidesetname[256];
          sprintf(sidesetname,"TEMPORARY_sideset_%-d",(*side_set_ids_glob)[i]);
          if (meshdim == 2)
            (*side_sets_glob)[i] = MSet_New(mesh,sidesetname,MEDGE);
          else if (meshdim == 3)
            (*side_sets_glob)[i] = MSet_New(mesh,sidesetname,MFACE);
        }
      }
      
      free(ssids_array_loc);
      free(ssids_array_glob);

    }
    else { /* No communicator given - treat it as a serial run */

      *num_side_set_glob = nsideset;
      *side_set_ids_glob = (int *) malloc(nsideset*sizeof(int));
      memcpy(*side_set_ids_glob,side_set_ids,nsideset*sizeof(int));

      *side_sets_glob = (MSet_ptr *) malloc(nsideset*sizeof(MSet_ptr));
      memcpy(*side_sets_glob,side_sets,nsideset*sizeof(MSet_ptr));

    }

#else

    *num_side_set_glob = nsideset;

    *side_set_ids_glob = (int *) malloc(nsideset*sizeof(int));
    memcpy(*side_set_ids_glob,side_set_ids,nsideset*sizeof(int));

    *side_sets_glob = (MSet_ptr *) malloc(nsideset*sizeof(MSet_ptr));
    memcpy(*side_sets_glob,side_sets,nsideset*sizeof(MSet_ptr));

#endif


  }
  

 
  void MESH_Get_Node_Set_Info(Mesh_ptr mesh, int enable_geometry_sets,
                              int *num_node_set_glob, 
			      MSet_ptr **node_sets_glob, 
                              int **node_set_ids_glob,
                              MSTK_Comm comm) {
    MVertex_ptr mv;
    int idx, idx2, i, j, k, found, nalloc;
    int nnodeset, nid, maxnum, maxnum1;
    int rank=0, numprocs=1;
    MSet_ptr mset;
    double rval;
    void *pval;

#ifdef MSTK_HAVE_MPI
    if (comm) {
      MPI_Comm_size(comm,&numprocs);
      MPI_Comm_rank(comm,&rank);
    }
#endif


    nnodeset = 0;
    nalloc = 10;
    MSet_ptr *node_sets = (MSet_ptr *) malloc(nalloc*sizeof(MSet_ptr));
    int *node_set_ids = (int *) malloc(nalloc*sizeof(int));



    idx = 0;
    while ((mset = MESH_Next_MSet(mesh,&idx))) {
      char mset_name[256];

      MSet_Name(mset,mset_name);
      int dim = MSet_EntDim(mset);
      
      if (dim != MVERTEX) continue; 
      if (strncmp(mset_name,"nodeset_",8) != 0) continue;
      
      sscanf(mset_name+8,"%d",&nid);
      
      if (nnodeset == nalloc) {
        nalloc *= 2;
        node_sets = (MSet_ptr *) realloc(node_sets,nalloc*sizeof(MSet_ptr));
        node_set_ids = (int *) realloc(node_set_ids,nalloc*sizeof(int));
      }
      
      node_set_ids[nnodeset] = nid;
      node_sets[nnodeset] = mset;
      nnodeset++;
    }
    
    if (enable_geometry_sets) { 
      /* If requested, create node sets formed by vertices classified
         on model faces, edges and vertices */

      idx = 0;
      while ((mv = MESH_Next_Vertex(mesh,&idx))) {
        if (MV_GEntDim(mv) >= 2) continue;
        
        int vowned;

#ifdef MSTK_HAVE_MPI
#ifdef MSTK_USE_MARKERS
        vowned = MEnt_IsMarked(mv,ownedmk);
#else
        MEnt_Get_AttVal(mv, ownedatt, &vowned, &rval, &pval);
#endif
        if (!vowned) continue; /* Vtx cnctd to only ghost elements */
#endif
        
        nid = MV_GEntDim(mv)*10000 + MV_GEntID(mv);
        
        found = 0;
        i = 0;
        while (!found && i < nnodeset) {
          if (node_set_ids[i] == nid) {
            found = 1;
            MSet_Add(node_sets[i],mv);
            break;
          }
          i++;
        }
        
        if (!found) {
          if (nnodeset == nalloc) {
            nalloc *= 2;
            node_sets = (MSet_ptr *) realloc(node_sets,nalloc*sizeof(MSet_ptr));
            node_set_ids = (int *) realloc(node_set_ids,nalloc*sizeof(int));
          }
          
          /* create a nodeset whose name starts with the string
             TEMPORARY so that we know to delete it when we are finished
             with it */
          
          char nodesetname[256];
          sprintf(nodesetname,"TEMPORARY_nodeset_%-d",nid);
          node_sets[nnodeset] = MSet_New(mesh,nodesetname,MVERTEX);
          MSet_Add(node_sets[nnodeset],mv);
          node_set_ids[nnodeset] = nid;
          nnodeset++;
        }
      }
    }

#ifdef MSTK_HAVE_MPI

    if (comm) {

      maxnum=0, maxnum1=0;
      MPI_Allreduce(&nnodeset,&maxnum,1,MPI_INT,MPI_MAX,comm);

      int *nsids_array_loc = (int *) calloc(maxnum,sizeof(int));

      for (i = 0; i < nnodeset; i++)
        nsids_array_loc[i] = node_set_ids[i];

      int *nsids_array_glob = (int *) calloc(maxnum*numprocs,sizeof(int));
    
      MPI_Gather(nsids_array_loc,maxnum,MPI_INT,nsids_array_glob,maxnum,
                 MPI_INT,0,comm);

      if (rank == 0) {

        /* Make a unique list of block IDs */

        maxnum1 = nnodeset;

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

        *num_node_set_glob = maxnum1;

      } /* if (rank == 0) */

      /* Tell everyone how many nodesets there really are */

      MPI_Bcast(num_node_set_glob,1,MPI_INT,0,comm);

      /* Send everyone the IDs of these nodesets */

      *node_set_ids_glob = (int *) malloc((*num_node_set_glob)*sizeof(int));
      if (rank == 0)
        for (i = 0; i < *num_node_set_glob; i++)
          (*node_set_ids_glob)[i] = nsids_array_glob[i];
  

      MPI_Bcast(*node_set_ids_glob,*num_node_set_glob,MPI_INT,0,comm);

      /* Populate the global nodeset data on each processor */

      *node_sets_glob = (MSet_ptr *) calloc((*num_node_set_glob),
                                            sizeof(MSet_ptr));

      for (i = 0; i < *num_node_set_glob; i++) {
        int found = 0;
        for (j = 0; j < nnodeset; j++) {
          if (node_set_ids[j] == (*node_set_ids_glob)[i]) {
            /* Global node set is present on this local partition */

            (*node_sets_glob)[i] = node_sets[j];
            found = 1;
            break;
          }
        }      
        if (!found) {
          /* Global node set is NOT present on this local partition */
          /* Put in an dummy nodeset - create a nodeset whose name
           starts with the string TEMPORARY so that we know to delete
           it when we are finished with it */

          char nodesetname[256];
          sprintf(nodesetname,"TEMPORARY_nodeset_%-d",(*node_set_ids_glob)[i]);
          (*node_sets_glob)[i] = MSet_New(mesh,nodesetname,MVERTEX); 
        }
      }

      free(nsids_array_loc);
      free(nsids_array_glob);

    }
    else { /* No communicator given - treat it as a serial run */

      *num_node_set_glob = nnodeset;
      
      *node_set_ids_glob = (int *) malloc(nnodeset*sizeof(int));
      memcpy(*node_set_ids_glob,node_set_ids,nnodeset*sizeof(int));

      *node_sets_glob = (MSet_ptr *) malloc(nnodeset*sizeof(MSet_ptr));
      memcpy(*node_sets_glob,node_sets,nnodeset*sizeof(MSet_ptr));
      
    }

#else

    *num_node_set_glob = nnodeset;

    *node_set_ids_glob = (int *) malloc(nnodeset*sizeof(int));
    memcpy(*node_set_ids_glob,node_set_ids,nnodeset*sizeof(int));

    *node_sets_glob = (MSet_ptr *) malloc(nnodeset*sizeof(MSet_ptr));
    memcpy(*node_sets_glob,node_sets,nnodeset*sizeof(MSet_ptr));

#endif


  }



  /* Collect element sets - these are done differently from node sets or 
     side sets - we convert explicitly defined sets of elements rather
     than use classification info to define the element sets */
 
  void MESH_Get_Element_Set_Info(Mesh_ptr mesh, int *num_element_set_glob, 
			      MSet_ptr **element_sets_glob, 
                              int **element_set_ids_glob,
                              MSTK_Comm comm) {
    MSet_ptr mset;
    int idx, i, j, k, nalloc;
    int maxnum, maxnum1, nelemset;
    int rank=0, numprocs=1;

#ifdef MSTK_HAVE_MPI
    if (comm) {
      MPI_Comm_size(comm,&numprocs);
      MPI_Comm_rank(comm,&rank);
    }
#endif


    nelemset = 0;
    nalloc = 10;
    MSet_ptr *element_sets = (MSet_ptr *) malloc(nalloc*sizeof(MSet_ptr));
    int *element_set_ids = (int *) malloc(nalloc*sizeof(int));

    int meshdim;
    if (MESH_Num_Regions(mesh))
      meshdim = 3;
    else if (MESH_Num_Faces(mesh))
      meshdim = 2;
    else if (MESH_Num_Edges(mesh))
      meshdim = 1;
    else
      meshdim = 0;

    int setid = 1;
    idx = 0;
    while ((mset = MESH_Next_MSet(mesh,&idx))) {
      int setdim = MSet_EntDim(mset);
      if (meshdim == setdim) {      
        char setname[256];
        MSet_Name(mset,setname);

        if (strncmp(setname,"elemset_",8) != 0) continue;

        if (nelemset == nalloc) {
          nalloc *= 2;
          element_set_ids = (int *) realloc(element_set_ids,nalloc*sizeof(int));
          element_sets = (MSet_ptr *) realloc(element_sets,nalloc*sizeof(MSet_ptr));
        }

        sscanf(setname+8,"%d",&(element_set_ids[nelemset]));
        element_sets[nelemset] = mset;
        setid++;
        nelemset++;
      }
    }
  

#ifdef MSTK_HAVE_MPI

    if (comm) {

      maxnum=0, maxnum1=0;
      MPI_Allreduce(&nelemset,&maxnum,1,MPI_INT,MPI_MAX,comm);

      int *esids_array_loc = (int *) calloc(maxnum,sizeof(int));

      for (i = 0; i < nelemset; i++)
        esids_array_loc[i] = element_set_ids[i];

      int *esids_array_glob = (int *) calloc(maxnum*numprocs,sizeof(int));
    
      MPI_Gather(esids_array_loc,maxnum,MPI_INT,esids_array_glob,maxnum,
                 MPI_INT,0,comm);

      if (rank == 0) {

        /* Make a unique list of element set IDs across all processors */

        maxnum1 = nelemset;

        for (i = 1; i < numprocs; i++) {

          int offset = maxnum*i;
          for (j = 0; j < maxnum; j++) {
            if (esids_array_glob[offset+j] == 0) continue;

            /* Is this element set ID in the list of global element
               set IDs being collected at the head of the
               esids_array_glob? */

            int found = 0;
            for (k = 0; k < maxnum1; k++) {
              if (esids_array_glob[offset+j] == esids_array_glob[k]) {
                found = 1;
                break;
              }
            }

            if (!found) {
            
              /* Found element set ID on processor i that is not in
                 the global list. Put it in. */

              esids_array_glob[maxnum1] = esids_array_glob[offset+j];
              maxnum1++;

            }
          }
        } /* for (i = 1; ....) */

        *num_element_set_glob = maxnum1;

      } /* if (rank == 0) */

      /* Tell everyone how many element sets there really are */

      MPI_Bcast(num_element_set_glob,1,MPI_INT,0,comm);

      /* Send everyone the IDs of these element sets */

      *element_set_ids_glob = (int *) malloc((*num_element_set_glob)*sizeof(int));
      if (rank == 0)
        for (i = 0; i < *num_element_set_glob; i++)
          (*element_set_ids_glob)[i] = esids_array_glob[i];
  

      MPI_Bcast(*element_set_ids_glob,*num_element_set_glob,MPI_INT,0,comm);

      /* Populate the global element set data on each processor */

      *element_sets_glob = (MSet_ptr *) calloc((*num_element_set_glob),
                                               sizeof(MSet_ptr));

      for (i = 0; i < *num_element_set_glob; i++) {
        int found = 0;
        for (j = 0; j < nelemset; j++) {
          if (element_set_ids[j] == (*element_set_ids_glob)[i]) {
            /* Global element set is present on this local partition */

            (*element_sets_glob)[i] = element_sets[j];
            found = 1;
            break;
          }
        }      
        if (!found) {
          /* Global element set is NOT present on this local partition */
          /* Put in an dummy set - create a element set whose name
           starts with the string TEMPORARY so that we know to delete
           it when we are finished with it */

          char elementsetname[256];
          sprintf(elementsetname,"TEMPORARY_elementset_%-d",
                  (*element_set_ids_glob)[i]);
          if (meshdim == 2)
            (*element_sets_glob)[i] = MSet_New(mesh,elementsetname,MFACE); 
          else
            (*element_sets_glob)[i] = MSet_New(mesh,elementsetname,MREGION); 
        }
      }

      free(esids_array_loc);
      free(esids_array_glob);

    }
    else { /* No communicator given - treat it as a serial run */

      *num_element_set_glob = nelemset;
      
      *element_set_ids_glob = (int *) malloc(nelemset*sizeof(int));
      memcpy(*element_set_ids_glob,element_set_ids,nelemset*sizeof(int));

      *element_sets_glob = (MSet_ptr *) malloc(nelemset*sizeof(MSet_ptr));
      memcpy(*element_sets_glob,element_sets,nelemset*sizeof(MSet_ptr));
      
    }

#else

    *num_element_set_glob = nelemset;

    *element_set_ids_glob = (int *) malloc(nelemset*sizeof(int));
    memcpy(*element_set_ids_glob,element_set_ids,nelemset*sizeof(int));

    *element_sets_glob = (MSet_ptr *) malloc(nelemset*sizeof(MSet_ptr));
    memcpy(*element_sets_glob,element_sets,nelemset*sizeof(MSet_ptr));

#endif


  }




  /* Discover the attributes on elements, nodes and sidesets */

  void MESH_Get_Attribute_Info(Mesh_ptr mesh, 
                               int *num_element_atts_glob, 
                               char ***element_att_names_glob, 
                               int *num_node_atts_glob, 
                               char ***node_att_names_glob, 
                               int *num_sideset_atts_glob, 
                               char ***sideset_att_names_glob, 
                               MSTK_Comm comm) {

    int MAXLEN=256;
    int i, j, k, idx, natt_loc, nalloc;
    int rank=0, numprocs=1;
    MAttrib_ptr att;
    int element_dim, face_dim, *att_dim_loc;
    char *att_names_loc; 

#ifdef MSTK_HAVE_MPI
    if (comm) {
      MPI_Comm_size(comm,&numprocs);
      MPI_Comm_rank(comm,&rank);
    }
#endif

    if (MESH_Num_Regions(mesh)) {
      element_dim = 3;
      face_dim = 2;
    }
    else if (MESH_Num_Faces(mesh)) {
      element_dim = 2;
      face_dim = 1;
    }
    else {
      element_dim = 1;
      face_dim = 0;  /* meaningless but here for completeness */
    }
    
    nalloc = MESH_Num_Attribs(mesh);

    // names arranged in 1-dimensional array

    att_names_loc = (char *) malloc(nalloc*MAXLEN*sizeof(char));
    att_dim_loc = (int *) malloc(nalloc*sizeof(int));
    natt_loc = 0;
 
    idx = 0;
    while ((att = MESH_Next_Attrib(mesh,&idx))) {
      
      // Can only export real valued (scalar or vector) attributes

      if (MAttrib_Get_Type(att) == INT || MAttrib_Get_Type(att) == POINTER) continue;
      
      MAttrib_Get_Name(att, &(att_names_loc[natt_loc*MAXLEN]));
      att_dim_loc[natt_loc] = MAttrib_Get_EntDim(att);
      ++natt_loc;

    }
      

    int natts_glob;
    char *att_names_glob=NULL; 
    int *att_dim_glob=NULL;

#ifdef MSTK_HAVE_MPI


    // Gather over all processors

    int maxnum=0, maxnum1=0;
    if (comm) {
      int *natt_all = (int *) malloc(numprocs*sizeof(int));

      MPI_Gather(&natt_loc,1,MPI_INT,natt_all,1,MPI_INT,0,comm);

      if (rank == 0) {
        maxnum = natt_loc;
        for (i = 1; i < numprocs; ++i) 
          maxnum = (natt_all[i] > maxnum) ? natt_all[i] : maxnum;
      }

      MPI_Bcast(&maxnum,1,MPI_INT,0,comm);

      if (maxnum) {

        att_names_glob = (char *) malloc(maxnum*numprocs*MAXLEN*sizeof(char));
      
        MPI_Gather(att_names_loc, maxnum*MAXLEN, MPI_CHAR, att_names_glob,
                   maxnum*MAXLEN, MPI_CHAR, 0, comm);
        
        
        att_dim_glob = (int *) malloc(maxnum*numprocs*sizeof(int));
        
        MPI_Gather(att_dim_loc, maxnum, MPI_INT, att_dim_glob, maxnum, MPI_INT, 
                   0, comm);
        
        if (rank == 0) {
          
          /* Make a unique list of att names and the dimension of the
             entities they apply to */
          
          maxnum1 = natt_loc;
          
          for (i = 1; i < numprocs; i++) {
            
            int offset = maxnum*i;
            for (j = 0; j < natt_all[i]; j++) {
              
              /* Is this attribute name in the list of attributes being
                 collected at the head of the global attribute list on P0? */
              
              int found = 0;
              for (k = 0; k < maxnum1; k++) {
                if (strcmp(&(att_names_glob[(offset+j)*MAXLEN]),
                           &(att_names_glob[k*MAXLEN])) == 0) {
                  found = 1;
                  break;
                }
              }
              
              if (!found) {
                
                /* Found element block ID on processor i that is not in
                   the global list. Put it in. */
                
                att_dim_glob[maxnum1] = att_dim_glob[offset+j];
                strcpy(&(att_names_glob[maxnum1*MAXLEN]),
                       &(att_names_glob[(offset+j)*MAXLEN]));
                ++maxnum1;
                
              }
            }
          } /* for (i = 1; ....) */
        
          natts_glob = maxnum1;
        
        } /* if (rank == 0) */
        
        
        /* Tell everyone how many attributes there really are */
        
        MPI_Bcast(&natts_glob,1,MPI_INT,0,comm);
        
        /* Send everyone the names of these attributes */
        
        MPI_Bcast(att_names_glob,natts_glob*MAXLEN,MPI_CHAR,0,comm);
        
        /* Send everyone the dimension of the entities these attributes are on */
        
        MPI_Bcast(att_dim_glob,natts_glob,MPI_INT,0,comm);
        
        free(natt_all);
      }
      else
        natts_glob = 0;
    }
    else { /* No MPI communicator given - treat it as a serial run */

      natts_glob = natt_loc;

      if (natt_loc) {
        att_names_glob = (char *) malloc(natt_loc*MAXLEN*sizeof(char));
        memcpy(att_names_glob,att_names_loc,natt_loc*MAXLEN*sizeof(char));
        
        att_dim_glob = (int *) malloc(natt_loc*sizeof(List_ptr));      
        memcpy(att_dim_glob,att_dim_loc,natt_loc*sizeof(int));
      }
    }

#else

    natts_glob = natt_loc;
    
    if (natt_loc) {
      att_names_glob = (char *) malloc(natt_loc*MAXLEN*sizeof(char));
      memcpy(att_names_glob,att_names_loc,natt_loc*MAXLEN*sizeof(char));
      
      att_dim_glob = (int *) malloc(natt_loc*sizeof(List_ptr));      
      memcpy(att_dim_glob,att_dim_loc,natt_loc*sizeof(int));
    }
#endif


    /* Now count up the number of element, node and sideset attributes
       to send out */

    *num_element_atts_glob = 0;
    *num_node_atts_glob = 0;
    *num_sideset_atts_glob = 0;
    
    for (i = 0; i < natts_glob; ++i) {
      if (att_dim_glob[i] == element_dim || att_dim_glob[i] == MALLTYPE)
        ++(*num_element_atts_glob);
      if (att_dim_glob[i] == face_dim || att_dim_glob[i] == MALLTYPE)
        ++(*num_sideset_atts_glob);
      if (att_dim_glob[i] == 0 || att_dim_glob[i] == MALLTYPE)
        ++(*num_node_atts_glob);      
    }

    if ((*num_element_atts_glob)) { 
      *element_att_names_glob = (char **) malloc((*num_element_atts_glob)*sizeof(char *));
      for (i = 0; i < *num_element_atts_glob; ++i)
        (*element_att_names_glob)[i] = (char *) malloc(MAXLEN*sizeof(char));
    }

    if ((*num_node_atts_glob)) {
      *node_att_names_glob = (char **) malloc((*num_node_atts_glob)*sizeof(char *));
      for (i = 0; i < *num_node_atts_glob; ++i)
        (*node_att_names_glob)[i] = (char *) malloc(MAXLEN*sizeof(char));
    }

    if ((*num_sideset_atts_glob)) {
      *sideset_att_names_glob = (char **) malloc((*num_sideset_atts_glob)*sizeof(char *));
      for (i = 0; i < *num_sideset_atts_glob; ++i)
        (*sideset_att_names_glob)[i] = (char *) malloc(MAXLEN*sizeof(char));
    }

    /* Copy into output variables. Point to note is that att_names_glob is a
       one-dimensional array containing all the attribute names while the 
       output arrays are arrays of string arrays (two dimensional) */

    int n1=0, n2=0, n3=0;
    for (i = 0; i < natts_glob; ++i) {
      if (att_dim_glob[i] == element_dim || att_dim_glob[i] == MALLTYPE) {
        strncpy((*element_att_names_glob)[n1],&(att_names_glob[i*MAXLEN]),MAXLEN);
        ++n1;
      }
      if (att_dim_glob[i] == 0 || att_dim_glob[i] == MALLTYPE) {
        strncpy((*node_att_names_glob)[n2],&(att_names_glob[i*MAXLEN]),MAXLEN);
        ++n2;
      }
      if (att_dim_glob[i] == face_dim || att_dim_glob[i] == MALLTYPE) {
        strncpy((*sideset_att_names_glob)[n3],&(att_names_glob[i*MAXLEN]),MAXLEN);
        ++n3;
      }
    }    

    if (nalloc) {
      free(att_names_loc);
      free(att_dim_loc);
    }
    if (natts_glob) {
      free(att_names_glob);
      free(att_dim_glob);
    }

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

