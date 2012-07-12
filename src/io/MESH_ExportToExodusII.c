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
				   char ***element_block_types);

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
     Global Parameter information:                                     
     num_node, num_elem, num_elem_blk, num_node_set, num_side_set  
     Quality assurance information(optional)                       
     The nodes(vertices) list and elements list  
     For nodes, there is only one attribute associate with them    
     opts[0] = 1 --- verbose (print stat information)
     opts[1] = 1 --- enable_set (sidesets, nodesets based on GEntID, GEntDim)
  */

  int MESH_ExportToExodusII(Mesh_ptr mesh, const char *filename, 
			    const int natt, const char **attnames, 
			    const int *opts) {

    int enable_set, verbose;
    int i, j, k, idx, idx2;
    int nv, ne, nf, nr;
    int num_element_block, num_side_set, num_node_set;
    int nvfblock, mkid, block_id, *connect, *nnpe;
    MVertex_ptr mv,  vertex;
    MEdge_ptr me;
    MFace_ptr mf;
    MRegion_ptr mr;
    RepType reptype;
    char MESH_rtype_str[5][3] = {"F1\0","F4\0","R1\0","R2\0","R4\0"};
    /*For pure 2D mesh, boundary_dim is 1, for surface mesh, it is 2 */
    int jv;
    double xyz_coord[3];
    int boundary_dim;
    reptype = MESH_RepType(mesh);
    int *element_block_ids, *side_set_ids, *node_set_ids;
    int num_face_block;
    char **element_block_types, block_name[256];
    List_ptr *element_blocks, *side_sets, *node_sets, face_block;
    List_ptr fverts, rverts;
    MAttrib_ptr fidatt;

    int mstk2exo_facemap[5][6]={{0,0,0,0,0,0},
				{4,1,2,3,0,0}, /* TET */
				{0,0,0,0,0,0}, /* PYRAMID, no support in exo*/
			        {4,5,1,2,3,0}, /* PRISM, must verify nums */
			        {5,6,1,2,3,4}};/* HEX */

#ifdef MSTK_HAVE_MPI
    int rank, numprocs;
    numprocs = MSTK_Comm_size();
    rank = MSTK_Comm_rank();
#endif    

    verbose = opts ? opts[0] : 0;
    enable_set = opts ? opts[1] : 1;

    if (verbose)
      fprintf(stdout,"\nMesh representation type is %s\n", 
	      MESH_rtype_str[reptype]);
    nv = MESH_Num_Vertices(mesh);
    ne = MESH_Num_Edges(mesh);
    nf = MESH_Num_Faces(mesh);
    nr = MESH_Num_Regions(mesh);
    
    if (nv == 0) {
      fprintf(stdout,"No vertices information \n");
      exit(-1);
    }
    
    int element_dim = -1, side_dim = -1;
    int num_element = 0, num_side = 0;
    if (nr) {
      element_dim = 3; side_dim = 2; boundary_dim = 2;
      num_element = nr; num_side = nf;
      if(verbose)
	fprintf(stdout,"\nThis is a 3D volume mesh with %d Vertices, %d Edges, %d Faces and %d Regions\n",\
		nv,ne,nf,nr);
    }
    else if(nf) {
      element_dim = 2; side_dim = 1; boundary_dim = 1;
      num_element = nf; num_side = ne;
      if(verbose)
      fprintf(stdout,"\nThis is a 3D surface mesh with %d Vertices, %d Edges and %d Faces\n", nv,ne,nf);
    }
    else if(ne) {
      element_dim = 1; side_dim = 0;
      num_element = ne; num_side = nv;
      if(verbose)
      fprintf(stdout,"\nThis is a 3D curve mesh with %d Vertices and %d Edges\n", nv,ne);
      exit(-1);
    }
    else if(nv) {
      element_dim = 0; side_dim = 0;
      num_element = nv; num_side = 0;
      if(verbose)
      fprintf(stdout,"\nThis is a 3D mesh with %d Vertices\n", nv);
      exit(-1);
    }
    else {
      fprintf(stdout,"\nThis is an empty mstk file\n");
      exit(-1);
    }  
    
    /* decide if it is 2d planar mesh or 3D surface mesh */

    if (element_dim == 2) {
      double xyz_coord0[3];
      vertex = MESH_Vertex(mesh,0);
      MV_Coords(vertex,xyz_coord0);

      idx = 0;
      while ((vertex = MESH_Next_Vertex(mesh,&idx))) {
	MV_Coords(vertex,xyz_coord);
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
    exoid = ex_create(filename,EX_CLOBBER,&CPU_word_size,&IO_word_size);
    if (exoid < 0) {
      fprintf(stderr, "after ex_create, error = %d\n", exoid);
      exit(-1);
    }


    
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


    MESH_Get_Element_Block_Info(mesh, &num_element_block, &element_blocks,
				 &element_block_ids, &element_block_types);
    
    if (verbose) {
      fprintf(stdout,"\nElemen1t block information:\n");
      fprintf(stdout,"total %d elements\n", num_element); 
      fprintf(stdout,"total %d element blocks:\n\n", num_element_block);
    }

    

    MESH_Get_Side_Set_Info(mesh, &num_side_set, &side_sets, &side_set_ids);


    MESH_Get_Node_Set_Info(mesh, &num_node_set, &node_sets, &node_set_ids);


    /* Put all faces of polyhedral elements in one block */

    face_block = List_New(10);
    num_face_block = 1;
    nvfblock = 0;
    mkid = MSTK_GetMarker();
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
	    if (!MEnt_IsMarked(rf,mkid)) {
	      List_Add(face_block,rf);
	      MEnt_Mark(rf,mkid);
	      MEnt_Set_AttVal(rf,fidatt,++j,0.0,NULL);
	      nvfblock += MF_Num_Vertices(rf);
	    }
	  }
	  List_Delete(rfaces);
	}

      }
    }

    List_Unmark(face_block,mkid);
    MSTK_FreeMarker(mkid);
   
    if (List_Num_Entries(face_block) == 0) {
      List_Delete(face_block);
      face_block = NULL;
      num_face_block = 0;
    }
    

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
      par.num_elem_blk = num_element_block;  
      par.num_node_sets = num_node_set;
      par.num_edge_sets = 0;
      par.num_face_sets = 0;
      par.num_side_sets = num_side_set;
      par.num_elem_sets = 0;
      par.num_node_maps = 0;
      par.num_edge_maps = 0;
      par.num_face_maps = 0;
      par.num_elem_maps = 0;
    }
    else {
      par.num_edge = 0;
      par.num_face = 0;
      par.num_face_blk = 0;
      par.num_edge_blk = 0;
      par.num_elem = nf;
      par.num_elem_blk = num_element_block;  
      par.num_node_sets = num_node_set;
      par.num_edge_sets = 0;
      par.num_face_sets = 0;
      par.num_side_sets = num_side_set;
      par.num_elem_sets = 0;
      par.num_node_maps = 0;
      par.num_edge_maps = 0;
      par.num_face_maps = 0;
      par.num_elem_maps = 0;
    }

    ex_put_init_ext(exoid,&par);
  
    if (verbose)
      fprintf(stdout,"Global parameters wrote into %s\n", filename);


    /* write coordinate values */

    double *xcoord, *ycoord, *zcoord, vxyz[3];
    xcoord = (double *) MSTK_malloc(nv*sizeof(double));
    ycoord = (double *) MSTK_malloc(nv*sizeof(double));
    zcoord = (double *) MSTK_malloc(nv*sizeof(double));
    for (jv = 0; jv < nv; jv++) {
      vertex = MESH_Vertex(mesh,jv);
      MV_Coords(vertex,vxyz);
      xcoord[jv] = vxyz[0];
      ycoord[jv] = vxyz[1];
      zcoord[jv] = vxyz[2];
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
    
    for (i = 0; i < num_element_block; i++) { 

      char block_name[256];

      block_id = element_block_ids[i];

      if (strncasecmp(element_block_types[i],"NFACED",6) == 0) {
	
	int nfrblock = 0;

	nnpe = (int *) MSTK_calloc(List_Num_Entries(element_blocks[i]),
				   sizeof(int));	
	idx = 0; j = 0;
	while ((mr = List_Next_Entry(element_blocks[i],&idx))) {	  
	  nnpe[j] = MR_Num_Faces(mr);
	  nfrblock += nnpe[j];
	  j++;
	}

	connect = (int *) MSTK_calloc(nfrblock,sizeof(int));
	idx = 0; j = 0;
	while ((mr = List_Next_Entry(element_blocks[i],&idx))) {
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
		     List_Num_Entries(element_blocks[i]),
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
      else if (strncasecmp(element_block_types[i],"NSIDED",6) == 0) {

	int nvfblock = 0;

	nnpe = (int *) MSTK_calloc(List_Num_Entries(element_blocks[i]),
				   sizeof(int));	
	idx = 0; j = 0;
	while ((mf = List_Next_Entry(element_blocks[i],&idx))) {	  
	  nnpe[j] = MF_Num_Vertices(mf);
	  nvfblock += nnpe[j];
	  j++;
	}

	connect = (int *) MSTK_calloc(nvfblock,sizeof(int));
	idx = 0; j = 0;
	while ((mf = List_Next_Entry(element_blocks[i],&idx))) {
	  List_ptr fverts = MF_Vertices(mf,1,0);
	  MVertex_ptr fv;
	  idx2 = 0;
	  while ((fv = List_Next_Entry(fverts,&idx2)))
	    connect[j++] = MV_ID(fv);
	  List_Delete(fverts);
	}


	sprintf(block_name,"POLYGON_BLOCK_%-d",block_id);
	
	ex_put_block(exoid, EX_ELEM_BLOCK, block_id, "NSIDED",
		     List_Num_Entries(element_blocks[i]),
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

	int nelem, nelnodes;

	sprintf(block_name,"BLOCK_%-d",block_id);
	nelem = List_Num_Entries(element_blocks[i]);
	
	if (MESH_Num_Regions(mesh)) {
	  if (strncasecmp(element_block_types[i],"TETRA",5) == 0) {
	    nelnodes = 4;
	    ex_put_elem_block(exoid, block_id, "TETRA", nelem, nelnodes, 1);
	  }
	  else if (strncasecmp(element_block_types[i],"WEDGE",5) == 0) {
	    nelnodes = 6;
	    ex_put_elem_block(exoid, block_id, "WEDGE", nelem, nelnodes, 1);
	  }
	  else if (strncasecmp(element_block_types[i],"HEX",3) == 0) {
	    nelnodes = 8;
	    ex_put_elem_block(exoid, block_id, "HEX", nelem, nelnodes, 1);
	  }
	  else
	    MSTK_Report("MESH_ExportToEXODUSII",
			"Element type unsupported by EXODUS II format",MSTK_FATAL);
	  
	  
	  connect = (int *) MSTK_malloc(nelnodes*nelem*sizeof(int));
	  
	  idx = 0; k = 0;
	  while ((mr = List_Next_Entry(element_blocks[i],&idx))) {
	    rverts = MR_Vertices(mr);
	    for (j = 0; j < nelnodes; j++) 
	      connect[nelnodes*k+j] = MV_ID(List_Entry(rverts,j));
	    List_Delete(rverts);	    
	    k++;
	  }
	  
	  ex_put_elem_conn(exoid, block_id, connect);

	  MSTK_free(connect);
	  
	}
	else if (MESH_Num_Faces(mesh)) {
	  if (strncasecmp(element_block_types[i],"TRIANGLE",8) == 0) {
	    nelnodes = 3;
	    ex_put_elem_block(exoid, block_id, "TRIANGLE", nelem, nelnodes, 1);
	  }
	  else if (strncasecmp(element_block_types[i],"QUAD",4) == 0) {
	    nelnodes = 4;
	    ex_put_elem_block(exoid, block_id, "QUAD", nelem, nelnodes, 1);
	  }
	  else
	    MSTK_Report("MESH_ExportToEXODUSII",
			"Element type unsupported by EXODUS II format",MSTK_FATAL);

	  connect = (int *) MSTK_malloc(nelnodes*nelem*sizeof(int));

	  idx = 0; k = 0;
	  while ((mf = List_Next_Entry(element_blocks[i],&idx))) {
	    fverts = MF_Vertices(mf,1,0);
	    for (j = 0; j < nelnodes; j++) 
	      connect[nelnodes*k+j] = MV_ID(List_Entry(fverts,j));	    
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

   
    for (i = 0; i < num_node_set; i++) {
      int nnodes = List_Num_Entries(node_sets[i]);
      ex_put_node_set_param(exoid, node_set_ids[i], nnodes, 0);

      int *node_list = (int *) MSTK_malloc(nnodes*sizeof(int));

      idx = 0; j = 0;
      while ((mv = List_Next_Entry(node_sets[i],&idx)))
	node_list[j++] = MV_ID(mv);

      ex_put_node_set(exoid, node_set_ids[i], node_list);

      MSTK_free(node_list);
    }

    if (verbose)
      fprintf(stdout,"Node sets information written into %s\n", filename);



    
    /* Write out the side set information */


    for (i = 0; i < num_side_set; i++) {

      int nsides = List_Num_Entries(side_sets[i]);
      ex_put_side_set_param(exoid, side_set_ids[i], nsides, 0);

      int *elem_list = (int *) MSTK_malloc(nsides*sizeof(int));
      int *side_list = (int *) MSTK_malloc(nsides*sizeof(int));

      if (nr) {
	idx = 0; j = 0;
	while ((mf = List_Next_Entry(side_sets[i],&idx))) {

	  List_ptr fregs = MF_Regions(mf);

	  if (!fregs || List_Num_Entries(fregs) == 0)
	    MSTK_Report("MESH_ExportToEXODUSII",
			"Standalone face with no regions in side set",
			MSTK_FATAL);
	  mr = List_Entry(fregs,0);

	  /* Since elements will be read back from the Exodus II file
	     element block by element block, we have to jump through
	     hoops to determine what the element number will be when
	     it is read */

	  int offset = 0;
	  int found = 0;
	  k = 0;
	  while (!found && k < num_element_block) {
	    int loc = List_Locate(element_blocks[k],mr);
	    if (loc == -1) {
	      /* This element block does not contain the element */
	      offset = List_Num_Entries(element_blocks[k]);
	      k++;
	    }
	    else {
	      elem_list[j] = offset+loc+1;
	      found = 1;
	    }
	  }
	  
          int lid = MF_LocalID_in_Region(mf,mr);
          MRType mrtype = MR_ElementType(mr);
          if (mrtype == TET || mrtype == PRISM || mrtype == HEX)
            side_list[j] = mstk2exo_facemap[mrtype][lid];
          else
            side_list[j] = lid+1;

	  List_Delete(fregs);
	  j++;
	}

	ex_put_side_set(exoid, side_set_ids[i], elem_list, side_list);

      }
      else {

	idx = 0; j = 0;
	while ((me = List_Next_Entry(side_sets[i],&idx))) {
	  List_ptr efaces = ME_Faces(me);

	  if (!efaces || List_Num_Entries(efaces) == 0)
	    MSTK_Report("MESH_ExportToEXODUSII",
			"Standalone edge with no faces in side set",
			MSTK_FATAL);
	  mf = List_Entry(efaces,0);

	  /* Since elements will be read back from the Exodus II file
	     element block by element block, we have to jump through
	     hoops to determine what the element number will be when
	     it is read */

	  int offset = 0;
	  int found = 0;
	  k = 0;
	  while (!found && k < num_element_block) {
	    int loc = List_Locate(element_blocks[k],mf);
	    if (loc == -1) {
	      /* This element block does not contain the element */
	      offset = List_Num_Entries(element_blocks[k]);
	      k++;
	    }
	    else {
	      elem_list[j] = offset+loc+1;
	      found = 1;
	    }
	  }
	  
          int lid = ME_LocalID_in_Face(me,mf);
          side_list[j] = lid+1;

	  List_Delete(efaces);
	  j++;
	}

	ex_put_side_set(exoid, side_set_ids[i], elem_list, side_list);

      }

      MSTK_free(elem_list);
      MSTK_free(side_list);

    }

    if (verbose)
      fprintf(stdout,"Side set information written into %s\n", filename);
    
    if (verbose)
      fprintf(stdout,"Node sets information written into %s\n", filename);







     
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
    return 1;

  }



  /* Categorize the elements into blocks based on their geometric
     entity ID and the element type */

  void MESH_Get_Element_Block_Info(Mesh_ptr mesh, int *num_element_block, List_ptr **element_blocks, int **element_block_ids, char ***element_block_types) {
    MEdge_ptr me;
    MFace_ptr mf;
    MRegion_ptr mr;
    int i, nb, nballoc, bid, found;

    nb = 0; nballoc=10;
    *element_block_ids = (int *) calloc(nballoc,sizeof(int));
    *element_block_types = (char **) malloc(nballoc*sizeof(char *));
    *element_blocks = (List_ptr *) calloc(nballoc,sizeof(List_ptr));

    int idx = 0;
    
    if (MESH_Num_Regions(mesh)) {
      while ((mr = MESH_Next_Region(mesh, &idx))) {
	int nrv = MR_Num_Vertices(mr);
	int nrf = MR_Num_Faces(mr);
        MRType mrtype = MR_ElementType(mr);

	if (mrtype == TET || mrtype == PRISM || mrtype == HEX)
	  bid = (MEnt_GEntID(mr)<<16) | nrv<<8 | nrf;
	else
	  bid = MEnt_GEntID(mr)<<16;
	
	found = 0; 
	i = 0;
	while (!found && i < nb) {
	  if (bid == (*element_block_ids)[i]) {
	    found = 1;
	    List_Add((*element_blocks)[i],mr);    
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

	  List_Add((*element_blocks)[nb],mr);
	  nb++;
	}
      }
    }
    else if (MESH_Num_Faces(mesh)) {

      while ((mf = MESH_Next_Face(mesh,&idx))) {

	int nfv = MF_Num_Vertices(mf);

	if (nfv == 3 || nfv == 4)
	  bid = (MEnt_GEntID(mf)<<16) | nfv;
	else
	  bid = MEnt_GEntID(mf)<<16;

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
	  if (nfv > 4)
	    strcpy((*element_block_types)[nb],"NSIDED");
	  else {
	    if (nfv == 4)
	      strcpy((*element_block_types)[nb],"QUAD");
	    else if (nfv == 3)
	      strcpy((*element_block_types)[nb],"TRIANGLE");
	  }
	  List_Add((*element_blocks)[nb],mf);
	  nb++;
	}
      }
    }
    else if (MESH_Num_Edges(mesh)) {
      while ((me = MESH_Next_Edge(mesh, &idx))) {
	bid = (MEnt_GEntID(me)<<16) | 2;

	found = 0; i = 0;
	while (!found && i < nb) {
	  if (bid == (*element_block_ids)[i]) {
	    found = 1;
	    List_Add((*element_blocks)[i],me);    
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
	  strcpy((*element_block_types)[nb],"BEAM");
	  List_Add((*element_blocks)[nb],me);
	  nb++;
	}
      }
    }

    *num_element_block = nb;
  }	


 
  void MESH_Get_Side_Set_Info(Mesh_ptr mesh, int *num_side_set, 
			      List_ptr **side_sets, int **side_set_ids) {
    MFace_ptr mf;
    MEdge_ptr me;
    int idx = 0, i = 0, found;
    int ns, nsalloc, nr, nf, ne, sid;


    ns = 0;
    nsalloc = 10;
    *side_sets = (List_ptr *) MSTK_malloc(nsalloc*sizeof(List_ptr));
    *side_set_ids = (int *) MSTK_malloc(nsalloc*sizeof(int));

    nr = MESH_Num_Regions(mesh);
    nf = MESH_Num_Faces(mesh);
    ne = MESH_Num_Edges(mesh);

    if (nr) {
      while ((mf = MESH_Next_Face(mesh,&idx))) {
	if (MF_GEntDim(mf) == 3) continue;

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
    int idx, idx2, idx3, idx4, i, found;
    int nn, nnalloc, nid, nid2, mkid1, mkid2;
    List_ptr fverts, fedges, vedges, efaces, elist, flist, vlist;


    nn = 0;
    nnalloc = 10;
    *node_sets = (List_ptr *) MSTK_malloc(nnalloc*sizeof(List_ptr));
    *node_set_ids = (int *) MSTK_malloc(nnalloc*sizeof(int));



    /* Node sets formed by vertices on geometric model vertices */

    idx = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx))) {
      if (MV_GEntDim(mv) != 0) continue;
      
      nid = MV_GEntID(mv)<<4;
      
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

      nid = ME_GEntID(me)<<4 | 1;

      /* Make a list of the mesh edges on this geometric model edge */
      
      idx2 = 0;
      while ((me2 = List_Next_Entry(elist,&idx2))) {
	for (i = 0; i < 2; i++) {
	  mv = ME_Vertex(me2,i);
	  vedges = MV_Edges(mv);
	  
	  idx3 = 0;
	  while ((ve = List_Next_Entry(vedges,&idx3))) {

	    nid2 = ME_GEntID(ve)<<4 | ME_GEntDim(ve);

	    if (nid != nid2) continue; /* Not on the same model edge */

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
	  if (!MEnt_IsMarked(mv,mkid2)) {
	    MEnt_Mark(mv,mkid2);
	    List_Add(vlist,mv);
	  }
	}
      }
      List_Unmark(vlist,mkid2);

      List_Unmark(elist,mkid2);
      List_Delete(elist);


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





    /* Node sets formed by vertices on the closure of geometric model faces */

    if (MESH_Num_Regions(mesh)) {

      idx = 0;
      while ((mf = MESH_Next_Face(mesh,&idx))) {
	if (MEnt_IsMarked(mf,mkid1)) continue;
	
	if (MF_GEntDim(mf) != 2) continue;

      
	/* Found a mesh face on an unprocessed geometric model face */
      
	flist = List_New(10);
	
	MEnt_Mark(mf,mkid1);
	MEnt_Mark(mf,mkid2);
	List_Add(flist,mf);
	
	nid = MF_GEntID(mf)<<4 | 2;
	
	
	/* Collect all mesh faces on this geometric model face */
	
	idx2 = 0;
	while ((mf2 = List_Next_Entry(flist,&idx2))) {
	  fedges = MF_Edges(mf2,1,0);
	  idx3 = 0;
	  while ((me = List_Next_Entry(fedges,&idx3))) {
	    efaces = ME_Faces(me);
	    
	    idx4 = 0;
	    while ((ef = List_Next_Entry(efaces,&idx4))) {
	      
	      nid2 = MF_GEntID(ef)<<4 | MF_GEntDim(ef);
	      
	      if (nid != nid2) continue; /* Not the same geometric model face */
	      
	      if (!MEnt_IsMarked(ef,mkid2)) {
		List_Add(flist,ef);
		MEnt_Mark(ef,mkid1);
		MEnt_Mark(ef,mkid2);
	      }
	    }
	    List_Delete(efaces);
	  }
	  List_Delete(fedges);
	}
	

	/* Collect the unique vertices of faces on this model face */
	
	vlist = List_New(10);
	idx2 = 0;
	while ((mf2 = List_Next_Entry(flist,&idx2))) {
	  fverts = MF_Vertices(mf2,1,0);
	  idx3 = 0;
	  while ((mv = List_Next_Entry(fverts,&idx3))) {
	    if (!MEnt_IsMarked(mv,mkid2)) {
	      MEnt_Mark(mv,mkid2);
	      List_Add(vlist,mv);
	    }
	  }
	  List_Delete(fverts);
	}
	List_Unmark(vlist,mkid2);
	
	List_Unmark(flist,mkid2);
	List_Delete(flist);
	
	
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
      while ((mf = MESH_Next_Face(mesh,&idx)))
	MEnt_Unmark(mf,mkid1);
    }


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
      List_ptr fedges = MF_Edges(me,1,0);
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

