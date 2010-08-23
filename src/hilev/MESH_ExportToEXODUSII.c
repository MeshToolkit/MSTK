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
  void MESH_Num_Element_Block(Mesh_ptr mesh, int *nevs,  int num_element, \
			      int element_dim, int *num_element_block, int enable_set);
  void MESH_Get_Element_Block_Info(Mesh_ptr mesh, int *nevs, int num_element,\
				   int* element_block_ids, int *num_element_per_block);
  /* this function collects side set information                                  */
  int MESH_Num_Face_Set(List_ptr esides, int *side_set_id, int enable_set);
  int MESH_Num_Edge_Set(List_ptr esides, int *side_set_id, int enable_set);
  void MESH_Get_Side_Set_Info(Mesh_ptr mesh, int num_side_set, int num_side, int side_dim, \
			      int *side_set_ids, int *num_side_per_set);
  /* this function get the local side number */
  int Get_Local_Face_Number(MFace_ptr mf, List_ptr mrvertices);
  int Get_Local_Edge_Number(MEdge_ptr me, List_ptr mfvertices);
  void itoa(int n, char s[]);
  void reverse(char s[]);
  int compare_GID(const void * a, const void * b);


  /*this defines the maximum number of side set per element block*/
#define MAX_SIDE_SET_IN_ONE_BLOCK 20

  /* Function to export MSTK mesh to EXODUSII format               
     First, decide the element type is region or face              
     Global Parameter information:                                     
     num_node, num_elem, num_elem_blk, num_node_set, num_side_set  
     Quality assurance information(optional)                       
     The nodes(vertices) list and elements list  
     For nodes, there is only one attribute associate with them    
     verbose: 1 to print stat information 
     enable_set: 1 to enable sets based on geometric entity ID
  */
  int MESH_ExportToEXODUSII(Mesh_ptr mesh, const char *filename, int enable_set, int verbose) {
    int i, j,k, idx;
    int nfe, nfv, nrv, nrf;
    int nv, ne, nf, nr;
    MVertex_ptr mv,  vertex;
    MEdge_ptr me;
    MFace_ptr mf;
    MRegion_ptr mr;
    MEntity_ptr ms;
    RepType reptype;
    char MESH_rtype_str[5][3] = {"F1\0","F4\0","R1\0","R2\0","R4\0"};
    int eid = 0;
    /*For pure 2D mesh, boundary_dim is 1, for surface mesh, it is 2 */
    int jv;
    double xyz_coord[3];
    int boundary_dim;
    reptype = MESH_RepType(mesh);
    
    if (verbose)
      fprintf(stdout,"\nMesh representation type is %s\n", MESH_rtype_str[reptype]);
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
    
    /* decide if it is 2d planary mesh or 3D surface mesh */
    if (element_dim == 2) {
      double xyz_coord0[3];
      vertex = MESH_Vertex(mesh,0);
      MV_Coords(vertex,xyz_coord0);
      for (jv = 1; jv < nv; jv++) {
	vertex = MESH_Vertex(mesh,jv);
	MV_Coords(vertex,xyz_coord);
	/* if z coordinate is not all the same, it is a surface mesh */
	if ((xyz_coord[2]!=xyz_coord0[2]))
	  boundary_dim = 2;
      }
      if (verbose) {
	if (boundary_dim == 1)
	  fprintf(stdout,"\nThis is a 2D planary mesh\n");
      }
    }


    /* create ExodusII file */
    int CPU_word_size=sizeof(double);
    int IO_word_size=sizeof(double);
    int exoid, err;
    exoid = ex_create(filename,EX_CLOBBER,&CPU_word_size,&IO_word_size);
    if (exoid < 0) {
      fprintf(stderr, "after ex_creat, error = %d\n", exoid);
      exit(-1);
    }
    
    /* Get element block information(based on element type, and geometric entity id) 
       element_block_id = geometric_entity_id << 6 | MR_Num_Vertices              */
    if(verbose) {
      if (enable_set)
	fprintf(stdout,"\nElement block and side set based on geometric id is enabled.\n");
      else
	fprintf(stdout,"\nElement block and side set based on geometric id is disabled.\n");
    }
    int num_element_block = 0;
    int *nevs = (int *) MSTK_malloc(num_element*sizeof(int));
    MESH_Num_Element_Block(mesh, nevs, num_element, element_dim, &num_element_block, enable_set);
    int *element_block_ids = (int *) MSTK_malloc((num_element_block)*sizeof(int));
    int *num_element_per_block = (int *) MSTK_malloc((num_element_block)*sizeof(int));
    MESH_Get_Element_Block_Info(mesh, nevs, num_element, element_block_ids, num_element_per_block);
    MSTK_free(nevs);
    
    if (verbose) {
      fprintf(stdout,"\nElemen1t block information:\n");
      fprintf(stdout,"total %d elements\n", num_element); 
      fprintf(stdout,"total %d element blocks:\n\n", num_element_block);
    }
    
    /* Get node set information */
    int num_side_set = 0;
    int *side_set_id = (int *)MSTK_malloc(num_element_block*MAX_SIDE_SET_IN_ONE_BLOCK*sizeof(int));
    int *num_side_set_in_block = (int *)MSTK_malloc(num_element_block*sizeof(int));
    for(i = 0; i < num_element_block; i++) {
      List_ptr esides=List_New(3*num_element_per_block[i]);
      int mksid = MSTK_GetMarker();
      if (element_dim == 3) {
	List_ptr mrfaces;
	idx = 0; 
	while ((mr = MESH_Next_Region(mesh,&idx))) {
	  nrv = MR_Num_Vertices(mr);
	  /* test if the region is in this block */
	  if (enable_set)
	    eid = (MEnt_GEntID(mr)<<6) | nrv;
	  else
	    eid = nrv;
	  if (eid == element_block_ids[i]) {
	    mrfaces = MR_Faces(mr);
	    nrf = List_Num_Entries(mrfaces);
	    for (j = 0; j < nrf; j++) {
	      mf = List_Entry(mrfaces,j);
	      /* if the face is not in the interior and not marked */
	      if(MEnt_GEntDim(mf)<=2 && (!MEnt_IsMarked(mf,mksid))) {
		List_Add(esides,mf);
		MEnt_Mark(mf,mksid);
	      }
	    }
	    List_Delete(mrfaces);
	  }
	}
	num_side_set_in_block[i] = MESH_Num_Face_Set(esides,&side_set_id[i*MAX_SIDE_SET_IN_ONE_BLOCK],\
						     enable_set);
	
	num_side_set = num_side_set + num_side_set_in_block[i];

	if(verbose)
	  fprintf(stdout,"Element block %d has %d nodes per region, %d regions, %d face sets \n",\
		  element_block_ids[i],element_block_ids[i]&31,num_element_per_block[i], \
		  num_side_set_in_block[i]);
      }
      if (element_dim == 2) {
	List_ptr mfedges;
	idx = 0;
	while ((mf = MESH_Next_Face(mesh,&idx))) {
	  nfv = MF_Num_Vertices(mf);
	  /* test if the region is in this block */
	  if (enable_set)
	    eid = (MEnt_GEntID(mf)<<6) | nfv;
	  else
	    eid = nfv;
	  if (eid == element_block_ids[i]) {
	    mfedges = MF_Edges(mf,1,0);
	    nfe = List_Num_Entries(mfedges);
	    for (j = 0; j<nfe; j++) {
	      me = List_Entry(mfedges,j);
	      /* if the face is not in the interior and not marked */
	      if(MEnt_GEntDim(me)<=boundary_dim && (!MEnt_IsMarked(me,mksid))) {
		List_Add(esides,me);
		MEnt_Mark(me,mksid);
	      }
	    }
	    List_Delete(mfedges);
	  }
	}
	num_side_set_in_block[i] = MESH_Num_Edge_Set(esides,&side_set_id[i*MAX_SIDE_SET_IN_ONE_BLOCK],\
						     enable_set);
	num_side_set = num_side_set + num_side_set_in_block[i];
	if(verbose)
	  fprintf(stdout,"Element block %d has %d nodes per face, %d faces, %d edge sets \n",\
		  element_block_ids[i],element_block_ids[i]&31, num_element_per_block[i], \
		  num_side_set_in_block[i]);
      }
      List_Unmark(esides,mksid);
      MSTK_FreeMarker(mksid);
      List_Delete(esides);
    }
	
    /* write global parameters */
    int num_node = nv;
    int num_node_set = num_side_set;
    err = ex_put_init (exoid, "MSTK to EXODUSII converter", 3,
		       num_node, num_element, num_element_block,
		       num_node_set, num_side_set);
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
    
    /* write element blocks, side sets and node sets information */
    int *num_node_per_element = MSTK_malloc(num_element_block*sizeof(int));
    int index_element = 0, index_side = 0;

    if (verbose) {
      fprintf(stdout,"\nSide sets and node sets information:\n");
      fprintf(stdout,"total %d side sets:\n",num_side_set);
    }
    
    for (i=0; i<num_element_block; i++) { 
      num_node_per_element[i] = element_block_ids[i] & 63;      
      /* evids store the vertex id of element 
	 seids store the element id of sides
	 slids store local number of side on the element
      */
      int *evids = (int *) MSTK_malloc(num_element_per_block[i]*num_node_per_element[i]*sizeof(int));
      int *seids = (int *) MSTK_malloc(num_element_per_block[i]*8*sizeof(int));
      int *slids = (int *) MSTK_malloc(num_element_per_block[i]*8*sizeof(int));
      /* 
	 esids store the unique sides pointers
	 evertices store the unique vertices pointers
       */
      List_ptr esides = List_New(100);
      int mksid = MSTK_GetMarker();
      char element_type[256]; 
      char buf[256];
      /* for volume mesh */
      if (element_dim == 3) {
	List_ptr mrfaces, mrvertices;
	idx = 0; index_element = 0; index_side = 0;
	while ((mr = MESH_Next_Region(mesh,&idx))) {
	  mrvertices = MR_Vertices(mr);
	  nrv = List_Num_Entries(mrvertices);
	  /* test if the region is in this block */
	  if (enable_set)
	    eid = (MEnt_GEntID(mr)<<6) | nrv;
	  else
	    eid = nrv;
	  if (eid == element_block_ids[i]) {
	    for (j = 0; j < nrv; j++) {
	      mv = List_Entry(mrvertices,j);
	      evids[index_element*num_node_per_element[i]+j] = MEnt_ID(mv);
	    }
	    index_element++;		   
	    mrfaces = MR_Faces(mr);
	    nrf = List_Num_Entries(mrfaces);
	    for (j = 0; j < nrf; j++) {
	      mf = List_Entry(mrfaces,j);
	      /* if the face is not in the interior and not marked */
	      if(MEnt_GEntDim(mf)<=2 && (!MEnt_IsMarked(mf,mksid))) {
		List_Add(esides,mf);
		MEnt_Mark(mf,mksid);
		/* if it is multimaterial mesh */
		seids[index_side] = 0;
		for (k = 0 ; k < i; k++)
		  seids[index_side] += num_element_per_block[k];
		seids[index_side] +=index_element;
		
		slids[index_side] = Get_Local_Face_Number(mf, mrvertices);		
		index_side++;
	      }
	    }
	    List_Delete(mrfaces);
	  }
	  List_Delete(mrvertices);
	}
	
	/* sort esides based on geometric ID */
	/* qsort(esides,index_side,sizeof(MFace_ptr),compare_GID);*/
       switch (num_node_per_element[i]) { 
       case 4 :
	 strcpy(buf,"TETRA");
	 break;
       case 6:
	 strcpy(buf,"WEDGE");
	 break;
       case 8:
	 strcpy(buf,"HEX");
	 break;
       default:
	 itoa(element_block_ids[i],buf);
	 strcat(buf,"_Element_Block");
       }
       strcpy(element_type,buf);
      } /*end if(element_dim == 3)*/
      /* for surface mesh */
      if (element_dim == 2) {
	List_ptr mfedges, mfvertices;
	idx = 0; index_element = 0; index_side = 0;
	while ((mf = MESH_Next_Face(mesh,&idx))) {
	  mfvertices = MF_Vertices(mf,1,0);
	  nfv = List_Num_Entries(mfvertices);
	  /* test if the region is in this block */
	  if (enable_set)
	    eid = (MEnt_GEntID(mf)<<6) | nfv;
	  else
	    eid = nfv;
	  if (eid == element_block_ids[i]) {
	    for (j = 0; j<nfv; j++) {
	      mv = List_Entry(mfvertices,j);
	      evids[index_element*num_node_per_element[i]+j] = MEnt_ID(mv);
	    }
	    index_element++;		   
	    mfedges = MF_Edges(mf,1,0);
	    nfe = List_Num_Entries(mfedges);
	    for (j = 0; j<nfe; j++) {
	      me = List_Entry(mfedges,j);
	      /* if the face is not in the interior and not marked */
	      if(MEnt_GEntDim(me)<=boundary_dim && (!MEnt_IsMarked(me,mksid))) {
		List_Add(esides,me);
		MEnt_Mark(me,mksid);
		seids[index_side] = 0;
		for (k = 0 ; k < i; k++)
		  seids[index_side] += num_element_per_block[k];
		seids[index_side] +=index_element;
		slids[index_side] = Get_Local_Edge_Number(me, mfvertices);
		index_side++;
	      }
	    }
	    List_Delete(mfedges);
	  }
	  List_Delete(mfvertices);
	}
	
       switch (num_node_per_element[i]) { 
       case 3 :
	 strcpy(buf,"TRIANGLE");
	 break;
       case 4:
	 strcpy(buf,"QUAD");
	 break;
       default:
	 itoa(element_block_ids[i],buf);
	 strcat(buf,"_Face_Block");
       }
       strcpy(element_type,buf);
      }
       /* write element block information and connectivity */
      int  num_attr_per_element = 0;
      err = ex_put_elem_block (exoid, element_block_ids[i], element_type, num_element_per_block[i],
			       num_node_per_element[i], num_attr_per_element);

      err = ex_put_elem_conn(exoid, element_block_ids[i], evids);
      MSTK_free(evids);

    
     
       /* write side and node set information*/
       int num_dist_fact_in_side_set = 0;
       int num_dist_fact_in_node_set = 0;
       int side_id;
       int nsv,jj;
       int mkvid;
       int local_index_side;
       int local_index_vertex;
       for (k = 0 ; k<num_side_set_in_block[i]; k++) {
	 mkvid = MSTK_GetMarker();
	 List_ptr msvertices;
	 List_ptr svertices=List_New(100);
	 int *local_nids = (int *)MSTK_malloc(num_element_per_block[i]*20*sizeof(int));
	 int *local_seids = (int *) MSTK_malloc(num_element_per_block[i]*8*sizeof(int));
	 int *local_slids = (int *) MSTK_malloc(num_element_per_block[i]*8*sizeof(int));
	 local_index_side = 0; local_index_vertex = 0;
	 for (j=0; j<index_side; j++) {
	   ms = List_Entry(esides,j);
	   if ((!enable_set) || (MEnt_GEntID(ms) == side_set_id[i*MAX_SIDE_SET_IN_ONE_BLOCK+k])) {
	     
	       local_seids[local_index_side] = seids[j];
	       local_slids[local_index_side++] = slids[j];
	     if (side_dim == 2)  {
	       msvertices = MF_Vertices((MFace_ptr)ms,1,0);
	       nsv = List_Num_Entries(msvertices);
	     }
	     else {
	       msvertices = List_New(2);
	       List_Add(msvertices,ME_Vertex((MEdge_ptr)ms,0));
	       List_Add(msvertices,ME_Vertex((MEdge_ptr)ms,1));
	       nsv = 2;
	     }
	     
	     for(jj = 0; jj<nsv; jj++) {
	       mv = List_Entry(msvertices,jj);
	       if (MEnt_GEntDim(mv)<=2 && !MEnt_IsMarked(mv,mkvid)) {
		 local_nids[local_index_vertex++] = MEnt_ID(mv);
		 List_Add(svertices,mv);
		 MEnt_Mark(mv,mkvid);
	       }
	     }
	     List_Delete(msvertices);
	   }
	 }
	 List_Unmark(svertices,mkvid);
	 MSTK_FreeMarker(mkvid);
	 List_Delete(svertices);
	 
	 side_id = (element_block_ids[i]<<8) | (k+1);
	 if (verbose)
	   fprintf(stdout,"side set %d on element block %d has %d sides and %d vertices\n",\
		   side_set_id[i*MAX_SIDE_SET_IN_ONE_BLOCK+k],element_block_ids[i],\
		   local_index_side,local_index_vertex);

	 if (local_index_side>0) {
	   err = ex_put_side_set_param(exoid,side_id,local_index_side,num_dist_fact_in_side_set);
	   err = ex_put_side_set(exoid,side_id,local_seids,local_slids);
	 }
	 if (local_index_vertex>0) {
	   err = ex_put_node_set_param(exoid,side_id,local_index_vertex,num_dist_fact_in_node_set);	
	   err = ex_put_node_set(exoid,side_id,local_nids);
	 }
	 if(err) {
	   fprintf(stdout, "error put side set and node side\n");
	   ex_close(exoid);
	   exit(-1);
	 }
	 MSTK_free(local_nids);
	 MSTK_free(local_seids);
	 MSTK_free(local_slids);
       }
      List_Unmark(esides,mksid);
      MSTK_FreeMarker(mksid);
      List_Delete(esides);      


      MSTK_free(seids);
      MSTK_free(slids);
    }/* end for(i=0;i<num_element_block */


    MSTK_free(side_set_id);
    MSTK_free(num_side_set_in_block);
    MSTK_free(element_block_ids);
    MSTK_free(num_element_per_block);
    MSTK_free(num_node_per_element);


    if (verbose)
      fprintf(stdout,"Side set information wrote into %s\n", filename);
    
    if (verbose)
      fprintf(stdout,"Node sets information wrote into %s\n", filename);

     
      /* write quality assurance information optional*/
  char *qa_record[1][4], date_str[256], time_str[256];
  qa_record[0][0] = "MSTK to ExodusII converter";
  qa_record[0][1] = "Duo Wang T-5 group LANL";
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


  int compare(const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
  }

  int compare_GID(const void * a, const void * b) {
    return (MEnt_GEntID((MEntity_ptr*)a) - MEnt_GEntID((MEntity_ptr*)b));
  }


  /* some functions used by the converter */
  void MESH_Num_Element_Block(Mesh_ptr mesh, int *nevs, int num_element, int element_dim, \
			      int *num_element_block, int enable_set) {
    MEdge_ptr me;
    MFace_ptr mf;
    MRegion_ptr mr;

    *num_element_block = 1;
    int idx = 0, i = 0;
    
    if (element_dim == 3) {
      while ((mr = MESH_Next_Region(mesh, &idx)))
	if (enable_set)
	  nevs[i++] = (MEnt_GEntID(mr)<<6) | MR_Num_Vertices(mr);
	else
	  nevs[i++] = MR_Num_Vertices(mr);
    }
    if (element_dim == 2) {
      while ((mf = MESH_Next_Face(mesh,&idx)))
	if (enable_set)
	  nevs[i++] = (MEnt_GEntID(mf)<<6) | MF_Num_Vertices(mf);
	else
	  nevs[i++] = MF_Num_Vertices(mf);
    }
    if (element_dim == 1) {
      while ((me = MESH_Next_Edge(mesh, &idx)))
	if (enable_set)
	  nevs[i++] = (MEnt_GEntID(me)<<6) | 2;
	else
	  nevs[i++] = 2;
	  
    }
    qsort(nevs,num_element,sizeof(int),compare);
    int prev_nev = nevs[0];
    /* first sweep counts the number of different node sets based on geometric entity id and dimension */
    for (i=0;i<num_element;i++) {
      if (nevs[i] != prev_nev) { 
	(*num_element_block)++;
	prev_nev = nevs[i];
      }
    }
  }	

  void MESH_Get_Element_Block_Info(Mesh_ptr mesh, int *nevs, int num_element, int* element_block_ids, \
				   int *num_element_per_block) {

    int i = 0;
    
    /* nevs is already sorted 
    qsort(nevs,num_element,sizeof(int),compare);
    */
    int prev_nev = nevs[0];
    int num_element_this_block = 0;
    /* if there is only one node set */
    element_block_ids[0] = nevs[0];
    num_element_per_block[0] = num_element;
    /* second sweep record node set id and number of nodes */    
    int k = 1;
    prev_nev = nevs[0];
    for (i=0;i<num_element;i++) {
      if (nevs[i] == prev_nev) 
	num_element_this_block++;
      else  {
	num_element_per_block[k-1] = num_element_this_block;
	k++;
	element_block_ids[k-1] = nevs[i];
	prev_nev = nevs[i];
	num_element_this_block = 1;
      }
    }
    num_element_per_block[k-1] = num_element_this_block;

  }
  

  void  MESH_Num_Side_Set(Mesh_ptr mesh, int num_side, int side_dim, int *num_side_set) {
    
    *num_side_set = 1;
    if (side_dim < 1) {
      MSTK_Report("MESH_Num_Side_Set", "No side information\n",WARN);
      return;
    }
    int *gids = (int *) MSTK_malloc(num_side*sizeof(int));
    MFace_ptr mf;
    MEdge_ptr me;
    int idx = 0, i = 0;
    if (side_dim == 1) {
    while ((me = MESH_Next_Edge(mesh,&idx)))
      gids[i++] = ((MEnt_GEntID(me)+1)<<3) | MEnt_GEntDim(me) ;
    }
    if (side_dim == 2) {
    while ((mf = MESH_Next_Face(mesh,&idx)))
      gids[i++] = ((MEnt_GEntID(mf)+1)<<3) | MEnt_GEntDim(mf) ;
    } 

    qsort(gids, num_side, sizeof(int), compare);
    int prev_gid = gids[0];
    /* first sweep counts the number of different node sets based on geometric entity id and dimension */
    for (i=0;i<num_side;i++) {
      if (gids[i] != prev_gid) { 
	(*num_side_set)++;
	prev_gid = gids[i];
      }
    }
    MSTK_free(gids);
  }

	

 
  void MESH_Get_Side_Set_Info(Mesh_ptr mesh, int num_side_set, int num_side, int side_dim, int *side_set_ids, int *num_side_per_set) {
    MFace_ptr mf;
    MEdge_ptr me;
    int idx = 0, i = 0;
    /* if num_side_set is 1, then quick return */
    if (num_side_set == 1) {
      if (side_dim == 2) {
	mf = MESH_Next_Face(mesh,&idx);
	side_set_ids[0] = ((MEnt_GEntID(mf)+1)<<3) | MEnt_GEntDim(mf) ;
      }
      else {
	me = MESH_Next_Edge(mesh, &idx);
	side_set_ids[0] = ((MEnt_GEntID(me)+1)<<3) | MEnt_GEntDim(me) ;
      }
      num_side_per_set[0] = num_side;
      return;
    }

    int *gids = (int *) MSTK_malloc(num_side*sizeof(int)); 
    
    if (side_dim == 2) {
      while ((mf = MESH_Next_Face(mesh,&idx)))
	gids[i++] = ((MEnt_GEntID(mf)+1)<<3) | MEnt_GEntDim(mf) ;
    }
    if (side_dim == 1) {
      while ((me = MESH_Next_Edge(mesh, &idx)))
	gids[i++] = ((MEnt_GEntID(me)+1)<<3) | MEnt_GEntDim(me) ;
    }
    qsort(gids,num_side,sizeof(int),compare);
    int prev_gid = gids[0];
    int num_side_this_set = 0;
    /* if there is only one node set */
    side_set_ids[0] = gids[0];
    num_side_per_set[0] = num_side;
    /* second sweep record node set id and number of nodes */    
    int k = 1;
    prev_gid = gids[0];
    for (i=0;i<num_side;i++) {
      if (gids[i] == prev_gid) 
	num_side_this_set++;
      else  {
	num_side_per_set[k-1] = num_side_this_set;
	k++;
	side_set_ids[k-1] = gids[i];
	prev_gid = gids[i];
	num_side_this_set = 1;
      }
    }
    num_side_per_set[k-1] = num_side_this_set;

    MSTK_free(gids);
  }
  

 void itoa(int n, char s[]) {
    int i, sign;
    sign = n;
    
    i = 0;
    do {
      s[i++] = abs(n % 10) + '0';
    } while ( n /= 10 );
    if (sign < 0)
      s[i++] = '-';
    
    s[i] = '\0';
    reverse(s);
  }
  void reverse(char s[]) {
    int c, i, j;
    for ( i = 0, j = strlen(s)-1; i < j; i++, j--) {
      c = s[i];
      s[i] = s[j];
      s[j] = c;
    }
  }  


  int Get_Local_Face_Number(MEdge_ptr mf, List_ptr mrvertices) {
    List_ptr mfvertices = MF_Vertices(mf,1,0);
    int i,j;
    MVertex_ptr mv;
    int num_node_per_element = List_Num_Entries(mrvertices);
    /* for tetra element */
    if (num_node_per_element == 4) {
      int v2l[4]={2,3,1,4};
      int svid[3], evid[4];
      for (i = 0; i < 3; i++) {
	mv = List_Entry(mfvertices,i);
	svid[i] = MEnt_ID(mv);
      }
      for (i = 0; i < 4; i++) {
	mv = List_Entry(mrvertices,i);
	evid[i] = MEnt_ID(mv);
      }
      for (i = 0; i < 4; i++) 
	for(j = 0; j < 3; j++) {
	  if (evid[i] == svid[j])
	    break;
	  if (j==2) {
	    return v2l[i];
	  }
	}
    }
    /*for hex element*/
    if (num_node_per_element == 8) {
      int svid[4], evid[8];
      int local_id[4];
      /*local_id stores the index of the face with respect to the region vertex index */
      for (i = 0; i < 4; i++) {
	mv = List_Entry(mfvertices,i);
	svid[i] = MEnt_ID(mv);
      }
      for (i = 0; i < 8; i++) {
	mv = List_Entry(mrvertices,i);
	evid[i] = MEnt_ID(mv);
      }
      /* for hex, need the first 3 local index of the face */
      int index = 0;
      for (i = 0; i < 8; i++) 
	for(j = 0; j < 4; j++) {
	  if (evid[i] == svid[j]) {
	    local_id[index++] = i+1;
	    break;
	  }
	  if(index>3)
	    break;
	}
      int sum = 0;
      for (i = 0; i < 3; i++)
	sum += local_id[i];
      switch (sum) {
      case 8:
	return 1;
      case 11:
	return 2;
      case 14:
	return 3;
      case 10:
	return 4;
      case 6:
	return 5;
      case 18:
	return 6;
      default:
	fprintf(stdout,"Error:no side index found on Hex\n");
	return 1;
      }
    }
    return 1;

  }
	    
  int Get_Local_Edge_Number(MEdge_ptr me, List_ptr mfvertices) {
    int i,j;
    MVertex_ptr mv;
    int num_node_per_element = List_Num_Entries(mfvertices);
    /* for tri element */
    if (num_node_per_element == 3) {
      /* see the cubit mannual appendix for 3d TRI side numbering */
      int v2l[3]={4,5,3};
      int svid[2], evid[3];
      for (i = 0; i < 2; i++) {
	mv = ME_Vertex(me,i);
	svid[i] = MEnt_ID(mv);
      }
      for (i = 0; i < 3; i++) {
	mv = List_Entry(mfvertices,i);
	evid[i] = MEnt_ID(mv);
      }
      for (i = 0; i < 3; i++) 
	for(j = 0; j < 2; j++) {
	  if (evid[i] == svid[j])
	    break;
	  if (j==1) {
	    return v2l[i];
	  }
	}
    }
    /*for quad element*/
    if (num_node_per_element == 4) {
      int svid[2], evid[4];
      for (i = 0; i < 2; i++) {
	mv = ME_Vertex(me,i);
	svid[i] = MEnt_ID(mv);
      }
      for (i = 0; i < 4; i++) {
	mv = List_Entry(mfvertices,i);
	evid[i] = MEnt_ID(mv);
      }
      int index = 0;
      int local_id[2];
      for (i = 0; i < 4; i++) 
	for(j = 0; j < 2; j++) {
	  if (evid[i] == svid[j]) {
	    local_id[index++] = i+1;
	    break;
	  }
	  if(index>2)
	    break;
	}
      /* Use the product of the two local id number to distinguishi*/
      int product = local_id[0]*local_id[1];
      switch (product) {
      case 2:
	return 1;
      case 6:
	return 2;
      case 12:
	return 3;
      case 4:
	return 4;
      default:
	fprintf(stdout,"Error:no side index found on Face\n");
	return 1;
      }
    } 
    return 1;
  }

  int MESH_Num_Face_Set(List_ptr esides, int *side_set_id, int enable_set) {
    int num_side = List_Num_Entries(esides);
    if (num_side == 0) 
      return 0;
    if (!enable_set) {
      side_set_id[0] = 1;
      return 1;
    }
    int num_side_set = 1;
    int *gids = (int *) MSTK_malloc(num_side*sizeof(int));
    MFace_ptr mf;
    int i = 0;
    for (i=0; i<num_side;i++) {
      mf = List_Entry(esides,i);
      gids[i] = MEnt_GEntID(mf);
    }    
    qsort(gids, num_side, sizeof(int), compare);
    int prev_gid = gids[0];
    /* first sweep counts the number of different node sets based on geometric entity id and dimension */
    side_set_id[0] = gids[0];
    for (i=1;i<num_side;i++) {
      if (gids[i] != prev_gid) { 
	side_set_id[num_side_set] = gids[i];
	num_side_set++;
	prev_gid = gids[i];
	if(num_side_set>=MAX_SIDE_SET_IN_ONE_BLOCK) {
	  MSTK_Report("MESH_Num_Edge_Side", "Maximum number of side set within one block is achieved\n",WARN);
	  break;
	}
      }
    }
    MSTK_free(gids);
    return num_side_set;
  }

  int MESH_Num_Edge_Set(List_ptr esides, int *side_set_id, int enable_set) {
    int num_side = List_Num_Entries(esides);
    if (num_side == 0) 
      return 0;
    if (!enable_set) {
      side_set_id[0] = 1;
      return 1;
    }
    int num_side_set = 1;
    int *gids = (int *) MSTK_malloc(num_side*sizeof(int));
    MEdge_ptr me;
    int i = 0;
    for (i=0; i<num_side;i++) {
      me = List_Entry(esides,i);
      gids[i] = MEnt_GEntID(me);
    }     
    qsort(gids, num_side, sizeof(int), compare);
    int prev_gid = gids[0];
    /* first sweep counts the number of different node sets based on geometric entity id and dimension */
    side_set_id[0] = gids[0];
    for (i=1;i<num_side;i++) {
      if (gids[i] != prev_gid) { 
	side_set_id[num_side_set] = gids[i];
	num_side_set++;
	prev_gid = gids[i];
	if(num_side_set>=MAX_SIDE_SET_IN_ONE_BLOCK) {
	  MSTK_Report("MESH_Num_Edge_Side", "Maximum number of side set within one block is achieved\n",WARN);
	  break;
	}
      }
    }
    MSTK_free(gids);
    return num_side_set;
  } 


#ifdef __cplusplus
}
#endif

