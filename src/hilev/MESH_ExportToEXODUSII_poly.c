
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"
#include "exodusII.h"

#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

  
  /* this function collects element block inforamtion based on element type       */

  int MESH_Num_Side_Set(Mesh_ptr mesh, int num_side, int side_dim, int enable_set);
  void MESH_Get_Side_Set_Info(Mesh_ptr mesh, int num_side_set, int num_side, int side_dim, int *side_set_ids, int *num_side_per_set);
  static  int compare(const void * a, const void * b);


  /* Function to export MSTK mesh to EXODUSII format               */

  /* First, decide the element type is region or face              */

  /* Global Parameter information:                                 */    
  /* num_node, num_elem, num_elem_blk, num_node_set, num_side_set  */
  
  /* Quality assurance information(optional)                       */

  /* The nodes(vertices) list and elements list  */
  /* For nodes, there is only one attribute associate with them    */

  /* notice: the set id of nodes and sides are (MEnt_GEntID+1)
     ID + 1 to avoid 0 id ExodusII does not permit node or side 
     set with 0 id
  */
int MESH_ExportToEXODUSII_poly(Mesh_ptr mesh, const char *filename, int enable_set, int verbose) {
  int i, j, idx;
  int  nfe, nfv,nev, nrf;
  int nv, ne, nf, nr;
  MVertex_ptr mv,  vertex;
  MEdge_ptr me;
  MFace_ptr mf, mrfaces;
  MRegion_ptr mr;
  List_ptr mfverts,  mrverts;
  RepType reptype;
  char MESH_rtype_str[5][3] = {"F1\0","F4\0","R1\0","R2\0","R4\0"};


  int enable_side_set = enable_set;
  int enable_node_set = enable_set;

  reptype = MESH_RepType(mesh);
  
  if (verbose)
  fprintf(stdout,"\nMesh representation type is %s\n", MESH_rtype_str[reptype]);
  nv = MESH_Num_Vertices(mesh);
  ne = MESH_Num_Edges(mesh);
  nf = MESH_Num_Faces(mesh);
  nr = MESH_Num_Regions(mesh);
 
  int element_dim = -1, side_dim = -1;
  int num_element = 0, num_side = 0;
  if (nr) {
      element_dim = 3; side_dim = 2;
      num_element = nr; num_side = nf;
      fprintf(stdout,"\nThis is a 3D volume mesh with %d Vertices, %d Edges, %d Faces and %d Regions\n", nv,ne,nf,nr);
    }
    else if(nf) {
      element_dim = 2; side_dim = 1;
      num_element = nf; num_side = ne;
      fprintf(stdout,"\nThis is a 3D surface mesh with %d Vertices, %d Edges and %d Faces\n", nv,ne,nf);
    }
    else if(ne) {
      element_dim = 1; side_dim = 0;
      num_element = ne; num_side = nv;
      fprintf(stdout,"\nThis is a 3D curve mesh with %d Vertices and %d Edges\n", nv,ne);
    }
    else if(nv) {
      element_dim = 0; side_dim = 0;
      num_element = nv; num_side = 0;
      fprintf(stdout,"\nThis is a 3D mesh with %d Vertices\n", nv);
    }
    else {
    fprintf(stdout,"\nThis is an empty mstk file\n");
    exit(-1);
    }  
  
  int num_node = nv;
  /* right now just for 3D volume mesh */
  if (side_dim != 2) {
    fprintf(stdout, "\n Stop: This is not a volume mesh\n");
    exit(-1);
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

 
  /* Get side set information(optional,based on geometric entity) */
  int num_side_set = 0;
  num_side_set = MESH_Num_Side_Set(mesh,num_side, side_dim, enable_set);
  if (verbose)
    fprintf(stdout,"\ntotal %d side set:\n", num_side_set);
  if (enable_side_set)
    fprintf(stdout,"\nSide set is enabled\n");
  else {
    fprintf(stdout,"\nSide set is disabled\n");
    num_side_set = 1;
  }
  int *side_set_ids = (int *) MSTK_malloc((num_side_set)*sizeof(int));
  int *num_side_per_set = (int *) MSTK_malloc((num_side_set)*sizeof(int));
  MESH_Get_Side_Set_Info(mesh, num_side_set, num_side, side_dim, side_set_ids, num_side_per_set);

  /* write global parameters */
  ex_init_params par;
  strcpy(par.title,filename);
  par.num_dim = 3;
  par.num_nodes = nv;
  par.num_edge = 0;
  par.num_face = nf;
  if (side_dim == 1) {
    par.num_edge_blk = num_side_set;
    par.num_face_blk = 0;
    par.num_elem = nf;
  }
  if (side_dim == 2) {
    par.num_face_blk = num_side_set;
    par.num_edge_blk = 0;
    par.num_elem = nr;
  }
  par.num_elem_blk = 1;  
  par.num_node_sets = num_side_set;
  par.num_edge_sets = 0;
  par.num_face_sets = 0;
  par.num_side_sets = 0;
  par.num_elem_sets = 0;
  par.num_node_maps = 0;
  par.num_edge_maps = 0;
  par.num_face_maps = 0;
  par.num_elem_maps = 0;

  ex_put_init_ext(exoid,&par);
  
  if (verbose)
    fprintf(stdout,"Globla parameters wrote into %s\n", filename);
  
  /* write coordinate values */
  double *xcoord, *ycoord, *zcoord, vxyz[3];
  int jv;
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
  if (err) {
    fprintf(stderr, "after ex_put_coord, error = %d\n", err);
    ex_close (exoid);
    exit(-1);
  }
  MSTK_free(xcoord);
  MSTK_free(ycoord);
  MSTK_free(zcoord);
  /* add coordinate names(optional) */
  char *coord_names[3] = {"xcoor","ycoor","zcoor"};
  err = ex_put_coord_names (exoid, coord_names);
  if (err) {
    fprintf(stderr, "after ex_put_coord_names, error = %d\n", err);
    ex_close (exoid);
    exit(-1);
    }
  if (verbose)
    fprintf(stdout,"Coordinate values wrote into %s\n", filename);


  /* Write face blocks and node sets */
  if (verbose)
  fprintf(stdout,"\nFace block information:\n");
  char buf[10];
  char face_block_name[256];
  int num_total_nodes = 0;
  int index_face = 0;
  int num_unique_node = 0;
  int mkid;
  for (i = 0; i < num_side_set; i++) {
    strcpy(face_block_name, "GFace_block_");
    itoa(side_set_ids[i],buf);
    strcat(face_block_name,buf);
    
    /* allocate maximum of num_face*20 nodes*/
    int *vids = (int *)MSTK_malloc(num_side_per_set[i]*20*sizeof(int));
    int *nvids = (int *)MSTK_malloc(num_side_per_set[i]*sizeof(int));

    int *unique_vids = (int *)MSTK_malloc(num_side_per_set[i]*20*sizeof(int));
  
    mkid = MSTK_GetMarker();
    idx = 0; num_total_nodes = 0; index_face = 0; num_unique_node = 0;
    while ((mf = MESH_Next_Face(mesh,&idx))) {
      if ((!enable_side_set) || (side_set_ids[i] == (MEnt_GEntID(mf)+1))) {
	nfv = MF_Num_Edges(mf);
	nvids[index_face++] = nfv;
	mfverts = MF_Vertices(mf,1,0);
	for (j = 0; j < nfv; j++) {
	  mv = List_Entry(mfverts,j);
	  vids[num_total_nodes++] = MEnt_ID(mv);
	  if (!MEnt_IsMarked(mv,mkid)) {
	    unique_vids[num_unique_node++] = MEnt_ID(mv);
	    MEnt_Mark(mv,mkid);
	  }
	}
	List_Delete(mfverts);
      }
    }
    idx = 0;
    /* Unmark vertices */
    while ((mv = MESH_Next_Vertex(mesh, &idx))) { 
      MEnt_Unmark(mv,mkid);
    }
    /* Release marker */		
    MSTK_FreeMarker(mkid);
    /* write face block */
    ex_put_block(exoid,EX_FACE_BLOCK,side_set_ids[i],"nsided",num_side_per_set[i],num_total_nodes,0,0,0);
    ex_put_name(exoid,EX_FACE_BLOCK,side_set_ids[i],face_block_name);
    ex_put_conn(exoid,EX_FACE_BLOCK,side_set_ids[i],vids,NULL,NULL);
    ex_put_entity_count_per_polyhedra(exoid,EX_FACE_BLOCK,side_set_ids[i],nvids);
    /* write node set */
    ex_put_node_set_param(exoid, side_set_ids[i],num_unique_node,0);  
    ex_put_node_set(exoid,side_set_ids[i],unique_vids);
    if (verbose) {
      fprintf(stdout, "Face block %s with geometric id %d has %d faces and %d nodes\n", face_block_name, side_set_ids[i]-1,num_side_per_set[i], num_unique_node);
    }
    MSTK_free(unique_vids);
    MSTK_free(vids);
    MSTK_free(nvids);

  }
  MSTK_free(side_set_ids);
  MSTK_free(num_side_per_set);
  /* write element block(only one) */
  int index_element = 0;
  int *fids = (int *) MSTK_malloc(num_side*30*sizeof(int));
  int *nfids = (int *) MSTK_malloc(num_side*sizeof(int));
  int num_total_faces = 0;
  idx = 0;
  while ((mr = MESH_Next_Region(mesh,&idx))) {
    nrf = MR_Num_Faces(mr);
    nfids[index_element++] = nrf;	
    mrfaces = MR_Faces(mr);
    for (j = 0; j < nrf; j++) {
      mf = List_Entry(mrfaces,j);
      fids[num_total_faces++] = MEnt_ID(mf);
    }
    List_Delete(mrfaces);
  }
  int element_block_id = 10;
  char element_block_name[256]="Global_Element_Block";
  ex_put_block(exoid, EX_ELEM_BLOCK, element_block_id, "nfaced", num_element, 0,0,num_total_faces,0);
  ex_put_name(exoid, EX_ELEM_BLOCK, element_block_id, element_block_name);
  ex_put_conn(exoid, EX_ELEM_BLOCK, element_block_id, NULL,NULL, fids);
  ex_put_entity_count_per_polyhedra(exoid, EX_ELEM_BLOCK, element_block_id, nfids);  
  MSTK_free(fids);
  MSTK_free(nfids);
  idx = 0;
  mr =MESH_Next_Region(mesh, &idx);
  if (verbose) {
    fprintf(stdout,"\nElement block information\n");
    fprintf(stdout, "Element block %s with geometric id %d and dimension %d has %d faces\n", element_block_name, MEnt_GEntID(mr), MEnt_GEntDim(mr), num_side);
  }
  
  
  if (verbose)
    fprintf(stdout,"Element blocks connectivity wrote into %s\n", filename);




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
  if (err) {
    fprintf(stderr, "after ex_put_qa, error = %d\n", err);
    ex_close (exoid);
    exit(-1);
  }
  if (verbose)
    fprintf(stdout,"\nQuality assurance information wrote into %s\n", filename);


  ex_close(exoid);
  return 1;
}


static  int compare(const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
  }



  /* some functions used by the converter */
  int  MESH_Num_Side_Set(Mesh_ptr mesh, int num_side, int side_dim, int enable_set) {
    
    int num_side_set = 1;
    if (side_dim < 1) {
      MSTK_Report("MESH_Num_Side_Set", "No side information\n",WARN);
      return;
    }
    if(!enable_set)
      return 1;


    int *gids = (int *) MSTK_malloc(num_side*sizeof(int));
    MFace_ptr mf;
    MEdge_ptr me;
    int idx = 0, i = 0;
    if (side_dim == 1) {
    while ((me = MESH_Next_Edge(mesh,&idx)))
      gids[i++] = MEnt_GEntID(me);
    }
    if (side_dim == 2) {
    while ((mf = MESH_Next_Face(mesh,&idx)))
      gids[i++] = MEnt_GEntID(mf);
    } 

    qsort(gids, num_side, sizeof(int), compare);
    int prev_gid = gids[0];
    /* first sweep counts the number of different node sets based on geometric entity id and dimension */
    for (i=0;i<num_side;i++) {
      if (gids[i] != prev_gid) { 
	num_side_set++;
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
	side_set_ids[0] = MEnt_GEntID(mf)+1;
      }
      else {
	me = MESH_Next_Edge(mesh, &idx);
	side_set_ids[0] = MEnt_GEntID(me)+1;
      }
      num_side_per_set[0] = num_side;
      return;
    }

    int *gids = (int *) MSTK_malloc(num_side*sizeof(int)); 
    
    if (side_dim == 2) {
      while ((mf = MESH_Next_Face(mesh,&idx)))
	gids[i++] = MEnt_GEntID(mf)+1;
    }
    if (side_dim == 1) {
      while ((me = MESH_Next_Edge(mesh, &idx)))
	gids[i++] = MEnt_GEntID(me)+1;
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


#ifdef __cplusplus
}
#endif

