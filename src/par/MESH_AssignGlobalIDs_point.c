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
     This function is a collective call
     It assigns each submesh the global IDs of vertices, edges, faces and regions

     It assumes the parallel neighboring information is already established

     If global IDs are already given, skip this function, call MESH_LabelPType()

     Author(s): Duo Wang, Rao Garimella
  */

int MESH_AssignGlobalIDs_point_Vertex(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);
int MESH_AssignGlobalIDs_point_Edge(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);
int MESH_AssignGlobalIDs_point_Face(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);
int MESH_AssignGlobalIDs_point_Region(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);


int MESH_AssignGlobalIDs_point(Mesh_ptr submesh, int rank, int num,  MPI_Comm comm) {
  int nf, nr;
  nf = MESH_Num_Faces(submesh);
  nr = MESH_Num_Regions(submesh);
  MESH_AssignGlobalIDs_point_Vertex(submesh, rank, num, comm);
  MESH_AssignGlobalIDs_point_Edge(submesh, rank, num, comm);
  if (nr)
    MESH_AssignGlobalIDs_point_Region(submesh, rank, num, comm);
  else if(nf) 
    MESH_AssignGlobalIDs_point_Face(submesh, rank, num, comm);
  else {
    MSTK_Report("MESH_AssignGlobalIDs()","only assign global id for volume or surface mesh",MSTK_ERROR);
    exit(-1);
  }
  return 1;
}


  /* 
     test if a vertex is on partition boundary, for surface mesh
     a vertex is on boundary if one of its incident edges has less than 2 neighboring faces
   */
static int vertex_on_face_boundary(MVertex_ptr mv) {
  int i, nve, ok = 0;
  MEdge_ptr me;
  List_ptr vedges = MV_Edges(mv);
  nve = List_Num_Entries(vedges);
  for(i = 0; i < nve; i++) {
    me = List_Entry(vedges,i);
    if(List_Num_Entries(ME_Faces(me)) <= 1)
      ok = 1;
  }
  List_Delete(vedges);
  return ok;
}

  /* 
     test if a vertex is on partition boundary, for volume mesh
     a vertex is on boundary if one of its incident faces has less than 2 neighboring regions
 */
static int vertex_on_region_boundary(MVertex_ptr mv) {
  int i, nrf, ok = 0;
  MFace_ptr mf;
  List_ptr vfaces = MV_Faces(mv);
  nrf = List_Num_Entries(vfaces);
  for(i = 0; i < nrf; i++) {
    mf = List_Entry(vfaces,i);
    if(List_Num_Entries(MF_Regions(mf)) <= 1)
      ok = 1;
  }
  List_Delete(vfaces);
  return ok;
}

 /* 
    Assign vertex global ID for both 2D and 3D meshes
    Also set PGHOST and POVERLAP, used for edge, face and region global ID 
    MESH_LabelPType_point() does not rely on PType set here, it will reset PType,
    so that user can call MESH_LabelPType_point() directly when global IDs are given
 */
     
int MESH_AssignGlobalIDs_point_Vertex(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  int i, j, nv, nbv, nov, ngv, nr, num_recv_procs, count;
  MVertex_ptr mv;
  List_ptr boundary_verts;
  MPI_Status status;
  RepType rtype;
  double coor[3], *loc;
  int max_nbv, index_nbv, iloc, global_id, mesh_info[10];
  int *global_mesh_info, *list_vertex, *recv_list_vertex, *mv_remote_info, *mv_ov_label;
  double *list_coor, *recv_list_coor;
  int (*func)(MVertex_ptr mv);                /* function pointer to check boundary vertex */
  for (i = 0; i < 10; i++) mesh_info[i] = 0;
  num_recv_procs = MESH_Num_GhostPrtns(submesh);

  rtype = MESH_RepType(submesh);
  nv = MESH_Num_Vertices(submesh);
  nr = MESH_Num_Regions(submesh);

  mesh_info[0] = rtype;
  mesh_info[1] = nv;
  mesh_info[4] = nr;

  if(nr) func = vertex_on_region_boundary;
  else   func = vertex_on_face_boundary;

  /* calculate number of boundary vertices */ 
  nbv = 0;  boundary_verts = List_New(10);
  for(i = 0; i < nv; i++) {
    mv = MESH_Vertex(submesh,i);
    if (func(mv)) {
      MV_Set_PType(mv,PBOUNDARY);
      List_Add(boundary_verts,mv);
      nbv++;
    }
  }
  mesh_info[5] = nbv;
  /* sort boundary vertices based on coordinate value, for binary search */
  List_Sort(boundary_verts,nbv,sizeof(MVertex_ptr),compareVertexCoor);

  global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  /* get largest number of boundary vertices of all the processors */
  max_nbv = 0;
  for(i = 0; i < num; i++)
    if(max_nbv < global_mesh_info[10*i+5])
      max_nbv = global_mesh_info[10*i+5];

  list_vertex = (int *)MSTK_malloc(nbv*sizeof(int));
  list_coor = (double *)MSTK_malloc(3*nbv*sizeof(double));
  recv_list_vertex = (int *)MSTK_malloc(max_nbv*sizeof(int));
  recv_list_coor = (double *)MSTK_malloc(3*max_nbv*sizeof(double));

  /* only coordinate values are sent */
  index_nbv = 0;
  for(i = 0; i < nbv; i++) {
    mv = List_Entry(boundary_verts,i);
    MV_Coords(mv,coor);
    list_coor[index_nbv*3] = coor[0];
    list_coor[index_nbv*3+1] = coor[1];
    list_coor[index_nbv*3+2] = coor[2];
    index_nbv++;
  }
  
  /* 
     used to store list id on incoming buffer of ghost vertex 
     No need to store processor id, already stored in master partition id
  */
  mv_remote_info = (int *)MSTK_malloc(nbv*sizeof(int));
  /* lable if a vertex is overlap */
  mv_ov_label = (int *)MSTK_malloc(max_nbv*sizeof(int));

  ngv = 0; nov = 0;
  /* label ghost and overlap vertex */
  for (i = 0; i < num; i++) {
    for(j = 0; j < max_nbv; j++) mv_ov_label[j] = 0;
    if(i == rank) continue;
    if(i < rank) {     
      if( MESH_Has_Ghosts_From_Prtn(submesh,i,MVERTEX) ) {
	MPI_Recv(recv_list_coor,3*max_nbv,MPI_DOUBLE,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_DOUBLE,&count);
	for(j = 0; j < nbv; j++) {
	  mv = List_Entry(boundary_verts,j);
	  if(MV_PType(mv) == PGHOST) continue;
	  MV_Coords(mv,coor);
	  loc = (double *)bsearch(&coor,
				  recv_list_coor,
				  count/3,
				  3*sizeof(double),
				  compareCoorDouble);
	  if(loc) {
	    iloc = (int)(loc - recv_list_coor)/3;
	    MV_Set_PType(mv,PGHOST);
	    MV_Set_MasterParID(mv,i);
	    ngv++;
	    mv_remote_info[j] = iloc;
	    mv_ov_label[iloc] = 1;
	  }
	}
	MPI_Send(mv_ov_label,count/3,MPI_INT,i,i,comm);
      }
    }
    if(i > rank) {     
      if( MESH_Has_Overlaps_On_Prtn(submesh,i,MVERTEX) ) {
	MPI_Send(list_coor,3*nbv,MPI_DOUBLE,i,i,comm);
	MPI_Recv(mv_ov_label,nbv,MPI_INT,i,rank,comm,&status);
	/* label overlap vertex */
	for(j = 0; j < nbv; j++) {
	  mv = List_Entry(boundary_verts,j);
	  if(mv_ov_label[j] && MV_PType(mv) != POVERLAP && MV_PType(mv) != PGHOST) {
	    nov++;
	    MV_Set_PType(mv,POVERLAP);
	    MV_Set_MasterParID(mv,rank);
	  }
	}
      }
    }
  }
  printf("num of ghost vertices %d, overlap vertices %d on rank %d\n", ngv, nov, rank);  
  mesh_info[9] = ngv;
  global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);
  
  /* Assign global ID for non ghost vertex */
  global_id = 1;
  for(i = 0; i < rank; i++) 
    global_id = global_id + global_mesh_info[10*i+1] - global_mesh_info[10*i+9];
  for(i = 0; i < nv; i++) {
    mv = MESH_Vertex(submesh,i);
    if (MV_PType(mv) == PGHOST) continue;
    MV_Set_GlobalID(mv,global_id++);
    MV_Set_MasterParID(mv,rank);
  }
  
  /* this time only global id are sent */
  index_nbv = 0;
  for(i = 0; i < nbv; i++) {
    mv = List_Entry(boundary_verts,i);
    list_vertex[index_nbv] = MV_GlobalID(mv);
    index_nbv++;
  }

  /* Assign ghost vertex global ID */
  for (i = 0; i < num; i++) {
    if(i == rank) continue;
    if(i < rank) {     
      if( MESH_Has_Ghosts_From_Prtn(submesh,i,MVERTEX) ) {
	MPI_Recv(recv_list_vertex,max_nbv,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	for(j = 0; j < nbv; j++) {
	  mv = List_Entry(boundary_verts,j);
	  if(MV_MasterParID(mv) != i) continue;
	  MV_Set_GlobalID(mv,recv_list_vertex[mv_remote_info[j]]);
	}
      }
    }
    if(i > rank) {     
      if( MESH_Has_Overlaps_On_Prtn(submesh,i,MVERTEX) )
	MPI_Send(list_vertex,nbv,MPI_INT,i,i,comm);
    }
  }

  List_Delete(boundary_verts);
  MSTK_free(global_mesh_info);
  MSTK_free(mv_remote_info);
  MSTK_free(mv_ov_label);
  MSTK_free(list_vertex);
  MSTK_free(list_coor);
  MSTK_free(recv_list_vertex);
  MSTK_free(recv_list_coor);
  return 1;
}


 /* 
    Assign edge global ID for both 2D and 3D meshes
    This function may produce more than necessary overlap edges
    Assume each processor knowns its overlap and ghost vertices 
 */
     
int MESH_AssignGlobalIDs_point_Edge(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  int i, j, k, nbe, noe, nge, ne, mesh_info[10], global_id, count;
  MVertex_ptr mv;
  MEdge_ptr me;
  List_ptr boundary_edges;
  MPI_Status status;
  int is_boundary;
  int *loc, edge_id[2],index_nbe, max_nbe, iloc;
  int *global_mesh_info, *list_edge, *recv_list_edge, *me_remote_info, *me_ov_label;
  for (i = 0; i < 10; i++) mesh_info[i] = 0;
  ne = MESH_Num_Edges(submesh);
  mesh_info[2] = ne;

  /* 
     collect 'boundary' edges
     if endpoints are either GHOST or OVERLAP, then it is a boundary edge 
  */
  nbe = 0;  boundary_edges = List_New(10);
  for(i = 0; i < ne; i++) {
    is_boundary = 1;       
    me = MESH_Edge(submesh,i);
    for( k = 0; k < 2; k++) {
      mv = ME_Vertex(me,k);
      if(MV_PType(mv)!=PGHOST && MV_PType(mv)!=POVERLAP)
	is_boundary = 0;
    }
    if(is_boundary) {
      ME_Set_PType(me,PBOUNDARY);
      List_Add(boundary_edges,me);
      nbe++;
    }
  }
  printf("num of boundary edges %d on rank %d\n", nbe, rank);  
  mesh_info[5] = nbe;

  List_Sort(boundary_edges,nbe,sizeof(MEdge_ptr),compareEdgeID);

  global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  /* get largest number of boundary vertices of all the processors */
  max_nbe = 0;
  for(i = 0; i < num; i++)
    if(max_nbe < global_mesh_info[10*i+5])
      max_nbe = global_mesh_info[10*i+5];

  list_edge = (int *)MSTK_malloc(4*nbe*sizeof(int));
  recv_list_edge = (int *)MSTK_malloc(4*max_nbe*sizeof(int));

  index_nbe = 0;
  for(i = 0; i < nbe; i++) {
    me = List_Entry(boundary_edges,i);
    list_edge[index_nbe] = MV_GlobalID(ME_Vertex(me,0));
    list_edge[index_nbe+1] = MV_GlobalID(ME_Vertex(me,1));
    list_edge[2*nbe+index_nbe] = ME_MasterParID(me);
    list_edge[2*nbe+index_nbe+1] = ME_GlobalID(me);
    index_nbe += 2;
  }
  
  /* 
     used to store list id on incoming buffer of ghost edge 
     No need to store processor id, already stored in master partition id
  */
  me_remote_info = (int *)MSTK_malloc(nbe*sizeof(int));
  /* lable if a edge is overlap */
  me_ov_label = (int *)MSTK_malloc(max_nbe*sizeof(int));
  nge = 0; noe = 0;
  /* label ghost and overlap edge */
  for (i = 0; i < num; i++) {
    for(j = 0; j < max_nbe; j++) me_ov_label[j] = 0;
    if(i == rank) continue;
    if(i < rank) {     
      if( MESH_Has_Ghosts_From_Prtn(submesh,i,MEDGE) ) {
	MPI_Recv(recv_list_edge,4*max_nbe,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	for(j = 0; j < nbe; j++) {
	  me = List_Entry(boundary_edges,j);
	  if(ME_GlobalID(me) > 0) continue;  /* if already assigned */
	  edge_id[0] = MV_GlobalID(ME_Vertex(me,0));
	  edge_id[1] = MV_GlobalID(ME_Vertex(me,1));
	  loc = (int *)bsearch(&edge_id,
			       recv_list_edge,
			       count/4,
			       2*sizeof(int),
			       compareEdgeINT);
	  if(loc) {
	    iloc = (int)(loc - recv_list_edge)/2;
	    ME_Set_PType(me,PGHOST);
	    ME_Set_MasterParID(me,i);
	    nge++;
	    me_remote_info[j] = iloc;
	    me_ov_label[iloc] = 1;
	  }
	}
	MPI_Send(me_ov_label,count/4,MPI_INT,i,i,comm);
      }
    }
    if(i > rank) {     
      if( MESH_Has_Overlaps_On_Prtn(submesh,i,MEDGE) ) {
	MPI_Send(list_edge,4*nbe,MPI_INT,i,i,comm);
	MPI_Recv(me_ov_label,nbe,MPI_INT,i,rank,comm,&status);
	/* label overlap edge */
	for(j = 0; j < nbe; j++) {
	  me = List_Entry(boundary_edges,j);
	  if(me_ov_label[j] && ME_PType(me) != POVERLAP && ME_PType(me) != PGHOST) {
	    noe++;
	    ME_Set_PType(me,POVERLAP);
	    ME_Set_MasterParID(me,rank);
	  }
	}
      }
    }
  }
  
  printf("num of ghost edges %d, overlap edges %d on rank %d\n", nge, noe, rank);  
  mesh_info[9] = nge;
  global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);
  
  /* Assign global ID for non ghost vertex */
  global_id = 1;
  for(i = 0; i < rank; i++) 
    global_id = global_id + global_mesh_info[10*i+2] - global_mesh_info[10*i+9];
  for(i = 0; i < ne; i++) {
    me = MESH_Edge(submesh,i);
    if (ME_PType(me) == PGHOST) continue;
    ME_Set_GlobalID(me,global_id++);
    ME_Set_MasterParID(me,rank);
  }

  /* this time only global id are sent */
  index_nbe = 0;
  for(i = 0; i < nbe; i++) {
    me = List_Entry(boundary_edges,i);
    list_edge[index_nbe] = ME_GlobalID(me);
    index_nbe++;
  }

  /* Assign ghost vertex global ID */
  for (i = 0; i < num; i++) {
    if(i == rank) continue;
    if(i < rank) {     
      if( MESH_Has_Ghosts_From_Prtn(submesh,i,MEDGE) ) {
	MPI_Recv(recv_list_edge,max_nbe,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	for(j = 0; j < nbe; j++) {
	  me = List_Entry(boundary_edges,j);
	  if(ME_MasterParID(me) != i) continue;
	  ME_Set_GlobalID(me,recv_list_edge[me_remote_info[j]]);
	}
      }
    }
    if(i > rank) {     
      if( MESH_Has_Overlaps_On_Prtn(submesh,i,MEDGE) )
	MPI_Send(list_edge,nbe,MPI_INT,i,i,comm);
    }
  }



  List_Delete(boundary_edges);
  MSTK_free(global_mesh_info);
  MSTK_free(me_remote_info);
  MSTK_free(me_ov_label);
  MSTK_free(list_edge);
  MSTK_free(recv_list_edge);
  return 1;
}

 /* 
    Assign face global ID for both 2D meshes
    Assume there are no overlapped faces
    Assume each processor knowns its overlap and ghost vertices 
 */

int MESH_AssignGlobalIDs_point_Face(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  int i, nf, global_id, mesh_info[10];
  MFace_ptr mf;
  int *global_mesh_info;
  for (i = 0; i < 10; i++) mesh_info[i] = 0;
  nf = MESH_Num_Faces(submesh);
  mesh_info[3] = nf;

  global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  /* calculate starting global id number for faces*/
  global_id = 1;
  for(i = 0; i < rank; i++) 
    global_id = global_id + global_mesh_info[10*i+3];
  for(i = 0; i < nf; i++) {
    mf = MESH_Face(submesh,i);
    MF_Set_PType(mf,PINTERIOR);
    MF_Set_GlobalID(mf,global_id++);
    MF_Set_MasterParID(mf,rank);
  }

  MSTK_free(global_mesh_info);
  return 1;
}


 /* 
    Assign region global ID for both 3D meshes
    Main work is to assign global ID for faces

    Assume there are no overlapped regions
    Assume each processor knowns its overlap and ghost vertices 
 */
     
int MESH_AssignGlobalIDs_point_Region(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  int i, j, k, nfv, nbf, nof, ngf, nf, nr, mesh_info[10], global_id, count;
  MVertex_ptr mv;
  MFace_ptr mf;
  MRegion_ptr mr;
  MPI_Status status;
  List_ptr boundary_faces, mfverts;
  RepType rtype;
  int *loc, face_id[MAXPV2+3],index_nbf, max_nbf, iloc, is_boundary;
  int *global_mesh_info, *list_face, *recv_list_face, *mf_remote_info, *mf_ov_label;
  for (i = 0; i < 10; i++) mesh_info[i] = 0;

  /* mesh_info store the mesh reptype, nv, ne, nf, nof */
  rtype = MESH_RepType(submesh);
  nf = MESH_Num_Faces(submesh);
  nr = MESH_Num_Regions(submesh);

  mesh_info[0] = rtype;
  mesh_info[3] = nf;
  mesh_info[4] = nr;

  /* 
     collect 'boundary' faces
     if all endpoints are either GHOST or OVERLAP, then it is a boundary face 
  */
  nbf = 0;  boundary_faces = List_New(10);
  for(i = 0; i < nf; i++) {
    is_boundary = 1;       
    mf = MESH_Face(submesh,i);
    mfverts = MF_Vertices(mf,1,0);
    nfv = List_Num_Entries(mfverts);
    for(j = 0; j < nfv; j++) {
      mv = List_Entry(mfverts,j);
      if(MV_PType(mv)!=PGHOST && MV_PType(mv)!=POVERLAP)
	is_boundary = 0;
    }
    if(is_boundary) {
      MF_Set_PType(mf,PBOUNDARY);
      List_Add(boundary_faces,mf);
      nbf++;
    }
  }

  mesh_info[5] = nbf;

  List_Sort(boundary_faces,nbf,sizeof(MFace_ptr),compareFaceID);

  global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  /* get largest number of boundary faces of all the processors */
  max_nbf = 0;
  for(i = 0; i < num; i++)
    if(max_nbf < global_mesh_info[10*i+5])
      max_nbf = global_mesh_info[10*i+5];

  list_face = (int *)MSTK_malloc((MAXPV2+3)*nbf*sizeof(int));
  recv_list_face = (int *)MSTK_malloc((MAXPV2+3)*max_nbf*sizeof(int));

  index_nbf = 0;
  for(i = 0; i < nbf; i++) {
    mf = List_Entry(boundary_faces,i);
    mfverts = MF_Vertices(mf,1,0);
    nfv = List_Num_Entries(mfverts);
    list_face[index_nbf] = nfv;
    for(j = 0; j < nfv; j++) 
      list_face[index_nbf+j+1] = MV_GlobalID(List_Entry(mfverts,j));
    list_face[index_nbf+nfv+1] = MF_MasterParID(mf);
    list_face[index_nbf+nfv+2] = MF_GlobalID(mf);
    index_nbf += MAXPV2+3;
    List_Delete(mfverts);
  }
  
  /* 
     used to store list id on incoming buffer of ghost edge 
     No need to store processor id, already stored in master partition id
  */
  mf_remote_info = (int *)MSTK_malloc(nbf*sizeof(int));
  /* lable if a edge is overlap */
  mf_ov_label = (int *)MSTK_malloc(max_nbf*sizeof(int));
  ngf = 0; nof = 0;
  /* label ghost and overlap edge */
  for (i = 0; i < num; i++) {
    for(j = 0; j < max_nbf; j++) mf_ov_label[j] = 0;
    if(i == rank) continue;
    if(i < rank) {     
      if( MESH_Has_Ghosts_From_Prtn(submesh,i,MFACE) ) {
	MPI_Recv(recv_list_face,(MAXPV2+3)*max_nbf,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	for(j = 0; j < nbf; j++) {
	  mf = List_Entry(boundary_faces,j);
	  if(MF_GlobalID(mf) > 0) continue;  /* if already assigned */
	  mfverts = MF_Vertices(mf,1,0);
	  nfv = List_Num_Entries(mfverts);
	  face_id[0] = nfv;                           
	  for(k = 0; k < nfv; k++) 
	    face_id[k+1] = MV_GlobalID(List_Entry(mfverts,k));
	  loc = (int *)bsearch(&face_id,
			       recv_list_face,
			       count/(MAXPV2+3),
			       (MAXPV2+3)*sizeof(int),
			       compareFaceINT);
	  if(loc) {
	    iloc = (int)(loc - recv_list_face)/(MAXPV2+3);
	    MF_Set_PType(mf,PGHOST);
	    MF_Set_MasterParID(mf,i);
	    ngf++;
	    mf_remote_info[j] = iloc;
	    mf_ov_label[iloc] = 1;
	  }
	  List_Delete(mfverts);
	}
	MPI_Send(mf_ov_label,count/(MAXPV2+3),MPI_INT,i,i,comm);
      }
    }
    if(i > rank) {     
      if( MESH_Has_Overlaps_On_Prtn(submesh,i,MFACE) ) {
	MPI_Send(list_face,(MAXPV2+3)*nbf,MPI_INT,i,i,comm);
	MPI_Recv(mf_ov_label,nbf,MPI_INT,i,rank,comm,&status);
	/* label overlap edge */
	for(j = 0; j < nbf; j++) {
	  mf = List_Entry(boundary_faces,j);
	  if(mf_ov_label[j] && MF_PType(mf) != POVERLAP && MF_PType(mf) != PGHOST) {
	    nof++;
	    MF_Set_PType(mf,POVERLAP);
	    MF_Set_MasterParID(mf,rank);
	  }
	}
      }
    }
  }
  
  printf("num of ghost faces %d, overlap faces %d on rank %d\n", ngf, nof, rank);  
  mesh_info[9] = ngf;
  global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);
  
  /* Assign global ID for non ghost vertex */
  global_id = 1;
  for(i = 0; i < rank; i++) 
    global_id = global_id + global_mesh_info[10*i+3] - global_mesh_info[10*i+9];
  for(i = 0; i < nf; i++) {
    mf = MESH_Face(submesh,i);
    if (MF_PType(mf) == PGHOST) continue;
    MF_Set_GlobalID(mf,global_id++);
    MF_Set_MasterParID(mf,rank);
  }

  /* this time only global id are sent */
  index_nbf = 0;
  for(i = 0; i < nbf; i++) {
    mf = List_Entry(boundary_faces,i);
    list_face[index_nbf] = MF_GlobalID(mf);
    index_nbf++;
  }

  /* Assign ghost vertex global ID */
  for (i = 0; i < num; i++) {
    if(i == rank) continue;
    if(i < rank) {     
      if( MESH_Has_Ghosts_From_Prtn(submesh,i,MFACE) ) {
	MPI_Recv(recv_list_face,max_nbf,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	for(j = 0; j < nbf; j++) {
	  mf = List_Entry(boundary_faces,j);
	  if(MF_MasterParID(mf) != i) continue;
	  MF_Set_GlobalID(mf,recv_list_face[mf_remote_info[j]]);
	}
      }
    }
    if(i > rank) {     
      if( MESH_Has_Overlaps_On_Prtn(submesh,i,MFACE) )
	MPI_Send(list_face,nbf,MPI_INT,i,i,comm);
    }
  }

  /* assign region global id */
  global_id = 1;
  for(i = 0; i < rank; i++) 
    global_id = global_id + global_mesh_info[10*i+4];
  for(i = 0; i < nr; i++) {
    mr = MESH_Region(submesh,i);
    MR_Set_PType(mr,PINTERIOR);
    MR_Set_GlobalID(mr,global_id++);
    MR_Set_MasterParID(mr,rank);
  }

  List_Delete(boundary_faces);
  MSTK_free(global_mesh_info);
  MSTK_free(mf_remote_info);
  MSTK_free(mf_ov_label);
  MSTK_free(list_face);
  MSTK_free(recv_list_face);

  return 1;
}


#ifdef __cplusplus
}
#endif

