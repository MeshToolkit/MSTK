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
  /*
  if (nr)
    MESH_AssignGlobalIDs_point_Region(submesh, rank, num, comm);
  else if(nf) 
    MESH_AssignGlobalIDs_point_Face(submesh, rank, num, comm);
  else {
    MSTK_Report("MESH_AssignGlobalIDs()","only assign global id for volume or surface mesh",MSTK_ERROR);
    exit(-1);
  }
  */
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
  int i, j, nv, nbv, nr, num_recv_procs, count;
  MVertex_ptr mv;
  List_ptr boundary_verts;
  MPI_Status status;
  RepType rtype;
  double coor[3], *loc;
  int max_nbv, index_nbv, iloc, num_ghost_verts, global_id, mesh_info[10];
  int *global_mesh_info, *list_vertex, *recv_list_vertex, *mv_remote_info;
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

  num_ghost_verts = 0;

  /* label ghost and overlap vertex */
  for (i = 0; i < num; i++) {
    if(i == rank) continue;
    if(i < rank) {     
      if( MESH_Has_Overlaps_On_Prtn(submesh,i,MVERTEX) ) {
	MPI_Recv(recv_list_coor,3*nv,MPI_DOUBLE,i,rank,comm,&status);
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
	    num_ghost_verts++;
	    mv_remote_info[j] = iloc;
	  }
	}
      }
    }
    if(i > rank) {     
      if( MESH_Has_Ghosts_From_Prtn(submesh,i,MVERTEX) )
	MPI_Send(list_coor,3*nbv,MPI_DOUBLE,i,i,comm);
    }
  }
  printf("num of ghost vertices on rank %d is %d\n",rank, num_ghost_verts);  
  mesh_info[9] = num_ghost_verts;
  global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);
  
  /* Assign global ID for non ghost vertex */
  global_id = 1;
  for(i = 0; i < rank; i++) 
    global_id = global_id + global_mesh_info[10*i+1] - global_mesh_info[10*i+9];
  for(i = 0; i < nv; i++) {
    mv = MESH_Vertex(submesh,i);
    if (MV_PType(mv) == PGHOST)
      continue;
    MV_Set_GlobalID(mv,global_id++);
    MV_Set_MasterParID(mv,rank);
  }
  
  /* this time only global id are sent */
  index_nbv = 0;
  for(i = 0; i < nbv; i++) {
    mv = List_Entry(boundary_verts,i);
    list_vertex[index_nbv] = MV_GlobalID(mv);
    index_nbv++;
    if(MV_PType(mv) != PGHOST) {
      MV_Set_PType(mv,POVERLAP);
      MV_Set_MasterParID(mv,rank);
    }
  }

  /* Assign ghost vertex global ID */
  for (i = 0; i < num; i++) {
    if(i == rank) continue;
    if(i < rank) {     
      if( MESH_Has_Overlaps_On_Prtn(submesh,i,MVERTEX) ) {
	MPI_Recv(recv_list_vertex,nv,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	for(j = 0; j < nbv; j++) {
	  mv = List_Entry(boundary_verts,j);
	  if(MV_MasterParID(mv) != i) continue;
	  MV_Set_GlobalID(mv,recv_list_vertex[mv_remote_info[j]]);
	}
      }
    }
    if(i > rank) {     
      if( MESH_Has_Ghosts_From_Prtn(submesh,i,MVERTEX) )
	MPI_Send(list_vertex,nbv,MPI_INT,i,i,comm);
    }
  }

  List_Delete(boundary_verts);
  MSTK_free(global_mesh_info);
  MSTK_free(mv_remote_info);

  MSTK_free(list_vertex);
  MSTK_free(list_coor);
  MSTK_free(recv_list_vertex);
  MSTK_free(recv_list_coor);
  return 1;
}


 /* 
    Assign edge global ID for both 2D and 3D meshes

    Assume each processor knowns its overlap and ghost vertices 
 */
     
int MESH_AssignGlobalIDs_point_Edge(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  int i, j, noe, nge, ne, mesh_info[10], global_id;
  MEdge_ptr me;
  List_ptr overlap_edges, ghost_edges;
  RepType rtype;
  int *loc, edge_id[2],index_noe, max_noe, iloc;
  int *global_mesh_info, *list_overlap_edge, *recv_list_edge;
  for (i = 0; i < 10; i++) mesh_info[i] = 0;

  rtype = MESH_RepType(submesh);
  ne = MESH_Num_Edges(submesh);

  mesh_info[0] = rtype;
  mesh_info[2] = ne;

  /* calculate number of GHOST and OVERLAP edges */ 
  noe = 0; nge = 0;
  overlap_edges = List_New(10);
  ghost_edges = List_New(10);
  for(i = 0; i < ne; i++) {
    me = MESH_Edge(submesh,i);
    /* if both ends are ghost vertices, then the edge is a ghost*/
    if(MV_PType(ME_Vertex(me,0)) == PGHOST && MV_PType(ME_Vertex(me,1)) == PGHOST) {
      List_Add(ghost_edges,me);
      ME_Set_PType(me,PGHOST);
      nge++;
    }

    /* 
       if either ends are overlap vertices, then the edge is a overlap
       this may produce more than necessary overlap edges 
    */
    if(MV_PType(ME_Vertex(me,0)) == POVERLAP || MV_PType(ME_Vertex(me,1)) == POVERLAP) {
      ME_Set_PType(me,POVERLAP);
      List_Add(overlap_edges,me);
      noe++;
    }
  }
  mesh_info[5] = nge;
  mesh_info[6] = noe;
  
  /* sort overlap edges based on coordinate value, for binary search */
  List_Sort(overlap_edges,noe,sizeof(MEdge_ptr),compareEdgeID);

  global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  /* Assign global id to non-ghost edges */
  global_id = 1;
  for(i = 0; i < rank; i++) 
    global_id = global_id + global_mesh_info[10*i+2] - global_mesh_info[10*i+5];
  for(i = 0; i < ne; i++) {
    me = MESH_Edge(submesh,i);
    if(ME_PType(me) != PGHOST) {
      ME_Set_GlobalID(me,global_id++);
      ME_Set_MasterParID(me,rank);
    }
  }

  max_noe = 0;
  for(i = 0; i < num; i++)
    if(max_noe < global_mesh_info[10*i+6])
      max_noe = global_mesh_info[10*i+6];

  list_overlap_edge = (int *)MSTK_malloc(max_noe*4*sizeof(int));

  recv_list_edge = (int *)MSTK_malloc(num*max_noe*4*sizeof(int));

  /* pack edge information to send  */
  index_noe = 0;
  for(i = 0; i < noe; i++) {
    me = List_Entry(overlap_edges,i);
    list_overlap_edge[index_noe] = MV_GlobalID(ME_Vertex(me,0));
    list_overlap_edge[index_noe+1] = MV_GlobalID(ME_Vertex(me,1));
    list_overlap_edge[2*noe+index_noe] = ME_MasterParID(me);
    list_overlap_edge[2*noe+index_noe+1] = ME_GlobalID(me);
    index_noe+=2;
  }

  MPI_Allgather(list_overlap_edge,4*max_noe,MPI_INT,recv_list_edge,4*max_noe,MPI_INT,comm);

  /* for processor other than 0 */
  if(rank > 0) {
    for(i = 0; i < nge; i++) {
      me = List_Entry(ghost_edges,i);
      edge_id[0] = MV_GlobalID(ME_Vertex(me,0));
      edge_id[1] = MV_GlobalID(ME_Vertex(me,1));
      for(j = 0; j < rank; j++) {
	loc = (int *)bsearch(&edge_id,
				&recv_list_edge[4*max_noe*j],
				global_mesh_info[10*j+6],
				2*sizeof(int),
				compareEdgeINT);
	/* if found the vertex on previous processors */
	if(loc) {
	  /* here the location iloc is relative to the beginning of the jth processor */
	  iloc = (int)(loc - &recv_list_edge[4*max_noe*j]);
	  ME_Set_MasterParID(me,recv_list_edge[4*max_noe*j+2*global_mesh_info[10*j+6]+iloc]);
	  ME_Set_GlobalID(me,recv_list_edge[4*max_noe*j+2*global_mesh_info[10*j+6]+iloc+1]);
	  break;
	}
      }
    }
  }

  List_Delete(overlap_edges);
  MSTK_free(global_mesh_info);

  MSTK_free(list_overlap_edge);
  MSTK_free(recv_list_edge);
  return 1;
}

 /* 
    Assign face global ID for both 2D meshes
    Assume there are no overlapped faces
    Assume each processor knowns its overlap and ghost vertices 
 */

int MESH_AssignGlobalIDs_point_Face(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  int i, nv, ne, nf, global_id, mesh_info[10];
  MFace_ptr mf;
  RepType rtype;
  int *global_mesh_info;
  for (i = 0; i < 10; i++) mesh_info[i] = 0;

  rtype = MESH_RepType(submesh);
  nv = MESH_Num_Vertices(submesh);
  ne = MESH_Num_Edges(submesh);
  nf = MESH_Num_Faces(submesh);

  mesh_info[0] = rtype;
  mesh_info[1] = nv;
  mesh_info[2] = ne;
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
  int i, j, k, nv, nfv, nof, ngf, ne, nf, nr, mesh_info[10], global_id;
  MVertex_ptr mv;
  MFace_ptr mf;
  MRegion_ptr mr;
  List_ptr overlap_faces, ghost_faces, mfverts;
  RepType rtype;
  int *loc, face_id[MAXPV2+3],index_nof, max_nof, iloc;
  int *global_mesh_info, *list_overlap_face, *recv_list_face;
  int is_ghost, is_overlap;
  for (i = 0; i < 10; i++) mesh_info[i] = 0;

  /* mesh_info store the mesh reptype, nv, ne, nf, nof */
  rtype = MESH_RepType(submesh);
  nv = MESH_Num_Vertices(submesh);
  ne = MESH_Num_Edges(submesh);
  nf = MESH_Num_Faces(submesh);
  nr = MESH_Num_Regions(submesh);

  mesh_info[0] = rtype;
  mesh_info[1] = nv;
  mesh_info[2] = ne;
  mesh_info[3] = nf;
  mesh_info[4] = nr;

  /* calculate number of GHOST and OVERLAP faces */ 
  nof = 0; ngf = 0;
  overlap_faces = List_New(10);
  ghost_faces = List_New(10);
  for(i = 0; i < nf; i++) {
    mf = MESH_Face(submesh,i);
    is_ghost = 1; is_overlap = 0;
    mfverts = MF_Vertices(mf,1,0);
    nfv = List_Num_Entries(mfverts);
    for(j = 0; j < nfv; j++) {
      mv = List_Entry(mfverts,j);
      if(MV_PType(mv) != PGHOST)  /* if all vertices are ghost, then the face is a ghost*/
	is_ghost = 0;
      if(MV_PType(mv) == POVERLAP) /*  if either vertex is a overlap, then the face is a overlap */
	is_overlap = 1;
    }
    if(is_ghost) {
      List_Add(ghost_faces,mf);
      MF_Set_PType(mf,PGHOST);
      ngf++;
    }

    if(is_overlap) {
      List_Add(overlap_faces,mf);
      MF_Set_PType(mf,POVERLAP);
      nof++;
    }
    List_Delete(mfverts);
  }
  mesh_info[5] = ngf;
  mesh_info[6] = nof;
  mesh_info[7] = nr;
  /* sort overlap faces based on coordinate value, for binary search */
  List_Sort(overlap_faces,nof,sizeof(MFace_ptr),compareFaceID);
  /* 
     gather submeshes information
  */
  global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  /* Assign global id to non-ghost faces */
  /* calculate starting global id number for faces */
  global_id = 1;
  for(i = 0; i < rank; i++) 
    global_id = global_id + global_mesh_info[3] - global_mesh_info[10*i+5];
  for(i = 0; i < nf; i++) {
    mf = MESH_Face(submesh,i);
    if(MF_PType(mf) != PGHOST) {
      MF_Set_GlobalID(mf,global_id++);
      MF_Set_MasterParID(mf,rank);
    }
  }


  /* get largest number of boundary vertices of all the processors */
  max_nof = 0;
  for(i = 0; i < num; i++)
    if(max_nof < global_mesh_info[10*i+6])
      max_nof = global_mesh_info[10*i+6];

  list_overlap_face = (int *)MSTK_malloc(max_nof*(MAXPV2+3)*sizeof(int));

  recv_list_face = (int *)MSTK_malloc(num*max_nof*(MAXPV2+3)*sizeof(int));

  /* pack face information to send  */
  index_nof = 0;
  for(i = 0; i < nof; i++) {
    mf = List_Entry(overlap_faces,i);
    mfverts = MF_Vertices(mf,1,0);
    nfv = List_Num_Entries(mfverts);
    list_overlap_face[index_nof] = nfv;
    for(j = 0; j < nfv; j++) 
      list_overlap_face[index_nof+j+1] = MV_GlobalID(List_Entry(mfverts,j));

    list_overlap_face[index_nof+nfv+1] = MF_MasterParID(mf);
    list_overlap_face[index_nof+nfv+2] = MF_GlobalID(mf);
    index_nof+=MAXPV2+3;
    List_Delete(mfverts);
  }

  /* gather overlap faces */
  MPI_Allgather(list_overlap_face,(MAXPV2+3)*max_nof,MPI_INT,recv_list_face,(MAXPV2+3)*max_nof,MPI_INT,comm);

  /* for processor other than 0 */
  if(rank > 0) {
    for(i = 0; i < ngf; i++) {
      mf = List_Entry(ghost_faces,i);
      mfverts = MF_Vertices(mf,1,0);
      nfv = List_Num_Entries(mfverts);
      face_id[0] = nfv;                           
      for(k = 0; k < nfv; k++) 
	face_id[k+1] = MV_GlobalID(List_Entry(mfverts,k));
      
      for(j = 0; j < rank; j++) {
	loc = (int *)bsearch(&face_id,
			     &recv_list_face[(MAXPV2+3)*max_nof*j],
			     global_mesh_info[10*j+6],
			     (MAXPV2+3)*sizeof(int),
			     compareFaceINT);
	/* if found the face on previous processors */
	if(loc) {
	  /* here the location iloc is relative to the beginning of the jth processor */
	  iloc = (int)(loc - &recv_list_face[(MAXPV2+3)*max_nof*j]);
	  MF_Set_MasterParID(mf,recv_list_face[(MAXPV2+3)*max_nof*j+iloc+nfv+1]);
	  MF_Set_GlobalID(mf,recv_list_face[(MAXPV2+3)*max_nof*j+iloc+nfv+2]);
	  break;
	}
      }
    }
  }

  /* assign region global id */
  global_id = 1;
  for(i = 0; i < rank; i++) 
    global_id = global_id + global_mesh_info[10*i+7];
  for(i = 0; i < nr; i++) {
    mr = MESH_Region(submesh,i);
    MR_Set_PType(mr,PINTERIOR);
    MR_Set_GlobalID(mr,global_id++);
    MR_Set_MasterParID(mr,rank);
  }

  List_Delete(overlap_faces);
  List_Delete(ghost_faces);
  MSTK_free(global_mesh_info);

  MSTK_free(list_overlap_face);
  MSTK_free(recv_list_face);
  return 1;
}


#ifdef __cplusplus
}
#endif

