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
     It assigns each submesh the global IDs of vertices and elements
     
     It also labels ghost and overlap entities
     If global IDs are already given, call MESH_LabelPType()

     Author(s): Duo Wang, Rao Garimella
  */

int MESH_AssignGlobalIDs_Vertex(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);
int MESH_AssignGlobalIDs_Edge(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);
int MESH_AssignGlobalIDs_Face(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);
int MESH_AssignGlobalIDs_Region(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);


int MESH_AssignGlobalIDs(Mesh_ptr submesh, int rank, int num,  MPI_Comm comm) {
  int nf, nr;
  RepType rtype;

  /* basic mesh information */

  rtype = MESH_RepType(submesh);
  nf = MESH_Num_Faces(submesh);
  nr = MESH_Num_Regions(submesh);
  /* 
     build geometric entity dimension, mark boundary vertices
     MESH_BuildVertexClassfn(submesh);
 */

  /* Assign Vertex Global ID */
  MESH_AssignGlobalIDs_Vertex(submesh, rank, num, comm);
  MESH_AssignGlobalIDs_Edge(submesh, rank, num, comm);
  if (nr)
    MESH_AssignGlobalIDs_Region(submesh, rank, num, comm);
  else if(nf) 
    MESH_AssignGlobalIDs_Face(submesh, rank, num, comm);
  else {
    MSTK_Report("MESH_AssignGlobalIDs()","only send volume or surface mesh",MSTK_ERROR);
    exit(-1);
  }
  return 1;
}



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

static int face_on_region_boundary(MFace_ptr mf) {
    int ok = 0;
    if( List_Num_Entries(MF_Regions(mf))<=1 )
      ok = 1;
    return ok;
}

 /* 
    Assign vertex global ID for both 2D and 3D meshes
    First build boundary list of vertices and broadcast them
    
 */
     
int MESH_AssignGlobalIDs_Vertex(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  int i, j, nv, nbv, ne, nf, nr, mesh_info[10];
  MVertex_ptr mv;
  MEdge_ptr me;
  MFace_ptr mf;
  List_ptr boundary_verts;
  RepType rtype;
  double coor[3], *loc;
  int index_nbv, max_nbv, iloc, num_ghost_verts, global_id;
  int *global_mesh_info, *list_boundary_vertex, *recv_list_vertex, *vertex_ov_label, *vertex_ov_global_id, *id_on_ov_list;
  double *list_boundary_coor, *recv_list_coor;
  int (*func)(MVertex_ptr mv);                /* function pointer to check boundary vertex */
  for (i = 0; i < 10; i++) mesh_info[i] = 0;

  /* mesh_info store the mesh reptype, nv, ne, nf, nbv */
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

  if(nr) func = vertex_on_region_boundary;
  else   func = vertex_on_face_boundary;
  /* calculate number of boundary vertices */ 
  nbv = 0;
  boundary_verts = List_New(10);
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
  /* 
     gather submeshes information
     right now we only need nv and nbv, and later num_ghost_verts, but we gather all mesh_info
  */
  global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  /* get largest number of boundary vertices of all the processors */
  max_nbv = 0;
  for(i = 0; i < num; i++)
    if(max_nbv < global_mesh_info[10*i+5])
      max_nbv = global_mesh_info[10*i+5];

  list_boundary_vertex = (int *)MSTK_malloc(max_nbv*sizeof(int));
  list_boundary_coor = (double *)MSTK_malloc(3*max_nbv*sizeof(double));

  recv_list_vertex = (int *)MSTK_malloc(num*max_nbv*sizeof(int));
  recv_list_coor = (double *)MSTK_malloc(3*num*max_nbv*sizeof(double));

  /* only local id and coordinate values are sent */
  index_nbv = 0;
  for(i = 0; i < nbv; i++) {
    mv = List_Entry(boundary_verts,i);
    list_boundary_vertex[index_nbv] = MV_ID(mv);
    MV_Coords(mv,coor);
    list_boundary_coor[index_nbv*3] = coor[0];
    list_boundary_coor[index_nbv*3+1] = coor[1];
    list_boundary_coor[index_nbv*3+2] = coor[2];
    index_nbv++;
  }

  /* gather boundary vertices */
  MPI_Allgather(list_boundary_vertex,max_nbv,MPI_INT,recv_list_vertex,max_nbv,MPI_INT,comm);
  MPI_Allgather(list_boundary_coor,3*max_nbv,MPI_DOUBLE,recv_list_coor,3*max_nbv,MPI_DOUBLE,comm);

  /* indicate if a vertex is overlapped */
  vertex_ov_label = (int *)MSTK_malloc(num*max_nbv*sizeof(int));

  /* 
     store the local boundary id on ov processor
     it is used to assign global id of local ghost vertices
     no need to store master partition id, MV_MasterParID(mv) is already assigned
  */
  id_on_ov_list = (int *)MSTK_malloc(max_nbv*sizeof(int));

  for (i = 0; i < num*max_nbv; i++)
    vertex_ov_label[i] = 0;
  num_ghost_verts = 0;
  /* for processor other than 0 */
  if(rank > 0) {
    for(i = 0; i < nbv; i++) {
      mv = List_Entry(boundary_verts,i);
      MV_Coords(mv,coor);
      /* check which previous processor has the same coordinate vertex */
      for(j = 0; j < rank; j++) {
	/* since each processor has sorted the boundary vertices, use binary search */
	loc = (double *)bsearch(&coor,
			     &recv_list_coor[3*max_nbv*j],
				global_mesh_info[10*j+5],
				3*sizeof(double),
				compareCoorDouble);
	/* if found the vertex on previous processors */
	if(loc) {
	  /* here the location iloc is relative to the beginning of the jth processor */
	  iloc = (int)(loc - &recv_list_coor[3*max_nbv*j])/3;
	  MV_Set_PType(mv,PGHOST);
	  MV_Set_MasterParID(mv,j);
	  
	  /*
	    printf("Assign: rank %d,boundary vertex local ID %d, found on processor %d, loc %d ",rank, MV_ID(mv), j, iloc);
	    printf("coor (%lf,%lf,%lf)\n",coor[0],coor[1],coor[2]);
	  */
	  num_ghost_verts++;
	  /* label the original vertex as overlapped */
	  vertex_ov_label[max_nbv*j+iloc] |= 1;
	  id_on_ov_list[i] = iloc;
	  /* if found on processor j, no need to test for j+1,j+2...*/
	  break;
	}
      }
    }
  }
  /* num of ghost verts */
  mesh_info[9] = num_ghost_verts;
  /* update ghost verts number */
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);
  /* since this is a OR reduction, we can use MPI_IN_PLACE, send buffer same as recv buffer */
  MPI_Allreduce(MPI_IN_PLACE,vertex_ov_label,num*max_nbv,MPI_INT,MPI_LOR,comm);    

  /* calculate starting global id number for vertices*/
  global_id = 1;
  for(i = 0; i < rank; i++) 
    global_id = global_id + global_mesh_info[10*i+1]-global_mesh_info[10*i+9];
  for(i = 0; i < nv; i++) {
    mv = MESH_Vertex(submesh,i);
    if (MV_PType(mv) == PGHOST)
      continue;
    MV_Set_GlobalID(mv,global_id++);
    MV_Set_MasterParID(mv,rank);
  }

      

  /* store overlapped vertices IDs and broadast*/
  vertex_ov_global_id = (int *)MSTK_malloc(num*max_nbv*sizeof(int));
  for(i = 0; i < num*max_nbv; i++) 
    vertex_ov_global_id[i] = 0;
  for(i = 0; i < nbv; i++) {
    if(vertex_ov_label[rank*max_nbv+i]) {
      mv = List_Entry(boundary_verts,i);
      MV_Set_PType(mv,POVERLAP);
      vertex_ov_global_id[rank*max_nbv+i] = MV_GlobalID(mv);
    }
  }

  MPI_Allreduce(MPI_IN_PLACE,vertex_ov_global_id,num*max_nbv,MPI_INT,MPI_MAX,comm);    

  for(i = 0; i < nbv; i++) {
    mv = List_Entry(boundary_verts,i);
    if(MV_PType(mv) == PGHOST) 
      MV_Set_GlobalID(mv,vertex_ov_global_id[MV_MasterParID(mv)*max_nbv+id_on_ov_list[i]]);
  }



  List_Delete(boundary_verts);
  MSTK_free(global_mesh_info);
  MSTK_free(vertex_ov_label);
  MSTK_free(vertex_ov_global_id);
  MSTK_free(id_on_ov_list);

  MSTK_free(list_boundary_vertex);
  MSTK_free(list_boundary_coor);
  MSTK_free(recv_list_vertex);
  MSTK_free(recv_list_coor);
  return 1;
}


 /* 
    as of now, each processor knowns its overlap and ghost vertices 

    
 */
     
int MESH_AssignGlobalIDs_Edge(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  int i, j, nv, noe, nge, ne, nf, nr, mesh_info[10], global_id;
  MVertex_ptr mv;
  MEdge_ptr me;
  MFace_ptr mf;
  List_ptr overlap_edges, ghost_edges;
  RepType rtype;
  double coor[3];
  int *loc, edge_id[2],index_noe, max_noe, iloc, num_ghost_verts;
  int *global_mesh_info, *list_overlap_edge, *recv_list_edge;
  double *list_boundary_coor, *recv_list_coor;
  printf("assign edge id on rank %d\n",rank);
  for (i = 0; i < 10; i++) mesh_info[i] = 0;

  /* mesh_info store the mesh reptype, nv, ne, nf, noe */
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
  /* 
     gather submeshes information
  */
  global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  /* Assign global id to non-ghost edges */
  /* calculate starting global id number for edges */
  global_id = 1;
  for(i = 0; i < rank; i++) 
    global_id = global_id + global_mesh_info[2] - global_mesh_info[10*i+5];
  for(i = 0; i < ne; i++) {
    me = MESH_Edge(submesh,i);
    if(ME_PType(me) != PGHOST) {
      ME_Set_GlobalID(me,global_id++);
      ME_Set_MasterParID(me,rank);
    }
  }


  /* get largest number of boundary vertices of all the processors */
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

  /* gather overlap edges */
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
	  /*
	    printf("Assign: rank %d,boundary vertex local ID %d, found on processor %d, loc %d ",rank, MV_ID(mv), j, iloc);
	    printf("coor (%lf,%lf,%lf)\n",coor[0],coor[1],coor[2]);
	  */
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

  /* right now assume there are no overlapped faces */
int MESH_AssignGlobalIDs_Face(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  int i, nv, ne, nf, global_id, mesh_info[10];
  MFace_ptr mf;
  RepType rtype;
  int *global_mesh_info;


  for (i = 0; i < 10; i++) mesh_info[i] = 0;

  /* mesh_info store the mesh reptype, nv, ne, nf, nbf */
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

int MESH_AssignGlobalIDs_Region_Face(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  int i, j, nv, nbf, ne, nf, nr, mesh_info[10];
  MVertex_ptr mv;
  MFace_ptr mf;
  List_ptr boundary_faces;
  RepType rtype;
  double coor[3], *loc;
  int index_nbf, max_nbf, iloc, num_ghost_verts, global_id;
  int *global_mesh_info, *list_boundary_vertex, *recv_list_vertex, *vertex_ov_label, *vertex_ov_global_id, *id_on_ov_list;
  double *list_boundary_coor, *recv_list_coor;

  for (i = 0; i < 10; i++) mesh_info[i] = 0;

  /* mesh_info store the mesh reptype, nv, ne, nf, nbf */
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

  /* calculate number of boundary vertices */ 
  nbf = 0;
  boundary_faces = List_New(10);
  for(i = 0; i < nf; i++) {
    mf = MESH_Face(submesh,i);
    if (face_on_region_boundary(mf)) {
      MF_Set_PType(mf,PBOUNDARY);
      List_Add(boundary_faces,mf);
      nbf++;
    }
  }
  mesh_info[5] = nbf;
  
  /* sort boundary vertices based on coordinate value, for binary search */
  List_Sort(boundary_faces,nbf,sizeof(MVertex_ptr),compareVertexCoor);
  /* 
     gather submeshes information
     right now we only need nv and nbf, and later num_ghost_verts, but we gather all mesh_info
  */
  global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  /* get largest number of boundary vertices of all the processors */
  max_nbf = 0;
  for(i = 0; i < num; i++)
    if(max_nbf < global_mesh_info[10*i+5])
      max_nbf = global_mesh_info[10*i+5];

  list_boundary_vertex = (int *)MSTK_malloc(max_nbf*sizeof(int));
  list_boundary_coor = (double *)MSTK_malloc(3*max_nbf*sizeof(double));

  recv_list_vertex = (int *)MSTK_malloc(num*max_nbf*sizeof(int));
  recv_list_coor = (double *)MSTK_malloc(3*num*max_nbf*sizeof(double));

  /* only local id and coordinate values are sent */
  index_nbf = 0;
  for(i = 0; i < nbf; i++) {
    mv = List_Entry(boundary_faces,i);
    list_boundary_vertex[index_nbf] = MV_ID(mv);
    MV_Coords(mv,coor);
    list_boundary_coor[index_nbf*3] = coor[0];
    list_boundary_coor[index_nbf*3+1] = coor[1];
    list_boundary_coor[index_nbf*3+2] = coor[2];
    index_nbf++;
  }

  /* gather boundary vertices */
  MPI_Allgather(list_boundary_vertex,max_nbf,MPI_INT,recv_list_vertex,max_nbf,MPI_INT,comm);
  MPI_Allgather(list_boundary_coor,3*max_nbf,MPI_DOUBLE,recv_list_coor,3*max_nbf,MPI_DOUBLE,comm);

  /* indicate if a vertex is overlapped */
  vertex_ov_label = (int *)MSTK_malloc(num*max_nbf*sizeof(int));

  /* 
     store the local boundary id on ov processor
     it is used to assign global id of local ghost vertices
     no need to store master partition id, MV_MasterParID(mv) is already assigned
  */
  id_on_ov_list = (int *)MSTK_malloc(max_nbf*sizeof(int));

  for (i = 0; i < num*max_nbf; i++)
    vertex_ov_label[i] = 0;
  num_ghost_verts = 0;
  /* for processor other than 0 */
  if(rank > 0) {
    for(i = 0; i < nbf; i++) {
      mv = List_Entry(boundary_faces,i);
      MV_Coords(mv,coor);
      /* check which previous processor has the same coordinate vertex */
      for(j = 0; j < rank; j++) {
	/* since each processor has sorted the boundary vertices, use binary search */
	loc = (double *)bsearch(&coor,
			     &recv_list_coor[3*max_nbf*j],
			     global_mesh_info[10*j+5],
			     3*sizeof(double),
			     compareCoorDouble);
	/* if found the vertex on previous processors */
	if(loc) {
	  /* here the location iloc is relative to the beginning of the jth processor */
	  iloc = (int)(loc - &recv_list_coor[3*max_nbf*j])/3;
	  MV_Set_PType(mv,PGHOST);
	  MV_Set_MasterParID(mv,j);
	  /*
	  printf("Assign: rank %d,boundary vertex local ID %d, found on processor %d, loc %d ",rank, MV_ID(mv), j, iloc);
	  printf("coor (%lf,%lf,%lf)\n",coor[0],coor[1],coor[2]);
	  */
	  num_ghost_verts++;
	  /* label the original vertex as overlapped */
	  vertex_ov_label[max_nbf*j+iloc] |= 1;
	  id_on_ov_list[i] = iloc;
	  /* if found on processor j, no need to test for j+1,j+2...*/
	  break;
	}
      }
    }

    /* num of ghost verts */
  
  }

  mesh_info[9] = num_ghost_verts;
  /* update ghost verts number */
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);
  /* since this is a OR reduction, we can use MPI_IN_PLACE, send buffer same as recv buffer */
  MPI_Allreduce(MPI_IN_PLACE,vertex_ov_label,num*max_nbf,MPI_INT,MPI_LOR,comm);    

  /* calculate starting global id number for vertices*/
  global_id = 1;
  for(i = 0; i < rank; i++) 
    global_id = global_id + global_mesh_info[10*i+1]-global_mesh_info[10*i+9];
  for(i = 0; i < nv; i++) {
    mv = MESH_Vertex(submesh,i);
    if (MV_PType(mv) == PGHOST)
      continue;
    MV_Set_GlobalID(mv,global_id++);
    MV_Set_MasterParID(mv,rank);
  }


  /* store overlapped vertices IDs and broadast*/
  vertex_ov_global_id = (int *)MSTK_malloc(num*max_nbf*sizeof(int));
  for(i = 0; i < num*max_nbf; i++) 
    vertex_ov_global_id[i] = 0;
  for(i = 0; i < nbf; i++) {
    if(vertex_ov_label[rank*max_nbf+i]) {
      mv = List_Entry(boundary_faces,i);
      MV_Set_PType(mv,POVERLAP);
      vertex_ov_global_id[rank*max_nbf+i] = MV_GlobalID(mv);
    }
  }
  MPI_Allreduce(MPI_IN_PLACE,vertex_ov_global_id,num*max_nbf,MPI_INT,MPI_MAX,comm);    
  for(i = 0; i < nbf; i++) {
    mv = List_Entry(boundary_faces,i);
    if(MV_PType(mv) == PGHOST) 
      MV_Set_GlobalID(mv,vertex_ov_global_id[MV_MasterParID(mv)*max_nbf+id_on_ov_list[i]]);
  }



  List_Delete(boundary_faces);
  MSTK_free(global_mesh_info);
  MSTK_free(vertex_ov_label);
  MSTK_free(vertex_ov_global_id);
  MSTK_free(id_on_ov_list);

  MSTK_free(list_boundary_vertex);
  MSTK_free(list_boundary_coor);
  MSTK_free(recv_list_vertex);
  MSTK_free(recv_list_coor);
  return 1;
}

  /* right now assume there are no overlapped regions */

int MESH_AssignGlobalIDs_Region(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  int i, nv, ne, nf, nr, global_id, mesh_info[10];
  MRegion_ptr mr;
  RepType rtype;
  int *global_mesh_info;


  for (i = 0; i < 10; i++) mesh_info[i] = 0;

  /* mesh_info store the mesh reptype, nv, ne, nf, nbf */
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

  global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  /* calculate starting global id number for faces*/
  global_id = 1;
  for(i = 0; i < rank; i++) 
    global_id = global_id + global_mesh_info[10*i+4];
  for(i = 0; i < nr; i++) {
    mr = MESH_Region(submesh,i);
    MR_Set_PType(mr,PINTERIOR);
    MR_Set_GlobalID(mr,global_id++);
    MR_Set_MasterParID(mr,rank);
  }

  MSTK_free(global_mesh_info);
  return 1;
}

#ifdef __cplusplus
}
#endif

