#include <stdio.h>
#include <stdlib.h>

#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* check is the parallel mesh is valid */
int MESH_Parallel_Check_Ghost(Mesh_ptr mesh, int rank);
int MESH_Parallel_Check_VertexGlobalID(Mesh_ptr mesh, int rank, int num, MSTK_Comm comm);
int MESH_Parallel_Check_EdgeGlobalID(Mesh_ptr mesh, int rank, int num, MSTK_Comm comm);
int MESH_Parallel_Check_FaceGlobalID(Mesh_ptr mesh, int rank, int num, MSTK_Comm comm);
int MESH_Parallel_Check_RegionGlobalID(Mesh_ptr mesh, int rank, int num, MSTK_Comm comm);
int MESH_Parallel_Check_GlobalID(Mesh_ptr mesh, int rank, int num, MSTK_Comm comm);


int MESH_Parallel_Check(Mesh_ptr mesh, MSTK_Comm comm) {
  int valid;

  int rank, num;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&num);

  if (num == 1) return 1;

#ifdef DEBUG
  fprintf(stderr,"Begin checking parallel information on submesh %d\n",rank);
#endif
  valid = MESH_Parallel_Check_Ghost(mesh,rank);
  if(valid)
    valid = MESH_Parallel_Check_GlobalID(mesh,rank,num,comm);
#ifdef DEBUG
  if (valid) fprintf(stderr,"Passed parallel checking on submesh %d\n",rank);
  else fprintf(stderr,"Failed parallel checking on submesh %d\n",rank);
#endif

  return valid;

}
 
 /* 
    Send ghost entities to their host processor and check if it exists and if 
    the information matches
 */

int MESH_Parallel_Check_GlobalID(Mesh_ptr mesh, int rank, int num, MSTK_Comm comm) {
  int valid = 1;
  valid = MESH_Parallel_Check_VertexGlobalID(mesh,rank,num,comm);
  valid &= MESH_Parallel_Check_EdgeGlobalID(mesh,rank,num,comm);
  valid &= MESH_Parallel_Check_FaceGlobalID(mesh,rank,num,comm);
  if(MESH_Num_Regions(mesh))
    valid &= MESH_Parallel_Check_RegionGlobalID(mesh,rank,num,comm);
  return valid;

}

int MESH_Parallel_Check_VertexGlobalID(Mesh_ptr mesh, int rank, int num, MSTK_Comm comm) {
  int i, j, idx, nv, nov, max_nv, index_mv, global_id, *loc, iloc, valid = 1;
  char mesg[256], funcname[32] = "MESH_Parallel_Check";
  MPI_Status status;
  MVertex_ptr mv;
  List_ptr ov_verts;
  double coor[3];
  int count, mesh_info[10], *global_mesh_info, *ov_list;
  int gid, gdim, ptype, master_id;
  for (i = 0; i < 10; i++) mesh_info[i] = 0;
  nv = MESH_Num_Vertices(mesh);

  mesh_info[0] = nv;

  /* First check if global IDs are continuous on each processor */
  int *vgidlist = (int *) malloc(nv*sizeof(int));
  int nown = 0;
  idx = 0; 
  while ((mv = MESH_Next_Vertex(mesh, &idx)))
    if (MV_PType(mv) != PGHOST)
      vgidlist[nown++] = MV_GlobalID(mv);
  qsort(vgidlist,nown,sizeof(int),compareINT);

  for (i = 1; i < nown; i++)
    if (vgidlist[i] != vgidlist[i-1]+1)
      MSTK_Report("MESH_Parallel_Check_VertexGlobalID",
                  "Vertex Global IDs are not contiguous",MSTK_WARN);
  if (vgidlist) free(vgidlist);

  /* Additional checks */
  /* collect overlap vertex list for fast checking */
  idx = 0; nov = 0; ov_verts = List_New(10);
  while ((mv = MESH_Next_Vertex(mesh, &idx)))  {
    if(MV_PType(mv) == POVERLAP)  {
      List_Add(ov_verts,mv);
      nov++;
    }
  }
  List_Sort(ov_verts,nov,sizeof(MVertex_ptr),compareGlobalID);

  /* first half stores the global id, second half are the local ids */
  ov_list = (int *) malloc(2*nov*sizeof(int));
  for(i = 0; i < nov; i++) {
    mv = List_Entry(ov_verts,i);
    ov_list[i] = MV_GlobalID(mv);
    ov_list[nov+i] = MV_ID(mv);
  }

  global_mesh_info = (int *)malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);
  max_nv = nv;
  for(i = 0; i < num; i++)
    if(max_nv < global_mesh_info[10*i])
      max_nv = global_mesh_info[10*i];
  
  /* collect data */
  int *list_vertex = (int *)malloc(3*max_nv*sizeof(int));
  double *list_coor = (double *)malloc(3*max_nv*sizeof(double));

  int *recv_list_vertex = (int *)malloc(3*max_nv*sizeof(int));
  double *recv_list_coor = (double *)malloc(3*max_nv*sizeof(double));

  /* similar as updateattr() and parallel_addghost(), order the send and recv process based on rank */
  for(i = 0; i < num; i++) {
    if(i == rank) continue;
    if(i < rank) {     
      if( MESH_Has_Overlaps_On_Prtn(mesh,i,MVERTEX) ) {
	MPI_Recv(recv_list_vertex,3*nv,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	//printf("rank %d receives %d verts from rank %d\n", rank,count/3, i);
	MPI_Recv(recv_list_coor,3*nv,MPI_DOUBLE,i,rank,comm,&status);
	for(j = 0; j < count/3;  j++) {
	  global_id = recv_list_vertex[3*j+2];
	  loc = (int *)bsearch(&global_id,
			       ov_list,
			       nov,
			       sizeof(int),
			       compareINT);
	  if(!loc) {              /* vertex not found */
	    valid = 0;
	    sprintf(mesg,"Global vertex %-d from processor %d is not on processor %d ", global_id, i, rank);
	    MSTK_Report(funcname,mesg,MSTK_ERROR);
	  }

	  if(loc) {                /* vertex found, but other information mismatch */
	    iloc = (int)(loc - ov_list);
	    mv = MESH_Vertex(mesh, ov_list[nov + iloc]-1);  /* get the vertex on current mesh */
	    gdim = (recv_list_vertex[3*j] & 7);
	    gid = (recv_list_vertex[3*j] >> 3);
	    ptype = (recv_list_vertex[3*j+1] & 3);
	    master_id = (recv_list_vertex[3*j+1] >> 2);
	    if(MV_GEntDim(mv) != gdim) {
	      valid = 0;
	      sprintf(mesg,"Global vertex %-d from processor %d and on processor %d GEntDim mismatch: %d vs %d ", global_id, i, rank, gdim, MV_GEntDim(mv));
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	    if(MV_GEntID(mv) != gid) {
	      valid = 0;
	      sprintf(mesg,"Global vertex %-d from processor %d and on processor %d GEntID mismatch: %d vs %d ", global_id, i, rank, gid, MV_GEntID(mv));
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }

	    if(ptype != PGHOST) {
	      valid = 0;
	      sprintf(mesg,"Global vertex %-d from processor %d is not a ghost", global_id, i);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	    
	    if(master_id != rank) {
	      valid = 0;
	      sprintf(mesg,"Global vertex %-d from processor %d is not on processor %d", global_id, i, rank);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }

	    MV_Coords(mv,coor);
	    if(compareCoorDouble(coor, &recv_list_coor[3*j]) != 0 ) {
	      valid = 0;
              // printf("The recv vertex %d coor: (%f,%f,%f)\n", global_id, recv_list_coor[3*j], recv_list_coor[3*j+1], recv_list_coor[3*j+2]);
	      sprintf(mesg,"Global vertex %-d from processor %d and on processor %d coordinate values mismatch", global_id, i, rank);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
            //	    if(!valid) 
            //	      printf("the mismatch vertex coor: (%f,%f,%f)\n", recv_list_coor[3*j], recv_list_coor[3*j+1], recv_list_coor[3*j+2]);
	  }
	}
      }
      if( MESH_Has_Ghosts_From_Prtn(mesh,i,MVERTEX) ) {
	idx = 0; index_mv = 0;
	while ((mv = MESH_Next_Vertex(mesh, &idx))) {
	  if(MV_MasterParID(mv) == i) {
	    list_vertex[3*index_mv] = (MV_GEntID(mv)<<3) | (MV_GEntDim(mv));
	    list_vertex[3*index_mv+1] = (MV_MasterParID(mv) <<2) | (MV_PType(mv));
	    list_vertex[3*index_mv+2] = MV_GlobalID(mv);
	    MV_Coords(mv,coor);
	    list_coor[index_mv*3] = coor[0];
	    list_coor[index_mv*3+1] = coor[1];
	    list_coor[index_mv*3+2] = coor[2];
	    index_mv++;
	  }
	}
	MPI_Send(list_vertex,3*index_mv,MPI_INT,i,i,comm);
	//	printf("rank %d sends %d verts to rank %d\n", rank,index_mv, i);
	MPI_Send(list_coor,3*index_mv,MPI_DOUBLE,i,i,comm);
	
      }
      
    }
    if(i > rank) {     
      if( MESH_Has_Ghosts_From_Prtn(mesh,i,MVERTEX) ) {
	idx = 0; index_mv = 0;
	while ((mv = MESH_Next_Vertex(mesh, &idx))) {
	  if(MV_MasterParID(mv) == i) {
	    list_vertex[3*index_mv] = (MV_GEntID(mv)<<3) | (MV_GEntDim(mv));
	    list_vertex[3*index_mv+1] = (MV_MasterParID(mv) <<2) | (MV_PType(mv));
	    list_vertex[3*index_mv+2] = MV_GlobalID(mv);
	    MV_Coords(mv,coor);
	    list_coor[index_mv*3] = coor[0];
	    list_coor[index_mv*3+1] = coor[1];
	    list_coor[index_mv*3+2] = coor[2];
	    index_mv++;
	  }
	}
	MPI_Send(list_vertex,3*index_mv,MPI_INT,i,i,comm);
	//	printf("rank %d sends %d verts to rank %d\n", rank,index_mv, i);
	MPI_Send(list_coor,3*index_mv,MPI_DOUBLE,i,i,comm);
	
      }
      if( MESH_Has_Overlaps_On_Prtn(mesh,i,MVERTEX) ) {
	MPI_Recv(recv_list_vertex,3*nv,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	//	printf("rank %d receives %d verts from rank %d\n", rank,count/3, i);
	MPI_Recv(recv_list_coor,3*nv,MPI_DOUBLE,i,rank,comm,&status);
	for(j = 0; j < count/3;  j++) {
	  global_id = recv_list_vertex[3*j+2];
	  loc = (int *)bsearch(&global_id,
			       ov_list,
			       nov,
			       sizeof(int),
			       compareINT);
	  if(!loc) {              /* vertex not found */
	    valid = 0;
	    sprintf(mesg,"Global vertex %-d from processor %d is not on processor %d" , global_id, i, rank);
	    printf("coor: (%f,%f,%f)\n", recv_list_coor[3*j], recv_list_coor[3*j+1], recv_list_coor[3*j+2]);
	    MSTK_Report(funcname,mesg,MSTK_ERROR);
	  }

	  if(loc) {                /* vertex found, but other information mismatch */
	    iloc = (int)(loc - ov_list);
	    mv = MESH_Vertex(mesh, ov_list[nov + iloc]-1);  /* get the vertex on current mesh */
	    gdim = (recv_list_vertex[3*j] & 7);
	    gid = (recv_list_vertex[3*j] >> 3);
	    ptype = (recv_list_vertex[3*j+1] & 3);
	    master_id = (recv_list_vertex[3*j+1] >> 2);
	    if(MV_GEntDim(mv) != gdim) {
	      valid = 0;
	      sprintf(mesg,"Global vertex %-d from processor %d and on processor %d GEntDim mismatch: %d vs %d ", global_id, i, rank, gdim, MV_GEntDim(mv));
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	    if(MV_GEntID(mv) != gid) {
	      valid = 0;
	      sprintf(mesg,"Global vertex %-d from processor %d and on processor %d GEntID mismatch: %d vs %d ", global_id, i, rank, gid, MV_GEntID(mv));
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }

	    if(ptype != PGHOST) {
	      valid = 0;
	      sprintf(mesg,"Global vertex %-d from processor %d is not a ghost", global_id, i);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	    
	    if(master_id != rank) {
	      valid = 0;
	      sprintf(mesg,"Global vertex %-d from processor %d is not on processor %d", global_id, i, rank);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }

	    MV_Coords(mv,coor);
	    if(compareCoorDouble(coor, &recv_list_coor[3*j]) != 0 ) {
	      valid = 0;
              //	      printf("the recv vertex %d coor: (%f,%f,%f)\n", global_id, recv_list_coor[3*j], recv_list_coor[3*j+1], recv_list_coor[3*j+2]);
	      sprintf(mesg,"Global vertex %-d from processor %d and on processor %d coordinate value mismatch", global_id, i, rank);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
            //	    if(!valid) 
            //	      printf("the mismatch vertex %d coor: (%f,%f,%f)\n", MV_GlobalID(mv), coor[0], coor[1], coor[2]);
	  }
	}
      }
    }
  }
  free(list_vertex);
  free(list_coor);
  free(recv_list_vertex);
  free(recv_list_coor);
  free(ov_list);
  free(global_mesh_info);
  List_Delete(ov_verts);
  return valid;
}


int MESH_Parallel_Check_EdgeGlobalID(Mesh_ptr mesh, int rank, int num, MSTK_Comm comm) {
  int i, j, idx, ne, noe, max_ne, index_me, global_id, *loc, iloc, valid = 1;
  char mesg[256], funcname[32] = "MESH_Parallel_Check";
  MPI_Status status;
  MEdge_ptr me;
  List_ptr ov_edges;
  int count, mesh_info[10], *global_mesh_info, *ov_list;
  int edge_int[2], gid, gdim, ptype, master_id;
  for (i = 0; i < 10; i++) mesh_info[i] = 0;
  ne = MESH_Num_Edges(mesh);

  mesh_info[0] = ne;

  /* First check if global IDs are continuous on each processor */
  int *egidlist = (int *) malloc(ne*sizeof(int));
  int nown = 0;
  idx = 0; 
  while ((me = MESH_Next_Edge(mesh, &idx)))
    if (ME_PType(me) != PGHOST)
      egidlist[nown++] = ME_GlobalID(me);
  qsort(egidlist,nown,sizeof(int),compareINT);

  for (i = 1; i < nown; i++)
    if (egidlist[i] != egidlist[i-1]+1)
      MSTK_Report("MESH_Parallel_Check_EdgeGlobalID",
                  "Edge Global IDs are not contiguous",MSTK_WARN);
  if (egidlist) free(egidlist);

  /* collect overlap edge list for fast checking */
  idx = 0; noe = 0; ov_edges = List_New(10);
  while ((me = MESH_Next_Edge(mesh, &idx)))  {
    if(ME_PType(me) == POVERLAP)  {
      List_Add(ov_edges,me);
      noe++;
    }
  }
  List_Sort(ov_edges,noe,sizeof(MEdge_ptr),compareGlobalID);

  /* first half stores the global id, second half are the local ids */
  ov_list = (int *) malloc(2*noe*sizeof(int));
  for(i = 0; i < noe; i++) {
    me = List_Entry(ov_edges,i);
    ov_list[i] = ME_GlobalID(me);
    ov_list[noe+i] = ME_ID(me);
  }

  global_mesh_info = (int *)malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  max_ne = ne;
  for(i = 0; i < num; i++)
    if(max_ne < global_mesh_info[10*i])
      max_ne = global_mesh_info[10*i];
  
  /* collect data */
  int *list_edge = (int *)malloc(5*max_ne*sizeof(int));
  int *recv_list_edge = (int *)malloc(5*max_ne*sizeof(int));

  for(i = 0; i < num; i++) {
    if(i == rank) continue;
    if(i < rank) {     
      if( MESH_Has_Overlaps_On_Prtn(mesh,i,MEDGE) ) {
	MPI_Recv(recv_list_edge,5*ne,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	//printf("rank %d receives %d edges from rank %d\n", rank,count/5, i);
	for(j = 0; j < count/5;  j++) {
	  global_id = recv_list_edge[5*j+4];
	  loc = (int *)bsearch(&global_id,
			       ov_list,
			       noe,
			       sizeof(int),
			       compareINT);
	  if(!loc) {
            //	    printf("the recv edge %d endpoints: (%d,%d)\n", global_id, recv_list_edge[5*j], recv_list_edge[5*j+1]);
	    sprintf(mesg,"Global edge %-d from processor %d not found on processor %d ", global_id, i, rank);
	    MSTK_Report(funcname,mesg,MSTK_ERROR);
	    valid = 0;
	  }
	  if(loc) {                /* vertex found, but other information mismatch */
	    iloc = (int)(loc - ov_list);
	    me = MESH_Edge(mesh, ov_list[noe + iloc]-1);  /* get the edge on current mesh */
	    gdim = (recv_list_edge[5*j+2] & 7);
	    gid = (recv_list_edge[5*j+2] >> 3);
	    ptype = (recv_list_edge[5*j+3] & 3);
	    master_id = (recv_list_edge[5*j+3] >> 2);
	    if(MV_GEntDim(me) != gdim) {
	      valid = 0;
	      sprintf(mesg,"Global edge %-d from processor %d and on processor %d GEntDim mismatch: %d vs %d ", global_id, i, rank, gdim, ME_GEntDim(me));
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	    if(ME_GEntID(me) != gid) {
	      valid = 0;
	      sprintf(mesg,"Global edge %-d from processor %d and on processor %d GEntID mismatch: %d vs %d ", global_id, i, rank, gid, ME_GEntID(me));
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }

	    if(ptype != PGHOST) {
	      valid = 0;
	      sprintf(mesg,"Global edge %-d from processor %d is not a ghost", global_id, i);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	    
	    if(master_id != rank) {
	      valid = 0;
	      sprintf(mesg,"Global edge %-d from processor %d is not on processor %d", global_id, i, rank);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }

	    edge_int[0] = recv_list_edge[5*j];
	    edge_int[1] = recv_list_edge[5*j+1];
	    if(compareEdgeINT(edge_int, &recv_list_edge[5*j]) != 0 ) {
	      valid = 0;
              //	      printf("the recv edge %d endpoints: (%d,%d)\n", global_id, recv_list_edge[5*j], recv_list_edge[5*j+1]);
	      sprintf(mesg,"Global edge %-d from processor %d and on processor %d endpoints mismatch", global_id, i, rank);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
            //	    if(!valid) 
            //	      printf("the mismatch edge %d endpoints: (%d,%d)\n", ME_GlobalID(me), edge_int[0], edge_int[1]);
	  }
	}
      }	
      if( MESH_Has_Ghosts_From_Prtn(mesh,i,MEDGE) ) {
	idx = 0; index_me = 0;
	while ((me = MESH_Next_Edge(mesh, &idx))) {
	  if(ME_MasterParID(me) == i) {
	    list_edge[5*index_me]   = MV_GlobalID(ME_Vertex(me,0));
	    list_edge[5*index_me+1] = MV_GlobalID(ME_Vertex(me,1));
	    list_edge[5*index_me+2] = (ME_GEntID(me)<<3) | (ME_GEntDim(me));
	    list_edge[5*index_me+3] = (ME_MasterParID(me) <<2) | (ME_PType(me));
	    list_edge[5*index_me+4] = ME_GlobalID(me);
	    index_me++;
	  }
	}
	MPI_Send(list_edge,5*index_me,MPI_INT,i,i,comm);
	//printf("rank %d sends %d edges to rank %d\n", rank,index_me, i);
      }
    }
    if(i > rank) {     
      if( MESH_Has_Ghosts_From_Prtn(mesh,i,MEDGE) ) {
	idx = 0; index_me = 0;
	while ((me = MESH_Next_Edge(mesh, &idx))) {
	  if(ME_MasterParID(me) == i) {
	    list_edge[5*index_me]   = MV_GlobalID(ME_Vertex(me,0));
	    list_edge[5*index_me+1] = MV_GlobalID(ME_Vertex(me,1));
	    list_edge[5*index_me+2] = (ME_GEntID(me)<<3) | (ME_GEntDim(me));
	    list_edge[5*index_me+3] = (ME_MasterParID(me) <<2) | (ME_PType(me));
	    list_edge[5*index_me+4] = ME_GlobalID(me);
	    index_me++;
	  }
	}
	MPI_Send(list_edge,5*index_me,MPI_INT,i,i,comm);
	//	printf("rank %d sends %d edges to rank %d\n", rank,index_me, i);
      }
      if( MESH_Has_Overlaps_On_Prtn(mesh,i,MEDGE) ) {
	MPI_Recv(recv_list_edge,5*ne,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	//printf("rank %d receives %d edges from rank %d\n", rank,count/5, i);
	for(j = 0; j < count/5;  j++) {
	  global_id = recv_list_edge[5*j+4];
	  loc = (int *)bsearch(&global_id,
			       ov_list,
			       noe,
			       sizeof(int),
			       compareINT);
	  if(!loc) {
            //	    printf("the recv edge %d endpoints: (%d,%d)\n", global_id, recv_list_edge[5*j], recv_list_edge[5*j+1]);
	    sprintf(mesg,"Global edge %-d from processor %d not found on processor %d ", global_id, i, rank);
	    MSTK_Report(funcname,mesg,MSTK_ERROR);
	    valid = 0;
	  }
	  if(loc) {                /* vertex found, but other information mismatch */
	    iloc = (int)(loc - ov_list);
	    me = MESH_Edge(mesh, ov_list[noe + iloc]-1);  /* get the edge on current mesh */
	    gdim = (recv_list_edge[5*j+2] & 7);
	    gid = (recv_list_edge[5*j+2] >> 3);
	    ptype = (recv_list_edge[5*j+3] & 3);
	    master_id = (recv_list_edge[5*j+3] >> 2);
	    if(ME_GEntDim(me) != gdim) {
	      valid = 0;
	      sprintf(mesg,"Global edge %-d from processor %d and on processor %d GEntDim mismatch: %d vs %d ", global_id, i, rank, gdim, ME_GEntDim(me));
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	    if(ME_GEntID(me) != gid) {
	      valid = 0;
	      sprintf(mesg,"Global edge %-d from processor %d and on processor %d GEntID mismatch: %d vs %d ", global_id, i, rank, gid, ME_GEntID(me));
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }

	    if(ptype != PGHOST) {
	      valid = 0;
	      sprintf(mesg,"Global edge %-d from processor %d is not a ghost", global_id, i);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	    
	    if(master_id != rank) {
	      valid = 0;
	      sprintf(mesg,"Global edge %-d from processor %d is not on processor %d", global_id, i, rank);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }

	    edge_int[0] = recv_list_edge[5*j];
	    edge_int[1] = recv_list_edge[5*j+1];
	    if(compareEdgeINT(edge_int, &recv_list_edge[5*j]) != 0 ) {
	      valid = 0;
              //	      printf("the recv edge %d endpoints: (%d,%d)\n", global_id, recv_list_edge[5*j], recv_list_edge[5*j+1]);
	      sprintf(mesg,"Global edge %-d from processor %d and on processor %d endpoints mismatch", global_id, i, rank);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
            //	    if(!valid) 
            //	      printf("the mismatch edge %d endpoints: (%d,%d)\n", ME_GlobalID(me), edge_int[0], edge_int[1]);
	  }
	}	
      }
    }
  }
  free(list_edge);
  free(recv_list_edge);
  free(ov_list);
  free(global_mesh_info);
  List_Delete(ov_edges);
  return valid;
}


int MESH_Parallel_Check_FaceGlobalID(Mesh_ptr mesh, int rank, int num, MSTK_Comm comm) {
  int i, j, idx, nf, nof, nfe, max_nf, index_mf, global_id, *loc, iloc, valid = 1;
  char mesg[256], funcname[32] = "MESH_Parallel_Check";
  MPI_Status status;
  MEdge_ptr me;
  MFace_ptr mf;
  List_ptr ov_faces, mfedges;
  int count, mesh_info[10], *global_mesh_info, *ov_list;
  int gid, gdim, ptype, master_id, dir;
  for (i = 0; i < 10; i++) mesh_info[i] = 0;
  nf = MESH_Num_Faces(mesh);

  mesh_info[0] = nf;

  /* First check if global IDs are continuous on each processor */
  int *fgidlist = (int *) malloc(nf*sizeof(int));
  int nown = 0;
  idx = 0; 
  while ((mf = MESH_Next_Face(mesh, &idx)))
    if (MF_PType(mf) != PGHOST)
      fgidlist[nown++] = MF_GlobalID(mf);
  qsort(fgidlist,nown,sizeof(int),compareINT);

  for (i = 1; i < nown; i++)
    if (fgidlist[i] != fgidlist[i-1]+1)
      MSTK_Report("MESH_Parallel_Check_FaceGlobalID",
                  "Face Global IDs are not contiguous",MSTK_WARN);
  if (fgidlist) free(fgidlist);

  /* collect overlap face list for fast checking */
  idx = 0; nof = 0; ov_faces = List_New(10);
  while ((mf = MESH_Next_Face(mesh, &idx)))  {
    if(MF_PType(mf) == POVERLAP)  {
      List_Add(ov_faces,mf);
      nof++;
    }
  }
  List_Sort(ov_faces,nof,sizeof(MFace_ptr),compareGlobalID);

  /* first half stores the global id, second half are the local ids */
  ov_list = (int *) malloc(2*nof*sizeof(int));
  for(i = 0; i < nof; i++) {
    mf = List_Entry(ov_faces,i);
    ov_list[i] = MF_GlobalID(mf);
    ov_list[nof+i] = MF_ID(mf);
  }

  global_mesh_info = (int *)malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  max_nf = nf;
  for(i = 0; i < num; i++)
    if(max_nf < global_mesh_info[10*i])
      max_nf = global_mesh_info[10*i];
  
  /* collect data */
  int *list_face = (int *)malloc((MAXPV3+4)*max_nf*sizeof(int));
  int *recv_list_face = (int *)malloc((MAXPV3+4)*max_nf*sizeof(int));
  for(i = 0; i < num; i++) {
    if(i == rank) continue;
    if(i < rank) {     
      if( MESH_Has_Overlaps_On_Prtn(mesh,i,MFACE) ) {
	MPI_Recv(recv_list_face,(MAXPV3+4)*nf,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	//printf("rank %d receives %d faces from rank %d\n", rank,count, i);
	for(index_mf = 0; index_mf < count;) {
	  nfe = recv_list_face[index_mf];
	  global_id = recv_list_face[index_mf+nfe+3];;
	  loc = (int *)bsearch(&global_id,
			       ov_list,
			       nof,
			       sizeof(int),
			       compareINT);
	  if(!loc) {
	    sprintf(mesg,"Global face %-d from processor %d not found on processor %d ", global_id, i, rank);
	    MSTK_Report(funcname,mesg,MSTK_ERROR);
	    valid = 0;
	  }
	  if(loc) {                /* vertex found, but other information mismatch */
	    iloc = (int)(loc - ov_list);
	    mf = MESH_Face(mesh, ov_list[nof + iloc]-1);  /* get the face on current mesh */
	    gdim = (recv_list_face[index_mf+nfe+1] & 7);
	    gid = (recv_list_face[index_mf+nfe+1] >> 3);
	    ptype = (recv_list_face[index_mf+nfe+2] & 3);
	    master_id = (recv_list_face[index_mf+nfe+2] >> 2);
	    if(MF_GEntDim(mf) != gdim) {
	      valid = 0;
	      sprintf(mesg,"Global face %-d from processor %d and on processor %d GEntDim mismatch: %d vs %d ", global_id, i, rank, gdim, MF_GEntDim(mf));
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	    if(MF_GEntID(mf) != gid) {
	      valid = 0;
	      sprintf(mesg,"Global face %-d from processor %d and on processor %d GEntID mismatch: %d vs %d ", global_id, i, rank, gid, MF_GEntID(mf));
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }

	    if(ptype != PGHOST) {
	      valid = 0;
	      sprintf(mesg,"Global face %-d from processor %d is not a ghost", global_id, i);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	    
	    if(master_id != rank) {
	      valid = 0;
	      sprintf(mesg,"Global face %-d from processor %d is not on processor %d", global_id, i, rank);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	  }
	  index_mf += (nfe + 4);
	}
      }
      if( MESH_Has_Ghosts_From_Prtn(mesh,i,MFACE) ) {
  	idx = 0; index_mf = 0;
	while ((mf = MESH_Next_Face(mesh, &idx))) {
	  if(MF_MasterParID(mf) == i) {
	    mfedges = MF_Edges(mf,1,0);
	    nfe = List_Num_Entries(mfedges);
	    list_face[index_mf] = nfe;
	    for(j = 0; j < nfe; j++) {
	      dir = MF_EdgeDir_i(mf,j) == 1 ? 1 : -1;
	      list_face[index_mf+j+1] = dir*ME_GlobalID(List_Entry(mfedges,j));
	    }
	    list_face[index_mf+nfe+1] = (MF_GEntID(mf)<<3) | (MF_GEntDim(mf));
	    list_face[index_mf+nfe+2] = (MF_MasterParID(mf)<<2) | (MF_PType(mf));
	    list_face[index_mf+nfe+3] = MF_GlobalID(mf);
	    index_mf += (nfe + 4);
	    List_Delete(mfedges);
	  }
	}
	MPI_Send(list_face,index_mf,MPI_INT,i,i,comm);
	//printf("rank %d sends %d faces to rank %d\n", rank,index_mf, i);
      }
    }
    if(i > rank) {     
      if( MESH_Has_Ghosts_From_Prtn(mesh,i,MFACE) ) {
  	idx = 0; index_mf = 0;
	while ((mf = MESH_Next_Face(mesh, &idx))) {
	  if(MF_MasterParID(mf) == i) {
	    mfedges = MF_Edges(mf,1,0);
	    nfe = List_Num_Entries(mfedges);
	    list_face[index_mf] = nfe;
	    for(j = 0; j < nfe; j++) {
	      dir = MF_EdgeDir_i(mf,j) == 1 ? 1 : -1;
	      list_face[index_mf+j+1] = dir*ME_GlobalID(List_Entry(mfedges,j));
	    }
	    list_face[index_mf+nfe+1] = (MF_GEntID(mf)<<3) | (MF_GEntDim(mf));
	    list_face[index_mf+nfe+2] = (MF_MasterParID(mf)<<2) | (MF_PType(mf));
	    list_face[index_mf+nfe+3] = MF_GlobalID(mf);
	    index_mf += (nfe + 4);
	    List_Delete(mfedges);
	  }
	}
	MPI_Send(list_face,index_mf,MPI_INT,i,i,comm);
	//printf("rank %d sends %d faces to rank %d\n", rank,index_mf, i);
      }

      if( MESH_Has_Overlaps_On_Prtn(mesh,i,MFACE) ) {
	MPI_Recv(recv_list_face,(MAXPV3+4)*nf,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	//printf("rank %d receives %d faces from rank %d\n", rank,count, i);
	for(index_mf = 0; index_mf < count;) {
	  nfe = recv_list_face[index_mf];
	  global_id = recv_list_face[index_mf+nfe+3];;
	  loc = (int *)bsearch(&global_id,
			       ov_list,
			       nof,
			       sizeof(int),
			       compareINT);
	  if(!loc) {
	    sprintf(mesg,"Global face %-d from processor %d not found on processor %d ", global_id, i, rank);
	    MSTK_Report(funcname,mesg,MSTK_ERROR);
	    valid = 0;
	  }
	  if(loc) {                /* vertex found, but other information mismatch */
	    iloc = (int)(loc - ov_list);
	    mf = MESH_Face(mesh, ov_list[nof + iloc]-1);  /* get the face on current mesh */
	    gdim = (recv_list_face[index_mf+nfe+1] & 7);
	    gid = (recv_list_face[index_mf+nfe+1] >> 3);
	    ptype = (recv_list_face[index_mf+nfe+2] & 3);
	    master_id = (recv_list_face[index_mf+nfe+2] >> 2);
	    if(MF_GEntDim(mf) != gdim) {
	      valid = 0;
	      sprintf(mesg,"Global face %-d from processor %d and on processor %d GEntDim mismatch: %d vs %d ", global_id, i, rank, gdim, MF_GEntDim(mf));
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	    if(MF_GEntID(mf) != gid) {
	      valid = 0;
	      sprintf(mesg,"Global face %-d from processor %d and on processor %d GEntID mismatch: %d vs %d ", global_id, i, rank, gid, MF_GEntID(mf));
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }

	    if(ptype != PGHOST) {
	      valid = 0;
	      sprintf(mesg,"Global face %-d from processor %d is not a ghost", global_id, i);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	    
	    if(master_id != rank) {
	      valid = 0;
	      sprintf(mesg,"Global face %-d from processor %d is not on processor %d", global_id, i, rank);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	  }
	  index_mf += (nfe + 4);
	}
      }
    }
  }
  free(list_face);
  free(recv_list_face);
  free(ov_list);
  free(global_mesh_info);
  List_Delete(ov_faces);
  return valid;
}

int MESH_Parallel_Check_RegionGlobalID(Mesh_ptr mesh, int rank, int num, MSTK_Comm comm) {
  int i, j, idx, nr, nor, nrf, max_nr, index_mr, global_id, *loc, iloc, valid = 1;
  char mesg[256], funcname[32] = "MESH_Parallel_Check";
  MPI_Status status;
  MFace_ptr mf;
  MRegion_ptr mr;
  List_ptr ov_regions, mrfaces;
  int count, mesh_info[10], *global_mesh_info, *ov_list;
  int gid, gdim, ptype, master_id, dir;
  for (i = 0; i < 10; i++) mesh_info[i] = 0;
  nr = MESH_Num_Regions(mesh);

  mesh_info[0] = nr;

  /* First check if global IDs are continuous on each processor */
  int *rgidlist = (int *) malloc(nr*sizeof(int));
  int nown = 0;
  idx = 0;
  while ((mr = MESH_Next_Region(mesh, &idx)))
    if (MR_PType(mr) != PGHOST)
      rgidlist[nown++] = MR_GlobalID(mr);
  qsort(rgidlist,nown,sizeof(int),compareINT);

  for (i = 1; i < nown; i++)
    if (rgidlist[i] != rgidlist[i-1]+1)
      MSTK_Report("MESH_Parallel_Check_RegionGlobalID",
                  "Region Global IDs are not contiguous",MSTK_WARN);
  if (rgidlist) free(rgidlist);

  /* collect overlap region list for fast checking */
  idx = 0; nor = 0; ov_regions = List_New(10);
  while ((mr = MESH_Next_Region(mesh, &idx)))  {
    if(MR_PType(mr) == POVERLAP)  {
      List_Add(ov_regions,mr);
      nor++;
    }
  }
  List_Sort(ov_regions,nor,sizeof(MRegion_ptr),compareGlobalID);

  /* first half stores the global id, second half are the local ids */
  ov_list = (int *) malloc(2*nor*sizeof(int));
  for(i = 0; i < nor; i++) {
    mr = List_Entry(ov_regions,i);
    ov_list[i] = MR_GlobalID(mr);
    ov_list[nor+i] = MR_ID(mr);
  }

  global_mesh_info = (int *)malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  max_nr = nr;
  for(i = 0; i < num; i++)
    if(max_nr < global_mesh_info[10*i])
      max_nr = global_mesh_info[10*i];
  
  /* collect data */
  int *list_region = (int *)malloc((MAXPF3+4)*max_nr*sizeof(int));
  int *recv_list_region = (int *)malloc((MAXPF3+4)*max_nr*sizeof(int));
  for(i = 0; i < num; i++) {
    if(i == rank) continue;
    if(i < rank) {     
      if( MESH_Has_Overlaps_On_Prtn(mesh,i,MREGION) ) {
	MPI_Recv(recv_list_region,(MAXPF3+4)*nr,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	//	printf("rank %d receives %d regions from rank %d\n", rank,count, i);
	for(index_mr = 0; index_mr < count;) {
	  nrf = recv_list_region[index_mr];
	  global_id = recv_list_region[index_mr+nrf+3];;
	  loc = (int *)bsearch(&global_id,
			       ov_list,
			       nor,
			       sizeof(int),
			       compareINT);
	  if(!loc) {
	    sprintf(mesg,"Global region %-d from processor %d not found on processor %d ", global_id, i, rank);
	    MSTK_Report(funcname,mesg,MSTK_ERROR);
	    valid = 0;
	  }
	  if(loc) {                /* vertex found, but other information mismatch */
	    iloc = (int)(loc - ov_list);
	    mr = MESH_Region(mesh, ov_list[nor + iloc]-1);  /* get the region on current mesh */
	    gdim = (recv_list_region[index_mr+nrf+1] & 7);
	    gid = (recv_list_region[index_mr+nrf+1] >> 3);
	    ptype = (recv_list_region[index_mr+nrf+2] & 3);
	    master_id = (recv_list_region[index_mr+nrf+2] >> 2);
	    if(MR_GEntDim(mr) != gdim) {
	      valid = 0;
	      sprintf(mesg,"Global region %-d from processor %d and on processor %d GEntDim mismatch: %d vs %d ", global_id, i, rank, gdim, MR_GEntDim(mr));
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	    if(MR_GEntID(mr) != gid) {
	      valid = 0;
	      sprintf(mesg,"Global region %-d from processor %d and on processor %d GEntID mismatch: %d vs %d ", global_id, i, rank, gid, MR_GEntID(mr));
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }

	    if(ptype != PGHOST) {
	      valid = 0;
	      sprintf(mesg,"Global region %-d from processor %d is not a ghost", global_id, i);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	    
	    if(master_id != rank) {
	      valid = 0;
	      sprintf(mesg,"Global region %-d from processor %d is not on processor %d", global_id, i, rank);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	  }
	  index_mr += (nrf + 4);
	}
      }
      if( MESH_Has_Ghosts_From_Prtn(mesh,i,MREGION) ) {
  	idx = 0; index_mr = 0;
	while ((mr = MESH_Next_Region(mesh, &idx))) {
	  if(MR_MasterParID(mr) == i) {
	    mrfaces = MR_Faces(mr);
	    nrf = List_Num_Entries(mrfaces);
	    list_region[index_mr] = nrf;
	    for(j = 0; j < nrf; j++) {
	      dir = MR_FaceDir_i(mr,j) == 1 ? 1 : -1;
	      list_region[index_mr+j+1] = dir*MF_GlobalID(List_Entry(mrfaces,j));
	    }
	    list_region[index_mr+nrf+1] = (MR_GEntID(mr)<<3) | (MR_GEntDim(mr));
	    list_region[index_mr+nrf+2] = (MR_MasterParID(mr)<<2) | (MR_PType(mr));
	    list_region[index_mr+nrf+3] = MR_GlobalID(mr);
	    index_mr += (nrf + 4);
	    List_Delete(mrfaces);
	  }
	}
	MPI_Send(list_region,index_mr,MPI_INT,i,i,comm);
	//printf("rank %d sends %d regions to rank %d\n", rank,index_mr, i);
      }
    }
    if(i > rank) {     
      if( MESH_Has_Ghosts_From_Prtn(mesh,i,MREGION) ) {
  	idx = 0; index_mr = 0;
	while ((mr = MESH_Next_Region(mesh, &idx))) {
	  if(MR_MasterParID(mr) == i) {
	    mrfaces = MR_Faces(mr);
	    nrf = List_Num_Entries(mrfaces);
	    list_region[index_mr] = nrf;
	    for(j = 0; j < nrf; j++) {
	      dir = MR_FaceDir_i(mr,j) == 1 ? 1 : -1;
	      list_region[index_mr+j+1] = dir*MF_GlobalID(List_Entry(mrfaces,j));
	    }
	    list_region[index_mr+nrf+1] = (MR_GEntID(mr)<<3) | (MR_GEntDim(mr));
	    list_region[index_mr+nrf+2] = (MR_MasterParID(mr)<<2) | (MR_PType(mr));
	    list_region[index_mr+nrf+3] = MR_GlobalID(mr);
	    index_mr += (nrf + 4);
	    List_Delete(mrfaces);
	  }
	}
	MPI_Send(list_region,index_mr,MPI_INT,i,i,comm);
	//printf("rank %d sends %d regions to rank %d\n", rank,index_mr, i);
      }
      if( MESH_Has_Overlaps_On_Prtn(mesh,i,MREGION) ) {
	MPI_Recv(recv_list_region,(MAXPF3+4)*nr,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	//printf("rank %d receives %d regions from rank %d\n", rank,count, i);
	for(index_mr = 0; index_mr < count;) {
	  nrf = recv_list_region[index_mr];
	  global_id = recv_list_region[index_mr+nrf+3];;
	  loc = (int *)bsearch(&global_id,
			       ov_list,
			       nor,
			       sizeof(int),
			       compareINT);
	  if(!loc) {
	    sprintf(mesg,"Global region %-d from processor %d not found on processor %d ", global_id, i, rank);
	    MSTK_Report(funcname,mesg,MSTK_ERROR);
	    valid = 0;
	  }
	  if(loc) {                /* vertex found, but other information mismatch */
	    iloc = (int)(loc - ov_list);
	    mr = MESH_Region(mesh, ov_list[nor + iloc]-1);  /* get the region on current mesh */
	    gdim = (recv_list_region[index_mr+nrf+1] & 7);
	    gid = (recv_list_region[index_mr+nrf+1] >> 3);
	    ptype = (recv_list_region[index_mr+nrf+2] & 3);
	    master_id = (recv_list_region[index_mr+nrf+2] >> 2);
	    if(MR_GEntDim(mr) != gdim) {
	      valid = 0;
	      sprintf(mesg,"Global region %-d from processor %d and on processor %d GEntDim mismatch: %d vs %d ", global_id, i, rank, gdim, MR_GEntDim(mr));
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	    if(MR_GEntID(mr) != gid) {
	      valid = 0;
	      sprintf(mesg,"Global region %-d from processor %d and on processor %d GEntID mismatch: %d vs %d ", global_id, i, rank, gid, MR_GEntID(mr));
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }

	    if(ptype != PGHOST) {
	      valid = 0;
	      sprintf(mesg,"Global region %-d from processor %d is not a ghost", global_id, i);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	    
	    if(master_id != rank) {
	      valid = 0;
	      sprintf(mesg,"Global region %-d from processor %d is not on processor %d", global_id, i, rank);
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	  }
	  index_mr += (nrf + 4);
	}
      }
    }
  }
  free(list_region);
  free(recv_list_region);
  free(ov_list);
  free(global_mesh_info);
  List_Delete(ov_regions);
  return valid;
}


  /* 
   check if every ghost entity has a master partition number that is different from current partition 
   check if other PType entity has the same master partition number as this parition
*/
int MESH_Parallel_Check_Ghost(Mesh_ptr mesh, int rank) {
  int nv, ne, nf, nr, idx, valid = 1;
  char mesg[256], funcname[32] = "MESH_Parallel_Check";
  MVertex_ptr mv;
  MEdge_ptr me;
  MEdge_ptr mf;
  MEdge_ptr mr;
  nv = MESH_Num_Vertices(mesh);
  ne = MESH_Num_Edges(mesh);
  nf = MESH_Num_Faces(mesh);
  nr = MESH_Num_Regions(mesh);
  idx = 0;
  /* check vertex */
  while ((mv = MESH_Next_Vertex(mesh, &idx))) {
    if(MV_PType(mv) == PGHOST) {
      if (MV_MasterParID(mv) == rank) {
	sprintf(mesg,"Ghost Vertex %-d has master partition id %d but it is on partition %d", MV_GlobalID(mv),MV_MasterParID(mv),rank);
	MSTK_Report(funcname,mesg,MSTK_ERROR);
	valid = 0;
      }
    }
    else  {
      if(MV_MasterParID(mv) != rank) {
	sprintf(mesg,"Non-Ghost Vertex %-d has master partition id %d but it is on partition %d", MV_GlobalID(mv), MV_MasterParID(mv),rank);
	MSTK_Report(funcname,mesg,MSTK_ERROR);
	valid = 0;
      }	
    }
  }

  /* check edge */
  idx = 0;
  while ((me = MESH_Next_Edge(mesh, &idx))) {
    if(ME_PType(me) == PGHOST) {
      if (ME_MasterParID(me) == rank) {
	sprintf(mesg,"Ghost Edge %-d has master partition id %d but it is on partition %d", ME_GlobalID(me), ME_MasterParID(me),rank);
        //	printf("the mismatch edge %d endpoints: (%d,%d)\n", ME_GlobalID(me), MV_GlobalID(ME_Vertex(me,0)), MV_GlobalID(ME_Vertex(me,1)));
	MSTK_Report(funcname,mesg,MSTK_ERROR);
	valid = 0;
      }
    }
    else  {
      if(ME_MasterParID(me) != rank) {
	sprintf(mesg,"Non-Ghost Edge %-d has master partition id %d but it is on partition %d", ME_GlobalID(me), ME_MasterParID(me),rank);
	MSTK_Report(funcname,mesg,MSTK_ERROR);
        //	printf("the mismatch edge %d endpoints: (%d,%d)\n", ME_GlobalID(me), MV_GlobalID(ME_Vertex(me,0)), MV_GlobalID(ME_Vertex(me,1)));
	valid = 0;
      }	
    }
  }

  /* check face */
  idx = 0;
  while ((mf = MESH_Next_Face(mesh, &idx))) {
    if(MF_PType(mf) == PGHOST) {
      if (MF_MasterParID(mf) == rank) {
	sprintf(mesg,"Ghost Face %-d has master partition id %d but it is on partition %d", MF_GlobalID(mf), MF_MasterParID(mf),rank);
	MSTK_Report(funcname,mesg,MSTK_ERROR);
	valid = 0;
      }
    }
    else  {
      if(MF_MasterParID(mf) != rank) {
	sprintf(mesg,"Non-Ghost Face %-d has master partition id %d but it is on partition %d", MF_GlobalID(mf), MF_MasterParID(mf),rank);
	MSTK_Report(funcname,mesg,MSTK_ERROR);
	valid = 0;
      }	
    }
  }
  if(!nr) return valid;
  /* check region */
  idx = 0;
  while ((mr = MESH_Next_Region(mesh, &idx))) {
    if(MR_PType(mr) == PGHOST) {
      if (MR_MasterParID(mr) == rank) {
	sprintf(mesg,"Ghost Region %-d has master partition id %d but it is on partition %d", MR_GlobalID(mr), MR_MasterParID(mr),rank);
	MSTK_Report(funcname,mesg,MSTK_ERROR);
	valid = 0;
      }
    }
    else  {
      if(MR_MasterParID(mr) != rank) {
	sprintf(mesg,"Non-Ghost Region %-d has master partition id %d but it is on partition %d", MR_GlobalID(mr), MR_MasterParID(mr),rank);
	MSTK_Report(funcname,mesg,MSTK_ERROR);
	valid = 0;
      }	
    }
  }
  return valid;
}
#ifdef __cplusplus
}
#endif
