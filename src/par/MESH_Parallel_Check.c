#include <stdio.h>
#include <stdlib.h>

#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Routine to check that the topological structure of the mesh is valid */
  /* Rao Garimella - 12/16/2010                                           */
int MESH_Parallel_Check_Ghost(Mesh_ptr mesh, int rank);
int MESH_Parallel_Check_VertexGlobalID(Mesh_ptr mesh, int rank, int num, MPI_Comm comm);
int MESH_Parallel_Check_EdgeGlobalID(Mesh_ptr mesh, int rank, int num, MPI_Comm comm);
int MESH_Parallel_Check_GlobalID(Mesh_ptr mesh, int rank, int num, MPI_Comm comm);
int MESH_Parallel_Check(Mesh_ptr mesh, int rank, int num, MPI_Comm comm) {
  int valid = 1;
  //valid = MESH_CheckTopo(mesh);
  valid = MESH_Parallel_Check_Ghost(mesh,rank);
  valid = MESH_Parallel_Check_GlobalID(mesh,rank,num,comm);
  return valid;

} /* int MESH_Parallel_Check */

int MESH_Parallel_Check_GlobalID(Mesh_ptr mesh, int rank, int num, MPI_Comm comm) {
  int valid = 1;
  valid = MESH_Parallel_Check_VertexGlobalID(mesh,rank,num,comm);
  valid = MESH_Parallel_Check_EdgeGlobalID(mesh,rank,num,comm);
  return 1;

}

int MESH_Parallel_Check_VertexGlobalID(Mesh_ptr mesh, int rank, int num, MPI_Comm comm) {
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

  /* collect overlap vertex list for fast checking */
  idx = 0; nov = 0; ov_verts = List_New(10);
  while(mv = MESH_Next_Vertex(mesh, &idx))  {
    if(MV_PType(mv) == POVERLAP)  {
      List_Add(ov_verts,mv);
      nov++;
    }
  }
  List_Sort(ov_verts,nov,sizeof(MVertex_ptr),compareGlobalID);

  /* first half stores the global id, second half are the local ids */
  ov_list = (int *) MSTK_malloc(2*nov*sizeof(int));
  for(i = 0; i < nov; i++) {
    mv = List_Entry(ov_verts,i);
    ov_list[i] = MV_GlobalID(mv);
    ov_list[nov+i] = MV_ID(mv);
  }

  global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);
  max_nv = nv;
  for(i = 0; i < num; i++)
    if(max_nv < global_mesh_info[10*i])
      max_nv = global_mesh_info[10*i];
  
  /* collect data */
  int *list_vertex = (int *)MSTK_malloc(3*max_nv*sizeof(int));
  double *list_coor = (double *)MSTK_malloc(3*max_nv*sizeof(double));

  int *recv_list_vertex = (int *)MSTK_malloc(3*max_nv*sizeof(int));
  double *recv_list_coor = (double *)MSTK_malloc(3*max_nv*sizeof(double));

  /* similar as updateattr() and parallel_addghost(), order the send and recv process based on rank */
  for(i = 0; i < num; i++) {
    if(i == rank) continue;
    if(i < rank) {     
      if( MESH_Has_Overlaps_On_Prtn(mesh,i,MVERTEX) ) {
	MPI_Recv(recv_list_vertex,3*nv,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	printf("rank %d receives %d verts from rank %d\n", rank,count/3, i);
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
	    mv = MESH_Vertex(mesh, ov_list[nov + iloc]);  /* get the vertex on current mesh */
	    gdim = (recv_list_vertex[3*j] & 7);
	    gid = (recv_list_vertex[3*j] >> 3);
	    ptype = (recv_list_vertex[3*j+1] & 3);
	    master_id = (recv_list_vertex[3*j+1] >> 2);
	    if(MV_GEntDim(mv) != gdim) {
	      valid = 0;
	      sprintf(mesg,"Global edge %-d from processor %d and on processor %d GEndDim mismatch: %d vs %d ", global_id, i, rank, gdim, MV_GEntDim(mv));
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	    if(MV_GEntID(mv) != gid) {
	      valid = 0;
	      sprintf(mesg,"Global edge %-d from processor %d and on processor %d GEndID mismatch: %d vs %d ", global_id, i, rank, gid, MV_GEntID(mv));
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

	    if(!valid) {
	      MV_Coords(mv,coor);
	      printf("the mismatch vertex coor: (%f,%f,%f)\n", recv_list_coor[3*j], recv_list_coor[3*j+1], recv_list_coor[3*j+2]);
	    }
	  }
	}
      }
      if( MESH_Has_Ghosts_From_Prtn(mesh,i,MVERTEX) ) {
	idx = 0; index_mv = 0;
	while(mv = MESH_Next_Vertex(mesh, &idx)) {
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
	printf("rank %d sends %d verts to rank %d\n", rank,index_mv, i);
	MPI_Send(list_coor,3*index_mv,MPI_DOUBLE,i,i,comm);
	
      }
      
    }
    if(i > rank) {     
      if( MESH_Has_Ghosts_From_Prtn(mesh,i,MVERTEX) ) {
	idx = 0; index_mv = 0;
	while(mv = MESH_Next_Vertex(mesh, &idx)) {
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
	printf("rank %d sends %d verts to rank %d\n", rank,index_mv, i);
	MPI_Send(list_coor,3*index_mv,MPI_DOUBLE,i,i,comm);
	
      }
      if( MESH_Has_Overlaps_On_Prtn(mesh,i,MVERTEX) ) {
	MPI_Recv(recv_list_vertex,3*nv,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	printf("rank %d receives %d verts from rank %d\n", rank,count/3, i);
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
	    mv = MESH_Vertex(mesh, ov_list[nov + iloc]);  /* get the vertex on current mesh */
	    gdim = (recv_list_vertex[3*j] & 7);
	    gid = (recv_list_vertex[3*j] >> 3);
	    ptype = (recv_list_vertex[3*j+1] & 3);
	    master_id = (recv_list_vertex[3*j+1] >> 2);
	    if(MV_GEntDim(mv) != gdim) {
	      valid = 0;
	      sprintf(mesg,"Global edge %-d from processor %d and on processor %d GEndDim mismatch: %d vs %d ", global_id, i, rank, gdim, MV_GEntDim(mv));
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	    }
	    if(MV_GEntID(mv) != gid) {
	      valid = 0;
	      sprintf(mesg,"Global edge %-d from processor %d and on processor %d GEndID mismatch: %d vs %d ", global_id, i, rank, gid, MV_GEntID(mv));
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

	    if(!valid) {
	      MV_Coords(mv,coor);
	      printf("the mismatch vertex coor: (%f,%f,%f)\n", recv_list_coor[3*j], recv_list_coor[3*j+1], recv_list_coor[3*j+2]);
	    }

	  }
	}
      }
    }
  }
  MSTK_free(list_vertex);
  MSTK_free(list_coor);
  MSTK_free(recv_list_vertex);
  MSTK_free(recv_list_coor);
  MSTK_free(ov_list);
  MSTK_free(global_mesh_info);
  List_Delete(ov_verts);
  return valid;
}


int MESH_Parallel_Check_EdgeGlobalID(Mesh_ptr mesh, int rank, int num, MPI_Comm comm) {
  int i, j, idx, ne, noe, max_ne, index_me, global_id, *loc, iloc, valid = 1;
  char mesg[256], funcname[32] = "MESH_Parallel_Check";
  MPI_Status status;
  MEdge_ptr me;
  List_ptr ov_edges;
  int count, mesh_info[10], *global_mesh_info, *ov_list;
  for (i = 0; i < 10; i++) mesh_info[i] = 0;
  ne = MESH_Num_Edges(mesh);

  mesh_info[0] = ne;

  /* collect overlap edge list for fast checking */
  idx = 0; noe = 0; ov_edges = List_New(10);
  while(me = MESH_Next_Edge(mesh, &idx))  {
    if(ME_PType(me) == POVERLAP)  {
      List_Add(ov_edges,me);
      noe++;
    }
  }
  List_Sort(ov_edges,noe,sizeof(MEdge_ptr),compareGlobalID);

  /* first half stores the global id, second half are the local ids */
  ov_list = (int *) MSTK_malloc(2*noe*sizeof(int));
  for(i = 0; i < noe; i++) {
    me = List_Entry(ov_edges,i);
    ov_list[i] = ME_GlobalID(me);
    ov_list[noe+i] = ME_ID(me);
  }

  global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  max_ne = ne;
  for(i = 0; i < num; i++)
    if(max_ne < global_mesh_info[10*i])
      max_ne = global_mesh_info[10*i];
  
  /* collect data */
  int *list_edge = (int *)MSTK_malloc(5*max_ne*sizeof(int));
  int *recv_list_edge = (int *)MSTK_malloc(5*max_ne*sizeof(int));

  for(i = 0; i < num; i++) {
    if(i == rank) continue;
    if(i < rank) {     
      if( MESH_Has_Overlaps_On_Prtn(mesh,i,MEDGE) ) {
	MPI_Recv(recv_list_edge,5*ne,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	printf("rank %d receives %d edges from rank %d\n", rank,count/5, i);
	for(j = 0; j < count/5;  j++) {
	  global_id = recv_list_edge[5*j+4];
	  loc = (int *)bsearch(&global_id,
			       ov_list,
			       noe,
			       sizeof(int),
			       compareINT);
	  if(!loc) {
	    sprintf(mesg,"Global edge %-d on processor %d is not on processor %d ", global_id, i, rank);
	    MSTK_Report(funcname,mesg,MSTK_ERROR);
	    valid = 0;
	  }
	  /*
	  if(loc) {
	    iloc = (int)(loc - ov_list);
	    ME_Set_PType(me,PGHOST);
	    ME_Set_MasterParID(me,j);
	  }
	  */
	}
      }	
      if( MESH_Has_Ghosts_From_Prtn(mesh,i,MEDGE) ) {
	idx = 0; index_me = 0;
	while(me = MESH_Next_Edge(mesh, &idx)) {
	  if(ME_MasterParID(me) == i) {
	    list_edge[5*index_me]   = MV_ID(ME_Vertex(me,0));
	    list_edge[5*index_me+1] = MV_ID(ME_Vertex(me,1));
	    list_edge[5*index_me+2] = (ME_GEntID(me)<<3) | (ME_GEntDim(me));
	    list_edge[5*index_me+3] = (ME_MasterParID(me) <<2) | (ME_PType(me));
	    list_edge[5*index_me+4] = ME_GlobalID(me);
	    index_me++;
	  }
	}
	MPI_Send(list_edge,5*index_me,MPI_INT,i,i,comm);
	printf("rank %d sends %d edges to rank %d\n", rank,index_me, i);
      }
    }
    if(i > rank) {     
      if( MESH_Has_Ghosts_From_Prtn(mesh,i,MEDGE) ) {
	idx = 0; index_me = 0;
	while(me = MESH_Next_Edge(mesh, &idx)) {
	  if(ME_MasterParID(me) == i) {
	    list_edge[5*index_me]   = MV_ID(ME_Vertex(me,0));
	    list_edge[5*index_me+1] = MV_ID(ME_Vertex(me,1));
	    list_edge[5*index_me+2] = (ME_GEntID(me)<<3) | (ME_GEntDim(me));
	    list_edge[5*index_me+3] = (ME_MasterParID(me) <<2) | (ME_PType(me));
	    list_edge[5*index_me+4] = ME_GlobalID(me);
	    index_me++;
	  }
	}
	MPI_Send(list_edge,5*index_me,MPI_INT,i,i,comm);
	printf("rank %d sends %d elements to rank %d\n", rank,index_me, i);
      }
      if( MESH_Has_Overlaps_On_Prtn(mesh,i,MEDGE) ) {
	MPI_Recv(recv_list_edge,5*ne,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	printf("rank %d receives %d edges from rank %d\n", rank,count/5, i);
	for(j = 0; j < count/5;  j++) {
	  global_id = recv_list_edge[5*j+4];
	  loc = (int *)bsearch(&global_id,
			       ov_list,
			       noe,
			       sizeof(int),
			       compareINT);
	  if(!loc) {
	    sprintf(mesg,"Global edge %-d on processor %d is not on processor %d ", global_id, i, rank);
	    MSTK_Report(funcname,mesg,MSTK_ERROR);
	    valid = 0;
	  }
	  /*
	  if(loc) {
	    iloc = (int)(loc - ov_list);
	    ME_Set_PType(me,PGHOST);
	    ME_Set_MasterParID(me,j);
	  }
	  */
	}
      }	
    }
  }
  MSTK_free(list_edge);
  MSTK_free(recv_list_edge);
  MSTK_free(ov_list);
  MSTK_free(global_mesh_info);
  List_Delete(ov_edges);
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
  printf("checking mesh %d \n", (int)mesh);
  /* check vertex */
  while(mv = MESH_Next_Vertex(mesh, &idx)) {
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
  while(me = MESH_Next_Edge(mesh, &idx)) {
    if(ME_PType(me) == PGHOST) {
      if (ME_MasterParID(me) == rank) {
	sprintf(mesg,"Ghost Edge %-d has master partition id %d but it is on partition %d", ME_GlobalID(me), ME_MasterParID(me),rank);
	MSTK_Report(funcname,mesg,MSTK_ERROR);
	valid = 0;
      }
    }
    else  {
      if(ME_MasterParID(me) != rank) {
	sprintf(mesg,"Non-Ghost Edge %-d has master partition id %d but it is on partition %d", ME_GlobalID(me), ME_MasterParID(me),rank);
	MSTK_Report(funcname,mesg,MSTK_ERROR);
	valid = 0;
      }	
    }
  }

  /* check face */
  idx = 0;
  while(mf = MESH_Next_Face(mesh, &idx)) {
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
  while(mr = MESH_Next_Region(mesh, &idx)) {
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
