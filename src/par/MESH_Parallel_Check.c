#include <stdio.h>
#include <stdlib.h>

#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif

  /* Routine to check that the topological structure of the mesh is valid */
  /* Rao Garimella - 12/16/2010                                           */
int MESH_Parallel_Check_Ghost(Mesh_ptr mesh, int rank);
int MESH_Parallel_Check_VertexGlobalID(Mesh_ptr mesh, int rank, int num, MPI_Comm comm);
int MESH_Parallel_Check_GlobalID(Mesh_ptr mesh, int rank, int num, MPI_Comm comm);
int MESH_Parallel_Check(Mesh_ptr mesh, int rank, int num, MPI_Comm comm) {
  int valid;
  //valid = MESH_CheckTopo(mesh);
  valid = MESH_Parallel_Check_Ghost(mesh,rank);
  valid = MESH_Parallel_Check_GlobalID(mesh,rank,num,comm);
  return valid;

} /* int MESH_Parallel_Check */

int MESH_Parallel_Check_GlobalID(Mesh_ptr mesh, int rank, int num, MPI_Comm comm) {
  MESH_Parallel_Check_VertexGlobalID(mesh,rank,num,comm);
  return 1;

}

int MESH_Parallel_Check_VertexGlobalID(Mesh_ptr mesh, int rank, int num, MPI_Comm comm) {
  int i, idx, nv, max_nv, index_mv;
  int nevs, nfes, nrfs, nfe, nrv, nrf, natt, nset, ncomp, dir;
  MPI_Status status;
  MVertex_ptr mv;
  MEdge_ptr me;
  MFace_ptr mf;
  MRegion_ptr mr;
  List_ptr mfedges, mrfaces, mrverts;
  double coor[3];
  int count, mesh_info[10], *global_mesh_info;
  for (i = 0; i < 10; i++) mesh_info[i] = 0;
  nv = MESH_Num_Vertices(mesh);
  mesh_info[0] = nv;
  global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);
  max_nv = nv;
  for(i = 0; i < num; i++)
    if(max_nv < global_mesh_info[10*i])
      max_nv = global_mesh_info[10*i];

  /* collect data */
  int *list_vertex = (int *)MSTK_malloc(3*nv*sizeof(int));
  double *list_coor = (double *)MSTK_malloc(3*nv*sizeof(double));

  int *recv_list_vertex = (int *)MSTK_malloc(3*max_nv*sizeof(int));
  double *recv_list_coor = (double *)MSTK_malloc(3*max_nv*sizeof(double));

  for(i = 0; i < num; i++) {
    if(i == rank) continue;
    if(i < rank) {     
      if( MESH_Has_Overlaps_On_Prtn(mesh,i,MVERTEX) ) {
	MPI_Recv(recv_list_vertex,3*nv,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	printf("rank %d receives %d elements from rank %d\n", rank,count/3, i);
	MPI_Recv(recv_list_coor,3*nv,MPI_DOUBLE,i,rank,comm,&status);

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
	printf("rank %d sends %d elements to rank %d\n", rank,index_mv, i);
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
	printf("rank %d sends %d elements to rank %d\n", rank,index_mv, i);
	MPI_Send(list_coor,3*index_mv,MPI_DOUBLE,i,i,comm);
	
      }
      if( MESH_Has_Overlaps_On_Prtn(mesh,i,MVERTEX) ) {
	MPI_Recv(recv_list_vertex,3*nv,MPI_INT,i,rank,comm,&status);
	MPI_Get_count(&status,MPI_INT,&count);
	printf("rank %d receives %d elements from rank %d\n", rank,count/3, i);
	MPI_Recv(recv_list_coor,3*nv,MPI_DOUBLE,i,rank,comm,&status);
      }
    }
  }
}
/* 
   check if every ghost entity has a master partition number that is different from current partition 
   check if other PType entity has the same master partition number as this parition
*/
int MESH_Parallel_Check_Ghost(Mesh_ptr mesh, int rank) {
  int nv, ne, nf, nr, idx, valid;
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
