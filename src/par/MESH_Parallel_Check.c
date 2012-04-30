#include <stdio.h>
#include <stdlib.h>

#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif

  /* Routine to check that the topological structure of the mesh is valid */
  /* Rao Garimella - 12/16/2010                                           */
int MESH_Parallel_Check_Ghost(Mesh_ptr mesh, int rank);
int MESH_Parallel_Check(Mesh_ptr mesh, int rank, int num, MPI_Comm comm) {
  int valid;
  //  valid = MESH_CheckTopo(mesh);
  valid = MESH_Parallel_Check_Ghost(mesh,rank);
  return valid;

} /* int MESH_Parallel_Check */


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
