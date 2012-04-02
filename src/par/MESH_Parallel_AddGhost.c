#define _H_Mesh_Private

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "Mesh.h"
#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif



  /* 
     This function is a collective call

     It send the OVERLAP entities to neighbor processors

     Assume OVERLAP list is properly built 

     Author(s): Duo Wang, Rao Garimella
  */

int MESH_Parallel_AddGhost_Face(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);
int MESH_Parallel_AddGhost_Region(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);


int MESH_Parallel_AddGhost(Mesh_ptr submesh, int rank, int num,  MPI_Comm comm) {
  int nf, nr;
  RepType rtype;

  /* basic mesh information */

  rtype = MESH_RepType(submesh);
  nf = MESH_Num_Faces(submesh);
  nr = MESH_Num_Regions(submesh);
  if (nr)
    MESH_Parallel_AddGhost_Region(submesh, rank, num, comm);
  else if(nf) 
    MESH_Parallel_AddGhost_Face(submesh, rank, num, comm);
  else {
    MSTK_Report("MESH_Parallel_AddGhost()","only send volume or surface mesh",MSTK_ERROR);
    exit(-1);
  }
  return 1;
}


  /* 
     Send 1-ring Faces to neighbor processors, and receive them 
     First update the parallel adjancy information, 
  */
int MESH_Parallel_AddGhost_Face(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  RepType rtype;
  int i, j, num_recv_procs;
  MVertex_ptr mv;
  MEdge_ptr me;
  List_ptr fedges, fverts;
  int adj;
  /* allocate meshes to receive from other processors */
  /*
  num_recv_procs = MESH_Num_GhostPrtns(submesh);
  printf(" number of recv_procs %d,on rank %d\n",rank, num_recv_procs);
  */
  /* 
     first update parallel adjancy information
     any two processor that has vertex connection now has face connection
  */


  for (i = 0; i < num; i++) {
    if(i == rank) continue;
    if( adj = MESH_Has_Ghosts_From_Prtn(submesh,i,MANYTYPE) ) {
      MESH_Flag_Has_Ghosts_From_Prtn(submesh,i,MFACE);
      MESH_Flag_Has_Overlaps_On_Prtn(submesh,i,MVERTEX);
      MESH_Flag_Has_Overlaps_On_Prtn(submesh,i,MFACE);
    }
    if( adj = MESH_Has_Overlaps_On_Prtn(submesh,i,MANYTYPE) ) {
      MESH_Flag_Has_Overlaps_On_Prtn(submesh,i,MFACE);
      MESH_Flag_Has_Ghosts_From_Prtn(submesh,i,MVERTEX);
      MESH_Flag_Has_Ghosts_From_Prtn(submesh,i,MFACE);
    }
  }
  /*
  for(i = 0; i < num; i++) {
    if (adj = MESH_Has_Overlaps_On_Prtn(submesh,i,MFACE))
      printf("after rank %d has overlap face on rank %d\n", rank,i);
    if (adj = MESH_Has_Ghosts_From_Prtn(submesh,i,MFACE))
      printf("after rank %d has ghost face from rank %d\n", rank,i);

  } 
  */
  /* 
     printf("before rank %d,  num of vertex %d, number of face %d\n",rank,MESH_Num_Vertices(submesh),MESH_Num_Faces(submesh));
     printf("before rank %d,  num of ov vertex %d, number of ov face %d\n",rank,MESH_Num_OverlapVertices(submesh),MESH_Num_OverlapFaces(submesh));
     send and receive overlap meshes 
  */

  for (i = 0; i < num; i++) {
    if(i == rank) continue;
    if(i < rank) {     
      if( MESH_Has_Ghosts_From_Prtn(submesh,i,MFACE) ) {
	MESH_RecvMesh_Elements(submesh,2,i,rank,comm);
	printf("rank %d receives elements from rank %d\n", rank,i);
      }
      if( MESH_Has_Overlaps_On_Prtn(submesh,i,MFACE) ) {
	MESH_SendMesh_Elements(submesh,MESH_OverlapFace_List(submesh),i,comm);
	printf("rank %d sends elements to rank %d\n", rank,i);
      }
    }
    if(i > rank) {     
      if( MESH_Has_Overlaps_On_Prtn(submesh,i,MFACE) ) {
	MESH_SendMesh_Elements(submesh,MESH_OverlapFace_List(submesh),i,comm);
	printf("rank %d sends elements to rank %d\n", rank,i);
      }
      if( MESH_Has_Ghosts_From_Prtn(submesh,i,MFACE) )
	MESH_RecvMesh_Elements(submesh,2,i,rank,comm); {
	printf("rank %d receives elements from rank %d\n", rank,i);
      }
    }
  }

  //  printf("after  rank %d, ov_mesh num of vertex %d, number of face %d\n",rank,MESH_Num_Vertices(submesh),MESH_Num_Faces(submesh));

  return 1;
}

  /* right now assume there are no overlapped regions */

int MESH_Parallel_AddGhost_Region(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
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
  for(i = 0; i < nf; i++) {
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

