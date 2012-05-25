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

     It builds connection information across processors based on vertex global id


     Author(s): Duo Wang, Rao Garimella
  */

int MESH_BuildConnection(Mesh_ptr submesh, int rank, int num,  MPI_Comm comm) {
  int i, j, nv, nr, nbv, mesh_info[10];
  MVertex_ptr mv;
  List_ptr boundary_verts;
  int *loc, *global_mesh_info, *list_boundary_vertex, *recv_list_vertex;
  int iloc,  global_id, max_nbv, index_nbv;

  for (i = 0; i < 10; i++) mesh_info[i] = 0;
  nv = MESH_Num_Vertices(submesh);
  nr = MESH_Num_Regions(submesh);
  mesh_info[1] = nv;

  /* calculate number of boundary vertices */ 
  nbv = 0;
  boundary_verts = List_New(10);
  for(i = 0; i < nv; i++) {
    mv = MESH_Vertex(submesh,i);
    if (MV_PType(mv)) {
      List_Add(boundary_verts,mv);
      nbv++;
    }
  }
  mesh_info[4] = nbv;
  
  /* sort boundary vertices based on global ID, for binary search */
  List_Sort(boundary_verts,nbv,sizeof(MVertex_ptr),compareGlobalID);

  global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  max_nbv = 0;
  for(i = 0; i < num; i++)
    if(max_nbv < global_mesh_info[10*i+4])
      max_nbv = global_mesh_info[10*i+4];

  list_boundary_vertex = (int *)MSTK_malloc(nbv*sizeof(int));
  recv_list_vertex = (int *)MSTK_malloc(num*max_nbv*sizeof(int));

  /* only global ID are sent */
  index_nbv = 0;
  for(i = 0; i < nbv; i++) {
    mv = List_Entry(boundary_verts,i);
    list_boundary_vertex[index_nbv] = MV_GlobalID(mv);
    index_nbv++;
  }

  /* gather boundary vertices */
  MPI_Allgather(list_boundary_vertex,max_nbv,MPI_INT,recv_list_vertex,max_nbv,MPI_INT,comm);

  for(i = 0; i < nbv; i++) {
    mv = List_Entry(boundary_verts,i);
    global_id = MV_GlobalID(mv);
    /* 
       check which processor has the vertex of the same global ID
       Different from assigning global id, we check all the other processors, from large to small
       by rank, whenever a vertex is found, mark it as neighbors, the masterparid is the smallest 
       rank processor among all the neighbors
    */
    for(j = num-1; j >= 0; j--) {
      if(j == rank) continue;
      loc = (int *)bsearch(&global_id,
			   &recv_list_vertex[max_nbv*j],
			   global_mesh_info[10*j+4],
			   sizeof(int),
			   compareINT);
      /* if found the vertex on previous processors */
      if(loc) {
	iloc = (int)(loc - &recv_list_vertex[max_nbv*j]);	
	/*
	if(j > rank) {
	  MESH_Flag_Has_Overlaps_On_Prtn(submesh,j,MVERTEX);
	}
	else {
	  MESH_Flag_Has_Ghosts_From_Prtn(submesh,j,MVERTEX);
	}
	*/	
	MESH_Flag_Has_Ghosts_From_Prtn(submesh,j,MVERTEX);
	MESH_Flag_Has_Ghosts_From_Prtn(submesh,j,MEDGE);
	//	MESH_Flag_Has_Overlaps_On_Prtn(submesh,j,MVERTEX);
	//MESH_Flag_Has_Overlaps_On_Prtn(submesh,j,MEDGE);
	if(nr) {
	  MESH_Flag_Has_Ghosts_From_Prtn(submesh,j,MFACE);
	  //MESH_Flag_Has_Overlaps_On_Prtn(submesh,j,MFACE);
	}

	/* if found on processor j, still need to check j-1, j-2 ...*/
      }
    }
  }

  MESH_Update_ParallelAdj(submesh, rank, num, MPI_COMM_WORLD);

  List_Delete(boundary_verts);
  MSTK_free(global_mesh_info);
  MSTK_free(list_boundary_vertex);
  MSTK_free(recv_list_vertex);
  return 1;
}
  
#ifdef __cplusplus
}
#endif

