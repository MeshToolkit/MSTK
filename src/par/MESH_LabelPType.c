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
     It assign vertex PType similar as MESH_AssignGlobalIDs(), 
     except that here global IDs are given

     It assigns all the elements with a POVERLAP vertex as POVERLAP element

     It also assign master partition ID

     Author(s): Duo Wang, Rao Garimella
  */

  int MESH_LabelPType_Vertex(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);
  int MESH_LabelPType_Edge(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);
  int MESH_LabelPType_Face(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);
  int MESH_LabelPType_Region(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);

  int MESH_LabelPType(Mesh_ptr submesh, int rank, int num,  MPI_Comm comm) {
  int nf, nr;
  RepType rtype;

   /* basic mesh information */
  rtype = MESH_RepType(submesh);
  nf = MESH_Num_Faces(submesh);
  nr = MESH_Num_Regions(submesh);
  /* first label vertex */
  MESH_LabelPType_Vertex(submesh, rank, num, comm);
  MESH_LabelPType_Edge(submesh, rank, num, comm);
  if (nr)
    MESH_LabelPType_Region(submesh, rank, num, comm);
  else if(nf) 
    MESH_LabelPType_Face(submesh, rank, num, comm);
  else {
    MSTK_Report("MESH_LabelPType()","only send volume or surface mesh",MSTK_ERROR);
    exit(-1);
  }
  return 1;
}

  /* 
     vertex PType information is the key for edge, face and region PType 
     
     inter process communication is required here in case global ID is already assigned 
  */
int MESH_LabelPType_Vertex(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  int i, j, nv, nbv, ne, nf, mesh_info[10];
  MVertex_ptr mv;
  List_ptr boundary_verts;
  RepType rtype;
  int *loc, *global_mesh_info, *list_boundary_vertex, *recv_list_vertex, *vertex_ov_label;
  int iloc,  global_id, max_nbv, index_nbv;

  for (i = 0; i < 10; i++) mesh_info[i] = 0;

  /* mesh_info store the mesh reptype, nv, ne, nf, nbv */
  rtype = MESH_RepType(submesh);
  nv = MESH_Num_Vertices(submesh);
  ne = MESH_Num_Edges(submesh);
  nf = MESH_Num_Faces(submesh);

  mesh_info[0] = rtype;
  mesh_info[1] = nv;
  mesh_info[2] = ne;
  mesh_info[3] = nf;

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

  /* 
     gather submeshes information
  */
  global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  /* get largest number of boundary vertices of all the processors */
  max_nbv = 0;
  for(i = 0; i < num; i++)
    if(max_nbv < global_mesh_info[10*i+4])
      max_nbv = global_mesh_info[10*i+4];

  list_boundary_vertex = (int *)MSTK_malloc(max_nbv*sizeof(int));
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

  /* indicate if a vertex is overlapped */
  vertex_ov_label = (int *)MSTK_malloc(num*max_nbv*sizeof(int));

  /* 
     store the local boundary id on ov processor
     it is used to assign global id of local ghost vertices
     do not need to store master partition id, MV_MasterParID(mv) is assigned
  */
  for (i = 0; i < num*max_nbv; i++)
    vertex_ov_label[i] = 0;
  /* for processor other than 0 */
  for(i = 0; i < nbv; i++) {
    mv = List_Entry(boundary_verts,i);
    global_id = MV_GlobalID(mv);
    /* 
       check which processor has the same coordinate vertex 
       Different from assigning global id, we check all the other processors, from large to small
       by rank, whenever a vertex is found, mark rank and j as neighbors, and the masterparid is still
       the smallest rank processor
    */
    for(j = num-1; j >= 0; j--) {
      if(j == rank) continue;
      /* since each processor has sorted the boundary vertices, use binary search */
      loc = (int *)bsearch(&global_id,
			   &recv_list_vertex[max_nbv*j],
			   global_mesh_info[10*j+4],
			   sizeof(int),
			   compareINT);
      /* if found the vertex on previous processors */
      if(loc) {
	/* here the location iloc is relative to the beginning of the jth processor */
	iloc = (int)(loc - &recv_list_vertex[max_nbv*j]);	
	MESH_Flag_Has_Ghosts_From_Prtn(submesh,j,MVERTEX);
	if(j < rank) {
	  MV_Set_PType(mv,PGHOST);
	  MV_Set_MasterParID(mv,j);
	} 
	/* label the original vertex as overlapped */
	vertex_ov_label[max_nbv*j+iloc] |= 1;
	/* if found on processor j, still need to check other processors*/
      }
    }
  }

  /* since this is a or reduction, we can use MPI_IN_PLACE, send buffer same as recv buffer */
  MPI_Allreduce(MPI_IN_PLACE,vertex_ov_label,num*max_nbv,MPI_INT,MPI_LOR,comm);    

  for(i = 0; i < nv; i++) {
    mv = MESH_Vertex(submesh,i);
    if (MV_PType(mv) == PGHOST)
      continue;
    MV_Set_MasterParID(mv,rank);
  }

  /* label OVERLAP vertices */
  for(i = 0; i < nbv; i++) {
    if(vertex_ov_label[rank*max_nbv+i]) {
      mv = List_Entry(boundary_verts,i);
      if(MV_PType(mv)!=PGHOST)
	MV_Set_PType(mv,POVERLAP);
    }
  }

  List_Delete(boundary_verts);
  MSTK_free(global_mesh_info);
  MSTK_free(vertex_ov_label);

  MSTK_free(list_boundary_vertex);
  MSTK_free(recv_list_vertex);
  return 1;
}
  
int MESH_LabelPType_Edge(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  int idx;
  MVertex_ptr mv1, mv2;
  MEdge_ptr me;
  idx = 0;
  while( (me = MESH_Next_Edge(submesh,&idx)) ) {
    if(ME_PType(me) != PGHOST)
      ME_Set_MasterParID(me,rank);
    /*
      if both ends are ghost vertices, then the edge is a ghost
      Do not set master partition id right now
    */
    mv1 = ME_Vertex(me,0); mv2 = ME_Vertex(me,1);  
    if(MV_PType(mv1) == PGHOST && MV_PType(mv2) == PGHOST) {
      ME_Set_PType(me,PGHOST);
    }
    
    /* 
       if either ends are overlap vertices, then the edge is a overlap
       this may produce more than necessary overlap edges 
    */
    if(MV_PType(mv1) == POVERLAP || MV_PType(mv2) == POVERLAP) {
      ME_Set_PType(me,POVERLAP);
    }
  }
  
  return 1;
}

  /* 
     Right now assume faces are not overlapped across processors
     Label the face that has OVERLAP or Ghost vertex as OVERLAP
     
  */
int MESH_LabelPType_Face(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  int i, idx;
  MVertex_ptr mv;
  MEdge_ptr me;
  MFace_ptr mf;
  List_ptr fverts, fedges;
  idx = 0;
  while( (mf = MESH_Next_Face(submesh,&idx)) ) {
    fverts = MF_Vertices(mf,1,0);
    for(i = 0; i < List_Num_Entries(fverts); i++) {
      mv = List_Entry(fverts,i);
      if( MV_PType(mv) == PGHOST || MV_PType(mv) == POVERLAP ) {
	MF_Set_PType(mf,POVERLAP);
	break;
      }
    }
    List_Delete(fverts);
  }
    /* if the face is marked, mark its edges and vertices */
  idx = 0;
  while( (mf = MESH_Next_Face(submesh,&idx)) ) {
    if(MF_PType(mf) != POVERLAP) continue;
    fverts = MF_Vertices(mf,1,0);
    for(i = 0; i < List_Num_Entries(fverts); i++) {
      mv = List_Entry(fverts,i);
      if(MV_PType(mv) != PGHOST)
	MV_Set_PType(mv,POVERLAP);
    }
    List_Delete(fverts);
    
    fedges = MF_Edges(mf,1,0);
    for(i = 0; i < List_Num_Entries(fedges); i++) {
      me = List_Entry(fedges,i);
      if(ME_PType(me) != PGHOST)
	ME_Set_PType(me,POVERLAP);
    }
    List_Delete(fedges);
  }
  return 1;
}

  /* 
     right now assume regions are not overlapped across processors
     Label the region that has OVERLAP vertex as OVERLAP
  */


int MESH_LabelPType_Region(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  int i, nfv, idx;
  MRegion_ptr mr;
  MFace_ptr mf;
  MVertex_ptr me;
  MVertex_ptr mv;
  List_ptr mfverts, rverts, redges, rfaces;
  int is_ghost, is_overlap;
  idx = 0;
  while( (mf = MESH_Next_Face(submesh,&idx)) ) {
    if(MF_PType(mf) != PGHOST)
      MF_Set_MasterParID(mf, rank);
    is_ghost = 1; is_overlap = 0;
    mfverts = MF_Vertices(mf,1,0);
    nfv = List_Num_Entries(mfverts);
    for(i = 0; i < nfv; i++) {
      mv = List_Entry(mfverts,i);
      if(MV_PType(mv) != PGHOST)  /* if all vertices are ghost, then the face is a ghost*/
	is_ghost = 0;
      if(MV_PType(mv) == POVERLAP) /*  if either vertex is a overlap, then the face is a overlap */
	is_overlap = 1;
    }
    if(is_ghost) 
      MF_Set_PType(mf,PGHOST);
    
    if(is_overlap) 
      MF_Set_PType(mf,POVERLAP);

    List_Delete(mfverts);
  }
  
  idx = 0;
  while( (mr = MESH_Next_Region(submesh,&idx)) ) {
    rverts = MR_Vertices(mr);
    for(i = 0; i < List_Num_Entries(rverts); i++) {
      mv = List_Entry(rverts,i);
      if(MV_PType(mv) == PGHOST || MV_PType(mv) == POVERLAP) {
	MR_Set_PType(mr,POVERLAP);
	break;
      }
    }
  }

  /* if the region is marked, mark its vertices, edges and faces */
  idx = 0;
  while( (mr = MESH_Next_Region(submesh,&idx)) ) {
    if(MV_PType(mr) != POVERLAP) continue;
    rverts = MR_Vertices(mr);
    for(i = 0; i < List_Num_Entries(rverts); i++) {
      mv = List_Entry(rverts,i);
      if(MV_PType(mv) != PGHOST)
	MV_Set_PType(mv,POVERLAP);
    }
    List_Delete(rverts);
    /* if the region is marked, mark its edges */
    redges = MR_Edges(mr);
    for(i = 0; i < List_Num_Entries(redges); i++) {
      me = List_Entry(redges,i);
      if(ME_PType(me) != PGHOST)
	ME_Set_PType(me,POVERLAP);
    }
    List_Delete(redges);
    /* if the region is marked, mark its faces */
    rfaces = MR_Faces(mr);
    for(i = 0; i < List_Num_Entries(rfaces); i++) {
      mf = List_Entry(rfaces,i);
      if(MF_PType(mf) != PGHOST)
	  MF_Set_PType(mf,POVERLAP);
    }
    List_Delete(rfaces);
  }
  return 1;
}


  
#ifdef __cplusplus
}
#endif

