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
     THIS IS AN INTERNAL MSTK CALL THAT SHOULD ***ONLY*** BE USED AS
     ONE STEP IN A SERIES OF WELL COORDINATED STEPS DURING
     DISTRIBUTING A MESH FROM ONE PROCESSOR TO MANY.

     Assign global IDs of vertices (by comparing coordinates, if not
     already given), and use global IDs of vertices to match up and
     assign global IDs of edges, faces and regions on each
     partition. Also assign the proper master partition ID and
     parallel type for entities on the parallel boundary.

     At this stage, there are no ghost layers and the partitions are
     not in their final state.

     **** IF THE MESH HAS DEGENERATE EDGES (COINCIDENT VERTICES) ON
     **** PARALLEL BOUNDARIES, CALLER MUST PROVIDE UNIQUE GIDs FOR THE
     **** VERTICES - OTHERWISE, THE WRONG GLOBAL ID MAY BE ASSIGNED TO
     **** GHOST VERTICES SINCE THE CODE HAS NO GOOD WAY OF KNOWING
     **** WHICH COINCIDENT VERTEX ON A REMOTE PROCESSOR A PARTICULAR
     **** GHOST VERTEX CORRESPONDS TO ***

     Author(s): Duo Wang, Rao Garimella
  */

  int MESH_AssignGlobalIDs_Vertex(Mesh_ptr submesh, int have_GIDs, MSTK_Comm comm);
  int MESH_AssignGlobalIDs_Edge(Mesh_ptr submesh, MSTK_Comm comm);
  int MESH_AssignGlobalIDs_Face(Mesh_ptr submesh, MSTK_Comm comm);
  int MESH_AssignGlobalIDs_Region(Mesh_ptr submesh, MSTK_Comm comm);


  int MESH_AssignGlobalIDs(Mesh_ptr submesh, int topodim, int have_GIDs, MSTK_Comm comm) {

  int nf, nr;
  RepType rtype;


  MESH_AssignGlobalIDs_Vertex(submesh, have_GIDs, comm);
  MESH_AssignGlobalIDs_Edge(submesh,comm);

  /* Use topological dimension of mesh to call appropriate function */
  /* Do not use the presence of mesh regions or faces since the mesh
     may be an empty mesh on some processors */

  if (topodim == 3)
    MESH_AssignGlobalIDs_Region(submesh,comm);
  else if (topodim == 2) 
    MESH_AssignGlobalIDs_Face(submesh,comm);
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
static int vertex_on_boundary2D(MVertex_ptr mv) {
  int i, nve, ok = 0;
  MEdge_ptr me;
  List_ptr vedges = MV_Edges(mv);
  nve = List_Num_Entries(vedges);
  for(i = 0; (i < nve) &!ok; i++) {
    List_ptr efaces;
    me = List_Entry(vedges,i);
    efaces = ME_Faces(me);
    if(List_Num_Entries(efaces) <= 1)
      ok = 1;
    List_Delete(efaces);
  }
  List_Delete(vedges);
  return ok;
}

  /* 
     test if a vertex is on partition boundary, for volume mesh
     a vertex is on boundary if one of its incident faces has less than 2 neighboring regions
 */
static int vertex_on_boundary3D(MVertex_ptr mv) {
  int i, nrf, ok = 0;
  MFace_ptr mf;
  List_ptr vfaces = MV_Faces(mv);
  nrf = List_Num_Entries(vfaces);
  for(i = 0; (i < nrf) && !ok; i++) {
    List_ptr fregs;
    mf = List_Entry(vfaces,i);
    fregs = MF_Regions(mf);
    if(List_Num_Entries(fregs) <= 1)
      ok = 1;
    List_Delete(fregs);
  }
  List_Delete(vfaces);
  return ok;
}

 /* 
    Assign vertex global ID for both 2D and 3D meshes
    Also set PGHOST and POVERLAP, used for edge, face and region global ID 
    MESH_LabelPType() does not rely on PType set here, it will reset PType,
    so that user can call MESH_LabelPType() directly when global IDs are given
 */
     
  int MESH_AssignGlobalIDs_Vertex(Mesh_ptr submesh, int have_GIDs, MSTK_Comm comm) {
  int i, j, nv, nbv, ne, nf, nr, mesh_info[10];
  MVertex_ptr mv;
  List_ptr boundary_verts;
  RepType rtype;
  int index_nbv, max_nbv, iloc, num_ghost_verts, global_id;
  int *global_mesh_info, *vertex_ov_label, *vertex_ov_global_id, *id_on_ov_list;

  int rank, num;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&num);

  for (i = 0; i < 10; i++) mesh_info[i] = 0;

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
  nbv = 0;  boundary_verts = List_New(10);
  if (nr) {
    for(i = 0; i < nv; i++) {
      mv = MESH_Vertex(submesh,i);
      if (vertex_on_boundary3D(mv)) {
        MV_Flag_OnParBoundary(mv);
        List_Add(boundary_verts,mv);
        nbv++;
      }
    }
  }
  else {
    for(i = 0; i < nv; i++) {
      mv = MESH_Vertex(submesh,i);
      if (vertex_on_boundary2D(mv)) {
        MV_Flag_OnParBoundary(mv);
        List_Add(boundary_verts,mv);
        nbv++;
      }
    }
  }
  mesh_info[5] = nbv;
  
  /* 
     gather submeshes information
     right now we only need nv and nbv, and later num_ghost_verts, but we gather all mesh_info
  */
  global_mesh_info = (int *)malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  /* get largest number of boundary vertices of all the processors */
  max_nbv = 0;
  for(i = 0; i < num; i++)
    if(max_nbv < global_mesh_info[10*i+5])
      max_nbv = global_mesh_info[10*i+5];

  if (have_GIDs) {
    int *list_boundary_vertex_gid = (int *)malloc(max_nbv*sizeof(int));

    int *recv_list_vertex_gid = (int *)malloc(num*max_nbv*sizeof(int));
    
    /* sort boundary vertices based on Global ID, for binary search */
    List_Sort(boundary_verts,nbv,sizeof(MVertex_ptr),compareGlobalID);

    /* only global ids are sent */
    index_nbv = 0;
    for(i = 0; i < nbv; i++) {
      mv = List_Entry(boundary_verts,i);
      list_boundary_vertex_gid[index_nbv] = MV_GlobalID(mv);
      index_nbv++;
    }
    
    MPI_Allgather(list_boundary_vertex_gid,max_nbv,MPI_INT,recv_list_vertex_gid,max_nbv,MPI_INT,comm);
    
    /* indicate if a vertex is overlapped */
    vertex_ov_label = (int *)malloc(num*max_nbv*sizeof(int));
    
    /* 
       store the local boundary id on ov processor
       it is used to assign global id of local ghost vertices
       no need to store master partition id, MV_MasterParID(mv) is already assigned
    */
    id_on_ov_list = (int *)malloc(max_nbv*sizeof(int));

    for (i = 0; i < num*max_nbv; i++)
      vertex_ov_label[i] = 0;
    num_ghost_verts = 0;
    /* for processor other than 0 */
    if(rank > 0) {
      for(i = 0; i < nbv; i++) {
        mv = List_Entry(boundary_verts,i);
        int gid = MV_GlobalID(mv);
        /* check which previous processor has a vertex with same Global ID*/
        for(j = 0; j < rank; j++) {
          /* since each processor has sorted the boundary vertices, use binary search */
          int *loc = (int *)bsearch(&gid,
                                    &recv_list_vertex_gid[max_nbv*j],
                                    global_mesh_info[10*j+5],
                                    sizeof(int),
                                    compareINT);
          /* if found the vertex on previous processors */
          if(loc) {
            /* here the location iloc is relative to the beginning of the jth processor */
            iloc = (int)(loc - &recv_list_vertex_gid[max_nbv*j]);
            MV_Set_PType(mv,PGHOST);
            MV_Set_MasterParID(mv,j);
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

    free(list_boundary_vertex_gid);
    free(recv_list_vertex_gid);
  }
  else {
    double coor[3];

    int *list_boundary_vertex = (int *)malloc(max_nbv*sizeof(int));
    double *list_boundary_coor = (double *)malloc(3*max_nbv*sizeof(double));

    int *recv_list_vertex = (int *)malloc(num*max_nbv*sizeof(int));
    double *recv_list_coor = (double *)malloc(3*num*max_nbv*sizeof(double));
    
    /* sort boundary vertices based on coordinate value, for binary search */
    List_Sort(boundary_verts,nbv,sizeof(MVertex_ptr),compareVertexCoor);

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
    
    MPI_Allgather(list_boundary_vertex,max_nbv,MPI_INT,recv_list_vertex,max_nbv,MPI_INT,comm);
    MPI_Allgather(list_boundary_coor,3*max_nbv,MPI_DOUBLE,recv_list_coor,3*max_nbv,MPI_DOUBLE,comm);
    
    /* indicate if a vertex is overlapped */
    vertex_ov_label = (int *)malloc(num*max_nbv*sizeof(int));
    
    /* 
       store the local boundary id on ov processor
       it is used to assign global id of local ghost vertices
       no need to store master partition id, MV_MasterParID(mv) is already assigned
    */
    id_on_ov_list = (int *)malloc(max_nbv*sizeof(int));

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
          double *loc = (double *)bsearch(&coor,
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

    free(list_boundary_coor);
    free(recv_list_coor);
    free(list_boundary_vertex);
    free(recv_list_vertex);
  }


  /* num of ghost verts */
  mesh_info[9] = num_ghost_verts;
  /* update ghost verts number */
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);
  /* since this is a OR reduction, we can use MPI_IN_PLACE, send buffer same as recv buffer */
  MPI_Allreduce(MPI_IN_PLACE,vertex_ov_label,num*max_nbv,MPI_INT,MPI_LOR,comm);    

  /* calculate starting global id number for vertices*/
  if (!have_GIDs) {
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
  }

      

  /* store overlapped vertices IDs and broadast */
  vertex_ov_global_id = (int *)malloc(num*max_nbv*sizeof(int));
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
  free(global_mesh_info);
  free(vertex_ov_label);
  free(vertex_ov_global_id);
  free(id_on_ov_list);
  
  return 1;
}


 /* 
    Assign edge global ID for both 2D and 3D meshes

    Assume each processor knowns its overlap and ghost vertices 
 */
     
  int MESH_AssignGlobalIDs_Edge(Mesh_ptr submesh, MSTK_Comm comm) {
  int i, j, k, nbe, noe, nge, ne, mesh_info[10], global_id;
  MVertex_ptr mv;
  MEdge_ptr me;
  List_ptr boundary_edges;
  int *loc, edge_id[2], max_nbe, index_nbe, iloc, is_boundary;
  int *global_mesh_info, *list_edge, *recv_list_edge, *edge_ov_label, *id_on_ov_list;
  
  int rank, num;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&num);

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
      ME_Flag_OnParBoundary(me);
      List_Add(boundary_edges,me);
      nbe++;
    }
  }
  mesh_info[6] = nbe;

  List_Sort(boundary_edges,nbe,sizeof(MEdge_ptr),compareEdgeID);

  global_mesh_info = (int *)malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  max_nbe = 0;
  for(i = 0; i < num; i++)
    if(max_nbe < global_mesh_info[10*i+6])
      max_nbe = global_mesh_info[10*i+6];

  list_edge = (int *)malloc(max_nbe*2*sizeof(int));
  recv_list_edge = (int *)malloc(num*max_nbe*2*sizeof(int));

  /* indicate if a edge is overlapped */
  edge_ov_label = (int *)malloc(num*max_nbe*sizeof(int));
  for (i = 0; i < num*max_nbe; i++)
    edge_ov_label[i] = 0;
  id_on_ov_list = (int *)malloc(max_nbe*sizeof(int));
  /* pack edge information to send  */
  index_nbe = 0;
  for(i = 0; i < nbe; i++) {
    me = List_Entry(boundary_edges,i);
    list_edge[index_nbe] = MV_GlobalID(ME_Vertex(me,0));
    list_edge[index_nbe+1] = MV_GlobalID(ME_Vertex(me,1));
    index_nbe += 2;
  }

  MPI_Allgather(list_edge,2*max_nbe,MPI_INT,recv_list_edge,2*max_nbe,MPI_INT,comm);

  nge = 0;
  /* for processor other than 0 */
  if(rank > 0) {
    for(i = 0; i < nbe; i++) {
      me = List_Entry(boundary_edges,i);
      if(ME_GlobalID(me) > 0) continue;  /* if already assigned */
      edge_id[0] = MV_GlobalID(ME_Vertex(me,0));
      edge_id[1] = MV_GlobalID(ME_Vertex(me,1));
      for(j = 0; j < rank; j++) {
	loc = (int *)bsearch(&edge_id,
			     &recv_list_edge[2*max_nbe*j],
			     global_mesh_info[10*j+6],
			     2*sizeof(int),
			     compareEdgeINT);
	if(loc) {
	  iloc = (int)(loc - &recv_list_edge[2*max_nbe*j])/2;
	  ME_Set_PType(me,PGHOST);
	  ME_Set_MasterParID(me,j);
	  edge_ov_label[max_nbe*j+iloc] |= 1;
	  id_on_ov_list[i] = iloc;
	  nge++;
	  break;
	}
      }
    }
  }

 
  /* num of ghost verts */
  mesh_info[9] = nge;
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);
  /* since this is a OR reduction, we can use MPI_IN_PLACE, send buffer same as recv buffer */
  MPI_Allreduce(MPI_IN_PLACE,edge_ov_label,num*max_nbe,MPI_INT,MPI_LOR,comm);    

  /* Assign global ID for non ghost edge */
  global_id = 1;
  for(i = 0; i < rank; i++) 
    global_id = global_id + global_mesh_info[10*i+2] - global_mesh_info[10*i+9];
  for(i = 0; i < ne; i++) {
    me = MESH_Edge(submesh,i);
    if (ME_PType(me) == PGHOST) continue;
    if (!ME_GlobalID(me)) ME_Set_GlobalID(me,global_id++);
    ME_Set_MasterParID(me,rank);
  }

  /* label OVERLAP edge */
  noe = 0;
  for(i = 0; i < nbe; i++) 
    if(edge_ov_label[rank*max_nbe+i]) {
      me = List_Entry(boundary_edges,i);
      ME_Set_PType(me,POVERLAP);
      noe++;
  }

  /* this time only global id are sent */
  for(i = 0; i < nbe; i++) {
    me = List_Entry(boundary_edges,i);
    list_edge[i] = ME_GlobalID(me);
  }

  MPI_Allgather(list_edge,max_nbe,MPI_INT,recv_list_edge,max_nbe,MPI_INT,comm);

  for(i = 0; i < nbe; i++) {
    me = List_Entry(boundary_edges,i);
    if(ME_PType(me)==PGHOST) {
      int gid = recv_list_edge[ME_MasterParID(me)*max_nbe+id_on_ov_list[i]];
#ifdef DEBUG
      if (ME_GlobalID(me) && ME_GlobalID(me) != gid)
	MSTK_Report("MESH_Assign_GlobalIDs_Edge",
		    "Ghost edge already has different global ID",
		    MSTK_WARN);
#endif
	ME_Set_GlobalID(me, gid);
    }
  }



  List_Delete(boundary_edges);
  free(global_mesh_info);
  free(edge_ov_label);
  free(id_on_ov_list);
  free(list_edge);
  free(recv_list_edge);
  return 1;
}

 /* 
    Assign face global ID for both 2D meshes
    Assume there are no overlapped faces
    Assume each processor knowns its overlap and ghost vertices 
 */

  int MESH_AssignGlobalIDs_Face(Mesh_ptr submesh, MSTK_Comm comm) {
  int i, nf, global_id, mesh_info[10];
  MFace_ptr mf;
  RepType rtype;
  int *global_mesh_info;

  int rank, num;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&num);

  for (i = 0; i < 10; i++) mesh_info[i] = 0;
  rtype = MESH_RepType(submesh);
  nf = MESH_Num_Faces(submesh);

  mesh_info[0] = rtype;
  mesh_info[3] = nf;

  global_mesh_info = (int *)malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  /* calculate starting global id number for faces*/
  global_id = 1;
  for(i = 0; i < rank; i++) 
    global_id = global_id + global_mesh_info[10*i+3];
  for(i = 0; i < nf; i++) {
    mf = MESH_Face(submesh,i);
    MF_Set_PType(mf,PINTERIOR);
    if (!MF_GlobalID(mf)) MF_Set_GlobalID(mf,global_id++);
    MF_Set_MasterParID(mf,rank);
  }

  free(global_mesh_info);
  return 1;
}

 /* 
    Assign region global ID for both 3D meshes
    Main work is to assign global ID for faces

    Assume there are no overlapped regions
    Assume each processor knowns its overlap and ghost vertices 
 */
     
  int MESH_AssignGlobalIDs_Region(Mesh_ptr submesh, MSTK_Comm comm) {
  int i, j, k, nfv, nbf, nof, ngf, nf, nr, mesh_info[10], global_id;
  MVertex_ptr mv;
  MFace_ptr mf;
  MRegion_ptr mr;
  List_ptr boundary_faces, mfverts;
  int *loc, face_id[MAXPV2+3],index_nbf, max_nbf, iloc, is_boundary;
  int *global_mesh_info, *list_face, *recv_list_face, *face_ov_label, *id_on_ov_list;

  int rank, num;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&num);

  for (i = 0; i < 10; i++) mesh_info[i] = 0;
  nf = MESH_Num_Faces(submesh);
  nr = MESH_Num_Regions(submesh);
  mesh_info[3] = nf;
  mesh_info[4] = nr;
  /* 
     collect 'boundary' faces
     if endpoints are either GHOST or OVERLAP, then it is a boundary face 
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
      MF_Flag_OnParBoundary(mf);
      List_Add(boundary_faces,mf);
      nbf++;
    }
    List_Delete(mfverts);
  }
  /* printf("num of boundary faces %d, on rank %d\n", nbf,rank); */
  mesh_info[6] = nbf;

  List_Sort(boundary_faces,nbf,sizeof(MFace_ptr),compareFaceID);

  global_mesh_info = (int *)malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  max_nbf = 0;
  for(i = 0; i < num; i++)
    if(max_nbf < global_mesh_info[10*i+6])
      max_nbf = global_mesh_info[10*i+6];

  list_face = (int *)malloc(max_nbf*(MAXPV2+1)*sizeof(int));
  recv_list_face = (int *)malloc(num*max_nbf*(MAXPV2+1)*sizeof(int));

  /* indicate if a face is overlapped */
  face_ov_label = (int *)malloc(num*max_nbf*sizeof(int));
  for (i = 0; i < num*max_nbf; i++)
    face_ov_label[i] = 0;
  id_on_ov_list = (int *)malloc(max_nbf*sizeof(int));
  /* pack face information to send  */
  index_nbf = 0;
  for(i = 0; i < nbf; i++) {
    mf = List_Entry(boundary_faces,i);
    mfverts = MF_Vertices(mf,1,0);
    nfv = List_Num_Entries(mfverts);
    list_face[index_nbf] = nfv;
    for(j = 0; j < nfv; j++) 
      list_face[index_nbf+j+1] = MV_GlobalID(List_Entry(mfverts,j));
    index_nbf += MAXPV2+1;
    List_Delete(mfverts);
  }

  MPI_Allgather(list_face,(MAXPV2+1)*max_nbf,MPI_INT,recv_list_face,(MAXPV2+1)*max_nbf,MPI_INT,comm);

  ngf = 0;
  /* for processor other than 0 */
  if(rank > 0) {
    for(i = 0; i < nbf; i++) {
      mf = List_Entry(boundary_faces,i);
      if(MF_GlobalID(mf) > 0) continue;  /* if already assigned */
      mfverts = MF_Vertices(mf,1,0);
      nfv = List_Num_Entries(mfverts);
      face_id[0] = nfv;                           
      for(k = 0; k < nfv; k++) 
	face_id[k+1] = MV_GlobalID(List_Entry(mfverts,k));
      List_Delete(mfverts);
      for(j = 0; j < rank; j++) {
	loc = (int *)bsearch(&face_id,
			     &recv_list_face[(MAXPV2+1)*max_nbf*j],
			     global_mesh_info[10*j+6],
			     (MAXPV2+1)*sizeof(int),
			     compareFaceINT);
	if(loc) {
	  iloc = (int)(loc - &recv_list_face[(MAXPV2+1)*max_nbf*j])/(MAXPV2+1);
	  MF_Set_PType(mf,PGHOST);
	  MF_Set_MasterParID(mf,j);
	  face_ov_label[max_nbf*j+iloc] |= 1;
	  id_on_ov_list[i] = iloc;
	  ngf++;
	  break;
	}
      }
    }
  }

 
  /* num of ghost verts */
  mesh_info[9] = ngf;
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);
  /* since this is a OR reduction, we can use MPI_IN_PLACE, send buffer same as recv buffer */
  MPI_Allreduce(MPI_IN_PLACE,face_ov_label,num*max_nbf,MPI_INT,MPI_LOR,comm);    

  /* Assign global ID for non ghost face */
  global_id = 1;
  for(i = 0; i < rank; i++) 
    global_id = global_id + global_mesh_info[10*i+3] - global_mesh_info[10*i+9];
  for(i = 0; i < nf; i++) {
    mf = MESH_Face(submesh,i);
    if (MF_PType(mf) == PGHOST) continue;
    if (!MF_GlobalID(mf)) MF_Set_GlobalID(mf,global_id++);
    MF_Set_MasterParID(mf,rank);
  }

  /* label OVERLAP face */
  nof = 0;
  for(i = 0; i < nbf; i++) 
    if(face_ov_label[rank*max_nbf+i]) {
      mf = List_Entry(boundary_faces,i);
      MF_Set_PType(mf,POVERLAP);
      nof++;
  }
  /* printf("num of ghost faces %d, overlap faces %d on rank %d\n", ngf, nof, rank);  */

  /* this time only global id are sent */
  for(i = 0; i < nbf; i++) {
    mf = List_Entry(boundary_faces,i);
    list_face[i] = MF_GlobalID(mf);
  }

  MPI_Allgather(list_face,max_nbf,MPI_INT,recv_list_face,max_nbf,MPI_INT,comm);

  for(i = 0; i < nbf; i++) {
    mf = List_Entry(boundary_faces,i);
    if(MF_PType(mf)==PGHOST) {
      int gid = recv_list_face[MF_MasterParID(mf)*max_nbf+id_on_ov_list[i]];
#ifdef DEBUG
      if (MF_GlobalID(mf) && MF_GlobalID(mf) != gid)
	MSTK_Report("MESH_AssignGlobalIDs_region",
		    "Ghost face already has different global ID", MSTK_WARN);
#endif
      MF_Set_GlobalID(mf, gid);
    }
  }

  /* assign region global id */
  global_id = 1;
  for(i = 0; i < rank; i++) 
    global_id = global_id + global_mesh_info[10*i+4];
  for(i = 0; i < nr; i++) {
    mr = MESH_Region(submesh,i);
    MR_Set_PType(mr,PINTERIOR);
    if (!MR_GlobalID(mr)) MR_Set_GlobalID(mr,global_id++);
    MR_Set_MasterParID(mr,rank);
  }


  List_Delete(boundary_faces);
  free(global_mesh_info);
  free(face_ov_label);
  free(id_on_ov_list);
  free(list_face);
  free(recv_list_face);


  return 1;
}


#ifdef __cplusplus
}
#endif

