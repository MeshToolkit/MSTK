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
     It label ghost and overlap entities based on global ID
     It also assign master partition ID

     Assume PBOUNDARY label is already assigned for processor boundary entities
     Author(s): Duo Wang, Rao Garimella
  */


  int MESH_Vol_LabelPType_FN(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);
  int MESH_Surf_LabelPType_FN(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);
  int MESH_Vol_LabelPType_R4(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);
  int MESH_Vol_LabelPType_R1R2(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);
  int MESH_Surf_LabelPType_R1R2R4(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);

  static int (*MESH_Vol_LabelPType_jmp[MSTK_MAXREP])(Mesh_ptr submesh, int rank, int num, 
						   MPI_Comm comm) = 
  {MESH_Vol_LabelPType_FN, MESH_Vol_LabelPType_FN, MESH_Vol_LabelPType_R1R2, 
   MESH_Vol_LabelPType_R1R2, MESH_Vol_LabelPType_R4};
  static int (*MESH_Surf_LabelPType_jmp[MSTK_MAXREP])(Mesh_ptr submesh, int rank, int num, 
						    MPI_Comm comm) = 
  {MESH_Surf_LabelPType_FN, MESH_Surf_LabelPType_FN, MESH_Surf_LabelPType_R1R2R4, 
   MESH_Surf_LabelPType_R1R2R4, MESH_Surf_LabelPType_R1R2R4};


int MESH_LabelPType(Mesh_ptr submesh, int rank, int num,  MPI_Comm comm) {
  int nf, nr;
  RepType rtype;

  /* basic mesh information */
  rtype = MESH_RepType(submesh);
  nf = MESH_Num_Faces(submesh);
  nr = MESH_Num_Regions(submesh);
  /* build geometric entity dimension, mark boundary vertices */
  MESH_BuildVertexClassfn(submesh);
  if (nr)
    (*MESH_Vol_LabelPType_jmp[rtype])(submesh, rank, num, comm);
  else if(nf) 
    (*MESH_Surf_LabelPType_jmp[rtype])(submesh, rank, num, comm);
  else {
    MSTK_Report("MESH_LabelPType()","only send volume or surface mesh",MSTK_ERROR);
    exit(-1);
  }
  return 1;
}

int MESH_Surf_LabelPType_FN(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  int i, j, nv, nbv, ne, nf, mesh_info[10], nevs, nfes, nfv, natt, nset, ncomp, dir;
  int nfe;
  MVertex_ptr mv;
  MEdge_ptr me;
  MFace_ptr mf;
  List_ptr boundary_verts, fverts, mfedges;
  RepType rtype;
  char attname[256], msetname[256];
  MType mtype;
  MAttType att_type;
  MAttrib_ptr attrib;
  MSet_ptr mset;
  int *list_attr, *list_mset;
  char *list_attr_names, *list_mset_names;
  int *loc;
  int iloc, num_ghost_verts, global_id;
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
    if (MV_PType(mv) == PBOUNDARY) {
      List_Add(boundary_verts,mv);
      nbv++;
    }
  }
  mesh_info[4] = nbv;
  
  /* sort boundary vertices based on global ID, for binary search */
  List_Sort(boundary_verts,nbv,sizeof(MVertex_ptr),compareGlobalID);
  for(i = 0; i < nbv; i++) {
    mv = List_Entry(boundary_verts,i);
  }


  /* 
     gather submeshes information
     right now we only need nv and nbv, and later num_ghost_verts, but we gather all mesh_info
  */
  int *global_mesh_info = (int *)MSTK_malloc(10*num*sizeof(int));
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);

  /* get largest number of boundary vertices of all the processors */
  int max_nbv = 0;
  for(i = 0; i < num; i++)
    if(max_nbv < global_mesh_info[10*i+4])
      max_nbv = global_mesh_info[10*i+4];

  int *list_boundary_vertex = (int *)MSTK_malloc(max_nbv*sizeof(int));
  int *recv_list_vertex = (int *)MSTK_malloc(num*max_nbv*sizeof(int));

  /* only global ID are sent */
  int index_nbv = 0;
  for(i = 0; i < nbv; i++) {
    mv = List_Entry(boundary_verts,i);
    list_boundary_vertex[index_nbv] = MV_GlobalID(mv);
    index_nbv++;
  }

  /* gather boundary vertices */
  MPI_Allgather(list_boundary_vertex,max_nbv,MPI_INT,recv_list_vertex,max_nbv,MPI_INT,comm);

  /* indicate if a vertex is overlapped */
  int *vertex_ov_label = (int *)MSTK_malloc(num*max_nbv*sizeof(int));

  /* store the local boundary id on ov processor
     it is used to assign global id of local ghost vertices
     do not need to store master partition id, MV_MasterParID(mv) is assigned
  */
  for (i = 0; i < num*max_nbv; i++)
    vertex_ov_label[i] = 0;
  num_ghost_verts = 0;
  /* for processor other than 0 */
  for(i = 0; i < nbv; i++) {
    mv = List_Entry(boundary_verts,i);
    global_id = MV_GlobalID(mv);
    /* check which processor has the same coordinate vertex */
    for(j = 0; j < num; j++) {
      if(j == rank)
	continue;
      /* since each processor has sorted the boundary vertices, use binary search */
      loc = (int *)bsearch(&global_id,
			   &recv_list_vertex[max_nbv*j],
			   global_mesh_info[10*j+4],
			   sizeof(int),
			   compareGlobalID);
      /* if found the vertex on previous processors */
      if(loc) {
	/* here the location iloc is relative to the beginning of the jth processor */
	iloc = (int)(loc - &recv_list_vertex[max_nbv*j]);
	MV_Set_PType(mv,PGHOST);
	MV_Set_MasterParID(mv,j);
	printf("rank %d,boundary vertex local ID %d, found on processor %d, loc %d\n",rank, MV_ID(mv), j, iloc);
	num_ghost_verts++;
	/* label the original vertex as overlapped */
	vertex_ov_label[max_nbv*j+iloc] |= 1;
	/* if found on processor j, no need to test other processors*/
	break;
      }
    }
  }
  /* num of ghost verts */
  mesh_info[5] = num_ghost_verts;
  printf("befor reduce rank %d:",rank);
  for (i = 0; i < num*max_nbv; i++)
    printf("%d",vertex_ov_label[i]);
  printf("\n");

  /* update ghost verts number */
  MPI_Allgather(mesh_info,10,MPI_INT,global_mesh_info,10,MPI_INT,comm);
  /* since this is a or reduction, we can use MPI_IN_PLACE, send buffer same as recv buffer */
  MPI_Allreduce(MPI_IN_PLACE,vertex_ov_label,num*max_nbv,MPI_INT,MPI_LOR,comm);    
  printf("after reduce rank %d:",rank);
  for (i = 0; i < num*max_nbv; i++) {
    printf("%d",vertex_ov_label[i]);
  }
  printf("\n");


  /* calculate starting global id number for vertices*/
  for(i = 0; i < nv; i++) {
    mv = MESH_Vertex(submesh,i);
    if (MV_PType(mv) == PGHOST)
      continue;
    MV_Set_MasterParID(mv,rank);
  }

  for(i = 0; i < nbv; i++) {
    if(vertex_ov_label[rank*max_nbv+i]) {
      mv = List_Entry(boundary_verts,i);
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


int MESH_Vol_LabelPType_FN(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  int i, j, nv, ne, nf, nr, mesh_info[10];
  int nevs, nfes, nrfs, nfe, nrv, nrf, natt, nset, ncomp, dir;
  MVertex_ptr mv;
  MEdge_ptr me;
  MFace_ptr mf;
  MRegion_ptr mr;
  List_ptr mfedges, mrfaces, mrverts;
  RepType rtype;
  char attname[256], msetname[256];
  MType mtype;
  MAttType att_type;
  MAttrib_ptr attrib;
  MSet_ptr mset;
  int *list_attr, *list_mset;
  char *list_attr_names, *list_mset_names;
  double coor[3];

  for (i = 0; i < 10; i++) mesh_info[i] = 0;

  /* mesh_info store the mesh reptype, nv, nf, nfvs and natt */
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
  /* collect data */
  int *list_vertex = (int *)MSTK_malloc(3*nv*sizeof(int));
  double *list_coor = (double *)MSTK_malloc(3*nv*sizeof(double));
  for(i = 0; i < nv; i++) {
    mv = MESH_Vertex(submesh,i);
    list_vertex[3*i] = (MV_GEntID(mv)<<3) | (MV_GEntDim(mv));
    list_vertex[3*i+1] = (MV_MasterParID(mv) <<2) | (MV_PType(mv));
    list_vertex[3*i+2] = MV_GlobalID(mv);
    MV_Coords(mv,coor);
    list_coor[i*3] = coor[0];
    list_coor[i*3+1] = coor[1];
    list_coor[i*3+2] = coor[2];
  }


 int *list_edge = (int *)MSTK_malloc(5*ne*sizeof(int));

  nevs = 0;
  /* Store the vertex ids, then the 3 auxilliary data fields */
  for(i = 0; i < ne; i++) {
    me = MESH_Edge(submesh,i);
    list_edge[nevs]   = MV_ID(ME_Vertex(me,0));
    list_edge[nevs+1] = MV_ID(ME_Vertex(me,1));
    list_edge[nevs+2] = (ME_GEntID(me)<<3) | (ME_GEntDim(me));
    list_edge[nevs+3] = (ME_MasterParID(me) <<2) | (ME_PType(me));
    list_edge[nevs+4] = ME_GlobalID(me);
    nevs += 5;
  }
  mesh_info[5] = nevs;


  int *list_face = (int *)MSTK_malloc((MAXPV2+4)*nf*sizeof(int));

  nfes = 0;
  /* first int store nfe, then the edge ids, then the 3 auxilliary data fields */
  for(i = 0; i < nf; i++) {
    mf = MESH_Face(submesh,i);
    mfedges = MF_Edges(mf,1,0);
    nfe = List_Num_Entries(mfedges);
    list_face[nfes] = nfe;
    for(j = 0; j < nfe; j++) {
      dir = MF_EdgeDir_i(mf,j) ? 1 : -1;
      list_face[nfes+j+1] = dir*ME_ID(List_Entry(mfedges,j));
    }
    list_face[nfes+nfe+1] = (MF_GEntID(mf)<<3) | (MF_GEntDim(mf));
    list_face[nfes+nfe+2] = (MF_MasterParID(mf) <<2) | (MF_PType(mf));
    list_face[nfes+nfe+3] = MF_GlobalID(mf);
    nfes += (nfe + 4);
    List_Delete(mfedges);
  }
  mesh_info[6] = nfes;



  int *list_region = (int *)MSTK_malloc((MAXPF3+4)*nr*sizeof(int));

  nrfs = 0;
  /* first store nrf, then the face ids, then the 3 auxilliary data fields */
  for(i = 0; i < nr; i++) {
    mr = MESH_Region(submesh,i);
    mrfaces = MR_Faces(mr);
    nrf = List_Num_Entries(mrfaces);
    list_region[nrfs] = nrf;
    for(j = 0; j < nrf; j++) {
      dir = MR_FaceDir_i(mr,j) == 1 ? 1 : -1;
      list_region[nrfs+j+1] = dir*MF_ID(List_Entry(mrfaces,j));
    }
    list_region[nrfs+nrf+1] = (MR_GEntID(mr)<<3) | (MR_GEntDim(mr));
    list_region[nrfs+nrf+2] = (MR_MasterParID(mr) <<2) | (MR_PType(mr));
    list_region[nrfs+nrf+3] = MR_GlobalID(mr);
    nrfs += (nrf + 4);
    List_Delete(mrfaces);
  }
  mesh_info[7] = nrfs;

  /* number of attr */
  natt = MESH_Num_Attribs(submesh);
  mesh_info[8] = natt;

  /* collect attrs */
  if(natt) {
    list_attr = (int *)MSTK_malloc((natt+1)*sizeof(int));
    list_attr_names = (char *)MSTK_malloc((natt+1)*256);
    for(i = 0; i < natt; i++) {
      attrib = MESH_Attrib(submesh,i);
      MAttrib_Get_Name(attrib,attname);
      att_type = MAttrib_Get_Type(attrib);
      ncomp = MAttrib_Get_NumComps(attrib);
      mtype = MAttrib_Get_EntDim(attrib);
      list_attr[i] = (ncomp << 6) | (mtype << 3) | (att_type);
      strcpy(&list_attr_names[i*256],attname);
    }
  }


  /* Mesh entity sets */

  nset = MESH_Num_MSets(submesh);
  mesh_info[9] = nset;

  if (nset) {
    list_mset = (int *) MSTK_malloc(nset*sizeof(int));
    list_mset_names = (char *) MSTK_malloc(nset*256*sizeof(char));

    for (i = 0; i < nset; i++) {
      mset = MESH_MSet(submesh,i);
      MSet_Name(mset,msetname);
      mtype = MSet_EntDim(mset);
      list_mset[i] = mtype;
      strcpy(&list_mset_names[i*256],msetname);
    }
  }



  /* send mesh_info */
  MPI_Send(mesh_info,10,MPI_INT,rank,rank,comm);

  /* send vertices */
  /* printf("%d vertices sent to rank %d\n",nv,rank); */
  MPI_Send(list_vertex,3*nv,MPI_INT,rank,rank,comm);
  MPI_Send(list_coor,3*nv,MPI_DOUBLE,rank,rank,comm);

  /* send edges */
  MPI_Send(list_edge,nevs,MPI_INT,rank,rank,comm);
  
  /* send faces */
  MPI_Send(list_face,nfes,MPI_INT,rank,rank,comm);

  /* send regions */
  /* printf("%d regions sent to rank %d\n",nr,rank); */
  MPI_Send(list_region,nrfs,MPI_INT,rank,rank,comm);

  /* send attr */
  /* printf("%d attr sent to rank %d\n",natt,rank); */
  if(natt) {
    MPI_Send(list_attr,natt,MPI_INT,rank,rank,comm);
    MPI_Send(list_attr_names,natt*256,MPI_CHAR,rank,rank,comm);
    MSTK_free(list_attr);
    MSTK_free(list_attr_names);
  }


  /* send sets */
  if (nset) {
    MPI_Send(list_mset,nset,MPI_INT,rank,rank,comm);
    MPI_Send(list_mset_names,nset*256,MPI_CHAR,rank,rank,comm);
    MSTK_free(list_mset);
    MSTK_free(list_mset_names);
  }



  MSTK_free(list_vertex);
  MSTK_free(list_region);
  MSTK_free(list_coor);
  return 1;
}




int MESH_Surf_LabelPType_R1R2R4(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  MSTK_Report("MESH_Surf_LabelPType_R1R2R4","Not implemented",MSTK_FATAL);
}


int MESH_Vol_LabelPType_R1R2(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  MSTK_Report("MESH_Vol_LabelPType_R4","Not implemented",MSTK_FATAL);
}


int MESH_Vol_LabelPType_R4(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  MSTK_Report("MESH_Vol_LabelPType_R4","Not implemented",MSTK_FATAL);
}





  
#ifdef __cplusplus
}
#endif

