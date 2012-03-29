#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif



  /* 
     This function send part of mesh to processor rank in communicator comm
     The list of elements to be sent is specified by list_send_elements
     
     attr list is sent but no attribute value is sent.
     call MESH_SendAttr() to send attribute values of elements

     Author(s): Duo Wang, Rao Garimella
  */


  int MESH_Vol_SendMesh_Elements_FN(Mesh_ptr mesh, List_ptr list_send_elements, int rank, MPI_Comm comm);
  int MESH_Surf_SendMesh_Elements_FN(Mesh_ptr mesh, List_ptr list_send_elements, int rank, MPI_Comm comm);
  int MESH_Vol_SendMesh_Elements_R4(Mesh_ptr mesh, List_ptr list_send_elements, int rank, MPI_Comm comm);
  int MESH_Vol_SendMesh_Elements_R1R2(Mesh_ptr mesh, List_ptr list_send_elements, int rank, MPI_Comm comm);
  int MESH_Surf_SendMesh_Elements_R1R2R4(Mesh_ptr mesh, List_ptr list_send_elements, int rank, MPI_Comm comm);

  static int (*MESH_Vol_SendMesh_Elements_jmp[MSTK_MAXREP])(Mesh_ptr mesh, List_ptr list_send_elements, int rank, 
						   MPI_Comm comm) = 
  {MESH_Vol_SendMesh_Elements_FN, MESH_Vol_SendMesh_Elements_FN, MESH_Vol_SendMesh_Elements_R1R2, 
   MESH_Vol_SendMesh_Elements_R1R2, MESH_Vol_SendMesh_Elements_R4};
  static int (*MESH_Surf_SendMesh_Elements_jmp[MSTK_MAXREP])(Mesh_ptr mesh, List_ptr list_send_elements, int rank, 
						    MPI_Comm comm) = 
  {MESH_Surf_SendMesh_Elements_FN, MESH_Surf_SendMesh_Elements_FN, MESH_Surf_SendMesh_Elements_R1R2R4, 
   MESH_Surf_SendMesh_Elements_R1R2R4, MESH_Surf_SendMesh_Elements_R1R2R4};


int MESH_SendMesh_Elements(Mesh_ptr mesh, List_ptr list_send_elements, int rank, MPI_Comm comm) {
  int nf, nr;
  RepType rtype;

  /* basic mesh information */
  rtype = MESH_RepType(mesh);
  nf = MESH_Num_Faces(mesh);
  nr = MESH_Num_Regions(mesh);
  if (nr)
    (*MESH_Vol_SendMesh_Elements_jmp[rtype])(mesh,list_send_elements,rank,comm);
  else if(nf) 
    (*MESH_Surf_SendMesh_Elements_jmp[rtype])(mesh,list_send_elements,rank,comm);
  else {
    MSTK_Report("MESH_SendMesh_Elements()","only send volume or surface mesh",MSTK_ERROR);
    exit(-1);
  }
  return 1;
}




  /* list_send_elements is a list of faces */


int MESH_Surf_SendMesh_Elements_FN(Mesh_ptr mesh, List_ptr list_send_elements, int rank, MPI_Comm comm) {
  int i, j, nv, nf, mesh_info[10], nfv,  nset;
  int nfvs, mkvid;
  MVertex_ptr mv;
  MFace_ptr mf;
  List_ptr mfverts, list_send_verts;
  RepType rtype;
  char msetname[256];
  MType mtype;
  MSet_ptr mset;
  int  *list_mset, *list_face, *list_vertex, *MV_to_list_ID;
  char  *list_mset_names;
  double coor[3], *list_coor;

  for (i = 0; i < 10; i++) mesh_info[i] = 0;
  nv = 0;
  nf = List_Num_Entries(list_send_elements);            /* number of faces to be sent */
  if (nf == 0)                                    
    return 1;


  mkvid = MSTK_GetMarker();

  list_face = (int *)MSTK_malloc((MAXPV2+4)*nf*sizeof(int));

  list_vertex = (int *)MSTK_malloc(3*nf*MAXPV2*sizeof(int));
  list_coor = (double *)MSTK_malloc(3*nf*MAXPV2*sizeof(double));

  /* stores the map from MV_ID to ID on list_send_verts */
  MV_to_list_ID = (int *)MSTK_malloc(MESH_Num_Vertices(mesh)*sizeof(int));


  /* 
     In the following code, we encode the face and vertex information at the same time
     The neighboring vertex IDs of a face are IDs on the list_send_verts, not MV_ID
  */

  list_send_verts = List_New(10);
  nfvs = 0;
  for (i = 0; i < nf; i++) {
    mf = List_Entry(list_send_elements,i);
    mfverts = MF_Vertices(mf,1,0);        /* add vertices */
    nfv = List_Num_Entries(mfverts);
    list_face[nfvs] = nfv;
    for (j = 0; j < nfv; j++) {
      mv = List_Entry(mfverts,j);
      if (!MEnt_IsMarked(mv,mkvid)) {
	List_Add(list_send_verts,mv);
	MEnt_Mark(mv,mkvid);

	list_vertex[3*nv] = (MV_GEntID(mv)<<3) | (MV_GEntDim(mv));  /* encode the vertex information */
	list_vertex[3*nv+1] = (MV_MasterParID(mv) <<2) | (MV_PType(mv));
	list_vertex[3*nv+2] = MV_GlobalID(mv);
	MV_Coords(mv,coor);
	list_coor[3*nv] = coor[0];
	list_coor[3*nv+1] = coor[1];
	list_coor[3*nv+2] = coor[2];

	MV_to_list_ID[MV_ID(mv)-1] = nv++;                         /* the receiver needs to know the vertex id on recv list */
      }
      list_face[nfvs+j+1] = MV_to_list_ID[MV_ID(mv)-1];            /* stores the vertex id on the lsit */
    }
    list_face[nfvs+nfv+1] = (MF_GEntID(mf)<<3) | (MF_GEntDim(mf));    /* encode the face information */
    list_face[nfvs+nfv+2] = (MF_MasterParID(mf)<<2) | (MF_PType(mf));
    list_face[nfvs+nfv+3] = MF_GlobalID(mf);
    nfvs += (nfv + 4);

    List_Delete(mfverts);    
  }
  List_Unmark(list_send_verts, mkvid);
  MSTK_FreeMarker(mkvid); 


  /* mesh_info store the mesh reptype, nv, nf, nfvs and natt */
  rtype = MESH_RepType(mesh);
  mesh_info[0] = rtype;
  mesh_info[1] = nv;

  mesh_info[3] = nf;

  mesh_info[5] = nfvs;

  /* Mesh entity sets */

  nset = MESH_Num_MSets(mesh);
  mesh_info[7] = nset;

  if (nset) {
    list_mset = (int *) MSTK_malloc(nset*sizeof(int));
    list_mset_names = (char *) MSTK_malloc(nset*256*sizeof(char));

    for (i = 0; i < nset; i++) {
      mset = MESH_MSet(mesh,i);
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

  /* send faces */
  /* printf("%d faces sent to rank %d\n",nf,rank); */
  MPI_Send(list_face,nfvs,MPI_INT,rank,rank,comm);
  
  /* send sets */
  if (nset) {
    MPI_Send(list_mset,nset,MPI_INT,rank,rank,comm);
    MPI_Send(list_mset_names,nset*256,MPI_CHAR,rank,rank,comm);
    MSTK_free(list_mset);
    MSTK_free(list_mset_names);
  }


  List_Delete(list_send_verts);

  MSTK_free(MV_to_list_ID);
  MSTK_free(list_vertex);
  MSTK_free(list_coor);
  MSTK_free(list_face);
  return 1;
}


int MESH_Vol_SendMesh_Elements_FN(Mesh_ptr mesh, List_ptr list_send_elements, int rank, MPI_Comm comm) {
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
  rtype = MESH_RepType(mesh);
  nv = MESH_Num_Vertices(mesh);
  ne = MESH_Num_Edges(mesh);
  nf = MESH_Num_Faces(mesh);
  nr = MESH_Num_Regions(mesh);
  mesh_info[0] = rtype;
  mesh_info[1] = nv;
  mesh_info[2] = ne;
  mesh_info[3] = nf;
  mesh_info[4] = nr;
  /* collect data */
  int *list_vertex = (int *)MSTK_malloc(3*nv*sizeof(int));
  double *list_coor = (double *)MSTK_malloc(3*nv*sizeof(double));
  for(i = 0; i < nv; i++) {
    mv = MESH_Vertex(mesh,i);
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
    me = MESH_Edge(mesh,i);
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
    mf = MESH_Face(mesh,i);
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
    mr = MESH_Region(mesh,i);
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
  natt = MESH_Num_Attribs(mesh);
  mesh_info[8] = natt;

  /* collect attrs */
  if(natt) {
    list_attr = (int *)MSTK_malloc((natt+1)*sizeof(int));
    list_attr_names = (char *)MSTK_malloc((natt+1)*256);
    for(i = 0; i < natt; i++) {
      attrib = MESH_Attrib(mesh,i);
      MAttrib_Get_Name(attrib,attname);
      att_type = MAttrib_Get_Type(attrib);
      ncomp = MAttrib_Get_NumComps(attrib);
      mtype = MAttrib_Get_EntDim(attrib);
      list_attr[i] = (ncomp << 6) | (mtype << 3) | (att_type);
      strcpy(&list_attr_names[i*256],attname);
    }
  }


  /* Mesh entity sets */

  nset = MESH_Num_MSets(mesh);
  mesh_info[9] = nset;

  if (nset) {
    list_mset = (int *) MSTK_malloc(nset*sizeof(int));
    list_mset_names = (char *) MSTK_malloc(nset*256*sizeof(char));

    for (i = 0; i < nset; i++) {
      mset = MESH_MSet(mesh,i);
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




int MESH_Surf_SendMesh_Elements_R1R2R4(Mesh_ptr mesh, List_ptr list_send_elements, int rank, MPI_Comm comm) {
  MSTK_Report("MESH_Surf_SendMesh_Elements_R1R2R4","Not implemented",MSTK_FATAL);
}


int MESH_Vol_SendMesh_Elements_R1R2(Mesh_ptr mesh, List_ptr list_send_elements, int rank, MPI_Comm comm) {
  MSTK_Report("MESH_Vol_SendMesh_Elements_R4","Not implemented",MSTK_FATAL);
}


int MESH_Vol_SendMesh_Elements_R4(Mesh_ptr mesh, List_ptr list_send_elements, int rank, MPI_Comm comm) {
  MSTK_Report("MESH_Vol_SendMesh_Elements_R4","Not implemented",MSTK_FATAL);
}





  
#ifdef __cplusplus
}
#endif

