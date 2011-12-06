#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif

  /* 
     This function receives mesh from processor rank in communicator comm
     attrib list is received, but no attrib values.
     call MESH_RecvAttr() to update entity attribute values
     send_rank: the rank of sending processor
     rank: the rank of receiving processor

     Author(s): Duo Wang, Rao Garimella
  */

  int MESH_Vol_RecvMesh_FN(Mesh_ptr mesh, int *mesh_info, int send_rank, int rank, MPI_Comm comm);
  int MESH_Surf_RecvMesh_FN(Mesh_ptr mesh, int *mesh_info, int send_rank, int rank, MPI_Comm comm);
  int MESH_Vol_RecvMesh_R4(Mesh_ptr mesh, int *mesh_info, int send_rank, int rank, MPI_Comm comm);
  int MESH_Vol_RecvMesh_R1R2(Mesh_ptr mesh, int *mesh_info, int send_rank, int rank, MPI_Comm comm);
  int MESH_Surf_RecvMesh_R1R2R4(Mesh_ptr mesh, int *mesh_info, int send_rank, int rank, MPI_Comm comm);

  static int (*MESH_Vol_RecvMesh_jmp[MSTK_MAXREP])(Mesh_ptr mesh, 
						   int *mesh_info, 
						   int send_rank, int rank, 
						   MPI_Comm comm) = 
  {MESH_Vol_RecvMesh_FN, MESH_Vol_RecvMesh_FN, MESH_Vol_RecvMesh_R1R2, 
   MESH_Vol_RecvMesh_R1R2, MESH_Vol_RecvMesh_R4};
  static int (*MESH_Surf_RecvMesh_jmp[MSTK_MAXREP])(Mesh_ptr mesh, 
						    int *mesh_info, 
						    int send_rank, int rank, 
						    MPI_Comm comm) = 
  {MESH_Surf_RecvMesh_FN, MESH_Surf_RecvMesh_FN, MESH_Surf_RecvMesh_R1R2R4, 
   MESH_Surf_RecvMesh_R1R2R4, MESH_Surf_RecvMesh_R1R2R4};



  int MESH_RecvMesh(Mesh_ptr mesh, int dim, int send_rank, int rank, MPI_Comm comm) {
    int mesh_info[10], count;
    MPI_Status status;
    RepType rtype;
    
    /* mesh_info store the mesh reptype, nv, nf, nfvs and natt */
    /* receive mesh_info */
    MPI_Recv(mesh_info,10,MPI_INT,send_rank,rank,comm,&status);
    MPI_Get_count(&status,MPI_INT,&count);
    
    rtype = mesh_info[0];
    
    if (dim == 3)
      (*MESH_Vol_RecvMesh_jmp[rtype])(mesh,mesh_info,send_rank,rank,comm);
    else if(dim == 2) 
      (*MESH_Surf_RecvMesh_jmp[rtype])(mesh,mesh_info,send_rank,rank,comm);
    else {
      MSTK_Report("MESH_RecvMesh()","only receive volume or surface mesh",MSTK_ERROR);
      exit(-1);
    }
    return 1;
  }
  





  int MESH_Surf_RecvMesh_FN(Mesh_ptr mesh, int *mesh_info, int send_rank, int rank, MPI_Comm comm) {
    int i, j, nevs, nevs_local, nfe, nfes, nfes_local, natt, nset, count, ncomp;
    int nv, ne, nf, *fedirs;
    MVertex_ptr *verts, *fverts;
    MEdge_ptr me, *edges, *fedges;
    MFace_ptr mf;
    MPI_Status status;
    RepType rtype;
    char attname[256], msetname[256];
    MType mtype;
    MAttType att_type;
    MAttrib_ptr attrib;
    MSet_ptr mset;
    int *list_attr, *list_mset;
    char *list_attr_names, *list_mset_names;
    double coor[3];


    /* basic mesh information */
    rtype = mesh_info[0];
    nv = mesh_info[1];
    ne = mesh_info[2];
    nf = mesh_info[3];
    nevs = mesh_info[4];
    nfes = mesh_info[5];
    natt = mesh_info[6];
    nset = mesh_info[7];
  
    MESH_SetRepType(mesh,rtype);

    /* allocate receive buffer */
    int *list_vertex = (int *) malloc(3*nv*sizeof(int));
    double *list_coor = (double *) malloc(3*nv*sizeof(double));

    /* receive vertex */
    MPI_Recv(list_vertex,3*nv,MPI_INT,send_rank,rank,comm,&status);
    MPI_Get_count(&status,MPI_INT,&count);
    /* printf("num of vertex received %d on rank %d\n",nv,rank); */
    MPI_Recv(list_coor,3*nv,MPI_DOUBLE,send_rank,rank,comm,&status);
    verts = (MVertex_ptr *) malloc(nv*sizeof(MVertex_ptr));
    for(i = 0; i < nv; i++) {
      verts[i] = MV_New(mesh);
      MV_Set_GEntDim(verts[i],(list_vertex[3*i] & 7));
      MV_Set_GEntID(verts[i],(list_vertex[3*i] >> 3));
      MV_Set_PType(verts[i],(list_vertex[3*i+1] & 3));
      MV_Set_MasterParID(verts[i],(list_vertex[3*i+1] >> 2));
      MV_Set_GlobalID(verts[i],list_vertex[3*i+2]);
      coor[0] = list_coor[i*3];
      coor[1] = list_coor[i*3+1];
      coor[2] = list_coor[i*3+2];
      MV_Set_Coords(verts[i],coor);
    }


    /* receive edge info */

    int *list_edge = (int *) malloc(nevs*sizeof(int));

    MPI_Recv(list_edge,nevs,MPI_INT,send_rank,rank,comm,&status);
    MPI_Get_count(&status,MPI_INT,&count);

    edges = (MEdge_ptr *) malloc(ne*sizeof(MEdge_ptr));

    nevs_local = 0;
    for (i = 0; i < ne; i++) {
      edges[i] = me = ME_New(mesh);
      ME_Set_GEntDim(me,(list_edge[nevs_local+2] & 7));
      ME_Set_GEntID(me,(list_edge[nevs_local+2] >> 3));
      ME_Set_PType(me,(list_edge[nevs_local+3] & 3));
      ME_Set_MasterParID(me,(list_edge[nevs_local+3] >> 2));
      ME_Set_GlobalID(me,list_edge[nevs_local+4]);

      ME_Set_Vertex(me,0,verts[list_edge[nevs_local]-1]);
      ME_Set_Vertex(me,1,verts[list_edge[nevs_local+1]-1]);
      nevs_local += 5;
    }
    if (verts) free(verts);


    /* receive face info */

    int *list_face = (int *) malloc(nfes*sizeof(int));

    MPI_Recv(list_face,nfes,MPI_INT,send_rank,rank,comm,&status);
    MPI_Get_count(&status,MPI_INT,&count);

    fedges = (MEdge_ptr *) malloc(MAXPV2*sizeof(MEdge_ptr));
    fedirs = (int *) malloc(MAXPV2*sizeof(int));
    nfes_local = 0;
    for (i = 0; i < nf; i++) {
      mf = MF_New(mesh);
      nfe = list_face[nfes_local];
      for (j = 0; j < nfe; j++) {
	fedges[j] = edges[abs(list_face[nfes_local+j+1])-1];
	fedirs[j] = list_face[nfes_local+j+1] > 0 ? 1 : 0;
      }
      MF_Set_GEntDim(mf,(list_face[nfes_local+nfe+1] & 7));
      MF_Set_GEntID(mf,(list_face[nfes_local+nfe+1] >> 3));
      MF_Set_PType(mf,(list_face[nfes_local+nfe+2] & 3));
      MF_Set_MasterParID(mf,(list_face[nfes_local+nfe+2] >> 2));
      MF_Set_GlobalID(mf,list_face[nfes_local+nfe+3]);

      MF_Set_Edges(mf,nfe,fedges,fedirs);
      nfes_local += (nfe + 4);
    }
    if (fedges) free(fedges);
    if (fedirs) free(fedirs);
    if (edges) free(edges);


    /* receive attr */
    if(natt) {
      list_attr = (int *) malloc(natt*sizeof(int));
      list_attr_names = (char *) malloc((natt)*256*sizeof(char));
      MPI_Recv(list_attr,natt,MPI_INT,send_rank,rank,comm,&status);
      MPI_Recv(list_attr_names,natt*256,MPI_CHAR,send_rank,rank,comm,&status);
      MPI_Get_count(&status,MPI_INT,&count);
      for(i = 0; i < natt; i++) {
	strcpy(attname,&list_attr_names[i*256]);
	/* see of the attrib exists */
	if(MESH_AttribByName(mesh,attname))
	  continue;
	att_type = list_attr[i] & 7;
	mtype = (list_attr[i] >> 3) & 7;
	ncomp = list_attr[i] >> 6;
	if(ncomp == 1)
	  attrib =  MAttrib_New(mesh, attname, att_type, mtype);
	else
	  attrib =  MAttrib_New(mesh, attname, att_type, mtype, ncomp);
	/* printf("attr %d with name %s received on rank %d\n",i, attname, rank); */
      }
      MSTK_free(list_attr);    
      MSTK_free(list_attr_names);
    }
    
    /* receive mesh entity sets */
    if(nset) {
      list_mset = (int *) malloc(nset*sizeof(int));
      list_mset_names = (char *) malloc((nset)*256*sizeof(char));
    
      MPI_Recv(list_mset,nset,MPI_INT,send_rank,rank,comm,&status);
      MPI_Recv(list_mset_names,nset*256,MPI_CHAR,send_rank,rank,comm,&status);
      MPI_Get_count(&status,MPI_INT,&count);

      for(i = 0; i < nset; i++) {
	strcpy(msetname,&list_mset_names[i*256]);
	/* see if the mset exists */
	if(MESH_MSetByName(mesh,msetname))
	  continue;
	mtype = list_mset[i];
	mset =  MSet_New(mesh, msetname, mtype);
      }
      MSTK_free(list_mset);    
      MSTK_free(list_mset_names);
    }


    MSTK_free(list_vertex);    
    MSTK_free(list_face);    
    MSTK_free(list_coor);
    return 1;
  }

  int MESH_Vol_RecvMesh_FN(Mesh_ptr mesh, int *mesh_info, int send_rank, int rank, MPI_Comm comm) {
    int i, j, nevs, nevs_local, nfes, nfes_local, nrv, nrfs, nrfs_local;
    int count, natt, nset, ncomp, nv, ne, nf, nr, nfe, nrf, *fedirs, *rfdirs;
    MVertex_ptr *verts, *rverts;
    MEdge_ptr me, *edges, *fedges;
    MFace_ptr mf, *faces, *rfaces;
    MRegion_ptr mr;
    MPI_Status status;
    RepType rtype;
    char attname[256], msetname[256];
    MType mtype;
    MAttType att_type;
    MAttrib_ptr attrib;
    MSet_ptr mset;
    int *list_attr, *list_mset;
    char *list_attr_names, *list_mset_names;
    double coor[3];
   

    /* mesh_info store the mesh reptype, nv, nr, nrvs and natt */
    MPI_Get_count(&status,MPI_INT,&count);

    /* basic mesh information */
    rtype = mesh_info[0];
    nv = mesh_info[1];
    ne = mesh_info[2];
    nf = mesh_info[3];
    nr = mesh_info[4];    
    nevs = mesh_info[5];
    nfes = mesh_info[6];
    nrfs = mesh_info[7];
    natt = mesh_info[8];
    nset = mesh_info[9];

    MESH_SetRepType(mesh,rtype);

    /* allocate receive buffer */
    int *list_vertex = (int *) malloc(3*nv*sizeof(int));
    double *list_coor = (double *) malloc(3*nv*sizeof(double));

    /* receive vertex */
    MPI_Recv(list_vertex,3*nv,MPI_INT,send_rank,rank,comm,&status);
    MPI_Get_count(&status,MPI_INT,&count);

    /* printf("num of vertex received %d on rank %d\n",nv,rank); */

    MPI_Recv(list_coor,3*nv,MPI_DOUBLE,send_rank,rank,comm,&status);
    verts = (MVertex_ptr *) malloc(nv*sizeof(MVertex_ptr));
    for(i = 0; i < nv; i++) {
      verts[i] = MV_New(mesh);
      MV_Set_GEntDim(verts[i],(list_vertex[3*i] & 7));
      MV_Set_GEntID(verts[i],(list_vertex[3*i] >> 3));
      MV_Set_PType(verts[i],(list_vertex[3*i+1] & 3));
      MV_Set_MasterParID(verts[i],(list_vertex[3*i+1] >> 2));
      MV_Set_GlobalID(verts[i],list_vertex[3*i+2]);
      coor[0] = list_coor[i*3];
      coor[1] = list_coor[i*3+1];
      coor[2] = list_coor[i*3+2];
      MV_Set_Coords(verts[i],coor);
    }


    /* receive edge info */

    int *list_edge = (int *) malloc(nevs*sizeof(int));

    MPI_Recv(list_edge,nevs,MPI_INT,send_rank,rank,comm,&status);
    MPI_Get_count(&status,MPI_INT,&count);

    edges = (MEdge_ptr *) malloc(ne*sizeof(MEdge_ptr));

    nevs_local = 0;
    for (i = 0; i < ne; i++) {
      edges[i] = me = ME_New(mesh);
      ME_Set_GEntDim(me,(list_edge[nevs_local+2] & 7));
      ME_Set_GEntID(me,(list_edge[nevs_local+2] >> 3));
      ME_Set_PType(me,(list_edge[nevs_local+3] & 3));
      ME_Set_MasterParID(me,(list_edge[nevs_local+3] >> 2));
      ME_Set_GlobalID(me,list_edge[nevs_local+4]);

      ME_Set_Vertex(me,0,verts[list_edge[nevs_local]-1]);
      ME_Set_Vertex(me,1,verts[list_edge[nevs_local+1]-1]);
      nevs_local += 5;
    }
    if (verts) free(verts);


    /* receive face info */

    int *list_face = (int *) malloc(nfes*sizeof(int));

    MPI_Recv(list_face,nfes,MPI_INT,send_rank,rank,comm,&status);
    MPI_Get_count(&status,MPI_INT,&count);

    fedges = (MEdge_ptr *) malloc(MAXPV2*sizeof(MEdge_ptr));
    fedirs = (int *) malloc(MAXPV2*sizeof(int));

    faces = (MFace_ptr *) malloc(nf*sizeof(MFace_ptr));

    nfes_local = 0;
    for (i = 0; i < nf; i++) {
      faces[i] = mf = MF_New(mesh);
      nfe = list_face[nfes_local];
      for (j = 0; j < nfe; j++) {
	fedges[j] = edges[abs(list_face[nfes_local+j+1])-1];
	fedirs[j] = list_face[nfes_local+j+1] > 0 ? 1 : 0;
      }
      MF_Set_GEntDim(mf,(list_face[nfes_local+nfe+1] & 7));
      MF_Set_GEntID(mf,(list_face[nfes_local+nfe+1] >> 3));
      MF_Set_PType(mf,(list_face[nfes_local+nfe+2] & 3));
      MF_Set_MasterParID(mf,(list_face[nfes_local+nfe+2] >> 2));
      MF_Set_GlobalID(mf,list_face[nfes_local+nfe+3]);

      MF_Set_Edges(mf,nfe,fedges,fedirs);
      nfes_local += (nfe + 4);
    }
    if (fedges) free(fedges);
    if (fedirs) free(fedirs);
    if (edges) free(edges);


    /* receive region info */

    int *list_region = (int *) malloc(nrfs*sizeof(int));

    MPI_Recv(list_region,nrfs,MPI_INT,send_rank,rank,comm,&status);
    MPI_Get_count(&status,MPI_INT,&count);

    rfaces = (MFace_ptr *) malloc(MAXPF3*sizeof(MFace_ptr));
    rfdirs = (int *) malloc(MAXPF3*sizeof(int));
    nrfs_local = 0;
    for (i = 0; i < nr; i++) {
      mr = MR_New(mesh);
      nrf = list_region[nrfs_local];
      for (j = 0; j < nrf; j++) {
	rfaces[j] = faces[abs(list_region[nrfs_local+j+1])-1];
	rfdirs[j] = list_region[nrfs_local+j+1] > 0 ? 1 : 0;
      }
      MR_Set_GEntDim(mr,(list_region[nrfs_local+nrf+1] & 7));
      MR_Set_GEntID(mr,(list_region[nrfs_local+nrf+1] >> 3));
      MR_Set_PType(mr,(list_region[nrfs_local+nrf+2] & 3));
      MR_Set_MasterParID(mr,(list_region[nrfs_local+nrf+2] >> 2));
      MR_Set_GlobalID(mr,list_region[nrfs_local+nrf+3]);

      MR_Set_Faces(mr,nrf,rfaces,rfdirs); // won't work for polyhedra
      nrfs_local += (nrf + 4);
    }
    if (rfaces) free(rfaces);
    if (rfdirs) free(rfdirs);
    if (faces) free(faces);

    /* receive attr */
    if(natt) {
      list_attr = (int *) malloc(natt*sizeof(int));
      list_attr_names = (char *) malloc((natt)*256*sizeof(char));
    
      MPI_Recv(list_attr,natt,MPI_INT,send_rank,rank,comm,&status);
      MPI_Recv(list_attr_names,natt*256,MPI_CHAR,send_rank,rank,comm,&status);
      MPI_Get_count(&status,MPI_INT,&count);

      for(i = 0; i < natt; i++) {
	strcpy(attname,&list_attr_names[i*256]);
	/* see of the attrib exists */
	if(MESH_AttribByName(mesh,attname))
	  continue;
	att_type = list_attr[i] & 7;
	mtype = (list_attr[i] >> 3) & 7;
	ncomp = list_attr[i] >> 6;
	if(ncomp == 1)
	  attrib =  MAttrib_New(mesh, attname, att_type, mtype);
	else
	  attrib =  MAttrib_New(mesh, attname, att_type, mtype, ncomp);
	/* printf("attr %d with name %s received on rank %d\n",i, attname, rank); */
      }
      MSTK_free(list_attr);    
      MSTK_free(list_attr_names);
    }


    /* receive mesh entity sets */
    if(nset) {
      list_mset = (int *) malloc(nset*sizeof(int));
      list_mset_names = (char *) malloc((nset)*256*sizeof(char));
    
      MPI_Recv(list_mset,nset,MPI_INT,send_rank,rank,comm,&status);
      MPI_Recv(list_mset_names,nset*256,MPI_CHAR,send_rank,rank,comm,&status);
      MPI_Get_count(&status,MPI_INT,&count);

      for(i = 0; i < nset; i++) {
	strcpy(msetname,&list_mset_names[i*256]);
	/* see if the mset exists */
	if(MESH_MSetByName(mesh,msetname))
	  continue;
	mtype = list_mset[i];
	mset =  MSet_New(mesh, msetname, mtype);
      }
      MSTK_free(list_mset);    
      MSTK_free(list_mset_names);
    }



    MSTK_free(list_vertex);    
    MSTK_free(list_region);    
    MSTK_free(list_coor);    
    return 1;
  }





  int MESH_Surf_RecvMesh_R1R2R4(Mesh_ptr mesh, int *mesh_info, int send_rank, int rank, MPI_Comm comm) {
    MSTK_Report("MESH_Surf_RecvMesh_R1R2R4","Not implemented",MSTK_FATAL);
  }

  int MESH_Vol_RecvMesh_R1R2(Mesh_ptr mesh, int *mesh_info, int send_rank, int rank, MPI_Comm comm) {
    MSTK_Report("MESH_Vol_RecvMesh_R1R2","Not implemented",MSTK_FATAL);
  }


  int MESH_Vol_RecvMesh_R4(Mesh_ptr mesh, int *mesh_info, int send_rank, int rank, MPI_Comm comm) {

    MSTK_Report("MESH_Vol_RecvMesh_R4","Not implemented",MSTK_FATAL);
  }





#ifdef __cplusplus
}
#endif

