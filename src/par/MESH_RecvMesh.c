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
     This function receives mesh from processor rank in communicator comm
     attrib list is received, but no attrib values.
     call MESH_RecvAttr() to update entity attribute values
     fromrank: the rank of sending processor
     rank: the rank of receiving processor

     Author(s): Duo Wang, Rao Garimella
  */

  int MESH_Vol_RecvMesh_FN(Mesh_ptr mesh, int *mesh_info, int fromrank, MSTK_Comm comm);
  int MESH_Surf_RecvMesh_FN(Mesh_ptr mesh, int *mesh_info, int fromrank, MSTK_Comm comm);
  int MESH_Vol_RecvMesh_R4(Mesh_ptr mesh, int *mesh_info, int fromrank, MSTK_Comm comm);
  int MESH_Vol_RecvMesh_R1R2(Mesh_ptr mesh, int *mesh_info, int fromrank, MSTK_Comm comm);
  int MESH_Surf_RecvMesh_R1R2R4(Mesh_ptr mesh, int *mesh_info, int fromrank, MSTK_Comm comm);

  static int (*MESH_Vol_RecvMesh_jmp[MSTK_MAXREP])(Mesh_ptr mesh, 
						   int *mesh_info, 
						   int fromrank, MSTK_Comm comm) =
  {MESH_Vol_RecvMesh_FN, MESH_Vol_RecvMesh_FN, MESH_Vol_RecvMesh_R1R2, 
   MESH_Vol_RecvMesh_R1R2, MESH_Vol_RecvMesh_R4};
  static int (*MESH_Surf_RecvMesh_jmp[MSTK_MAXREP])(Mesh_ptr mesh, 
						    int *mesh_info, 
						    int fromrank, MSTK_Comm comm) =
  {MESH_Surf_RecvMesh_FN, MESH_Surf_RecvMesh_FN, MESH_Surf_RecvMesh_R1R2R4, 
   MESH_Surf_RecvMesh_R1R2R4, MESH_Surf_RecvMesh_R1R2R4};



  int MESH_RecvMesh(Mesh_ptr mesh, int dim, int fromrank, MSTK_Comm comm) {
    int mesh_info[9], count;
    RepType rtype;
    MPI_Request request;
    MPI_Status status;
    char mesg[256], errorstr[256];
    int len, errcode;

    int rank;
    MPI_Comm_rank(comm,&rank);

    /* mesh_info store the mesh reptype, nv, nf, nfvs and natt */
    /* receive mesh_info */
    errcode = MPI_Recv(mesh_info,9,MPI_INT,fromrank,rank,comm,&status);
    if (errcode != MPI_SUCCESS) {
      MPI_Error_string(errcode,errorstr,&len);
      sprintf(mesg,"MPI_Recv failed with error %s",errorstr);
      MSTK_Report("MESH_RecvMesh",mesg,MSTK_FATAL);
    }
    /* MPI_Get_count(&status,MPI_INT,&count); */
    
    rtype = mesh_info[0];
    
    if (dim == 3)
      (*MESH_Vol_RecvMesh_jmp[rtype])(mesh,mesh_info,fromrank,comm);
    else if(dim == 2) 
      (*MESH_Surf_RecvMesh_jmp[rtype])(mesh,mesh_info,fromrank,comm);
    else {
      MSTK_Report("MESH_RecvMesh()","only receive volume or surface mesh",MSTK_ERROR);
      exit(-1);
    }
    return 1;
  }
  





  int MESH_Surf_RecvMesh_FN(Mesh_ptr mesh, int *mesh_info, int fromrank, MSTK_Comm comm) {
    int i, j, nevs, nfe, nfes, natt, nset, count, ncomp;
    int nv, ne, nf, nr, nrfs, *fedirs;
    int maxnfe, maxnrf;
    MVertex_ptr *verts, *fverts;
    MEdge_ptr me, *edges, *fedges;
    MFace_ptr mf;
    RepType rtype;
    char attname[256], msetname[256];
    MType mtype;
    MAttType att_type;
    MAttrib_ptr attrib;
    MSet_ptr mset;
    int *list_attr, *list_mset;
    char *list_attr_names, *list_mset_names;
    double coor[3];
    MPI_Request request[15];
    /*    MPI_Status status[15]; */
    MPI_Status status;
    char mesg[256], errorstr[256];
    int errcode, len, nreq=0;

    int rank;
    MPI_Comm_rank(comm,&rank);

    /* basic mesh information */
    rtype = mesh_info[0];
    nv = mesh_info[1];
    ne = mesh_info[2];
    nf = mesh_info[3];
    nr = mesh_info[4]; /* should be 0 since its a surface mesh */
    maxnfe = mesh_info[5];
    maxnrf = mesh_info[6]; /* should be 0 since its a surface mesh */
    natt = mesh_info[7];
    nset = mesh_info[8];
    
    MESH_SetRepType(mesh,rtype);

    /* copied from MESH_SendMesh */

    nevs = (2+3)*ne;
    nfes = (1 + maxnfe + 3)*nf; 


    /* allocate receive buffer */
    int *list_vertex = (int *) malloc(3*nv*sizeof(int));
    double *list_coor = (double *) malloc(3*nv*sizeof(double));

    /* receive vertex info */
    errcode = MPI_Irecv(list_vertex,3*nv,MPI_INT,fromrank,rank,comm,
			&(request[nreq++]));
    if (errcode != MPI_SUCCESS)
      MSTK_Report("MESH_Surf_RecvMesh_FN","Trouble receiving mesh",
		  MSTK_FATAL);
    
    errcode = MPI_Irecv(list_coor,3*nv,MPI_DOUBLE,fromrank,rank,comm,
			&(request[nreq++]));
    if (errcode != MPI_SUCCESS)
      MSTK_Report("MESH_Surf_RecvMesh_FN","Trouble receiving mesh",
		  MSTK_FATAL);
    

    /* receive edge-vertex info */

    int *list_edge = (int *) malloc(nevs*sizeof(int));

    errcode = MPI_Irecv(list_edge,nevs,MPI_INT,fromrank,rank,comm,
                        &(request[nreq++]));
    if (errcode != MPI_SUCCESS)
      MSTK_Report("MESH_Surf_RecvMesh_FN","Trouble receiving mesh",
		  MSTK_FATAL);

    /* receive face-edge info */

    int *list_face = (int *) malloc(nfes*sizeof(int));

    errcode = MPI_Irecv(list_face,nfes,MPI_INT,fromrank,rank,comm,
                        &(request[nreq++]));
    if (errcode != MPI_SUCCESS)
      MSTK_Report("MESH_Surf_RecvMesh_FN","Trouble receiving mesh",
		  MSTK_FATAL);


    if (natt) {
      list_attr = (int *) malloc(natt*sizeof(int));
      list_attr_names = (char *) malloc((natt)*256*sizeof(char));

      errcode = MPI_Irecv(list_attr,natt,MPI_INT,fromrank,rank,comm,
			  &(request[nreq++]));
      

      errcode = MPI_Irecv(list_attr_names,natt*256,MPI_CHAR,fromrank,rank,comm,
			  &(request[nreq++]));
    }

    if (nset) {
      list_mset = (int *) malloc(nset*sizeof(int));
      list_mset_names = (char *) malloc((nset)*256*sizeof(char));
      
      errcode = MPI_Irecv(list_mset,nset,MPI_INT,fromrank,rank,comm,
			  &(request[nreq++]));

      errcode = MPI_Irecv(list_mset_names,nset*256,MPI_CHAR,fromrank,rank,comm,
			  &(request[nreq++]));
    }


    errcode = MPI_Waitall(nreq,request,MPI_STATUSES_IGNORE);
    if (errcode != MPI_SUCCESS)
      MSTK_Report("MESH_Surf_RecvMesh_FN","Trouble receiving mesh",
		  MSTK_FATAL);

    nreq = 0;  // reset requests

    
    /* build vertices */


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


    /* build edges */

    edges = (MEdge_ptr *) malloc(ne*sizeof(MEdge_ptr));

    nevs = 0;
    for (i = 0; i < ne; i++) {
      edges[i] = me = ME_New(mesh);
      ME_Set_GEntDim(me,(list_edge[nevs+2] & 7));
      ME_Set_GEntID(me,(list_edge[nevs+2] >> 3));
      ME_Set_PType(me,(list_edge[nevs+3] & 3));
      ME_Set_MasterParID(me,(list_edge[nevs+3] >> 2));
      ME_Set_GlobalID(me,list_edge[nevs+4]);

      ME_Set_Vertex(me,0,verts[list_edge[nevs]-1]);
      ME_Set_Vertex(me,1,verts[list_edge[nevs+1]-1]);
      nevs += 5;
    }
    if (verts) free(verts);

    
    /* build faces */

    fedges = (MEdge_ptr *) malloc(maxnfe*sizeof(MEdge_ptr));
    fedirs = (int *) malloc(maxnfe*sizeof(int));
    nfes = 0;
    for (i = 0; i < nf; i++) {
      mf = MF_New(mesh);
      nfe = list_face[nfes];
      for (j = 0; j < nfe; j++) {
	fedges[j] = edges[abs(list_face[nfes+j+1])-1];
	fedirs[j] = list_face[nfes+j+1] > 0 ? 1 : 0;
      }
      MF_Set_GEntDim(mf,(list_face[nfes+nfe+1] & 7));
      MF_Set_GEntID(mf,(list_face[nfes+nfe+1] >> 3));
      MF_Set_PType(mf,(list_face[nfes+nfe+2] & 3));
      MF_Set_MasterParID(mf,(list_face[nfes+nfe+2] >> 2));
      MF_Set_GlobalID(mf,list_face[nfes+nfe+3]);

      MF_Set_Edges(mf,nfe,fedges,fedirs);
      nfes += (nfe + 4);
    }
    if (fedges) free(fedges);
    if (fedirs) free(fedirs);
    if (edges) free(edges);


    free(list_vertex);    
    free(list_face);    
    free(list_coor);


    /* receive and build attributes */

    for(i = 0; i < natt; i++) {
      strcpy(attname,&list_attr_names[i*256]);
      
      if (MESH_AttribByName(mesh,attname)) 	/* see if the attrib exists */
        continue;
      
      att_type = list_attr[i] & 7;
      mtype = (list_attr[i] >> 3) & 7;
      ncomp = list_attr[i] >> 6;
      if(ncomp == 1)
        attrib =  MAttrib_New(mesh, attname, att_type, mtype);
      else
        attrib =  MAttrib_New(mesh, attname, att_type, mtype, ncomp);
    }
    
    if (natt) {
      free(list_attr);    
      free(list_attr_names);
    }
    

    /* receive mesh entity sets */
   
    for(i = 0; i < nset; i++) {
      strcpy(msetname,&list_mset_names[i*256]);
      
      if(MESH_MSetByName(mesh,msetname)) /* see if the mset exists */
        continue;
      mtype = list_mset[i];
      mset =  MSet_New(mesh, msetname, mtype);
    }


    if (nset) {
      free(list_mset);    
      free(list_mset_names);
    }

    return 1;
  }


  int MESH_Vol_RecvMesh_FN(Mesh_ptr mesh, int *mesh_info, int fromrank, MPI_Comm comm) {
    int i, j, nevs, nfes, nrv, nrfs, maxnfe, maxnrf;
    int count, natt, nset, ncomp, nv, ne, nf, nr, nfe, nrf, *fedirs, *rfdirs;
    MVertex_ptr *verts, *rverts;
    MEdge_ptr me, *edges, *fedges;
    MFace_ptr mf, *faces, *rfaces;
    MRegion_ptr mr;
    RepType rtype;
    char attname[256], msetname[256];
    MType mtype;
    MAttType att_type;
    MAttrib_ptr attrib;
    MSet_ptr mset;
    int *list_attr, *list_mset;
    char *list_attr_names, *list_mset_names;
    double coor[3];
    MPI_Request request[15];
    /* MPI_Status status[15]; */
    MPI_Status status;
    char mesg[256], errorstr[256];
    int errcode, len, nreq=0;
   
    int rank;
    MPI_Comm_rank(comm,&rank);

    /* mesh_info store the mesh reptype, nv, nr, nrvs and natt */

    /* basic mesh information */
    rtype = mesh_info[0];
    nv = mesh_info[1];
    ne = mesh_info[2];
    nf = mesh_info[3];
    nr = mesh_info[4];    
    maxnfe = mesh_info[5];
    maxnrf = mesh_info[6];
    natt = mesh_info[7];
    nset = mesh_info[8];
    
    MESH_SetRepType(mesh,rtype);

    /* copied from MESH_SendMesh */

    nevs = (2+3)*ne;
    nfes = (1 + maxnfe + 3)*nf; 
    nrfs = (1 + maxnrf + 3)*nr;

    /* allocate receive buffer */
    int *list_vertex = (int *) malloc(3*nv*sizeof(int));
    double *list_coor = (double *) malloc(3*nv*sizeof(double));

    /* receive vertex info */

    errcode = MPI_Irecv(list_vertex,3*nv,MPI_INT,fromrank,rank,comm,
			&(request[nreq++]));
    if (errcode != MPI_SUCCESS)
      MSTK_Report("MESH_Vol_RecvMesh_FN","Trouble receiving mesh",
		  MSTK_FATAL);
    

    errcode = MPI_Irecv(list_coor,3*nv,MPI_DOUBLE,fromrank,rank,comm,
		       &(request[nreq++]));
    if (errcode != MPI_SUCCESS)
      MSTK_Report("MESH_Vol_RecvMesh_FN","Trouble receiving mesh",
		  MSTK_FATAL);


    /* receive edge-vertex info */

    int *list_edge = (int *) malloc(nevs*sizeof(int));

    errcode = MPI_Irecv(list_edge,nevs,MPI_INT,fromrank,rank,comm,
		       &(request[nreq++]));
    if (errcode != MPI_SUCCESS)
      MSTK_Report("MESH_Vol_RecvMesh_FN","Trouble receiving mesh",
		  MSTK_FATAL);
    

    /* receive face-edge info */
    
    int *list_face = (int *) malloc(nfes*sizeof(int));

    errcode = MPI_Irecv(list_face,nfes,MPI_INT,fromrank,rank,comm,
		       &(request[nreq++]));
    if (errcode != MPI_SUCCESS)
      MSTK_Report("MESH_Vol_RecvMesh_FN","Trouble receiving mesh",
		  MSTK_FATAL);
    


    /* receive region-face info */
    
    int *list_region = (int *) malloc(nrfs*sizeof(int));
    
    errcode = MPI_Irecv(list_region,nrfs,MPI_INT,fromrank,rank,comm,
		       &(request[nreq++]));
    if (errcode != MPI_SUCCESS)
      MSTK_Report("MESH_Vol_RecvMesh_FN","Trouble receiving mesh",
		  MSTK_FATAL);
    
    
    /* receive and build attributes */
    
    if (natt) {
      list_attr = (int *) malloc(natt*sizeof(int));
      list_attr_names = (char *) malloc((natt)*256*sizeof(char));
    
      errcode = MPI_Irecv(list_attr,natt,MPI_INT,fromrank,rank,comm,
			  &(request[nreq++]));
      if (errcode != MPI_SUCCESS)
	MSTK_Report("MESH_Vol_RecvMesh_FN","Trouble receiving mesh",
		    MSTK_FATAL);

      
      errcode = MPI_Irecv(list_attr_names,natt*256,MPI_CHAR,fromrank,rank,comm,
                          &(request[nreq++]));
      if (errcode != MPI_SUCCESS)
	MSTK_Report("MESH_Vol_RecvMesh_FN","Trouble receiving mesh",
		    MSTK_FATAL);
      
    }
    

    /* receive mesh entity sets */

    if (nset) {
      list_mset = (int *) malloc(nset*sizeof(int));
      list_mset_names = (char *) malloc((nset)*256*sizeof(char));
    
      errcode = MPI_Irecv(list_mset,nset,MPI_INT,fromrank,rank,comm,
			  &(request[nreq++]));
          if (errcode != MPI_SUCCESS)
      MSTK_Report("MESH_Vol_RecvMesh_FN","Trouble receiving mesh",
		  MSTK_FATAL);


      errcode = MPI_Irecv(list_mset_names,nset*256,MPI_CHAR,fromrank,rank,comm,
			  &(request[nreq++]));
      if (errcode != MPI_SUCCESS)
	MSTK_Report("MESH_Vol_RecvMesh_FN","Trouble receiving mesh",
		    MSTK_FATAL);
    }


    errcode = MPI_Waitall(nreq,request,MPI_STATUSES_IGNORE);
    if (errcode != MPI_SUCCESS)
      MSTK_Report("MESH_Vol_RecvMesh_FN","Trouble receiving mesh",MSTK_FATAL);
    nreq = 0;  // reset requests



    /* build vertices */

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


    /* build edges */

    edges = (MEdge_ptr *) malloc(ne*sizeof(MEdge_ptr));

    nevs = 0;
    for (i = 0; i < ne; i++) {
      edges[i] = me = ME_New(mesh);
      ME_Set_GEntDim(me,(list_edge[nevs+2] & 7));
      ME_Set_GEntID(me,(list_edge[nevs+2] >> 3));
      ME_Set_PType(me,(list_edge[nevs+3] & 3));
      ME_Set_MasterParID(me,(list_edge[nevs+3] >> 2));
      ME_Set_GlobalID(me,list_edge[nevs+4]);

      ME_Set_Vertex(me,0,verts[list_edge[nevs]-1]);
      ME_Set_Vertex(me,1,verts[list_edge[nevs+1]-1]);
      nevs += 5;
    }
    if (verts) free(verts);



    /* build faces */

    fedges = (MEdge_ptr *) malloc(maxnfe*sizeof(MEdge_ptr));
    fedirs = (int *) malloc(maxnfe*sizeof(int));

    faces = (MFace_ptr *) malloc(nf*sizeof(MFace_ptr));

    nfes = 0;
    for (i = 0; i < nf; i++) {
      faces[i] = mf = MF_New(mesh);
      nfe = list_face[nfes];
      for (j = 0; j < nfe; j++) {
	fedges[j] = edges[abs(list_face[nfes+j+1])-1];
	fedirs[j] = list_face[nfes+j+1] > 0 ? 1 : 0;
      }
      MF_Set_GEntDim(mf,(list_face[nfes+nfe+1] & 7));
      MF_Set_GEntID(mf,(list_face[nfes+nfe+1] >> 3));
      MF_Set_PType(mf,(list_face[nfes+nfe+2] & 3));
      MF_Set_MasterParID(mf,(list_face[nfes+nfe+2] >> 2));
      MF_Set_GlobalID(mf,list_face[nfes+nfe+3]);

      MF_Set_Edges(mf,nfe,fedges,fedirs);
      nfes += (nfe + 4);
    }
    if (fedges) free(fedges);
    if (fedirs) free(fedirs);
    if (edges) free(edges);


    /* build regions */

    rfaces = (MFace_ptr *) malloc(maxnrf*sizeof(MFace_ptr));
    rfdirs = (int *) malloc(maxnrf*sizeof(int));
    nrfs = 0;
    for (i = 0; i < nr; i++) {
      mr = MR_New(mesh);
      nrf = list_region[nrfs];
      for (j = 0; j < nrf; j++) {
	rfaces[j] = faces[abs(list_region[nrfs+j+1])-1];
	rfdirs[j] = list_region[nrfs+j+1] > 0 ? 1 : 0;
      }
      MR_Set_GEntDim(mr,(list_region[nrfs+nrf+1] & 7));
      MR_Set_GEntID(mr,(list_region[nrfs+nrf+1] >> 3));
      MR_Set_PType(mr,(list_region[nrfs+nrf+2] & 3));
      MR_Set_MasterParID(mr,(list_region[nrfs+nrf+2] >> 2));
      MR_Set_GlobalID(mr,list_region[nrfs+nrf+3]);

      MR_Set_Faces(mr,nrf,rfaces,rfdirs); // won't work for polyhedra? why not?
      nrfs += (nrf + 4);
    }
    if (rfaces) free(rfaces);
    if (rfdirs) free(rfdirs);
    if (faces) free(faces);


    free(list_vertex);   
    free(list_face);
    free(list_edge);
    free(list_region);    
    free(list_coor);    


    /* build attributes */

    for(i = 0; i < natt; i++) {
      strcpy(attname,&list_attr_names[i*256]);
      
      if(MESH_AttribByName(mesh,attname)) /* see of the attrib exists */
        continue;
      
      att_type = list_attr[i] & 7;
      mtype = (list_attr[i] >> 3) & 7;
      ncomp = list_attr[i] >> 6;
      if(ncomp == 1)
        attrib =  MAttrib_New(mesh, attname, att_type, mtype);
      else
        attrib =  MAttrib_New(mesh, attname, att_type, mtype, ncomp);
    }

    if (natt) {
      free(list_attr);    
      free(list_attr_names);
    }


    /* build sets */
    
    for(i = 0; i < nset; i++) {
      strcpy(msetname,&list_mset_names[i*256]);
      
      if (MESH_MSetByName(mesh,msetname)) /* see if the mset exists */
        continue;
      mtype = list_mset[i];
      mset =  MSet_New(mesh, msetname, mtype);
    }
    
    if (nset) {
      free(list_mset);    
      free(list_mset_names);
    }
    
    return 1;
  }





  int MESH_Surf_RecvMesh_R1R2R4(Mesh_ptr mesh, int *mesh_info, int fromrank, MSTK_Comm comm) {
    MSTK_Report("MESH_Surf_RecvMesh_R1R2R4","Not implemented",MSTK_FATAL);
    return 0;
  }

  int MESH_Vol_RecvMesh_R1R2(Mesh_ptr mesh, int *mesh_info, int fromrank, MSTK_Comm comm) {
    MSTK_Report("MESH_Vol_RecvMesh_R1R2","Not implemented",MSTK_FATAL);
    return 0;
  }


  int MESH_Vol_RecvMesh_R4(Mesh_ptr mesh, int *mesh_info, int fromrank, MSTK_Comm comm) {
    MSTK_Report("MESH_Vol_RecvMesh_R4","Not implemented",MSTK_FATAL);
    return 0;
  }





#ifdef __cplusplus
}
#endif

