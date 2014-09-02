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
     This function send mesh to processor torank in communicator comm
     attr list is sent but no attribute value is sent.
     call MESH_SendAttr() to send attribute values of entities

     Author(s): Duo Wang, Rao Garimella
  */


  int MESH_Vol_SendMesh_FN(Mesh_ptr mesh, int torank, MSTK_Comm comm, 
                           int *numreq, int *maxreq, MPI_Request **requests,
                           int *numptrs2free, int *maxptrs2free,
                           void ***ptrs2free);
  int MESH_Surf_SendMesh_FN(Mesh_ptr mesh, int torank, MSTK_Comm comm, 
                            int *numreq, int *maxreq, MPI_Request **requests,
                            int *numptrs2free, int *maxptrs2free,
                            void ***ptrs2free);
  int MESH_Vol_SendMesh_R4(Mesh_ptr mesh, int torank, MSTK_Comm comm, 
                           int *numreq, int *maxreq, MPI_Request **requests,
                           int *numptrs2free, int *maxptrs2free,
                           void ***ptrs2free);
  int MESH_Vol_SendMesh_R1R2(Mesh_ptr mesh, int torank, MSTK_Comm comm, 
                             int *numreq, int *maxreq, MPI_Request **requests,
                             int *numptrs2free, int *maxptrs2free,
                             void ***ptrs2free);
  int MESH_Surf_SendMesh_R1R2R4(Mesh_ptr mesh, int torank, MSTK_Comm comm, 
                                int *numreq, int *maxreq, 
                                MPI_Request **requests,
                                int *numptrs2free, int *maxptrs2free,
                                void ***ptrs2free);

  static int (*MESH_Vol_SendMesh_jmp[MSTK_MAXREP])(Mesh_ptr mesh, int torank, 
                                                   MSTK_Comm comm, int *numreq,
                                                   int *maxreq, 
                                                   MPI_Request **requests,
                                                   int *numptrs2free, 
                                                   int *maxptrs2free,
                                                   void ***ptrs2free) = 
  {MESH_Vol_SendMesh_FN, MESH_Vol_SendMesh_FN, MESH_Vol_SendMesh_R1R2, 
   MESH_Vol_SendMesh_R1R2, MESH_Vol_SendMesh_R4};
  static int (*MESH_Surf_SendMesh_jmp[MSTK_MAXREP])(Mesh_ptr mesh, int torank, 
                                                    MSTK_Comm comm, int *numreq,
                                                    int *maxreq, 
                                                    MPI_Request **requests,
                                                    int *numptrs2free, 
                                                    int *maxptrs2free,
                                                    void ***ptrs2free) =
  {MESH_Surf_SendMesh_FN, MESH_Surf_SendMesh_FN, MESH_Surf_SendMesh_R1R2R4, 
   MESH_Surf_SendMesh_R1R2R4, MESH_Surf_SendMesh_R1R2R4};


  int MESH_SendMesh(Mesh_ptr mesh, int torank, MSTK_Comm comm, int *numreq, 
                    int *maxreq, MPI_Request **requests,
                    int *numptrs2free, int *maxptrs2free, void ***ptrs2free) {
    int nf, nr;
    RepType rtype;

    /* basic mesh information */
    rtype = MESH_RepType(mesh);
    nf = MESH_Num_Faces(mesh);
    nr = MESH_Num_Regions(mesh);
    if (nr)
      (*MESH_Vol_SendMesh_jmp[rtype])(mesh,torank,comm,numreq,maxreq,requests,
                                      numptrs2free,maxptrs2free,ptrs2free);
    else if(nf) 
      (*MESH_Surf_SendMesh_jmp[rtype])(mesh,torank,comm,numreq,maxreq,requests,
                                       numptrs2free,maxptrs2free,ptrs2free);
    else {
      MSTK_Report("MESH_SendMesh()","only send volume or surface mesh",MSTK_ERROR);
      exit(-1);
    }

    return 1;
  }







  int MESH_Surf_SendMesh_FN(Mesh_ptr mesh, int torank, MSTK_Comm comm, 
                            int *numreq, int *maxreq, MPI_Request **requests,
                            int *numptrs2free, int *maxptrs2free,
                            void ***ptrs2free) {
    int i, j, nv, ne, nf, nevs, nfes, nfv, natt, nset, ncomp, dir;
    int nfe;
    int mesh_info[10];
    MVertex_ptr mv;
    MEdge_ptr me;
    MFace_ptr mf;
    List_ptr fverts, mfedges;
    RepType rtype;
    char attname[256], msetname[256];
    MType mtype;
    MAttType att_type;
    MAttrib_ptr attrib;
    MSet_ptr mset;
    int *list_attr=NULL, *list_mset=NULL;
    char *list_attr_names=NULL, *list_mset_names=NULL;
    double coor[3];
    MPI_Request mpirequest;
  
    natt = MESH_Num_Attribs(mesh);
    nset = MESH_Num_MSets(mesh);

    if (requests == NULL)
      MSTK_Report("MESH_Surf_SendMesh","MPI requests array is NULL",MSTK_FATAL);
    
    if (*maxreq == 0) {
      *maxreq = 25;
      *requests = (MPI_Request *) malloc(*maxreq*sizeof(MPI_Request));
      *numreq = 0;
    }
    else if (*maxreq < (*numreq) + 11 + 2*natt + 2*nset) {
      *maxreq = 2*(*maxreq) + 11 + 2*natt + 2*nset;
      *requests = (MPI_Request *) realloc(*requests,*maxreq*sizeof(MPI_Request));
    }
  

    /* mesh_info store the mesh reptype, nv, nf, nfvs and natt */

    rtype = MESH_RepType(mesh);
    nv = MESH_Num_Vertices(mesh);
    ne = MESH_Num_Edges(mesh);
    nf = MESH_Num_Faces(mesh);

    /* some other known quantitites - 5 items per edge (2 for verts
       and 3 for extra data), MAXPV2+4 items per face (1 for number of
       edges, MAXPV2 for edge indices, anad 3 for extra data) */

    nevs = (2+3)*ne;  
    nfes = (1 + MAXPV2 + 3)*nf;
    
    mesh_info[0] = rtype;
    mesh_info[1] = nv;
    mesh_info[2] = ne;
    mesh_info[3] = nf;
    mesh_info[4] = 0;
    mesh_info[5] = nevs;
    mesh_info[6] = nfes;
    mesh_info[7] = 0;
    mesh_info[8] = natt;
    mesh_info[9] = nset;


    /* Send some global mesh info */

    MPI_Isend(mesh_info,10,MPI_INT,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;


    /* Now send out detailed vertex info */

    int *list_vertex = (int *) malloc(3*nv*sizeof(int));
    double *list_coor = (double *) malloc(3*nv*sizeof(double));

    /* Store the 3 auxilliary data fields */
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

    /* send vertices */
    MPI_Isend(list_vertex,3*nv,MPI_INT,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;;
    MPI_Isend(list_coor,3*nv,MPI_DOUBLE,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;;



    /* Reserve nevs spots for each edge */

    int *list_edge = (int *) malloc(nevs*sizeof(int));

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

    /* send detailed edge info */

    MPI_Isend(list_edge,nevs,MPI_INT,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;;



    /* Reserve nfes spots for each face */

    int *list_face = (int *) malloc((MAXPV2+4)*nf*sizeof(int));

    nfes = 0;

    /* first store nfe, then the edge ids, then the 3 auxilliary data fields */
    for(i = 0; i < nf; i++) {
      mf = MESH_Face(mesh,i);
      mfedges = MF_Edges(mf,1,0);
      nfe = List_Num_Entries(mfedges);
      list_face[nfes] = nfe;
      for(j = 0; j < nfe; j++) {
        dir = MF_EdgeDir_i(mf,j) == 1 ? 1 : -1;
        list_face[nfes+j+1] = dir*ME_ID(List_Entry(mfedges,j));
      }
      list_face[nfes+nfe+1] = (MF_GEntID(mf)<<3) | (MF_GEntDim(mf));
      list_face[nfes+nfe+2] = (MF_MasterParID(mf)<<2) | (MF_PType(mf));
      list_face[nfes+nfe+3] = MF_GlobalID(mf);
      nfes += (nfe + 4);
      List_Delete(mfedges);
    }

    /* send detailed face info */
    MPI_Isend(list_face,nfes,MPI_INT,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;;
  


    /* collect attrs */
    if(natt) {
      list_attr = (int *) malloc((natt)*sizeof(int));
      list_attr_names = (char *) malloc((natt)*256*sizeof(char));
    
      for(i = 0; i < natt; i++) {
        attrib = MESH_Attrib(mesh,i);
        MAttrib_Get_Name(attrib,attname);
        att_type = MAttrib_Get_Type(attrib);
        ncomp = MAttrib_Get_NumComps(attrib);
        mtype = MAttrib_Get_EntDim(attrib);
        list_attr[i] = (ncomp << 6) | (mtype << 3) | (att_type);
        strcpy(&list_attr_names[i*256],attname);
      }

      /* send attr */
      /* printf("%d attrs sent to torank %d\n",natt,torank); */
      MPI_Isend(list_attr,natt,MPI_INT,torank,torank,comm,&mpirequest);
      (*requests)[*numreq] = mpirequest;
      (*numreq)++;;
      MPI_Isend(list_attr_names,natt*256,MPI_CHAR,torank,torank,comm,&mpirequest);
      (*requests)[*numreq] = mpirequest;
      (*numreq)++;
    }

    /* Mesh entity sets */

    if (nset) {
      list_mset = (int *) malloc(nset*sizeof(int));
      list_mset_names = (char *) malloc(nset*256*sizeof(char));

      for (i = 0; i < nset; i++) {
        mset = MESH_MSet(mesh,i);
        MSet_Name(mset,msetname);
        mtype = MSet_EntDim(mset);
        list_mset[i] = mtype;
        strcpy(&list_mset_names[i*256],msetname);
      }

      MPI_Isend(list_mset,nset,MPI_INT,torank,torank,comm,&mpirequest);
      (*requests)[*numreq] = mpirequest;
      (*numreq)++;;
      MPI_Isend(list_mset_names,nset*256,MPI_CHAR,torank,torank,comm,&mpirequest);
      (*requests)[*numreq] = mpirequest;
      (*numreq)++;;
    }

    int nptrs = 5;
    if (natt) nptrs += 2;
    if (nset) nptrs += 2;

    if (*maxptrs2free == 0) {
      *maxptrs2free = 25;
      *ptrs2free = (void **) malloc(*maxptrs2free*sizeof(void *));
      *numptrs2free = 0;
    }
    else if (*maxptrs2free < (*numptrs2free) + nptrs) {
      *maxptrs2free *= 2;
      *ptrs2free = (void **) realloc(*ptrs2free,(*maxptrs2free)*sizeof(void *));
    }

    (*ptrs2free)[(*numptrs2free)++] = list_vertex;
    (*ptrs2free)[(*numptrs2free)++] = list_coor;  
    (*ptrs2free)[(*numptrs2free)++] = list_edge;
    (*ptrs2free)[(*numptrs2free)++] = list_face;
    if (natt) {
      (*ptrs2free)[(*numptrs2free)++] = list_attr;
      (*ptrs2free)[(*numptrs2free)++] = list_attr_names;
    }
    if (nset) {
      (*ptrs2free)[(*numptrs2free)++] = list_mset;
      (*ptrs2free)[(*numptrs2free)++] = list_mset_names;
    }

    return 1;
  }


  int MESH_Vol_SendMesh_FN(Mesh_ptr mesh, int torank, MSTK_Comm comm,
                           int *numreq, int *maxreq, MPI_Request **requests,
                           int *numptrs2free, int *maxptrs2free,
                           void ***ptrs2free) {
    int i, j, nv, ne, nf, nr;
    int nevs, nfes, nrfs, nfe, nrv, nrf, natt, nset, ncomp, dir;
    int mesh_info[10];
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
    int *list_attr=NULL, *list_mset=NULL;
    char *list_attr_names=NULL, *list_mset_names=NULL;
    double coor[3];
    MPI_Request mpirequest;

    natt = MESH_Num_Attribs(mesh);
    nset = MESH_Num_MSets(mesh);

    if (requests == NULL)
      MSTK_Report("MESH_Surf_SendMesh","MPI requests array is NULL",MSTK_FATAL);

    if (*maxreq == 0) {
      *maxreq = 25;
      *requests = (MPI_Request *) malloc(*maxreq*sizeof(MPI_Request));
      *numreq = 0;
    }
    else if (*maxreq < (*numreq) + 13 + 2*natt + 2*nset) {
      *maxreq = 2*(*maxreq) + 11 + 2*natt + 2*nset;
      *requests = (MPI_Request *) realloc(*requests,*maxreq*sizeof(MPI_Request));
    }
  

    /* mesh_info store the mesh reptype, nv, ne, nf, nr and natt, nset */

    rtype = MESH_RepType(mesh);
    nv = MESH_Num_Vertices(mesh);
    ne = MESH_Num_Edges(mesh);
    nf = MESH_Num_Faces(mesh);
    nr = MESH_Num_Regions(mesh);

    /* some other known quantitites - 5 items per edge (2 for verts
       and 3 for extra data), MAXPV2+4 items per face (1 for number of
       edges, MAXPV2 for edge indices, anad 3 for extra data),
       MAXPF3+4 items per region (1 for number of faces, MAXPF3 for
       face indices and 3 for extra data */

    nevs = (2+3)*ne;  
    nfes = (1 + MAXPV2 + 3)*nf;
    nrfs = (1 + MAXPF3 + 3)*nr;
    

    mesh_info[0] = rtype;
    mesh_info[1] = nv;
    mesh_info[2] = ne;
    mesh_info[3] = nf;
    mesh_info[4] = nr;
    mesh_info[5] = nevs;
    mesh_info[6] = nfes;
    mesh_info[7] = nrfs;
    mesh_info[8] = natt;
    mesh_info[9] = nset;

    /* send mesh_info */

    MPI_Isend(mesh_info,10,MPI_INT,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;;


    /* collect data */

    int *list_vertex = (int *) malloc(3*nv*sizeof(int));
    double *list_coor = (double *) malloc(3*nv*sizeof(double));
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

    /* send vertices */

    MPI_Isend(list_vertex,3*nv,MPI_INT,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;
    MPI_Isend(list_coor,3*nv,MPI_DOUBLE,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;



    int *list_edge = (int *) malloc(5*ne*sizeof(int));

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

    /* send detailed edge info */

    MPI_Isend(list_edge,nevs,MPI_INT,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;
  


    int *list_face = (int *) malloc(nfes*sizeof(int));

    nfes = 0;

    /* first store nfe, then the edge ids, then the 3 auxilliary data fields */

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


    /* send detailed face info */

    MPI_Isend(list_face,nfes,MPI_INT,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;


    int *list_region = (int *) malloc(nrfs*sizeof(int));

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

    /* send detailed region info */

    MPI_Isend(list_region,nrfs,MPI_INT,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;




    /* collect attrs */

    if (natt) {
      list_attr = (int *) malloc((natt+1)*sizeof(int));
      list_attr_names = (char *) malloc((natt+1)*256*sizeof(char));
      for(i = 0; i < natt; i++) {
        attrib = MESH_Attrib(mesh,i);
        MAttrib_Get_Name(attrib,attname);
        att_type = MAttrib_Get_Type(attrib);
        ncomp = MAttrib_Get_NumComps(attrib);
        mtype = MAttrib_Get_EntDim(attrib);
        list_attr[i] = (ncomp << 6) | (mtype << 3) | (att_type);
        strcpy(&list_attr_names[i*256],attname);
      }


      /* send attr */
      MPI_Isend(list_attr,natt,MPI_INT,torank,torank,comm,&mpirequest);
      (*requests)[*numreq] = mpirequest;
      (*numreq)++;
      MPI_Isend(list_attr_names,natt*256,MPI_CHAR,torank,torank,comm,&mpirequest);
      (*requests)[*numreq] = mpirequest;
      (*numreq)++;
    }

    /* Mesh entity sets */

    if (nset) {
      list_mset = (int *) malloc(nset*sizeof(int));
      list_mset_names = (char *) malloc(nset*256*sizeof(char));

      for (i = 0; i < nset; i++) {
        mset = MESH_MSet(mesh,i);
        MSet_Name(mset,msetname);
        mtype = MSet_EntDim(mset);
        list_mset[i] = mtype;
        strcpy(&list_mset_names[i*256],msetname);
      }

      MPI_Isend(list_mset,nset,MPI_INT,torank,torank,comm,&mpirequest);
      (*requests)[*numreq] = mpirequest;
      (*numreq)++;
      MPI_Isend(list_mset_names,nset*256,MPI_CHAR,torank,torank,comm,
                &mpirequest);
      (*requests)[*numreq] = mpirequest;
      (*numreq)++;
    }

    /* collect allocated memory so it can be freed in a higher level
       routine after MPI_Waitall or MPI_Test has ensured that the send
       has been completed */

    if (ptrs2free == NULL) 
      MSTK_Report("MESH_Surf_SendMesh_FN","ptrs2free array is NULL",MSTK_FATAL);

    int nptrs = 6;
    if (natt) nptrs += 2;
    if (nset) nptrs += 2;

    if (*maxptrs2free == 0) {
      *maxptrs2free = 25;
      *ptrs2free = (void **) malloc(*maxptrs2free*sizeof(void *));
      *numptrs2free = 0;
    }
    else if (*maxptrs2free < (*numptrs2free) + nptrs) {
      *maxptrs2free = 2*(*maxptrs2free) + nptrs;
      *ptrs2free = (void **) realloc(*ptrs2free,(*maxptrs2free)*sizeof(void *));
    }

    (*ptrs2free)[(*numptrs2free)++] = list_vertex;
    (*ptrs2free)[(*numptrs2free)++] = list_coor;  
    (*ptrs2free)[(*numptrs2free)++] = list_edge;
    (*ptrs2free)[(*numptrs2free)++] = list_face;
    (*ptrs2free)[(*numptrs2free)++] = list_region;
    if (natt) {
      (*ptrs2free)[(*numptrs2free)++] = list_attr;
      (*ptrs2free)[(*numptrs2free)++] = list_attr_names;
    }
    if (nset) {
      (*ptrs2free)[(*numptrs2free)++] = list_mset;
      (*ptrs2free)[(*numptrs2free)++] = list_mset_names;
    }


    return 1;
  }




  int MESH_Surf_SendMesh_R1R2R4(Mesh_ptr mesh, int torank, MSTK_Comm comm,
                                int *numreq, int *maxreq, 
                                MPI_Request **requests,
                                int *numptrs2free, int *maxptrs2free,
                                void ***ptrs2free) {
    MSTK_Report("MESH_Surf_SendMesh_R1R2R4","Not implemented",MSTK_FATAL);
    return 0;
  }


  int MESH_Vol_SendMesh_R1R2(Mesh_ptr mesh, int torank, MSTK_Comm comm,
                             int *numreq, int *maxreq, 
                             MPI_Request **requests,
                             int *numptrs2free, int *maxptrs2free,
                             void ***ptrs2free) {
    MSTK_Report("MESH_Vol_SendMesh_R4","Not implemented",MSTK_FATAL);
    return 0;
  }


  int MESH_Vol_SendMesh_R4(Mesh_ptr mesh, int torank, MSTK_Comm comm,
                           int *numreq, int *maxreq, 
                           MPI_Request **requests,
                           int *numptrs2free, int *maxptrs2free,
                           void ***ptrs2free) {
    MSTK_Report("MESH_Vol_SendMesh_R4","Not implemented",MSTK_FATAL);
    return 0;
  }





  
#ifdef __cplusplus
}
#endif

