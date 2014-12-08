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


  int MESH_Send_NonVertexEntities_FN(Mesh_ptr mesh, int torank, MSTK_Comm comm, 
                                     int *numreq, int *maxreq, 
                                     MPI_Request **requests,
                                     int *numptrs2free, int *maxptrs2free,
                                     void ***ptrs2free);

  int MESH_Send_NonVertexEntities_R1R2(Mesh_ptr mesh, int torank, MSTK_Comm comm, 
                                       int *numreq, int *maxreq, 
                                       MPI_Request **requests,
                                       int *numptrs2free, int *maxptrs2free,
                                       void ***ptrs2free);

  int MESH_Send_NonVertexEntities_R4(Mesh_ptr mesh, int torank, MSTK_Comm comm, 
                                     int *numreq, int *maxreq, 
                                     MPI_Request **requests,
                                     int *numptrs2free, int *maxptrs2free,
                                     void ***ptrs2free);

  static int (*MESH_Send_NonVertexEntities_jmp[MSTK_MAXREP])(Mesh_ptr mesh, 
                                                             int torank, 
                                                             MSTK_Comm comm, 
                                                             int *numreq,
                                                             int *maxreq, 
                                                             MPI_Request **requests,
                                                             int *numptrs2free, 
                                                             int *maxptrs2free,
                                                             void ***ptrs2free) = 
  {MESH_Send_NonVertexEntities_FN, MESH_Send_NonVertexEntities_FN, 
   MESH_Send_NonVertexEntities_R1R2, MESH_Send_NonVertexEntities_R1R2, 
   MESH_Send_NonVertexEntities_R4};

  int MESH_Send_NonVertexEntities(Mesh_ptr mesh, int fromrank, MSTK_Comm comm,
                                  int *numreq, int *maxreq, 
                                  MPI_Request **requests,
                                  int *numptrs2free, int *maxptrs2free,
                                  void ***ptrs2free) {
    RepType rtype = MESH_RepType(mesh);
    return MESH_Send_NonVertexEntities_jmp[rtype](mesh,fromrank,comm,
                                                  numreq,maxreq,requests,
                                                  numptrs2free,maxptrs2free,
                                                  ptrs2free);
  }


  int MESH_Send_NonVertexEntities_FN(Mesh_ptr mesh, int torank, MSTK_Comm comm,
                           int *numreq, int *maxreq, MPI_Request **requests,
                           int *numptrs2free, int *maxptrs2free,
                           void ***ptrs2free) {
    int i, j, nv, ne, nf, nr;
    int nevs, nfes, nrfs, nfe, nrv, nrf, dir;
    int maxnfe, maxnrf;
    int *mesh_info;
    int *list_edge=NULL, *list_face=NULL, *list_region=NULL;
    MVertex_ptr mv;
    MEdge_ptr me;
    MFace_ptr mf;
    MRegion_ptr mr;
    List_ptr mfedges, mrfaces, mrverts;
    RepType rtype;
    double coor[3];
    MPI_Request mpirequest;

    if (requests == NULL)
      MSTK_Report("MESH_Surf_SendMesh","MPI requests array is NULL",MSTK_FATAL);

    if (*maxreq == 0) {
      *maxreq = 25;
      *requests = (MPI_Request *) malloc(*maxreq*sizeof(MPI_Request));
      *numreq = 0;
    }
    else if (*maxreq < (*numreq) + 13) {
      *maxreq = 2*(*maxreq) + 11;
      *requests = (MPI_Request *) realloc(*requests,*maxreq*sizeof(MPI_Request));
    }
  

    ne = MESH_Num_Edges(mesh);
    nf = MESH_Num_Faces(mesh);
    nr = MESH_Num_Regions(mesh);

    /* some other known quantitites - 5 items per edge (2 for verts
       and 3 for extra data), maxnfe+4 items per face (1 for number of
       edges, maxnfe for edge indices, anad 3 for extra data),
       maxnrf+4 items per region (1 for number of faces, maxnrf for
       face indices and 3 for extra data */


    maxnfe = 0;
    for (i = 0; i < nf; i++) {
      mf = MESH_Face(mesh,i);
      nfe = MF_Num_Edges(mf);
      if (nfe > maxnfe)
        maxnfe = nfe;
    }

    maxnrf = 0;
    for (i = 0; i < nr; i++) {
      mr = MESH_Region(mesh,i);
      nrf = MR_Num_Faces(mr);
      if (nrf > maxnrf)
        maxnrf = nrf;
    }

    // The amount of extra info we are sending and their meaning is obviously
    // known on the receiving side too. So nevs, nfes and nrfs can be 
    // calculated without us sending it


    nevs = (2+3)*ne;    
    nfes = (1 + maxnfe + 3)*nf;
    nrfs = (1 + maxnrf + 3)*nr;
    
    /* Reserve nevs spots for each edge */

    list_edge = (int *) malloc(5*ne*sizeof(int));

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
  

    /* Reserve nfes spots for each face */

    list_face = (int *) malloc(nfes*sizeof(int));

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

    
    if (nr) {

      list_region = (int *) malloc(nrfs*sizeof(int));
      
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
      
    }
      

    /* collect allocated memory so it can be freed in a higher level
       routine after MPI_Waitall or MPI_Test has ensured that the send
       has been completed */

    if (ptrs2free == NULL) 
      MSTK_Report("MESH_Surf_SendMesh_FN","ptrs2free array is NULL",MSTK_FATAL);

    int nptrs = 3;

    if (*maxptrs2free == 0) {
      *maxptrs2free = 25;
      *ptrs2free = (void **) malloc(*maxptrs2free*sizeof(void *));
      *numptrs2free = 0;
    }
    else if (*maxptrs2free < (*numptrs2free) + nptrs) {
      *maxptrs2free = 2*(*maxptrs2free) + nptrs;
      *ptrs2free = (void **) realloc(*ptrs2free,(*maxptrs2free)*sizeof(void *));
    }

    if (ne)
      (*ptrs2free)[(*numptrs2free)++] = list_edge;
    if (nf)
      (*ptrs2free)[(*numptrs2free)++] = list_face;
    if (nr)
      (*ptrs2free)[(*numptrs2free)++] = list_region;

    return 1;
  }


  int MESH_Send_NonVertexEntities_R1R2(Mesh_ptr mesh, int torank, 
                                       MSTK_Comm comm,
                                       int *numreq, int *maxreq, 
                                       MPI_Request **requests,
                                       int *numptrs2free, int *maxptrs2free,
                                       void ***ptrs2free) {
    MSTK_Report("MESH_Send_NonVertexEntities_R1R2","Not implemented",MSTK_FATAL);
    return 0;
  }


  int MESH_Send_NonVertexEntities_R4(Mesh_ptr mesh, int torank, 
                                     MSTK_Comm comm,
                                     int *numreq, int *maxreq, 
                                     MPI_Request **requests,
                                     int *numptrs2free, int *maxptrs2free,
                                     void ***ptrs2free) {
    MSTK_Report("MESH_Send_NonVertexEntities_R4","Not implemented",MSTK_FATAL);
    return 0;
  }


  
#ifdef __cplusplus
}
#endif

