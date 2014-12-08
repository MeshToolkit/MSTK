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



  int MESH_Send_VertexCoords(Mesh_ptr mesh, int torank, MSTK_Comm comm, 
                         int *numreq, int *maxreq, MPI_Request **requests,
                         int *numptrs2free, int *maxptrs2free,
                         void ***ptrs2free) {
    int i, j, nv;
    MVertex_ptr mv;
    double coor[3];
    MPI_Request mpirequest;
  
    if (requests == NULL)
      MSTK_Report("MESH_Surf_SendMesh","MPI requests array is NULL",MSTK_FATAL);
    
    if (*maxreq == 0) {
      *maxreq = 25;
      *requests = (MPI_Request *) malloc(*maxreq*sizeof(MPI_Request));
      *numreq = 0;
    }
    else if (*maxreq < (*numreq) + 11) {
      *maxreq = 2*(*maxreq) + 11;
      *requests = (MPI_Request *) realloc(*requests,*maxreq*sizeof(MPI_Request));
    }

  
    /* Now send out detailed vertex info */

    nv = MESH_Num_Vertices(mesh);
    double *list_coor = (double *) malloc(3*nv*sizeof(double));

    /* Store the 3 auxilliary data fields */
    for(i = 0; i < nv; i++) {
      mv = MESH_Vertex(mesh,i);
      MV_Coords(mv,list_coor+i*3);
    }

    MPI_Isend(list_coor,3*nv,MPI_DOUBLE,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;;


    int nptrs = 1;

    if (*maxptrs2free == 0) {
      *maxptrs2free = 25;
      *ptrs2free = (void **) malloc(*maxptrs2free*sizeof(void *));
      *numptrs2free = 0;
    }
    else if (*maxptrs2free < (*numptrs2free) + nptrs) {
      *maxptrs2free = 2*(*maxptrs2free) + nptrs;
      *ptrs2free = (void **) realloc(*ptrs2free,(*maxptrs2free)*sizeof(void *));
    }

    (*ptrs2free)[(*numptrs2free)++] = list_coor;  

    return 1;
  }
  
#ifdef __cplusplus
}
#endif

