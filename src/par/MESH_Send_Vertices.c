/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

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



  int MESH_Send_Vertices(Mesh_ptr mesh, int torank, MSTK_Comm comm, 
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
    int *list_vertex = (int *) malloc(3*nv*sizeof(int));

    /* Store the 3 auxilliary data fields - would be nice if we didn't
     * have to send the data in such an error prone way (with bit
     * shifting) that requires knowledge of the internal structure of
     * MEntity */
    
    for(i = 0; i < nv; i++) {
      mv = MESH_Vertex(mesh,i);
      list_vertex[3*i] = (MV_GEntID(mv)<<3) | (MV_GEntDim(mv));
      list_vertex[3*i+1] = (MV_MasterParID(mv) <<3) | MV_OnParBoundary(mv)<<2 | (MV_PType(mv));
      list_vertex[3*i+2] = MV_GlobalID(mv);
    }

    /* send vertices */
    MPI_Isend(list_vertex,3*nv,MPI_INT,torank,torank,comm,&mpirequest);
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

    (*ptrs2free)[(*numptrs2free)++] = list_vertex;

    return 1;
  }
  
#ifdef __cplusplus
}
#endif

