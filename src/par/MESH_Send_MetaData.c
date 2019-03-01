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
     This function send mesh meta data to processor torank in communicator comm

     Author(s): Rao Garimella
  */


  int MESH_Send_MetaData(Mesh_ptr mesh, int torank, MSTK_Comm comm,
                         int *numreq, int *maxreq, MPI_Request **requests,
                         int *numptrs2free, int *maxptrs2free,
                         void ***ptrs2free) {
    int i, j, nv, ne, nf, nr;
    int nevs, nfes, nrfs, nfe, nrv, nrf, dir;
    int maxnfe, maxnrf;
    int *mesh_info;
    MVertex_ptr mv;
    MEdge_ptr me;
    MFace_ptr mf;
    MRegion_ptr mr;
    List_ptr mfedges, mrfaces, mrverts;
    RepType rtype;
    double coor[3];
    MPI_Request mpirequest;
    char funcname[256] = "MESH_Send_MetaData";

    if (requests == NULL)
      MSTK_Report(funcname,"MPI requests array is NULL",MSTK_FATAL);

    if (*maxreq == 0) {
      *maxreq = 25;
      *requests = (MPI_Request *) malloc(*maxreq*sizeof(MPI_Request));
      *numreq = 0;
    }
    else if (*maxreq < (*numreq) + 1) {
      *maxreq = 2*(*maxreq) + 11;
      *requests = (MPI_Request *) realloc(*requests,*maxreq*sizeof(MPI_Request));
    }
  

    /* mesh_info store the mesh reptype, nv, ne, nf, nr */

    rtype = MESH_RepType(mesh);
    nv = MESH_Num_Vertices(mesh);
    ne = MESH_Num_Edges(mesh);
    nf = MESH_Num_Faces(mesh);
    nr = MESH_Num_Regions(mesh);

    mesh_info = (int *) malloc(5*sizeof(int));

    mesh_info[0] = rtype;
    mesh_info[1] = nv;
    mesh_info[2] = ne;
    mesh_info[3] = nf;
    mesh_info[4] = nr;

    /* send mesh_info */

    MPI_Isend(mesh_info,5,MPI_INT,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;;

    /* collect allocated memory so it can be freed in a higher level
       routine after MPI_Waitall or MPI_Test has ensured that the send
       has been completed */

    if (ptrs2free == NULL) 
      MSTK_Report(funcname,"ptrs2free array is NULL",MSTK_FATAL);

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

    (*ptrs2free)[(*numptrs2free)++] = mesh_info;

    return 1;
  }


#ifdef __cplusplus
}
#endif

