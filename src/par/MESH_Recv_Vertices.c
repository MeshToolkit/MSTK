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
     This function receives mesh vertices from processor 'fromrank'

     Author(s): Rao Garimella
  */


  int MESH_Recv_Vertices(Mesh_ptr mesh, int fromrank, int nvertices, 
                         MSTK_Comm comm) {
    int i;
    MVertex_ptr v;
    MPI_Status status;
    MPI_Request request;
    char mesg[256], errorstr[256], funcname[256]="MESH_Recv_Vertices";
    int errcode, len, nreq=0;

    int rank;
    MPI_Comm_rank(comm,&rank);

    /* allocate receive buffer */
    int *list_vertex = (int *) malloc(3*nvertices*sizeof(int));

    /* receive vertex info */
    errcode = MPI_Irecv(list_vertex,3*nvertices,MPI_INT,fromrank,rank,comm,
                        &request);
    if (errcode != MPI_SUCCESS)
      MSTK_Report(funcname,"Trouble receiving mesh vertex info",MSTK_FATAL);
    

    /* Create the vertices while waiting for the message to complete */

    for (i = 0; i < nvertices; i++)
      v = MV_New(mesh);


    errcode = MPI_Wait(&request,MPI_STATUS_IGNORE);
    if (errcode != MPI_SUCCESS)
      MSTK_Report(funcname,"Trouble receiving mesh vertex info",MSTK_FATAL);    

    for(i = 0; i < nvertices; i++) {
      v = MESH_Vertex(mesh,i);
      int gentdim = list_vertex[3*i] & 7; /* first 3 bits; 7 is 0...00111 */
      int gentid = list_vertex[3*i] >> 3; /* All but the first 3 bits */   
      MV_Set_GEntDim(v,gentdim);
      MV_Set_GEntID(v,gentid);

      int ptype = list_vertex[3*i+1] & 3; /* first 2 bits; 3 is 0...00011 */  
      int on_par_bdry = list_vertex[3*i+1] & 4; /* 3rd bit; 4 is 0...00100 */ 
      int masterparid = list_vertex[3*i+1] >> 3; /* All but the first 3 bits */ 
      MV_Set_PType(v,ptype);
      if (on_par_bdry)
        MV_Flag_OnParBoundary(v);
      MV_Set_MasterParID(v,masterparid);

      MV_Set_GlobalID(v,list_vertex[3*i+2]);
    }

    free(list_vertex);    

    return 1;
  }


#ifdef __cplusplus
}
#endif

