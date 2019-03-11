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
     This function receives mesh from processor rank in communicator comm
     attrib list is received, but no attrib values.
     call MESH_RecvAttr() to update entity attribute values
     fromrank: the rank of sending processor
     rank: the rank of receiving processor

     Author(s): Duo Wang, Rao Garimella
  */

 
  int MESH_Recv_MetaData(Mesh_ptr mesh, int fromrank, RepType *rtype,
                             int *nv, int *ne, int *nf, int *nr,
                             MSTK_Comm comm) {
    int mesh_info[5], count;
    MPI_Request request;
    MPI_Status status;
    char mesg[256], errorstr[256];
    int len, errcode;

    int rank;
    MPI_Comm_rank(comm,&rank);

    /* mesh_info store the mesh reptype, nv, nf, nfvs */
    /* receive mesh_info */

    errcode = MPI_Recv(mesh_info,5,MPI_INT,fromrank,rank,comm,&status);
    if (errcode != MPI_SUCCESS) {
      MPI_Error_string(errcode,errorstr,&len);
      sprintf(mesg,"MPI_Recv failed with error %s",errorstr);
      MSTK_Report("MESH_Recv_MetaData",mesg,MSTK_FATAL);
    }
    
    *rtype = mesh_info[0];
    MESH_SetRepType(mesh,*rtype);
    
    *nv = mesh_info[1];
    *ne = mesh_info[2];
    *nf = mesh_info[3];
    *nr = mesh_info[4];
    
    return 1;
  }
  

#ifdef __cplusplus
}
#endif

