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


  int MESH_Recv_VertexCoords(Mesh_ptr mesh, int fromrank, int nvertices, 
                             MSTK_Comm comm) {
    int i;
    MVertex_ptr v;
    double coor[3];
    MPI_Status status;
    MPI_Request vrequest[2];
    char mesg[256], errorstr[256], funcname[256]="MESH_Recv_VertexCoords";
    int errcode, len, nreq=0;

    int rank;
    MPI_Comm_rank(comm,&rank);

    /* allocate receive buffer */
    double *list_coor = (double *) malloc(3*nvertices*sizeof(double));

    errcode = MPI_Recv(list_coor,3*nvertices,MPI_DOUBLE,fromrank,rank,comm,
                        &status);
    if (errcode != MPI_SUCCESS)
      MSTK_Report(funcname,"Trouble receiving mesh coordinate info",MSTK_FATAL);
    

    for(i = 0; i < nvertices; i++) {
      v = MESH_Vertex(mesh,i);
      MV_Set_Coords(v,list_coor+3*i);
    }

    free(list_coor);

    return 1;
  }


#ifdef __cplusplus
}
#endif

