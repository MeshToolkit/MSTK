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

     Author(s): Rao Garimella
                Duo Wang
  */


  int MESH_RecvMesh(Mesh_ptr mesh, int fromrank, int with_attr, 
                    MSTK_Comm comm) {
    int i, j, a, m;
    int nv, ne, nf, nr;
    RepType rtype;
    MAttrib_ptr attrib;
    MSet_ptr mset;
    MPI_Status status;
    char mesg[256], errorstr[256], funcname[256]="MESH_RecvMesh";
    int errcode, len, nreq=0;
   
    int rank;
    MPI_Comm_rank(comm,&rank);

    
    /* Receive meta data about mesh */

    MESH_Recv_MetaData(mesh,fromrank,&rtype,&nv,&ne,&nf,&nr,comm);
    
    /* Receive mesh vertex info */

    MESH_Recv_Vertices(mesh,fromrank,nv,comm);


    /* Receive mesh vertex coordinate info */

    MESH_Recv_VertexCoords(mesh,fromrank,nv,comm);


    /* Receive higher dimensional entities */

    MESH_Recv_NonVertexEntities(mesh,fromrank,ne,nf,nr,comm);


    if (with_attr) {

      /* Receive attribute meta data and create attributes */
      
      MESH_Recv_AttributeMetaData(mesh,fromrank,comm);
      
      /* Receive each attribute - assumption is this processor is
         receiving attributes in the same order that the sending
         processor put them out */
      
      int natt_local = MESH_Num_Attribs(mesh);
      
      for (a = 0; a < natt_local; a++) {
        attrib = MESH_Attrib(mesh,a);
        MESH_Recv_Attribute(mesh,attrib,fromrank,comm);
      }
      
      /* Receive mesh set meta data and create mesh sets */
      
      MESH_Recv_MSetMetaData(mesh,fromrank,comm);
      
      /* Receive each mesh set - assumption is this processor is
         receiving mesh sets in the same order that the sending
         processor put them out */
      
      int nset_local = MESH_Num_MSets(mesh);
      
      for (m = 0; m < nset_local; m++) {
        mset = MESH_MSet(mesh,m);
        MESH_Recv_MSet(mesh,mset,fromrank,comm);
      }

    } /* if (with_attr) */

    return 1;
  }


#ifdef __cplusplus
}
#endif

