#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Receive one mesh set - the assumption is that this mesh set has
     been created with the right type already based on some meta data
     and that the sending processor has sent this particular mset
     THERE IS NO CHECKING WHETHER THE MESH SET BEING RECEIVED IS THE
     CORRECT ONE

     Authors: Rao Garimella
  */
              

  int MESH_Recv_MSet(Mesh_ptr mesh, MSet_ptr mset, int fromrank, 
                     MSTK_Comm comm) {
    int i, idx, gid, count;
    int num, nent;
    int *list_info, *list_value_int;
    MType mtype, enttype;
    MEntity_ptr ment;
    MVertex_ptr vtx;
    MEdge_ptr edge;
    MFace_ptr face;
    MRegion_ptr rgn;
    MPI_Status status;
    MPI_Request request;
    int result;
    
    int rank;
    MPI_Comm_rank(comm,&rank);
    
    MESH_Enable_GlobalIDSearch(mesh); /* no harm in calling repeatedly */
    
    /* get set properties */
    mtype = MSet_EntDim(mset);
    
    /* receive info */
    list_info = (int *)malloc(sizeof(int));
    result = MPI_Recv(list_info,1,MPI_INT,fromrank,rank,comm,&status); 
    if (result != MPI_SUCCESS)
      MSTK_Report("MESH_RecvMSet","Error receiving mesh set info",MSTK_FATAL);
    
    nent = list_info[0];
    
    if (!nent) {
      free(list_info);
      return 1;
    }
    
    num = 2*nent;
    list_value_int = (int *)malloc(num*sizeof(int));
    
    result = MPI_Recv(list_value_int,num,MPI_INT,fromrank,rank,comm,&status);
    if (result != MPI_SUCCESS)
      MSTK_Report("MESH_RecvMSet","Error receiving mesh set info",MSTK_FATAL);
    
    
    for (i = 0; i < nent; i++) {
      enttype = list_value_int[2*i];
      gid = list_value_int[2*i+1];
      
      ment = 0;
      idx = 0;
      switch (enttype) {
      case MVERTEX:
        ment = MESH_VertexFromGlobalID(mesh,gid);
        break;
      case MEDGE:
        ment = MESH_EdgeFromGlobalID(mesh,gid);
        break;
      case MFACE:
        ment = MESH_FaceFromGlobalID(mesh,gid);
        break;
      case MREGION:
        ment = MESH_RegionFromGlobalID(mesh,gid);
        break;
      case MALLTYPE:
        MSTK_Report("MESH_RecvMSet","An entity cannot be of type MALLTYPE",MSTK_ERROR);
        break;
      default:
        MSTK_Report("MESH_RecvMSet","Unrecognized entity type",MSTK_ERROR);
        break;
      }
      
      if (!ment) {
        MSTK_Report("MESH_RecvMSet","Cannot find entity with this global ID",MSTK_ERROR);
        continue;
      }
      
      MSet_Add(mset,ment);
    }
    
    free(list_info);
    free(list_value_int);
    
    return 1;
  }
  
#ifdef __cplusplus
}
#endif



