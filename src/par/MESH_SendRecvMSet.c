#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* 
     this function send entity sets

     Used only when distributing mesh
     caller should be the root processor, and recv processor call MESH_RecvAttr()
     
     set_name: the name of the attribute 
  */

  int MESH_SendMSet(Mesh_ptr mesh, const char *mset_name, int torank, MSTK_Comm comm) {
  int i, idx;
  int num, nent;
  int *list_info, *list_value_int;

  MType mtype;
  MEntity_ptr ment;
  MSet_ptr mset;


  mset = MESH_MSetByName(mesh,mset_name);
  if(!mset) {
    MSTK_Report("MESH_SendMSet","Mesh has no such entity set",MSTK_ERROR);
    return 0;
  }

  mtype = MSet_EntDim(mset);
  nent = MSet_Num_Entries(mset);

  list_info = (int *)MSTK_malloc(sizeof(int));
  list_info[0] = nent;

  MPI_Send(list_info,1,MPI_INT,torank,torank,comm);

  if (!nent) {
    MSTK_free(list_info);
    return 1;
  }

  /* attribute index and global id */ 
  num = 2*nent;
  list_value_int = (int *)MSTK_malloc(num*sizeof(int));

  idx = 0; i = 0;
  while ((ment = MSet_Next_Entry(mset,&idx))) {
    list_value_int[i++] = MEnt_Dim(ment);
    list_value_int[i++] = MEnt_GlobalID(ment);
  }
  

  /* send info */
  MPI_Send(list_value_int,num,MPI_INT,torank,torank,comm);

  /* release the send buffer */
  MSTK_free(list_info);
  MSTK_free(list_value_int);
  
  return 1;
}



  /* 
     This function receives an entity set
     Called by slave processor.

     set_name: name of the entity set
     fromrank: rank of the send processor
     rank: rank of the receiving processor

  */

  int MESH_RecvMSet(Mesh_ptr mesh, const char *mset_name, int fromrank, MSTK_Comm comm) {
  int i, idx, gid, count;
  int num, nent;
  int *list_info, *list_value_int;
  MType mtype, enttype;
  MEntity_ptr ment;
  MVertex_ptr vtx;
  MEdge_ptr edge;
  MFace_ptr face;
  MRegion_ptr rgn;
  MSet_ptr mset;
  MPI_Status status;


  int rank;
  MPI_Comm_rank(comm,&rank);

  MESH_Enable_GlobalIDSearch(mesh); /* no harm in calling repeatedly */
  
  mset = MESH_MSetByName(mesh,mset_name);

  if (!mset) {
    MSTK_Report("MESH_RecvMSet","Mesh contains no such set",MSTK_ERROR);
    return 0;
  }

  /* get set properties */
  mtype = MSet_EntDim(mset);

  /* receive info */
  list_info = (int *)MSTK_malloc(sizeof(int));
  MPI_Recv(list_info,1,MPI_INT,fromrank,rank,comm,&status);
  MPI_Get_count(&status,MPI_INT,&count);

  assert(1==count);

  nent = list_info[0];

  if (!nent) {
    MSTK_free(list_info);
    return 1;
  }

  num = 2*nent;
  list_value_int = (int *)MSTK_malloc(num*sizeof(int));

  MPI_Recv(list_value_int,num,MPI_INT,fromrank,rank,comm, &status);
  MPI_Get_count(&status,MPI_INT,&count);
  assert(num==count);


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
  
  MSTK_free(list_info);
  MSTK_free(list_value_int);

  return 1;
}
  
#ifdef __cplusplus
}
#endif



