/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "MSTK.h"

#ifdef __cplusplus
extern "C" {
#endif


  /* 
     Receive one attribute - the assumption is that this mesh set has
     been created with the right type already based on some meta data
     and that the sending processor has sent this particular mset
     THERE IS NO CHECKING WHETHER THE MESH SET BEING RECEIVED IS THE
     CORRECT ONE

     Author: Rao Garimella
  */

  int MESH_Recv_Attribute(Mesh_ptr mesh, MAttrib_ptr attrib, int fromrank, 
                          MPI_Comm comm) {
  int j, k, count;
  int num, ncomp, ival=0;
  double rval=0.0;
  double *rval_arr=NULL;
  int *list_info=NULL, *list_value_int=NULL;
  double *list_value_double=NULL;
  MType mtype;
  MAttType att_type;
  MEntity_ptr ment=NULL;
  MPI_Request request;
  MPI_Status status;
  int result;
  
  int rank;
  MPI_Comm_rank(comm,&rank);


  /* get attribute properties */

  att_type = MAttrib_Get_Type(attrib);
  ncomp = MAttrib_Get_NumComps(attrib);
  mtype = MAttrib_Get_EntDim(attrib);

  /* attribute entity type */

  switch (mtype) {
  case MVERTEX:
    num = MESH_Num_Vertices(mesh);
    break;
  case MEDGE:
    num = MESH_Num_Edges(mesh);
    break;
  case MFACE:
    num = MESH_Num_Faces(mesh);
    break;
  case MREGION:
    num = MESH_Num_Regions(mesh);
    break;
  case MALLTYPE:
    num = 0;
#ifdef DEBUG2
    MSTK_Report("MESH_RecvAttr()","Cannot receive attributes on entity type MALLTYPE",MSTK_WARN);
#endif
    return 0;
  default:
    num = 0;
    MSTK_Report("MESH_RecvAttr()","Invalid entity type MALLTYPE",MSTK_WARN);
    return 0;
  }

  /* receive entity global IDs */

  list_info = (int *) malloc((num)*sizeof(int));
  result = MPI_Recv(list_info,num,MPI_INT,fromrank,rank,comm,&status);
  if (result != MPI_SUCCESS)
    MSTK_Report("MESH_RecvAttr","Trouble with receiving attributes",MSTK_FATAL);


  /* receive values */

  if (att_type == INT) {
    list_value_int = (int *) malloc((num)*ncomp*sizeof(int));
    result = MPI_Recv(list_value_int,(num)*ncomp,MPI_INT,fromrank,rank,comm, &status);
    if (result != MPI_SUCCESS)
      MSTK_Report("MESH_RecvAttr","Trouble with receiving attributes",MSTK_FATAL);    
  }
  else {
    list_value_double = (double *) malloc(num*ncomp*sizeof(double));
    result = MPI_Recv(list_value_double,num*ncomp,MPI_DOUBLE,fromrank,rank,comm, &status);
    if (result != MPI_SUCCESS)
      MSTK_Report("MESH_RecvAttr","Trouble with receiving attributes",MSTK_FATAL);    
  }

  /* associate attributes with entities */  
  for(j = 0; j < num; j++) {
    switch (mtype) {
    case MVERTEX:
      ment = MESH_Vertex(mesh,j);
      break;
    case MEDGE:
      ment = MESH_Edge(mesh,j);
      break;
    case MFACE:
      ment = MESH_Face(mesh,j);
      break;
    case MREGION:
      ment = MESH_Region(mesh,j);
      break;
    default:
      MSTK_Report("MESH_RecvAttr()","Invalid entity type",MSTK_FATAL);
    }

    assert(MEnt_GlobalID(ment) == list_info[j]);
    
    int nullval = 1;
    if (att_type == INT) {
      ival = list_value_int[j];
      if (ival != 0) nullval = 0;
    }
    else if (att_type == DOUBLE) {
      rval = list_value_double[j];
      if (rval != 0.0) nullval = 0;
    }
    else {
      if (ncomp > 1) {
        nullval = 0;  /* always add vectors */
        rval_arr = (double *) malloc(ncomp*sizeof(double));
        for(k = 0; k < ncomp; k++)
          rval_arr[k] = list_value_double[j*ncomp+k];
      }
    }
    if (nullval) continue; /* Don't waste time initializing to zero */

    MEnt_Set_AttVal(ment,attrib,ival,rval,(void*) rval_arr);
  }
  
   free(list_info);
   if (att_type == INT) 
     free(list_value_int);
   else if (att_type == DOUBLE)
     free(list_value_double);

  return 1;
}
  

#ifdef __cplusplus
}
#endif

