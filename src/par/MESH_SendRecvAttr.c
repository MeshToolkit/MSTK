#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "MSTK.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* 
     this function send attributes

     Used only when distributing mesh, for communication, use MESH_UpdateAttr()
     caller should be the root processor, and recv processor call MESH_RecvAttr()
     
     attr_name: the name of the attribute 

     Author(s): Rao Garimella
  */
int MESH_SendAttr(Mesh_ptr mesh, const char *attr_name, int torank) {
  int j, k;
  int num, ncomp, ival;
  double rval;
  void *pval;
  double *rval_arr;
  int *list_info;

  MType mtype;
  MAttType att_type;
  MEntity_ptr ment;
  MAttrib_ptr attrib;

  MPI_Comm comm = MSTK_Comm();

  attrib = MESH_AttribByName(mesh,attr_name);
  /* if there is no such attribute */
  if(!attrib) {
    MSTK_Report("MESH_SendAttr","There is no such attribute",MSTK_ERROR);
    printf("attribute name %s\n",attr_name);
    return 0;
  }
  /* get attribute properties */
  att_type = MAttrib_Get_Type(attrib);
  ncomp = MAttrib_Get_NumComps(attrib);
  mtype = MAttrib_Get_EntDim(attrib);
  
  /* attribute entity type, used before ghost list established, so no ghost */
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
  default:
    num = 0;
    MSTK_Report("MESH_SendAttr()","Cannot send attributes on entity type MALLTYPE",MSTK_WARN);
    return 0;
  }
    
  /* attribute index and global id */ 
  list_info = (int *)MSTK_malloc(num*sizeof(int));
  /* attribute values */
  int *list_value_int = (int *)MSTK_malloc(num*ncomp*sizeof(int));
  double *list_value_double = (double *)MSTK_malloc(num*ncomp*sizeof(double));
  
  /* collect data */
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
      MSTK_Report("MESH_SendAttr()","Invalid entity type",MSTK_WARN);
      return 0;
    }
    
    MEnt_Get_AttVal(ment,attrib,&ival,&rval,&pval);

    list_info[j] = MEnt_GlobalID(ment);
    if (att_type == INT)
      list_value_int[j] = ival;
    else {
      if(ncomp == 1)
	list_value_double[j] = rval;
      if(ncomp > 1) {
	rval_arr = (double *)pval;
	for(k = 0; k < ncomp; k++)
	  list_value_double[ncomp*j+k] = rval_arr[k];
      }
    }
  }
  
  /* send info */
  MPI_Send(list_info,num,MPI_INT,torank,torank,comm);
  /* send value */
  if (att_type == INT)
    MPI_Send(list_value_int,num*ncomp,MPI_INT,torank,torank,comm);
  else
    MPI_Send(list_value_double,num*ncomp,MPI_DOUBLE,torank,torank,comm);

  /* release the send buffer */
  MSTK_free(list_info);
  MSTK_free(list_value_int);
  MSTK_free(list_value_double);
  
  return 1;
}


  /* 
     this function receive attribute
     called by slave processor.
     attr_name: name of the attribute
     fromrank: rank of the send processor
     rank: rank of the receiving processor
     
  */

int MESH_RecvAttr(Mesh_ptr mesh, const char *attr_name, int fromrank) {
  int j, k, count;
  int num, ncomp, ival=0;
  double rval=0.0;
  double *rval_arr=NULL;
  int *list_info, *list_value_int;
  double *list_value_double;
  MType mtype;
  MAttType att_type;
  MEntity_ptr ment=NULL;
  MAttrib_ptr attrib;
  MPI_Status status;
  
  MPI_Comm comm = MSTK_Comm();
  int rank = MSTK_Comm_rank();

  attrib = MESH_AttribByName(mesh,attr_name);

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
    MSTK_Report("MESH_SendAttr()","Cannot receive attributes on entity type MALLTYPE",MSTK_WARN);
    return 0;
  default:
    num = 0;
    MSTK_Report("MESH_SendAttr()","Invalid entity type MALLTYPE",MSTK_WARN);
    return 0;
  }

  /* receive info */
  list_info = (int *)MSTK_malloc((num)*sizeof(int));
  MPI_Recv(list_info,num,MPI_INT,fromrank,rank,comm,&status);
  MPI_Get_count(&status,MPI_INT,&count);

  assert((num)==count);
  /* printf("received %d attributes of attribute index %d in MESH_RecvAttr() from rank %d on rank %d\n",count,attr_index,fromrank,rank); */

  list_value_int = (int *)MSTK_malloc((num)*ncomp*sizeof(int));
  list_value_double = (double *)MSTK_malloc(num*ncomp*sizeof(double));
  /* reveive value */
  if (att_type == INT) {
    MPI_Recv(list_value_int,(num)*ncomp,MPI_INT,fromrank,rank,comm, &status);
    MPI_Get_count(&status,MPI_INT,&count);
    assert(num*ncomp==count);
  }
  else {
    MPI_Recv(list_value_double,num*ncomp,MPI_DOUBLE,fromrank,rank,comm, &status);
    MPI_Get_count(&status,MPI_DOUBLE,&count);
    assert(ncomp*num==count);
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
    if (att_type == INT)
      ival = list_value_int[j];
    else {
      if(ncomp == 1) 
	rval = list_value_double[j];
      if(ncomp > 1)
	{
	  rval_arr = (double *)MSTK_malloc(ncomp*sizeof(double));
	  for(k = 0; k < ncomp; k++)
	    rval_arr[k] = list_value_double[j*ncomp+k];
	}	
    }
    MEnt_Set_AttVal(ment,attrib,ival,rval,(void*) rval_arr);
  }
  
  MSTK_free(list_info);
  MSTK_free(list_value_int);
  MSTK_free(list_value_double);

  return 1;
}
  

#ifdef __cplusplus
}
#endif

