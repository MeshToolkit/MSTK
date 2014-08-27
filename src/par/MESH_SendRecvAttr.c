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
  
  int MESH_SendAttr(Mesh_ptr mesh, const char *attr_name, int torank, 
                    MSTK_Comm comm, int *numreq, int *maxreq, 
                    MPI_Request **requests,
                    int *numptrs2free, int *maxptrs2free,
                    void ***ptrs2free) {
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
    MPI_Request mpirequest;

    if (requests == NULL)
      MSTK_Report("MSTK_SendMSet","Invalid MPI request buffer",MSTK_FATAL);
  
    if (*maxreq == 0) {
      *maxreq = 25;
      *requests = (MPI_Request *) malloc(*maxreq*sizeof(MPI_Request));
      *numreq = 0;
    }
    else if (*maxreq < (*numreq) + 2) {
      *maxreq *= 2;
      *requests = (MPI_Request *) realloc(*requests,*maxreq*sizeof(MPI_Request));
    }
  
    attrib = MESH_AttribByName(mesh,attr_name);
    
    if (!attrib) { /* if there is no such attribute */
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
#ifdef DEBUG2
      MSTK_Report("MESH_SendAttr()","Cannot send attributes on entity type MALLTYPE",MSTK_WARN);
#endif
      return 0;
    }
    
    /* attribute index and global id */ 
    list_info = (int *) malloc(num*sizeof(int));

    /* attribute values */
    int *list_value_int = NULL;
    double *list_value_double = NULL;
    if (att_type == INT)
      list_value_int = (int *) malloc(num*ncomp*sizeof(int));
    else
      list_value_double = (double *) malloc(num*ncomp*sizeof(double));
  
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
    MPI_Isend(list_info,num,MPI_INT,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;

    /* send value */
    if (att_type == INT) {
      MPI_Isend(list_value_int,num*ncomp,MPI_INT,torank,torank,comm,&mpirequest);
      (*requests)[*numreq] = mpirequest;
    }
    else {
      MPI_Isend(list_value_double,num*ncomp,MPI_DOUBLE,torank,torank,comm,
                &mpirequest);
      (*requests)[*numreq] = mpirequest;
    }
    (*numreq)++;

    /* release the send buffer */
    int nptrs = 2;

    if (*maxptrs2free == 0) {
      *maxptrs2free = 25;
      *ptrs2free = (void **) malloc(*maxptrs2free*sizeof(void *));
      *numptrs2free = 0;
    }
    else if (*maxptrs2free < (*numptrs2free) + 3) {
      *maxptrs2free *= 2;
      *ptrs2free = (void **) realloc(*ptrs2free,(*maxptrs2free)*sizeof(void *));
    }

    (*ptrs2free)[(*numptrs2free)++] = list_info;
    if (att_type == INT) 
      (*ptrs2free)[(*numptrs2free)++] = list_value_int;  
    else
      (*ptrs2free)[(*numptrs2free)++] = list_value_double;

    return 1;
  }


  /* 
     this function receive attribute
     called by slave processor.
     attr_name: name of the attribute
     fromrank: rank of the send processor
     rank: rank of the receiving processor
     
  */

  int MESH_RecvAttr(Mesh_ptr mesh, const char *attr_name, int fromrank, MPI_Comm comm) {
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
  MPI_Request request;
  MPI_Status status;
  int result;
  
  int rank;
  MPI_Comm_rank(comm,&rank);

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
#ifdef DEBUG2
    MSTK_Report("MESH_SendAttr()","Cannot receive attributes on entity type MALLTYPE",MSTK_WARN);
#endif
    return 0;
  default:
    num = 0;
    MSTK_Report("MESH_SendAttr()","Invalid entity type MALLTYPE",MSTK_WARN);
    return 0;
  }

  /* receive info */
  list_info = (int *) malloc((num)*sizeof(int));
  result = MPI_Recv(list_info,num,MPI_INT,fromrank,rank,comm,&status);
  if (result != MPI_SUCCESS)
    MSTK_Report("MESH_RecvAttr","Trouble with receiving attributes",MSTK_FATAL);

  /* printf("received %d attributes of attribute index %d in MESH_RecvAttr() from rank %d on rank %d\n",count,attr_index,fromrank,rank); */

  list_value_int = (int *) malloc((num)*ncomp*sizeof(int));
  list_value_double = (double *) malloc(num*ncomp*sizeof(double));
  /* reveive value */
  if (att_type == INT) {
    result = MPI_Recv(list_value_int,(num)*ncomp,MPI_INT,fromrank,rank,comm, &status);
    if (result != MPI_SUCCESS)
      MSTK_Report("MESH_RecvAttr","Trouble with receiving attributes",MSTK_FATAL);    
  }
  else {
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
      if(ncomp > 1)
	{
          for (k = 0; k < ncomp; k++)
            if (list_value_double[j*ncomp+k] != 0) nullval = 0;
          if (!nullval) {
            rval_arr = (double *) malloc(ncomp*sizeof(double));
            for(k = 0; k < ncomp; k++)
              rval_arr[k] = list_value_double[j*ncomp+k];
          }
        }
    }
    if (nullval) continue; /* Don't waste time initializing to zero */

    MEnt_Set_AttVal(ment,attrib,ival,rval,(void*) rval_arr);
  }
  
   free(list_info);
   free(list_value_int);
   free(list_value_double);

  return 1;
}
  

#ifdef __cplusplus
}
#endif

