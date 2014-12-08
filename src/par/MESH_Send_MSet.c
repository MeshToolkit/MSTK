#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Send one mesh set to a particular processor 

     Author: Rao Garimella
  */
  

  int MESH_Send_MSet(Mesh_ptr mesh, MSet_ptr mset, int torank, 
                     MSTK_Comm comm, int *numreq, int *maxreq,
                     MPI_Request **requests, int *numptrs2free, 
                     int *maxptrs2free, void ***ptrs2free) {
    int i, idx;
    int num, nent;
    int *list_info, *list_value_int;
    MType mtype;
    MEntity_ptr ment;
    MPI_Request mpirequest;

    if (requests == NULL)
      MSTK_Report("MSTK_SendMSet","Invalid MPI request buffer",MSTK_FATAL);
  
    if (*maxreq == 0) {
      *maxreq = 10;
      *requests = (MPI_Request *) malloc(*maxreq*sizeof(MPI_Request));
      *numreq = 0;
    }
    else if (*maxreq < (*numreq) + 2) {
      *maxreq *= 2;
      *requests = (MPI_Request *) realloc(*requests,*maxreq*sizeof(MPI_Request));
    }
  
    mtype = MSet_EntDim(mset);
    nent = MSet_Num_Entries(mset);

    list_info = (int *)malloc(sizeof(int));
    list_info[0] = nent;

    MPI_Isend(list_info,1,MPI_INT,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;

    if (!nent) {
      free(list_info);
      return 1;
    }

    /* attribute index and global id */ 
    num = 2*nent;
    list_value_int = (int *)malloc(num*sizeof(int));

    idx = 0; i = 0;
    while ((ment = MSet_Next_Entry(mset,&idx))) {
      list_value_int[i++] = MEnt_Dim(ment);
      list_value_int[i++] = MEnt_GlobalID(ment);
    }
  

    /* send info */
    MPI_Isend(list_value_int,num,MPI_INT,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;

    /* track the buffers used for sending */

    if (*maxptrs2free == 0) {
      *maxptrs2free = 25;
      *ptrs2free = (void **) malloc(*maxptrs2free*sizeof(void *));
      *numptrs2free = 0;
    }
    else if (*maxptrs2free < (*numptrs2free) + 2) {
      *maxptrs2free = 2*(*maxptrs2free) + 2 ;
      *ptrs2free = (void **) realloc(*ptrs2free,(*maxptrs2free)*sizeof(void *));
    }

    (*ptrs2free)[(*numptrs2free)++] = list_info;
    (*ptrs2free)[(*numptrs2free)++] = list_value_int;  
  
    return 1;
  }


  
#ifdef __cplusplus
}
#endif



