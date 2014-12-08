#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Send a entity sets of a mesh to a processor 'torank' 

     Authors: Rao Garimella
  */

  int MESH_Send_MSetMetaData(Mesh_ptr mesh, int torank, MSTK_Comm comm,
                             int *numreq, int *maxreq, MPI_Request **requests,
                             int *numptrs2free, int *maxptrs2free, 
                             void ***ptrs2free) {
    int i, nset, mtype;
    char msetname[256];
    MSet_ptr mset;
    int *list_mset_num, *list_mset_types;
    char *list_mset_names;
    MPI_Request mpirequest;
    
    if (requests == NULL)
      MSTK_Report("MSTK_SendMSets","Invalid MPI request buffer",MSTK_FATAL);
    
    if (*maxreq == 0) {
      *maxreq = 25;
      *requests = (MPI_Request *) malloc(*maxreq*sizeof(MPI_Request));
      *numreq = 0;
    }
    else if (*maxreq < (*numreq) + 3) {
      *maxreq = 2*(*maxreq) + 3;
      *requests = (MPI_Request *) realloc(*requests,*maxreq*sizeof(MPI_Request));
    }
    
    if (ptrs2free == NULL)
      MSTK_Report("MSTK_SendMSets","Invalid ptrs2free buffer",MSTK_FATAL);

    if (*maxptrs2free == 0) {
      *maxptrs2free = 25;
      *ptrs2free = (void **) malloc(*maxptrs2free*sizeof(void *));
      *numptrs2free = 0;
    }
    else if (*maxptrs2free < (*numptrs2free) + 3) {
      *maxptrs2free *= 2*(*maxptrs2free) + 3;
      *ptrs2free = (void **) realloc(*ptrs2free,(*maxptrs2free)*sizeof(void *));
    }
    
    
    /* Send info about how many mesh sets there are, what kind and
       what there names are in a packed fashion */

    nset = MESH_Num_MSets(mesh);
    list_mset_num = (int *) malloc(sizeof(int));
    list_mset_num[0] = nset;

    MPI_Isend(list_mset_num,1,MPI_INT,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;


    (*ptrs2free)[(*numptrs2free)++] = list_mset_num;

    if (!nset) return 1;



    list_mset_types = (int *) malloc(nset*sizeof(int));
    list_mset_names = (char *) malloc(nset*256*sizeof(char));
    
    for(i = 0; i < nset; i++) {
      mset = MESH_MSet(mesh,i);
      MSet_Name(mset,msetname);
      mtype = MSet_EntDim(mset);
      list_mset_types[i] = mtype;
      strcpy(&(list_mset_names[i*256]),msetname);
    }

    MPI_Isend(list_mset_types,nset,MPI_INT,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;
    
    MPI_Isend(list_mset_names,nset*256,MPI_CHAR,torank,torank,comm,
              &mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;



    (*ptrs2free)[(*numptrs2free)++] = list_mset_types;
    (*ptrs2free)[(*numptrs2free)++] = list_mset_names;

    return 1;
  }


#ifdef __cplusplus
}
#endif

