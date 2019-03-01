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

  /* Send a list of all mesh attributes including their name, type, type of
     entity they live on etc. to a processor 'torank'

     Authors: Rao Garimella
  */

  int MESH_Send_AttributeMetaData(Mesh_ptr mesh, int torank, MSTK_Comm comm,
                                  int *numreq, int *maxreq, 
                                  MPI_Request **requests,
                                  int *numptrs2free, int *maxptrs2free, 
                                  void ***ptrs2free) {
    int i, natt, mtype, att_type, ncomp;
    char attname[256];
    MAttrib_ptr attrib;
    int *list_attr_num, *list_attr_info;
    char *list_attr_names;
    MPI_Request mpirequest;
    
    if (requests == NULL)
      MSTK_Report("MSTK_SendMesh","Invalid MPI request buffer",MSTK_FATAL);
    
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
      MSTK_Report("MSTK_SendMesh","Invalid ptrs2free buffer",MSTK_FATAL);
    
    int nptrs = 3;
    if (*maxptrs2free == 0) {
      *maxptrs2free = 25;
      *ptrs2free = (void **) malloc(*maxptrs2free*sizeof(void *));
      *numptrs2free = 0;
    }
    else if (*maxptrs2free < (*numptrs2free) + nptrs) {
      *maxptrs2free = 2*(*maxptrs2free) + 3;
      *ptrs2free = (void **) realloc(*ptrs2free,(*maxptrs2free)*sizeof(void *));
    }    

    /* Send info about how many attributes there are, what kind and
       what there names are in a packed fashion */

    natt = MESH_Num_Attribs(mesh);

    list_attr_num = (int *) malloc(sizeof(int));
    list_attr_info = (int *) malloc(natt*sizeof(int));
    list_attr_names = (char *) malloc(natt*256*sizeof(char));
    
    int nattsend = 0;
    for(i = 0; i < natt; i++) {
      attrib = MESH_Attrib(mesh,i);
      MAttrib_Get_Name(attrib,attname);
      att_type = MAttrib_Get_Type(attrib);
      if (att_type == POINTER) continue;
      
      ncomp = MAttrib_Get_NumComps(attrib);
      mtype = MAttrib_Get_EntDim(attrib);
      list_attr_info[nattsend] = (ncomp << 6) | (mtype << 3) | (att_type);
      strcpy(&list_attr_names[nattsend*256],attname);
      nattsend++;
    }

    list_attr_num[0] = nattsend;
    MPI_Isend(list_attr_num,1,MPI_INT,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;

    (*ptrs2free)[(*numptrs2free)++] = list_attr_num;

    if (!nattsend) {
      free(list_attr_info);
      free(list_attr_names);
      return 1;
    }
    
    MPI_Isend(list_attr_info,nattsend,MPI_INT,torank,torank,comm,&mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;
    
    MPI_Isend(list_attr_names,nattsend*256,MPI_CHAR,torank,torank,comm,
              &mpirequest);
    (*requests)[*numreq] = mpirequest;
    (*numreq)++;


    (*ptrs2free)[(*numptrs2free)++] = list_attr_info;
    (*ptrs2free)[(*numptrs2free)++] = list_attr_names;

    return 1;
  }


#ifdef __cplusplus
}
#endif

