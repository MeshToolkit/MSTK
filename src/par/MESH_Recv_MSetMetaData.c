#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif


  /* Receive meta data about entity sets (number of sets, their names,
     types of entities) of a mesh from distributing processor

     Authors: Rao Garimella
  */

  int MESH_Recv_MSetMetaData(Mesh_ptr mesh, int fromrank, MSTK_Comm comm) {
    int i, nset, errcode, rank, len;
    int mtype;   
    char msetname[256], mesg[256], errorstr[256];
    MSet_ptr mset;
    char funcname[256] = "MESH_RecvMSets";
    int *list_mset_types;
    char *list_mset_names;
    MPI_Status status;
    
    MPI_Comm_rank(comm,&rank);

    /* receive number of mesh sets and meta data about them */

    errcode = MPI_Recv(&nset,1,MPI_INT,fromrank,rank,comm,&status);
    if (errcode != MPI_SUCCESS) {
      MPI_Error_string(errcode,errorstr,&len);
      sprintf(mesg,"MPI_Recv failed with error %s",errorstr);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }

    if (!nset) return 1;

    list_mset_types = (int *) malloc(nset*sizeof(int));
    
    errcode = MPI_Recv(list_mset_types,nset,MPI_INT,fromrank,rank,comm,&status);
    if (errcode != MPI_SUCCESS) {
      MPI_Error_string(errcode,errorstr,&len);
      sprintf(mesg,"MPI_Recv failed with error %s",errorstr);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }


    list_mset_names = (char *) malloc(nset*256*sizeof(char));
    errcode = MPI_Recv(list_mset_names,nset*256,MPI_CHAR,fromrank,rank,comm,
                        &status);
    if (errcode != MPI_SUCCESS) {
      MPI_Error_string(errcode,errorstr,&len);
      sprintf(mesg,"MPI_Recv failed with error %s",errorstr);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }
    

    for(i = 0; i < nset; i++) {
      strcpy(msetname,&list_mset_names[i*256]);      
      if (MESH_MSetByName(mesh,msetname)) {
        sprintf(mesg,"Mesh set by the name %s already exists on partition %d",
                msetname,rank);
        MSTK_Report(funcname,mesg,MSTK_FATAL);
      }
      
      mtype = list_mset_types[i];
      mset =  MSet_New(mesh, msetname, mtype);
    }
    
    if (nset) {
      free(list_mset_types);    
      free(list_mset_names);
    }
    
    return 1;
  }
   

#ifdef __cplusplus
}
#endif

