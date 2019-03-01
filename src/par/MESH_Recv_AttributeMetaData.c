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


  /* Receive a list of mesh attributes and their defining meta data
     (name, type, type of entity they live on, etc.) from distributing
     processor

     Authors: Rao Garimella
  */

  int MESH_Recv_AttributeMetaData(Mesh_ptr mesh, int fromrank, MSTK_Comm comm) {
    int i, natt, ncomp, errcode, rank, len;
    int mtype, att_type;
    char attname[256], mesg[256], errorstr[256];
    int *list_attr_info;
    char *list_attr_names;
    MAttrib_ptr attrib;
    char funcname[256] = "MESH_Recv_AttributeMetaData";
    MPI_Status status;

    MPI_Comm_rank(comm,&rank);

    /* receive number of attributes and meta data about them */

    errcode = MPI_Recv(&natt,1,MPI_INT,fromrank,rank,comm,&status);
    if (errcode != MPI_SUCCESS) {
      MPI_Error_string(errcode,errorstr,&len);
      sprintf(mesg,"MPI_Recv failed with error %s",errorstr);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }

    if (!natt) return 1;

    list_attr_info = (int *) malloc(natt*sizeof(int));
    
    errcode = MPI_Recv(list_attr_info,natt,MPI_INT,fromrank,rank,comm,&status);
    if (errcode != MPI_SUCCESS) {
      MPI_Error_string(errcode,errorstr,&len);
      sprintf(mesg,"MPI_Recv failed with error %s",errorstr);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }


    list_attr_names = (char *) malloc(natt*256*sizeof(char));
    errcode = MPI_Recv(list_attr_names,natt*256,MPI_CHAR,fromrank,rank,comm,
                        &status);
    if (errcode != MPI_SUCCESS) {
      MPI_Error_string(errcode,errorstr,&len);
      sprintf(mesg,"MPI_Recv failed with error %s",errorstr);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }
    
    for(i = 0; i < natt; i++) {
      strcpy(attname,&list_attr_names[i*256]);
      
      if (MESH_AttribByName(mesh,attname)) {
        sprintf(mesg,
                "Mesh attribute by the name %s already exists on partition %d",
                attname,rank);
        MSTK_Report(funcname,mesg,MSTK_FATAL);
      }
      
      att_type = list_attr_info[i] & 7;
      mtype = (list_attr_info[i] >> 3) & 7;
      ncomp = list_attr_info[i] >> 6;
      if(ncomp == 1)
        attrib =  MAttrib_New(mesh, attname, att_type, mtype);
      else
        attrib =  MAttrib_New(mesh, attname, att_type, mtype, ncomp);
    }
    
    if (natt) {
      free(list_attr_info);    
      free(list_attr_names);
    }
    
    return 1;
  }
    


#ifdef __cplusplus
}
#endif

