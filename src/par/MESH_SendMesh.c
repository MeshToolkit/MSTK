#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif



  /* 
     This function sends a mesh to processor torank in communicator comm

     Author(s): Rao Garimella
                Duo Wang
  */



  int MESH_SendMesh(Mesh_ptr mesh, int torank, int with_attr, MSTK_Comm comm, 
                    int *numreq, int *maxreq, MPI_Request **requests,
                    int *numptrs2free, int *maxptrs2free, void ***ptrs2free) {
    int i, j, a, m, p, rank;
    MAttrib_ptr attrib;
    MSet_ptr mset;
    RepType rtype;
    MPI_Request mpirequest;
    int maxpendreq = 0;

    if (requests == NULL)
      MSTK_Report("MESH_Surf_SendMesh","MPI requests array is NULL",MSTK_FATAL);

    if (*maxreq == 0) {
      *maxreq = 25;
      *requests = (MPI_Request *) malloc(*maxreq*sizeof(MPI_Request));
      *numreq = 0;
    }
    else if (*maxreq < (*numreq) + 13) {
      *maxreq = 2*(*maxreq) + 11;
      *requests = (MPI_Request *) realloc(*requests,*maxreq*sizeof(MPI_Request));
    }


    MESH_Send_MetaData(mesh, torank, comm,
                       numreq, maxreq, requests,
                       numptrs2free, maxptrs2free, ptrs2free);
        
    /* check if we buffered too many requests - if so, we wait
       until all the data is sent out; if not, we continue. One
       can control how frequently we do a blocking wait for the
       send requests by adjusting maxpendreq */
    
    if (*numreq > maxpendreq) {
      if (MPI_Waitall(*numreq,*requests,MPI_STATUSES_IGNORE) != MPI_SUCCESS)
        MSTK_Report("MSTK_Mesh_Distribute","Could not send mesh",MSTK_FATAL);
      else {
        *numreq = 0;
        for (p = 0; p < *numptrs2free; ++p) free((*ptrs2free)[p]);
        *numptrs2free = 0;
      }
    }	



    /* Send Mesh Vertices */
    
    MESH_Send_Vertices(mesh, torank, comm,
                       numreq, maxreq, requests,
                       numptrs2free, maxptrs2free, ptrs2free);
    
    
    if (*numreq > maxpendreq) {
      if (MPI_Waitall(*numreq,*requests,MPI_STATUSES_IGNORE) != MPI_SUCCESS)
        MSTK_Report("MSTK_Mesh_Distribute","Could not send mesh",MSTK_FATAL);
      else {
        *numreq = 0;
        for (p = 0; p < *numptrs2free; ++p) free((*ptrs2free)[p]);
        *numptrs2free = 0;
      }
    }	



    /* Send Mesh Vertex Coordinates */
    
    MESH_Send_VertexCoords(mesh, torank, comm,
                           numreq, maxreq, requests,
                           numptrs2free, maxptrs2free, ptrs2free);
    
    
    if (*numreq > maxpendreq) {
      if (MPI_Waitall(*numreq,*requests,MPI_STATUSES_IGNORE) != MPI_SUCCESS)
        MSTK_Report("MSTK_Mesh_Distribute","Could not send mesh",MSTK_FATAL);
      else {
        *numreq = 0;
        for (p = 0; p < *numptrs2free; ++p) free((*ptrs2free)[p]);
        *numptrs2free = 0;
      }
    }	


    /* Send higher dimensional mesh entities  */

    MESH_Send_NonVertexEntities(mesh, torank, comm,
                                numreq, maxreq, requests,
                                numptrs2free, maxptrs2free, ptrs2free);
  
    if (*numreq > maxpendreq) {
      if (MPI_Waitall(*numreq,*requests,MPI_STATUSES_IGNORE) != MPI_SUCCESS)
        MSTK_Report("MSTK_Mesh_Distribute","Could not send mesh",MSTK_FATAL);
      else {
        *numreq = 0;
        for (p = 0; p < *numptrs2free; ++p) free((*ptrs2free)[p]);
        *numptrs2free = 0;
      }
    }	


    if (with_attr) {

      /* Send Attribute meta data */
      
      MESH_Send_AttributeMetaData(mesh, torank, comm,
                                  numreq, maxreq, requests,
                                  numptrs2free, maxptrs2free, ptrs2free);
      
      if (*numreq > maxpendreq) {
        if (MPI_Waitall(*numreq,*requests,MPI_STATUSES_IGNORE) != MPI_SUCCESS)
          MSTK_Report("MSTK_Mesh_Distribute","Could not send mesh",MSTK_FATAL);
        else {
          *numreq = 0;
          for (p = 0; p < *numptrs2free; ++p) free((*ptrs2free)[p]);
          *numptrs2free = 0;
        }
      }	
      
      
      /* Send Attributes */
      int natt = MESH_Num_Attribs(mesh);
      for (a = 0; a < natt; a++) {
        
        attrib = MESH_Attrib(mesh,a);
        if (MAttrib_Get_Type(attrib) == POINTER) continue;
        
        MESH_Send_Attribute(mesh, attrib, torank, comm,
                            numreq, maxreq, requests,
                            numptrs2free, maxptrs2free, ptrs2free);
        
        
        if (*numreq > maxpendreq) {
          if (MPI_Waitall(*numreq,*requests,MPI_STATUSES_IGNORE) != MPI_SUCCESS)
            MSTK_Report("MSTK_Mesh_Distribute","Could not send mesh",MSTK_FATAL);
          else {
            *numreq = 0;
            for (p = 0; p < *numptrs2free; ++p) free((*ptrs2free)[p]);
            *numptrs2free = 0;
          }
        }	
      }
      
      
      /* Send Mesh Set Meta Data */
      
      MESH_Send_MSetMetaData(mesh, torank, comm,
                             numreq, maxreq, requests,
                             numptrs2free, maxptrs2free, ptrs2free);
      
      
      if (*numreq > maxpendreq) {
        if (MPI_Waitall(*numreq,*requests,MPI_STATUSES_IGNORE) != MPI_SUCCESS)
          MSTK_Report("MSTK_Mesh_Distribute","Could not send mesh",MSTK_FATAL);
        else {
          *numreq = 0;
          for (p = 0; p < *numptrs2free; ++p) free((*ptrs2free)[p]);
          *numptrs2free = 0;
        }
      }	
      
      
      /* Send Mesh Sets */
      
      int nset = MESH_Num_MSets(mesh);
      for (m = 0; m < nset; m++) {
        
        mset = MESH_MSet(mesh,m);
        
        MESH_Send_MSet(mesh, mset, torank, comm,
                       numreq, maxreq, requests,
                       numptrs2free, maxptrs2free, ptrs2free);
        
        if (*numreq > maxpendreq) {
          if (MPI_Waitall(*numreq,*requests,MPI_STATUSES_IGNORE) != MPI_SUCCESS)
            MSTK_Report("MSTK_Mesh_Distribute","Could not send mesh",MSTK_FATAL);
          else {
            *numreq = 0;
            for (p = 0; p < *numptrs2free; ++p) free((*ptrs2free)[p]);
            *numptrs2free = 0;
          }
        }	
      }
      
      if (*numreq) {
        if (MPI_Waitall(*numreq,*requests,MPI_STATUSES_IGNORE) != MPI_SUCCESS)
          MSTK_Report("MSTK_Mesh_Distribute","Could not send mesh",MSTK_FATAL);
        else {
          *numreq = 0;
          for (p = 0; p < *numptrs2free; ++p) free((*ptrs2free)[p]);
          *numptrs2free = 0;
        }
      }

    } /* if (with_attr) */
  

    return 1;
  }


#ifdef __cplusplus
}
#endif

