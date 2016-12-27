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
     Partition a mesh into as many submeshes as requested and distribute them

     Author(s): Duo Wang, Rao Garimella
  */


  /* Partition a given mesh into 'num' submeshes, adding a 'ring'
     layers of ghost elements around each partition. If 'with_attr' is
     1, attributes from the mesh are copied onto the submeshes. This
     routine does not send the meshes to other partitions */

  int MESH_Partition_and_Send(Mesh_ptr parentmesh,  int num, int *part, 
                              int *toranks, int ring, int with_attr, 
                              int del_inmesh, MSTK_Comm comm,
                              Mesh_ptr *mysubmesh) {
    int i, j, a, m, idx, ival, rank, numprocs, atttype;
    double rval;
    MAttrib_ptr attrib, g2latt, l2gatt;
    MSet_ptr mset;
    MVertex_ptr mv;
    MEdge_ptr me;
    MFace_ptr mf;
    MRegion_ptr mr;
    List_ptr g2llist;
    int numreq=0, maxreq=25, numptrs2free=0, maxptrs2free=25;
    MPI_Request *requests=NULL;
    void **ptrs2free = NULL;
    int maxpendreq = 200;
    int p, n;
    int torank;
    int maxparts = 2;
    int beginparts, endparts;

    int to_alloc = maxparts;
    if (num < maxparts) to_alloc = num;
    Mesh_ptr *submeshes = (Mesh_ptr *) malloc(to_alloc*sizeof(Mesh_ptr));
    
    char funcname[256] = "MESH_Partition_And_Send";


    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&numprocs);

    /* Check will go away once we let multiple processors do the distribution */
    if (numprocs != num)
      MSTK_Report(funcname,
                  "Number of partitions not equal to number of processors",
                  MSTK_FATAL);

    /* Create an attribute to keep track of the connections from
       entities of the global mesh to entities of submeshes */
    /* This attribute value will be populated with lists in
       MESH_Partition along with another attribute called Local2Global
       and used in MESH_BuildPBoundary and MESH_AddGhost. They are not
       needed subsequently */
    /* NOTE: PROBABLY SHOULD CHANGE IT TO PARENT2CHILD AND CHILD2PARENT */

    g2latt = MAttrib_New(parentmesh,"Global2Local",POINTER,MALLTYPE);

    
    beginparts = 0;
    while (beginparts < num) {
      endparts = beginparts + maxparts;
      if (endparts > num) endparts = num;

      /* Split the mesh into 'maxparts' submeshes */
      MESH_Partition_Some(parentmesh, num, part, beginparts, endparts, submeshes);

      for (i = beginparts; i < endparts; i++) {
        /* Tag entities as being in the partition interior or on the
           partition boundary */
        MESH_BuildPBoundary(parentmesh,submeshes[i-beginparts]);

        /* Add ghost layers */
        MESH_AddGhost(parentmesh,submeshes[i-beginparts],i,ring);
      }

      /* Send/receive mesh */
      requests = (MPI_Request *) malloc(maxreq*sizeof(MPI_Request));
      ptrs2free = (void **) malloc(maxptrs2free*sizeof(void *));

      /* Send Mesh Meta Data */
      for (n = beginparts; n < endparts; n++) {
        torank = toranks[n];
        if (torank == rank) continue;
        
        MESH_Send_MetaData(submeshes[torank-beginparts], torank, comm,
                           &numreq, &maxreq, &requests,
                           &numptrs2free, &maxptrs2free, &ptrs2free);
        
        /* check if we buffered too many requests - if so, we wait
           until all the data is sent out; if not, we continue. One
           can control how frequently we do a blocking wait for the
           send requests by adjusting maxpendreq */
        if (numreq > maxpendreq) {
          if (MPI_Waitall(numreq,requests,MPI_STATUSES_IGNORE) != MPI_SUCCESS)
            MSTK_Report("MSTK_Mesh_Distribute","Could not send mesh",MSTK_FATAL);
          else {
            numreq = 0;
            for (p = 0; p < numptrs2free; ++p) free(ptrs2free[p]);
            numptrs2free = 0;
          }
        }
      }

      /* Send Mesh Vertices */
      for (n = beginparts; n < endparts; n++) {
        torank = toranks[n];
        if (torank == rank) continue;

        MESH_Send_Vertices(submeshes[torank-beginparts], torank, comm,
                           &numreq, &maxreq, &requests,
                           &numptrs2free, &maxptrs2free, &ptrs2free);
        
        if (numreq > maxpendreq) {
          if (MPI_Waitall(numreq,requests,MPI_STATUSES_IGNORE) != MPI_SUCCESS)
            MSTK_Report("MSTK_Mesh_Distribute","Could not send mesh",MSTK_FATAL);
          else {
            numreq = 0;
            for (p = 0; p < numptrs2free; ++p) free(ptrs2free[p]);
            numptrs2free = 0;
          }
        }	
      }

      /* Send Mesh Vertices */
      for (n = beginparts; n < endparts; n++) {
        torank = toranks[n];
        if (torank == rank) continue;

        MESH_Send_VertexCoords(submeshes[torank-beginparts], torank, comm,
                &numreq, &maxreq, &requests,
                &numptrs2free, &maxptrs2free, &ptrs2free);
        
        if (numreq > maxpendreq) {
          if (MPI_Waitall(numreq,requests,MPI_STATUSES_IGNORE) != MPI_SUCCESS)
            MSTK_Report("MSTK_Mesh_Distribute","Could not send mesh",MSTK_FATAL);
          else {
            numreq = 0;
            for (p = 0; p < numptrs2free; ++p) free(ptrs2free[p]);
            numptrs2free = 0;
          }
        }	
      }

      /* Send higher dimensional mesh entities  */
      for (n = beginparts; n < endparts; n++) {
        torank = toranks[n];
        if (torank == rank) continue;

        MESH_Send_NonVertexEntities(submeshes[torank-beginparts], torank, comm,
                &numreq, &maxreq, &requests,
                &numptrs2free, &maxptrs2free, &ptrs2free);
        
        if (numreq > maxpendreq) {
          if (MPI_Waitall(numreq,requests,MPI_STATUSES_IGNORE) != MPI_SUCCESS)
            MSTK_Report("MSTK_Mesh_Distribute","Could not send mesh",MSTK_FATAL);
          else {
            numreq = 0;
            for (p = 0; p < numptrs2free; ++p) free(ptrs2free[p]);
            numptrs2free = 0;
          }
        }	
      }


      /* If requested, send attributes and mesh sets to partitions */
      if (with_attr) {

        /* First collect attribute information and copy to submeshes */
        int natt_global = MESH_Num_Attribs(parentmesh);
        char (*attnames)[256] = 
            (char (*)[256]) malloc(natt_global*sizeof(char [256]));

        for (a = 0; a < natt_global; a++) {
          attrib = MESH_Attrib(parentmesh,a);          

          MAttrib_Get_Name(attrib,attnames[a]);
          if (attrib == g2latt) continue;

          atttype = MAttrib_Get_Type(attrib);
          if (atttype == POINTER) continue;
          
          MESH_CopyAttr(parentmesh,beginparts-endparts,submeshes,attnames[a]);
        }        

        /* Send Attribute meta data */

        for (n = beginparts; n < endparts; n++) {
          torank = toranks[n];
          if (torank == rank) continue;
          
          MESH_Send_AttributeMetaData(submeshes[torank-beginparts], torank, comm,
                  &numreq, &maxreq, &requests,
                  &numptrs2free, &maxptrs2free, &ptrs2free);
          
          if (numreq > maxpendreq) {
            if (MPI_Waitall(numreq,requests,MPI_STATUSES_IGNORE) != MPI_SUCCESS)
              MSTK_Report("MSTK_Mesh_Distribute","Could not send mesh",MSTK_FATAL);
            else {
              numreq = 0;
              for (p = 0; p < numptrs2free; ++p) free(ptrs2free[p]);
              numptrs2free = 0;
            }
          }	
        }


        /* Send each attribute to the various processors */

        for (a = 0; a < natt_global; a++) {
          
          for (n = beginparts; n < endparts; n++) {
            torank = toranks[n];
            if (torank == rank) continue;
            
            attrib = MESH_AttribByName(submeshes[torank-beginparts],attnames[a]);
            if (!attrib) continue; /* this attribute does not exist on this  */
            /* processor - right now its not possible */
            /* but it might be in the future          */
            
            if (MAttrib_Get_Type(attrib) == POINTER) continue;
            
            MESH_Send_Attribute(submeshes[torank-beginparts], attrib, torank, comm,
                    &numreq, &maxreq, &requests,
                    &numptrs2free, &maxptrs2free, &ptrs2free);
            
            if (numreq > maxpendreq) {
              if (MPI_Waitall(numreq,requests,MPI_STATUSES_IGNORE) != MPI_SUCCESS)
                MSTK_Report("MSTK_Mesh_Distribute","Could not send mesh",MSTK_FATAL);
              else {
                numreq = 0;
                for (p = 0; p < numptrs2free; ++p) free(ptrs2free[p]);
                numptrs2free = 0;
              }
            }	
          }
        }
          
     
        
        /* First collect the mesh set information and copy into submeshes */

        int nset_global = MESH_Num_MSets(parentmesh);
        char (*msetnames)[256] = 
            (char (*)[256]) malloc(nset_global*sizeof(char [256]));

        for (m = 0; m < nset_global; m++) {
          mset = MESH_MSet(parentmesh,m);
          MSet_Name(mset,msetnames[m]);
          MESH_CopySet(parentmesh,endparts-beginparts,submeshes,mset);
        }
        
        /* Send Mesh Set Meta Data */
        for (n = beginparts; n < endparts; n++) {
          torank = toranks[n];
          if (torank == rank) continue;
          
          MESH_Send_MSetMetaData(submeshes[torank-beginparts], torank, comm,
                  &numreq, &maxreq, &requests,
                  &numptrs2free, &maxptrs2free, &ptrs2free);
          
          /* check if we buffered too many requests - if so, we wait
             until all the data is sent out; if not, we continue. One
             can control how frequently we do a blocking wait for the
             send requests by adjusting maxpendreq */
          
          if (numreq > maxpendreq) {
            if (MPI_Waitall(numreq,requests,MPI_STATUSES_IGNORE) != MPI_SUCCESS)
              MSTK_Report("MSTK_Mesh_Distribute","Could not send mesh",MSTK_FATAL);
            else {
              numreq = 0;
              for (p = 0; p < numptrs2free; ++p) free(ptrs2free[p]);
              numptrs2free = 0;
            }
          }	
        }
        
        
        /* Send Mesh Sets */
        
        for (m = 0; m < nset_global; m++) {
          
          for (n = beginparts; n < endparts; n++) {
            torank = toranks[n];
            if (torank == rank) continue;
            
            mset = MESH_MSetByName(submeshes[torank-beginparts],msetnames[m]);
            if (!mset) continue; /* this mset does not exist on this processor */
            
            MESH_Send_MSet(submeshes[torank-beginparts], mset, torank, comm,
                           &numreq, &maxreq, &requests,
                           &numptrs2free, &maxptrs2free, &ptrs2free);
            
            if (numreq > maxpendreq) {
              if (MPI_Waitall(numreq,requests,MPI_STATUSES_IGNORE) != MPI_SUCCESS)
                MSTK_Report("MSTK_Mesh_Distribute","Could not send mesh",MSTK_FATAL);
              else {
                numreq = 0;
                for (p = 0; p < numptrs2free; ++p) free(ptrs2free[p]);
                numptrs2free = 0;
              }
            }	
          }
          
        }
      }


      if (*mysubmesh == NULL)
        *mysubmesh = submeshes[rank];
      else if ((beginparts <= rank) && (rank < endparts))
        MESH_Copy(submeshes[rank],*mysubmesh,1,1);


      /* Final flush of all requests */

      if (numreq) {
        if (MPI_Waitall(numreq,requests,MPI_STATUSES_IGNORE) != MPI_SUCCESS)
          MSTK_Report("MSTK_Mesh_Distribute","Could not send mesh",MSTK_FATAL);
        else {
          numreq = 0;
          for (p = 0; p < numptrs2free; ++p) free(ptrs2free[p]);
          numptrs2free = 0;
        }
      }
      
      if (maxptrs2free) free(ptrs2free);
      if (maxreq) free(requests);


      /* Cleanup */
      /* Delete temporary submeshes */
      if (*mysubmesh != submeshes[rank])
        MESH_Delete(submeshes[rank]);
      else
        submeshes[rank] = NULL;

      for (n = 1; n < beginparts - endparts; ++n) {
        MESH_Delete(submeshes[n]);
        submeshes[n] = NULL;
      }

      /* increment and loop to cover all parts */
      beginparts = endparts;
    }

    /* cleanup */
    free(submeshes);


    /* Delete the lists associated with g2latt attribute but don't
       remove attribute itself from each of these entities - it will
       get deleted when the mesh gets deleted */

    idx = 0;
    while ((mv = MESH_Next_Vertex(parentmesh,&idx))) {
      MEnt_Get_AttVal(mv,g2latt,&ival,&rval,&g2llist);
      if (g2llist) List_Delete(g2llist);
    }
	
    idx = 0;
    while ((me = MESH_Next_Edge(parentmesh,&idx))) {
      MEnt_Get_AttVal(me,g2latt,&ival,&rval,&g2llist);
      if (g2llist) List_Delete(g2llist);
    }
	
    idx = 0;
    while ((mf = MESH_Next_Face(parentmesh,&idx))) {
      MEnt_Get_AttVal(mf,g2latt,&ival,&rval,&g2llist);
      if (g2llist) List_Delete(g2llist);
    }
	
    idx = 0;
    while ((mr = MESH_Next_Region(parentmesh,&idx))) {
      MEnt_Get_AttVal(mr,g2latt,&ival,&rval,&g2llist);
      if (g2llist) List_Delete(g2llist);
    }	
	
    if (del_inmesh) MESH_Delete(parentmesh);
    return 1;
  }

#ifdef __cplusplus
}
#endif

