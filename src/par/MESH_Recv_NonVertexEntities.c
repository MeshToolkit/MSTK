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
     This function receives mesh entities other than vertices 
     from processor rank in communicator comm

     Author(s): Rao Garimella
  */

  int MESH_Recv_NonVertexEntities_FN(Mesh_ptr mesh, int fromrank, int ne,
                                     int nf, int nr, MSTK_Comm comm);
  int MESH_Recv_NonVertexEntities_R1R2(Mesh_ptr mesh, int fromrank, int ne,
                                       int nf, int nr, MSTK_Comm comm);
  int MESH_Recv_NonVertexEntities_R4(Mesh_ptr mesh, int fromrank, int ne,
                                     int nf, int nr, MSTK_Comm comm);

  static int (*MESH_Recv_NonVertexEntities_jmp[MSTK_MAXREP])(Mesh_ptr mesh, 
                                                             int fromrank, 
                                                             int ne, int nf,
                                                             int nr,
                                                             MSTK_Comm comm) =
  {MESH_Recv_NonVertexEntities_FN, MESH_Recv_NonVertexEntities_FN, 
   MESH_Recv_NonVertexEntities_R1R2, MESH_Recv_NonVertexEntities_R1R2, 
   MESH_Recv_NonVertexEntities_R4};

  int MESH_Recv_NonVertexEntities(Mesh_ptr mesh, int fromrank, int ne, int nf, 
                                  int nr, MSTK_Comm comm) {
    RepType rtype = MESH_RepType(mesh);
    return MESH_Recv_NonVertexEntities_jmp[rtype](mesh,fromrank,ne,nf,nr,comm);
  }


  int MESH_Recv_NonVertexEntities_FN(Mesh_ptr mesh, 
                                     int fromrank, 
                                     int ne, int nf, int nr, 
                                     MSTK_Comm comm) {
    int i, j, nevs, nfes, nrv, nrfs, maxnfe, maxnrf;
    int nfe, nrf, *fedirs=NULL, *rfdirs=NULL;
    int *list_edge=NULL, *list_face=NULL, *list_region=NULL;
    MVertex_ptr *verts, *rverts;
    MEdge_ptr me, *edges=NULL, *fedges=NULL;
    MFace_ptr mf, *faces=NULL, *rfaces=NULL;
    MRegion_ptr mr;
    RepType rtype;
    double coor[3];
    MPI_Status status;
    MPI_Request vrequest[2], erequest, frequest, rrequest;
    char mesg[256], errorstr[256];
    char funcname[256]="MESH_Recv_NonVertexEntities_FN";
    int errcode, len, nreq=0;
   
    int rank;
    MPI_Comm_rank(comm,&rank);

    /* Sanity check */
    
    if (!MESH_Num_Vertices(mesh))
      MSTK_Report(funcname,"Must receive mesh vertices first",MSTK_FATAL);

    if (!ne & !nf & !nr)
      MSTK_Report(funcname,"No entities to receive",MSTK_FATAL);

    
    maxnfe = MAXPV2;
    maxnrf = MAXPF3;

    /* Copied from MSTK_Send_MeshMetaData */

    nevs = (2+3)*ne;
    nfes = (1 + maxnfe + 3)*nf; 
    nrfs = (1 + maxnrf + 3)*nr;

    if (ne) {
      /* receive edge-vertex info */
      
      list_edge = (int *) malloc(nevs*sizeof(int));
      
      errcode = MPI_Irecv(list_edge,nevs,MPI_INT,fromrank,rank,comm,&erequest);
      if (errcode != MPI_SUCCESS)
        MSTK_Report(funcname,"Trouble receiving mesh edges",MSTK_FATAL);
    }
    

    if (nf) {
      /* receive face-edge info */
      
      list_face = (int *) malloc(nfes*sizeof(int));
      
      errcode = MPI_Irecv(list_face,nfes,MPI_INT,fromrank,rank,comm,&frequest);
      if (errcode != MPI_SUCCESS)
        MSTK_Report(funcname,"Trouble receiving mesh faces",MSTK_FATAL);
    }


    if (nr) {
      /* receive region-face info */
      
      list_region = (int *) malloc(nrfs*sizeof(int));
      
      errcode = MPI_Irecv(list_region,nrfs,MPI_INT,fromrank,rank,comm,&rrequest);
      if (errcode != MPI_SUCCESS)
        MSTK_Report(funcname,"Trouble receiving mesh regions",MSTK_FATAL);
    }


    /* while waiting for messages to complete, create the entitites */

    if (ne) {
      edges = (MEdge_ptr *) malloc(ne*sizeof(MEdge_ptr));
      for (i = 0; i < ne; i++) 
        edges[i] = ME_New(mesh);
    }
    if (nf) {
      faces = (MFace_ptr *) malloc(nf*sizeof(MFace_ptr));
      for (i = 0; i < nf; i++)
        faces[i] = MF_New(mesh);
    }
    if (nr) {
      for (i = 0; i < nr; i++)
        mr = MR_New(mesh);
    }


    /* Now fill in the data for the entities */

    if (ne) {
      /* if edge info is received, build edges */

      errcode = MPI_Wait(&erequest,MPI_STATUS_IGNORE);
      if (errcode != MPI_SUCCESS)
        MSTK_Report(funcname,"Trouble receiving edge info",MSTK_FATAL);
      
      nevs = 0;
      for (i = 0; i < ne; i++) {
        me = edges[i];
        int gentdim = list_edge[nevs+2] & 7; /* first 3 bits; 7 is 0...00111 */
        int gentid = list_edge[nevs+2] >> 3; /* All but the first 3 bits */
        ME_Set_GEntDim(me,gentdim);
        ME_Set_GEntID(me,gentid);

        int ptype = list_edge[nevs+3] & 3; /* first 2 bits; 3 is 0...00011 */
        int on_par_bdry = list_edge[nevs+3] & 4; /* 3rd bit; 4 is 0...00100 */
        int masterparid = list_edge[nevs+3] >> 3; /* All but the first 3 bits */
        ME_Set_PType(me,ptype);
        if (on_par_bdry)
          ME_Flag_OnParBoundary(me);
        ME_Set_MasterParID(me,masterparid);

        ME_Set_GlobalID(me,list_edge[nevs+4]);
        
        int vid0 = list_edge[nevs];
        int vid1 = list_edge[nevs+1];
        ME_Set_Vertex(me,0,MESH_VertexFromID(mesh,vid0));
        ME_Set_Vertex(me,1,MESH_VertexFromID(mesh,vid1));
        nevs += 5;
      }
    }



    if (nf) {
      /* if face info is received, build faces */
      
      errcode = MPI_Wait(&frequest,MPI_STATUS_IGNORE);
      if (errcode != MPI_SUCCESS)
        MSTK_Report(funcname,"Trouble receiving face info",MSTK_FATAL);
      
      fedges = (MEdge_ptr *) malloc(maxnfe*sizeof(MEdge_ptr));
      fedirs = (int *) malloc(maxnfe*sizeof(int));
      
      nfes = 0;
      for (i = 0; i < nf; i++) {
        mf = faces[i];
        nfe = list_face[nfes];
        for (j = 0; j < nfe; j++) {
          fedges[j] = edges[abs(list_face[nfes+j+1])-1];
          fedirs[j] = list_face[nfes+j+1] > 0 ? 1 : 0;
        }

        int gentdim = list_face[nfes+nfe+1] & 7; /* first 3 bits; 7 is 0...00111 */
        int gentid = list_face[nfes+nfe+1] >> 3; /* All but the first 3 bits */        
        MF_Set_GEntDim(mf,gentdim);
        MF_Set_GEntID(mf,gentid);

        int ptype = list_face[nfes+nfe+2] & 3; /* first 2 bits; 3 is 0...00011 */
        int on_par_bdry = list_face[nfes+nfe+2] & 4; /* 3rd bit; 4 is 0...00100 */
        int masterparid = list_face[nfes+nfe+2] >> 3; /* All but the first 3 bits */
        MF_Set_PType(mf,ptype);
        if (on_par_bdry)
          MF_Flag_OnParBoundary(mf);
        MF_Set_MasterParID(mf,masterparid);

        MF_Set_GlobalID(mf,list_face[nfes+nfe+3]);
        
        MF_Set_Edges(mf,nfe,fedges,fedirs);
        nfes += (nfe + 4);
      }
    }


    if (nr) {
      /* if region info is received, build regions */
      
      errcode = MPI_Wait(&rrequest,MPI_STATUS_IGNORE);
      if (errcode != MPI_SUCCESS)
        MSTK_Report(funcname,"Trouble receiving region info",MSTK_FATAL);
      
      
      rfaces = (MFace_ptr *) malloc(maxnrf*sizeof(MFace_ptr));
      rfdirs = (int *) malloc(maxnrf*sizeof(int));
      nrfs = 0;
      for (i = 0; i < nr; i++) {
        mr = MESH_Region(mesh,i);
        nrf = list_region[nrfs];
        for (j = 0; j < nrf; j++) {
          rfaces[j] = faces[abs(list_region[nrfs+j+1])-1];
          rfdirs[j] = list_region[nrfs+j+1] > 0 ? 1 : 0;
        }


        int gentdim = list_region[nrfs+nrf+1] & 7; /* first 3 bits; 7 is 0...00111 */
        int gentid = list_region[nrfs+nrf+1] >> 3; /* All but the first 3 bits */        
        MR_Set_GEntDim(mr,gentdim);
        MR_Set_GEntID(mr,gentid);

        int ptype = list_region[nrfs+nrf+2] & 3; /* first 2 bits; 3 is 0...00011 */
        /* int on_par_bdry = list_region[nrfs+nrf+2] & 4; */ /* 3rd bit should always be 0 for regions */
        int masterparid = list_region[nrfs+nrf+2] >> 3; /* All but the first 3 bits */
        MR_Set_PType(mr,ptype);
        MR_Set_MasterParID(mr,masterparid);

        MR_Set_GlobalID(mr,list_region[nrfs+nrf+3]);
        
        MR_Set_Faces(mr,nrf,rfaces,rfdirs); 
        nrfs += (nrf + 4);
      }
    }


    /* clean up */

    if (edges) free(edges);
    if (fedges) free(fedges);
    if (fedirs) free(fedirs);
    if (faces) free(faces);
    if (rfaces) free(rfaces);
    if (rfdirs) free(rfdirs);


    if (list_edge) free(list_edge);
    if (list_face) free(list_face);
    if (list_region) free(list_region);    

    return 1;
  }





  int MESH_Recv_NonVertexEntities_R1R2(Mesh_ptr mesh, int fromrank, int ne,
                                       int nf, int nr, MSTK_Comm comm) {
    MSTK_Report("MESH_Surf_RecvMesh_R1R2","Not implemented",MSTK_FATAL);
    return 0;
  }



  int MESH_Recv_NonVertexEntities_R4(Mesh_ptr mesh, int fromrank, int ne,
                                       int nf, int nr, MSTK_Comm comm) {
    MSTK_Report("MESH_Surf_RecvMesh_R4","Not implemented",MSTK_FATAL);
    return 0;
  }




#ifdef __cplusplus
}
#endif

