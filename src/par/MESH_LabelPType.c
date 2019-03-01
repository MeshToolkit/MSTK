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



  /* 
     This function label 1-ring boundary layer
     
     It assigns all the elements with a POVERLAP or PGHOST vertex as POVERLAP

     Author(s): Duo Wang, Rao Garimella
  */

  int MESH_LabelPType_Face(Mesh_ptr submesh, MSTK_Comm comm);
  int MESH_LabelPType_Region(Mesh_ptr submesh, MSTK_Comm comm);

  int MESH_LabelPType(Mesh_ptr submesh, int topodim, MSTK_Comm comm) {
  if (topodim == 3)
    MESH_LabelPType_Region(submesh, comm);
  else if (topodim == 2) 
    MESH_LabelPType_Face(submesh, comm);
  else {
    MSTK_Report("MESH_LabelPType()","only send volume or surface mesh",MSTK_ERROR);
    exit(-1);
  }
  return 1;
}

  int MESH_LabelPType_Face(Mesh_ptr submesh, MSTK_Comm comm) {
  int i, idx;
  MVertex_ptr mv;
  MEdge_ptr me;
  MFace_ptr mf;
  List_ptr fverts, fedges;

  int rank, size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  /* label 1-ring boundary face as POVERLAP */
  idx = 0;
  while( (mf = MESH_Next_Face(submesh,&idx)) ) {
    fverts = MF_Vertices(mf,1,0);
    for(i = 0; i < List_Num_Entries(fverts); i++) {
      mv = List_Entry(fverts,i);
      if( MV_PType(mv) == PGHOST || MV_PType(mv) == POVERLAP ) {
	MF_Set_PType(mf,POVERLAP);
	break;
      }
    }
    List_Delete(fverts);
  }
  /* 
     if the face is marked, mark its edges and vertices
     ghost should be still ghost
  */
  idx = 0;
  while( (mf = MESH_Next_Face(submesh,&idx)) ) {
    if(MF_PType(mf) != POVERLAP) continue;
    fverts = MF_Vertices(mf,1,0);
    for(i = 0; i < List_Num_Entries(fverts); i++) {
      mv = List_Entry(fverts,i);
      if(MV_PType(mv) != PGHOST)
	MV_Set_PType(mv,POVERLAP);
    }
    List_Delete(fverts);
    
    fedges = MF_Edges(mf,1,0);
    for(i = 0; i < List_Num_Entries(fedges); i++) {
      me = List_Entry(fedges,i);
      if(ME_PType(me) != PGHOST)
	ME_Set_PType(me,POVERLAP);
    }
    List_Delete(fedges);
  }
  return 1;
}

  /* 
     For 3D mesh only, label 1-ring inner neighbor
     Assume regions are not overlapped across processors before this     
     Label the region that has OVERLAP vertex as OVERLAP
  */

  int MESH_LabelPType_Region(Mesh_ptr submesh, MSTK_Comm comm) {
  int i, idx;
  MRegion_ptr mr;
  MFace_ptr mf;
  MVertex_ptr me;
  MVertex_ptr mv;
  List_ptr rverts, redges, rfaces;
  
  int rank, size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  /* label 1-ring boundary region as POVERLAP */
  idx = 0;
  while( (mr = MESH_Next_Region(submesh,&idx)) ) {
    rverts = MR_Vertices(mr);
    for(i = 0; i < List_Num_Entries(rverts); i++) {
      mv = List_Entry(rverts,i);
      if(MV_PType(mv) == PGHOST || MV_PType(mv) == POVERLAP) {
	MR_Set_PType(mr,POVERLAP);
	break;
      }
    }
    List_Delete(rverts);
  }

  /* if the region is marked, mark its vertices, edges and faces */
  idx = 0;
  while( (mr = MESH_Next_Region(submesh,&idx)) ) {
    if(MV_PType(mr) != POVERLAP) continue;
    rverts = MR_Vertices(mr);
    for(i = 0; i < List_Num_Entries(rverts); i++) {
      mv = List_Entry(rverts,i);
      if(MV_PType(mv) != PGHOST)
	MV_Set_PType(mv,POVERLAP);
    }
    List_Delete(rverts);

    redges = MR_Edges(mr);
    for(i = 0; i < List_Num_Entries(redges); i++) {
      me = List_Entry(redges,i);
      if(ME_PType(me) != PGHOST)
	ME_Set_PType(me,POVERLAP);
    }
    List_Delete(redges);

    rfaces = MR_Faces(mr);
    for(i = 0; i < List_Num_Entries(rfaces); i++) {
      mf = List_Entry(rfaces,i);
      if(MF_PType(mf) != PGHOST)
	  MF_Set_PType(mf,POVERLAP);
    }
    List_Delete(rfaces);
  }
  return 1;
}


  
#ifdef __cplusplus
}
#endif

