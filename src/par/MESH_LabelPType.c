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

  int MESH_LabelPType_Face(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);
  int MESH_LabelPType_Region(Mesh_ptr submesh, int rank, int num, MPI_Comm comm);

  int MESH_LabelPType(Mesh_ptr submesh, int rank, int num,  MPI_Comm comm) {
  int nf, nr;
  nf = MESH_Num_Faces(submesh);
  nr = MESH_Num_Regions(submesh);
  if (nr)
    MESH_LabelPType_Region(submesh, rank, num, comm);
  else if(nf) 
    MESH_LabelPType_Face(submesh, rank, num, comm);
  else {
    MSTK_Report("MESH_LabelPType()","only send volume or surface mesh",MSTK_ERROR);
    exit(-1);
  }
  return 1;
}

int MESH_LabelPType_Face(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  int i, idx;
  MVertex_ptr mv;
  MEdge_ptr me;
  MFace_ptr mf;
  List_ptr fverts, fedges;

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

int MESH_LabelPType_Region(Mesh_ptr submesh, int rank, int num, MPI_Comm comm) {
  int i, nfv, idx;
  MRegion_ptr mr;
  MFace_ptr mf;
  MVertex_ptr me;
  MVertex_ptr mv;
  List_ptr mfverts, rverts, redges, rfaces;
  int is_ghost, is_overlap;
  idx = 0;

  /* first label faces */
  while( (mf = MESH_Next_Face(submesh,&idx)) ) {
    if(MF_PType(mf) != PGHOST)
      MF_Set_MasterParID(mf, rank);
    is_ghost = 1; is_overlap = 1;
    mfverts = MF_Vertices(mf,1,0);
    nfv = List_Num_Entries(mfverts);
    for(i = 0; i < nfv; i++) {
      mv = List_Entry(mfverts,i);
      if(MV_PType(mv) != PGHOST) {  /* if all vertices are ghost, then the face is a ghost*/ 
	is_ghost = 0;
	if(MV_PType(mv) != POVERLAP) /*  if all vertices are ghost but at least one of them is overlap */
	  is_overlap = 0;
      }
    }
    if(is_ghost) 
      MF_Set_PType(mf,PGHOST);
    else if(is_overlap) 
      MF_Set_PType(mf,POVERLAP);
    List_Delete(mfverts);
  }
  
  /* label overlap region */
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

