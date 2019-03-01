/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdlib.h>
#include <stdio.h>
#include "MSTK.h"


MFace_ptr MFs_Join(MFace_ptr f0, MFace_ptr f1, MEdge_ptr e) {
  int i, k0, k1, nfe0, nfe1, *fedir0, *fedir1, *fedir3, gdim, gid;
  MEdge_ptr *fe0, *fe1, *fe3;
  MFace_ptr  nuface;
  MRegion_ptr fr00=NULL, fr01=NULL, fr10=NULL, fr11=NULL;
  Mesh_ptr   mesh;
  List_ptr   fedges, fregs;

  mesh = MF_Mesh(f0);
  gid = MF_GEntID(f0);
  gdim = MF_GEntDim(f0);

  if (mesh != MF_Mesh(f1)) {
    MSTK_Report("MFs_Join","Faces not from same mesh",MSTK_ERROR);
    return 0;
  }
  else if (gid != MF_GEntID(f0) || gdim != MF_GEntDim(f0)) {
    MSTK_Report("MFs_Join","Faces not from same geometric entity",MSTK_ERROR);
    return 0;
  }

  fregs = MF_Regions(f0);
  if (fregs) {
    int nfr = List_Num_Entries(fregs);
    fr00 = (nfr > 0) ? List_Entry(fregs,0) : NULL;
    fr01 = (nfr > 1) ? List_Entry(fregs,1) : NULL;
    List_Delete(fregs);
  }

  fregs = MF_Regions(f1);
  if (fregs) {
    int nfr = List_Num_Entries(fregs);
    fr10 = (nfr > 0) ? List_Entry(fregs,0) : NULL;
    fr11 = (nfr > 1) ? List_Entry(fregs,1) : NULL;
    List_Delete(fregs);
  }

  if (fr00 != fr10 || fr01 != fr11) {
    MSTK_Report("MFs_Join","Faces have incompatible regions connected to them",
                MSTK_ERROR);
    return 0;
  }


  nfe0 = MF_Num_Edges(f0);
  nfe1 = MF_Num_Edges(f1);

  fe0 = (MEdge_ptr *) malloc(nfe0*sizeof(MEdge_ptr));
  fedir0 = (int *) malloc(nfe0*sizeof(int));
  fe1 = (MEdge_ptr *) malloc(nfe1*sizeof(MEdge_ptr));
  fedir1 = (int *) malloc(nfe1*sizeof(int));
  fe3 = (MEdge_ptr *) malloc((nfe0+nfe1-2)*sizeof(MEdge_ptr));
  fedir3 = (int *) malloc((nfe0+nfe1-2)*sizeof(int));

  fedges = MF_Edges(f0,1,0);
  for (i = 0, k0 = -1; i < nfe0; i++) {
    fe0[i] = List_Entry(fedges,i);
    fedir0[i] = MF_EdgeDir_i(f0,i);
    if (fe0[i] == e)
      k0 = i;
  }
  List_Delete(fedges);

  fedges = MF_Edges(f1,1,0);
  for (i = 0, k1 = -1; i < nfe1; i++) {
    fe1[i] = List_Entry(fedges,i);
    fedir1[i] = MF_EdgeDir_i(f1,i);
    if (fe1[i] == e)
      k1 = i;
  }
  List_Delete(fedges);

  if (k0 == -1 || k1 == -1) {
    MSTK_Report("MFs_Join","Cannot find edge in face",MSTK_ERROR);
    return 0;
  }

  for (i = 0; i < nfe0-1; i++) {
    fe3[i] = fe0[(k0+1+i)%nfe0];
    fedir3[i] = fedir0[(k0+1+i)%nfe0];
  }

  for (i = 0; i < nfe1-1; i++) {
    fe3[(nfe0-1)+i] = fe1[(k1+1+i)%nfe1];
    fedir3[(nfe0-1)+i] = fedir1[(k1+1+i)%nfe1];
  }

  
  nuface = MF_New(mesh);
  MF_Set_Edges(nuface,(nfe0+nfe1-2),fe3,fedir3);

  MF_Set_GEntDim(nuface,gdim);
  MF_Set_GEntID(nuface,gid);

  MFace_ptr oldf[2];
  oldf[0] = f0; oldf[1] = f1;
  if (fr00) {
    int rfdir = MR_FaceDir(fr00,f0);
    MR_Replace_Faces(fr00,2,oldf,1,&nuface,&rfdir);
  }
  if (fr01) {
    int rfdir = MR_FaceDir(fr01,f0);
    MR_Replace_Faces(fr01,2,oldf,1,&nuface,&rfdir);
  }

  MF_Delete(f0,0);
  MF_Delete(f1,0);
  ME_Delete(e,0);

  free(fe0); free(fedir0);
  free(fe1); free(fedir1);
  free(fe3); free(fedir3);

  return nuface;
}
