#include <stdio.h>
#include <stdlib.h>

#include "MSTK.h"

/* Routine to split a face along a line connecting two given vertices */
/* A new edge is created and two new faces are created in place of the
   old face. The two new faces are incorporated into connected regions */

#ifdef __cplusplus
extern "C" {
#endif

MEdge_ptr MF_Split(MFace_ptr fsplit, MVertex_ptr v0, MVertex_ptr v1) {
  Mesh_ptr mesh;
  MVertex_ptr vnext;
  MEdge_ptr enew, *fedges0, *fedges1, fe;
  MFace_ptr fnew[2];
  MRegion_ptr fr;
  int gid, gdim, i, idx, found, done;
  int nfe, n0, n1, *fedirs, *fedirs0, *fedirs1;;
  List_ptr fregions, felist;

  enew = MVs_CommonEdge(v0,v1);
  if (enew)
    return enew;

  gdim = MF_GEntDim(fsplit);
  gid = MF_GEntID(fsplit);
 
  /* Collect information */

  mesh = MF_Mesh(fsplit);

  
  felist = MF_Edges(fsplit,1,0);
  nfe = List_Num_Entries(felist);
  fedirs = (int *) malloc(nfe*sizeof(int));
  for (i = 0; i < nfe; i++)
    fedirs[i] = MF_EdgeDir_i(fsplit,i);

  
  /* Setup for making the left and right faces after splitting */

  fedges0 = (MEdge_ptr *) malloc(nfe*sizeof(MEdge_ptr));
  fedirs0 = (int *) malloc(nfe*sizeof(int));

  fedges1 = (MEdge_ptr *) malloc(nfe*sizeof(MEdge_ptr));
  fedirs1 = (int *) malloc(nfe*sizeof(int));

  /* Create the splitting edge */

  enew = ME_New(mesh);
  ME_Set_Vertex(enew,0,v0);
  ME_Set_Vertex(enew,1,v1);
  ME_Set_GEntDim(enew,gdim);
  ME_Set_GEntID(enew,gid);


  /* collect edges and dirs for right face */
  found = 0;
  for (i = 0; i < nfe; i++) {
    fe = List_Entry(felist,i);
    if (ME_Vertex(fe,!fedirs[i]) == v0) {
      found = 1;
      break;
    }
  }

  if (!found)
    MSTK_Report("MF_Split","Cannot find first input vertex in face to be split",
                MSTK_FATAL);


  n0 = 0;
  done = 0;
  while (!done) {
    fedges0[n0] = (MEdge_ptr) List_Entry(felist,i%nfe);
    fedirs0[n0] = fedirs[i%nfe];
    vnext = ME_Vertex(fedges0[n0],fedirs0[n0]);
    n0++;

    if (vnext == v1)
      done = 1;
    else if (vnext == v0) 
      MSTK_Report("MF_Split",
                  "Could not find second input vertex in face to be split",
                  MSTK_FATAL);
    else
      i++;
  }

  fedges0[n0] = enew;
  fedirs0[n0] = 0;
  n0++;

  fnew[0] = MF_New(mesh);
  MF_Set_Edges(fnew[0],n0,fedges0,fedirs0);
  MF_Set_GEntDim(fnew[0],gdim);
  MF_Set_GEntID(fnew[0],gid);


  /* Construct the left face */
  
  n1 = 0;
  i = i+1;
  done = 0;
  while (!done) {
    fedges1[n1] = (MEdge_ptr) List_Entry(felist,i%nfe);
    fedirs1[n1] = fedirs[i%nfe];
    vnext = ME_Vertex(fedges1[n1],fedirs1[n1]);
    n1++;

    if (vnext == v0)
      done = 1;
    else if (vnext == v1)
      MSTK_Report("MF_Split",
                  "Could not find first input vertex in face to be split",
                  MSTK_FATAL);
    else
      i++;
  }

  fedges1[n1] = enew;
  fedirs1[n1] = 1;
  n1++;

  fnew[1] = MF_New(mesh);
  MF_Set_Edges(fnew[1],n1,fedges1,fedirs1);
  MF_Set_GEntDim(fnew[1],gdim);
  MF_Set_GEntID(fnew[1],gid);


  free(fedges0);
  free(fedirs0);
  free(fedges1);
  free(fedirs1);

  /* Update face lists of of connected regions */

  fregions = MF_Regions(fsplit);

  idx = 0;
  while ((fr = List_Next_Entry(fregions,&idx))) {
    int fdir[2];
    fdir[0] = fdir[1] = MR_FaceDir(fr,fsplit);
    MR_Replace_Faces(fr,1,&fsplit,2,fnew,fdir);
  }

  List_Delete(fregions);

  /* Delete the split face */

  MF_Delete(fsplit,1);


  return enew;
}

#ifdef __cplusplus
}
#endif
