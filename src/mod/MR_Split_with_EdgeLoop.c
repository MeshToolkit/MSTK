/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>

#include "MSTK.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Routine to split a region along a loop of edges */
/* A new face is created and two new regions are created in place of the
   old region */

/* MISSING A UNIT TEST */

MFace_ptr MR_Split_with_EdgeLoop(MRegion_ptr rsplit, int nfe, MEdge_ptr *fedges) {
  Mesh_ptr mesh;
  int i, j, idx1, idx2, fedirs[MAXPV2], rfdir_adj, curdir, edir, edir_adj;;
  int gid, mkid, nrf1, nrf2, rfdirs1[MAXPF3], rfdirs2[MAXPF3];
  MEdge_ptr fe;
  MFace_ptr fnew, eface, rfarray1[MAXPF3], rfarray2[MAXPF3], curface;
  MRegion_ptr rnew[2];
  List_ptr felist, efaces;

#ifdef DEBUG
  List_ptr redges;
#endif

  gid = MR_GEntID(rsplit);

#ifdef DEBUG
  /* check to make sure we got meaningful input */
  redges = MR_Edges(rsplit);
  for (i = 0; i < nfe; i++)
    if (!List_Contains(redges,fedges[i]))
      MSTK_Report("MR_Split","Input edges are not part of the region boundary",
                  MSTK_FATAL);
  List_Delete(redges);
#endif

  mesh = MR_Mesh(rsplit);


  /* Fix a set of directions for the edges */

  fedirs[0] = 1;
  for (i = 1; i < nfe; i++) {
    MVertex_ptr vprev, v0, v1;

    vprev = ME_Vertex(fedges[i-1],fedirs[i-1]);
    v0 = ME_Vertex(fedges[i],0);
    v1 = ME_Vertex(fedges[i],1);
    if (vprev == v0)
      fedirs[i] = 1;
    else if (vprev == v1)
      fedirs[i] = 0;
    else
      MSTK_Report("MR_Split","Input edges do not form a loop as listed",
                  MSTK_FATAL);
  }

  /* Create the splitting face */

  fnew = MF_New(mesh);
  MF_Set_GEntDim(fnew,3);
  MF_Set_GEntID(fnew,gid);
  MF_Set_Edges(fnew,nfe,fedges,fedirs);


  /* Collect info for the first region */

  List_ptr processed_faces = List_New(0);

  rfarray1[0] = fnew;
  rfdirs1[0] = 1; 
  nrf1 = 1;
  List_Add(processed_faces,rfarray1[0]);

  i = 0;
  while (i < nrf1) {
    curface = rfarray1[i];
    curdir = rfdirs1[i];
    i++;

    /* Get adjacent faces in region of current face and if they are
       not already in the new region face list (not marked), then add
       them */
 
    felist = MF_Edges(curface,1,0);
    idx1 = 0; j = 0;
    while ((fe = List_Next_Entry(felist,&idx1))) {
      edir = MF_EdgeDir_i(curface,j);      
      j++;

      efaces = ME_Faces(fe);
      if (curface != fnew && List_Contains(efaces,fnew)) {
        /* we have come back to the starting or splitting face - don't
           go across this edge */
        List_Delete(efaces);
        continue;
      }

      /* Add an adjacent unprocessed face of the region to the list of 
       faces for the new region */
      idx2 = 0; 
      while ((eface = List_Next_Entry(efaces,&idx2))) { 
        if (eface == curface) continue;
        if (List_Contains(processed_faces,eface)) continue;
        if (!MR_UsesEntity(rsplit,eface,MFACE)) continue; /* does not belong to region */

        edir_adj = MF_EdgeDir(eface,fe);
        rfdir_adj = MR_FaceDir(rsplit,eface);

        /* add adjacent face based on the check that if two adjacent faces of
           region are used by the region in the same sense, then their common
           edge should be used by the two faces in opposite senses (or the opposite
           of both the conditions should be true) */

        if ((edir != edir_adj && curdir == rfdir_adj) ||
            (edir == edir_adj && curdir != rfdir_adj)) {
          rfarray1[nrf1] = eface;
          rfdirs1[nrf1] = rfdir_adj;
          List_Add(processed_faces,rfarray1[nrf1]);
          nrf1++;
          break;
        }
      }
      List_Delete(efaces);
    }
    List_Delete(felist);
  }


  /* collect info for the second region */

  rfarray2[0] = fnew;
  rfdirs2[0] = !rfdirs1[0];
  nrf2 = 1;
  List_Add(processed_faces,rfarray2[0]);
  
  i = 0;
  while (i < nrf2) {
    curface = rfarray2[i];
    curdir = rfdirs2[i];
    i++;

    /* Get adjacent faces in region of current face and if they are
       not already in the new region face list (not marked), then add
       them */
 
    felist = MF_Edges(curface,1,0);
    idx1 = 0; j = 0;
    while ((fe = List_Next_Entry(felist,&idx1))) {
      edir = MF_EdgeDir_i(curface,j);      
      j++;

      efaces = ME_Faces(fe);
      if (curface != fnew && List_Contains(efaces,fnew)) {
        /* we have come back to the starting or splitting face - don't
           go across this edge */
        List_Delete(efaces);
        continue;
      }

      
      /* Add an adjacent unprocessed face of the region to the list of 
       faces for the new region */
      idx2 = 0; 
      while ((eface = List_Next_Entry(efaces,&idx2))) { 
        if (eface == curface) continue;
        if (List_Contains(processed_faces,eface)) continue;
        if (!MR_UsesEntity(rsplit,eface,MFACE)) continue; /* does not belong to region */

        edir_adj = MF_EdgeDir(eface,fe);
        rfdir_adj = MR_FaceDir(rsplit,eface);

        /* add adjacent face based on the check that if two adjacent faces of
           region are used by the region in the same sense, then their common
           edge should be used by the two faces in opposite senses (or the opposite
           of both the conditions should be true) */

        if ((edir != edir_adj && curdir == rfdir_adj) ||
            (edir == edir_adj && curdir != rfdir_adj)) {
          rfarray2[nrf2] = eface;
          rfdirs2[nrf2] = rfdir_adj;
          List_Add(processed_faces,rfarray2[nrf2]);
          nrf2++;
          break;
        }
      }
      List_Delete(efaces);
    }
    List_Delete(felist);
  }

  /* Delete the original region */

  MR_Delete(rsplit,0);

  /* Make the two new regions */

  rnew[0] = MR_New(mesh);
  MR_Set_GEntDim(rnew[0],3);
  MR_Set_GEntID(rnew[0],gid);
  MR_Set_Faces(rnew[0],nrf1,rfarray1,rfdirs1);

  rnew[1] = MR_New(mesh);
  MR_Set_GEntDim(rnew[1],3);
  MR_Set_GEntID(rnew[1],gid);
  MR_Set_Faces(rnew[1],nrf2,rfarray2,rfdirs2);

  List_Delete(processed_faces);

  return fnew;
}

#ifdef __cplusplus
}
#endif
