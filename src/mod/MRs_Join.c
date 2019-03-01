/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

/* Join two regions sharing a common face.  The faces of the second
   region (excluding the common face) are added to the first
   region. The common face and the second region are deleted and the
   first region is returned */

#include <stdlib.h>
#include <stdio.h>
#include "MSTK.h"


MRegion_ptr MRs_Join(MRegion_ptr r1, MRegion_ptr r2, MFace_ptr f) {
  int i, j, nrf1, nrf2, gdim, gid, *rfdir2;
  MFace_ptr *rf2, fcmn=f;
  Mesh_ptr   mesh;
  List_ptr   rfaces2;


  mesh = MF_Mesh(r1);
  gid = MF_GEntID(r1);

  if (mesh != MF_Mesh(r2)) {
    MSTK_Report("MRs_Join","Regions not from same mesh",MSTK_ERROR);
    return 0;
  }
  else if (gid != MR_GEntID(r2)) {
    MSTK_Report("MRs_Join","Regions not from same geometric entity",MSTK_ERROR);
    return 0;
  }


  rfaces2 = MR_Faces(r2);
  nrf2 = List_Num_Entries(rfaces2);
  
  if (fcmn) {
    if (!MR_UsesEntity(r1,fcmn,MFACE)) {
      MSTK_Report("MRs_Join","Cannot find common face in region",MSTK_ERROR);
      return 0;
    }
  }
  else { /* find the common face */

    List_ptr rfaces1 = MR_Faces(r1);

    int idx = 0;
    MFace_ptr rf;
    while ((rf = List_Next_Entry(rfaces2,&idx))) {
      if (List_Contains(rfaces1,rf)) {
        fcmn = rf;
        break;
      }
    }

    List_Delete(rfaces1);

  }

  rf2 = (MFace_ptr) malloc(nrf2*sizeof(MFace_ptr));
  rfdir2 = (int *) malloc(nrf2*sizeof(int));

  int found;
  for (i = 0, j = 0, found = 0; i < nrf2; i++) {
    MFace_ptr rface = List_Entry(rfaces2,i);
    if (rface == fcmn) 
      found = 1;
    else {
      rf2[j] = rface;
      rfdir2[j] = MR_FaceDir_i(r2,i);
      j++;
    }
  }
  List_Delete(rfaces2);

  if (!found) {
    MSTK_Report("MRs_Join","Cannot find common face in region",MSTK_ERROR);
    return 0;
  }

  MR_Delete(r2,0);

  MR_Replace_Faces(r1,1,&fcmn,nrf2-1,rf2,rfdir2);
  
  MF_Delete(fcmn,0);

  free(rf2); free(rfdir2);

  return r1;
}
