/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MSTK.h"

/* Routine to split an edge at a given location - two new edges are
   created in place of this edge and the faces connected to the old
   edge are updated to reflect the new topology. The split edge is
   deleted */

#ifdef __cplusplus
extern "C" {
#endif

MVertex_ptr ME_Split(MEdge_ptr esplit, double *xyz) {
  MVertex_ptr vnew, v[2];
  MEdge_ptr enew;
  List_ptr elist;
  double vxyz[3];

#ifdef DEBUG
  v[0] = ME_Vertex(esplit,0);
  v[1] = ME_Vertex(esplit,1);
#endif
  
  if (xyz) {
    memcpy(vxyz,xyz,3*sizeof(double));
    elist = ME_MultiSplit(esplit,1,&vxyz);
  }
  else
    elist = ME_MultiSplit(esplit,1,NULL);

  enew = List_Entry(elist,0);
  vnew = ME_Vertex(enew,1);

#ifdef DEBUG
  if (vnew == v[0] || vnew == v[1]) {
    MSTK_Report("ME_Split","Unexpected failure to retrieve split vertex",
                MSTK_FATAL);
  }
#endif

  List_Delete(elist);

  return vnew;
}

#ifdef __cplusplus
}
#endif
