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

void MF_Coords(MFace_ptr f, int *n, double (*xyz)[3]) {
  int i;
  MVertex_ptr fv;
  List_ptr fverts;

  fverts = MF_Vertices(f,1,0);
  *n = List_Num_Entries(fverts);

  for (i = 0; i < *n; i++) {
    fv = List_Entry(fverts,i);
    MV_Coords(fv,xyz[i]);
  }
  List_Delete(fverts);
}

#ifdef __cplusplus
}
#endif
