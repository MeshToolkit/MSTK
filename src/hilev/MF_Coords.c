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

  fverts = MF_Vertices(f,1);
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
