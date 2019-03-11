/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#define _H_MEdge_Private

#include <math.h>
#include "MEdge.h"
#include "MEdge_jmp.h"
#include "MSTK_private.h"


#ifdef __cplusplus
extern "C" {
#endif

  double ME_LenSqr(MEdge_ptr e) {
    double xyz0[3],xyz1[3],len2;

    MV_Coords(e->vertex[0],xyz0);
    MV_Coords(e->vertex[1],xyz1);

    len2 = ((xyz1[0]-xyz0[0])*(xyz1[0]-xyz0[0])+
	    (xyz1[1]-xyz0[1])*(xyz1[1]-xyz0[1])+
	    (xyz1[2]-xyz0[2])*(xyz1[2]-xyz0[2]));
    
    return len2;
  }

  double ME_Len(MEdge_ptr e) {
    return sqrt(ME_LenSqr(e));
  }

  void ME_Vec(MEdge_ptr e, double *evec) {
    double xyz0[3],xyz1[3];
    int i;

    MV_Coords(e->vertex[0],xyz0);
    MV_Coords(e->vertex[1],xyz1);

    for (i = 0; i < 3; i++)
      evec[i] = xyz1[i]-xyz0[i];
  }

#ifdef __cplusplus
}
#endif

