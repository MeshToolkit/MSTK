#define _H_MEdge_Private

#include <math.h>
#include "MEdge.h"
#include "MEdge_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"


#ifdef __cplusplus
extern "C" {
#endif

  double ME_LenSqr(MEdge_ptr e) {
    double xyz0[3],xyz1[3],len2;

    MV_Coords(e->vertex[0],xyz0);
    MV_Coords(e->vertex[1],xyz1);

    len2 = ((xyz1[0]-xyz0[0])*(xyz1[0]-xyz0[0])+
	    (xyz1[0]-xyz0[0])*(xyz1[0]-xyz0[0])+
	    (xyz1[0]-xyz0[0])*(xyz1[0]-xyz0[0]));
    
    return len2;
  }

  double ME_Len(MEdge_ptr e) {
    return sqrt(ME_LenSqr(e));
  }

  void ME_Vec(MEdge_ptr e, double *evec) {
    double xyz0[3],xyz1[3];

    MV_Coords(e->vertex[0],xyz0);
    MV_Coords(e->vertex[1],xyz1);

    VDiff3(xyz1,xyz0,evec);
  }

#ifdef __cplusplus
}
#endif

