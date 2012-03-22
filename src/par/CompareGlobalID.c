#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  int compareGlobalID(MEntity_ptr *a, MEntity_ptr *b) {
    if ( MEnt_GlobalID(*a) > MEnt_GlobalID(*b) ) 
      return 1;
    else if ( MEnt_GlobalID(*a) < MEnt_GlobalID(*b) ) 
      return -1;
    else
      return 0;
  }


  int compareVertexCoor(MVertex_ptr *a, MVertex_ptr *b) {
    double coor1[3], coor2[3];
    double tol = 1e-8;
    int i;
    MV_Coords(*a,coor1);
    MV_Coords(*b,coor2);
    for(i = 0; i < 3; i++) {
      if ( (coor1[i] - coor2[i]) > tol )
	return 1;
      if ( (coor2[i] - coor1[i]) > tol )
	return -1;
    }
    return 0;
  }

#ifdef __cplusplus
}
#endif
