#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif



  int compareINT(const void *a, const void *b) {
    if ( *(int*)a > *(int*)b )  return 1;
    else if ( *(int*)a < *(int*)b )  return -1;
    else return 0;
  }

  int compareGlobalID(const void *a, const void *b) {
     if ( MEnt_GlobalID(*(MEntity_ptr*)a) > MEnt_GlobalID(*(MEntity_ptr*)b) ) return 1 ;
     else if ( MEnt_GlobalID(*(MEntity_ptr*)a) < MEnt_GlobalID(*(MEntity_ptr*)b) ) return -1 ;
     else return 0;
  }


  int compareCoorDouble(const void * a, const void * b) {
    double tol = 1e-8;
    int i;
    double *coor1 = (double *)a;
    double *coor2 = (double *)b;
    for(i = 0; i < 3; i++) {
      if ( (coor1[i] - coor2[i]) > tol )
	return 1;
      if ( (coor2[i] - coor1[i]) > tol )
	return -1;
    }
    return 0;
  }
  
  int compareVertexCoor(const void *a, const void *b) {
    double coor1[3], coor2[3];
    double tol = 1e-8;
    int i;
    MV_Coords(*(MVertex_ptr*)a,coor1);
    MV_Coords(*(MVertex_ptr*)b,coor2);
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
