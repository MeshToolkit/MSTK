#include <stdio.h>
#include <math.h>
#include "MSTK.h"


/* Condition number of polygon corners */

#ifdef __cplusplus
extern "C" {
#endif

void MF_CondNums(MFace_ptr f, double *conds) {
  double xyz[MAXPV2][3], L10_sqr, L20_sqr, L21_sqr, A;
  double xyz1[3][3], eps=1.0E-14;
  int i, k, nfv;
  MVertex_ptr vtx[MAXPV2];
  Set_ptr fverts;

  fverts = MF_Vertices(f,1,0);
  nfv = Set_Num_Entries(fverts);
  for (i = 0; i < nfv; i++) {
    vtx[i] = Set_Entry(fverts,i);
    MV_Coords(vtx[i],xyz[i]);
  }
  Set_Delete(fverts);
    

  for (i = 0; i < nfv; i++) {
    for (k = 0; k < 3; k++)
      xyz1[0][k] = xyz[i][k];
    for (k = 0; k < 3; k++)
      xyz1[1][k] = xyz[(i+1)%nfv][k];
    for (k = 0; k < 3; k++)
      xyz1[2][k] = xyz[(i+nfv-1)%nfv][k];

    L10_sqr = ((xyz1[1][0]-xyz1[0][0])*(xyz1[1][0]-xyz1[0][0]) + 
	       (xyz1[1][1]-xyz1[0][1])*(xyz1[1][1]-xyz1[0][1]) + 
	       (xyz1[1][2]-xyz1[0][2])*(xyz1[1][2]-xyz1[0][2]));
    
    if (L10_sqr < eps*eps)
      L10_sqr = eps*eps;
      
    L20_sqr = ((xyz1[2][0]-xyz1[0][0])*(xyz1[2][0]-xyz1[0][0]) + 
	       (xyz1[2][1]-xyz1[0][1])*(xyz1[2][1]-xyz1[0][1]) + 
	       (xyz1[2][2]-xyz1[0][2])*(xyz1[2][2]-xyz1[0][2]));
    
    if (L20_sqr < eps*eps)
      L20_sqr = eps*eps;
      
    L21_sqr = ((xyz1[2][0]-xyz1[1][0])*(xyz1[2][0]-xyz1[1][0]) + 
	       (xyz1[2][1]-xyz1[1][1])*(xyz1[2][1]-xyz1[1][1]) + 
	       (xyz1[2][2]-xyz1[1][2])*(xyz1[2][2]-xyz1[1][2]));
    
    if (L21_sqr < eps*eps)
      L21_sqr = eps*eps;
    
    /* modified version of Herron's formula sqrt(s*(s-a)*(s-b)*(s-c)) */
    A = 0.25*sqrt(4*L10_sqr*L20_sqr - (L10_sqr+L20_sqr-L21_sqr)*(L10_sqr+L20_sqr-L21_sqr));
    
    conds[i] = (L10_sqr+L20_sqr)/(2*A);
  }

  return;
}

#ifdef __cplusplus
}
#endif

