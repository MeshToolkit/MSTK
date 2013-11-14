#include <stdio.h>
#include <math.h>
#include "MSTK_VecFuncs.h"

#ifdef __cplusplus
extern "C" {
#endif

void MSTK_VDiff3(double *a, double *b, double *c) {
  c[0] = a[0]-b[0]; c[1] = a[1]-b[1]; c[2] = a[2]-b[2];
}

void MSTK_VSum3(double *a, double *b, double *c) {
  c[0] = a[0]+b[0]; c[1] = a[1]+b[1]; c[2] = a[2]+b[2];
}

void MSTK_VScale3(double *a, double s) {
  a[0] *= s; a[1] *= s; a[2] *= s;
}

double MSTK_VLen3(double *a) {
  return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

double MSTK_VLenSqr3(double *a) {
  return (a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

void MSTK_VNormalize3(double *a) {
  double len = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);

  if (len <= 1.0e-28) {
    return;
#ifdef DEBUG
    fprintf(stderr,"Zero Length vector\n");
#endif
  }

  a[0] /= len; a[1] /= len; a[2] /= len;
}

double MSTK_VDot3(double *a, double *b) {
  return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}

void MSTK_VCross3(double *a, double *b, double *x) {
  x[0] = a[1]*b[2]-a[2]*b[1];
  x[1] = a[2]*b[0]-a[0]*b[2];
  x[2] = a[0]*b[1]-a[1]*b[0];
}

void MSTK_VCopy3(double *a, double *b) {
  a[0] = b[0]; a[1] = b[1]; a[2] = b[2];
}

void MSTK_VNeg3(double *a) {
  a[0] = -a[0]; a[1] = -a[1]; a[2] = -a[2];
}

#ifdef __cplusplus
}
#endif
