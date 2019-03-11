/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#ifndef _H_MSTK_VecFuncs
#define _H_MSTK_VecFuncs


/* Vector functions for internal use in MSTK */

#ifdef __cplusplus
extern "C" {
#endif

void MSTK_VDiff3(double *a, double *b, double *c);
void MSTK_VSum3(double *a, double *b, double *c);
void MSTK_VScale3(double *a, double s);
double MSTK_VLen3(double *a);
double MSTK_VLenSqr3(double *a);
void MSTK_VNormalize3(double *a);
double MSTK_VDot3(double *a, double *b);
void MSTK_VCross3(double *a, double *b, double *x);
void MSTK_VCopy3(double *a, double *b);
void MSTK_VNeg3(double *a);

#ifdef __cplusplus
}
#endif

#endif
