/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#ifndef _H_MKSTRUC
#define _H_MKSTRUC

/* function to transform coordinates */
/* assume that xyz is 3 dimensional  */
/* 'boundary' = 1 indicates that point is on the domain boundary */
/*            = 0 indicates that point is in domain interior     */

void transform_xyz(double *xyz, int boundary);




#endif
