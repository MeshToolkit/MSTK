/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdlib.h>
#include <stdio.h>
#include "MFace_jmp.h"
#include "MSTK_private.h"


MFace_ptr MFs_Merge_R1R2(MFace_ptr f1, MFace_ptr f2, int topoflag) {

  /* Faces are not explicit entities in reduced representations -
     nothing to do */

  MF_Delete(f2,0);
  return f1;
}
