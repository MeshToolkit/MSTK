/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdlib.h>
#include <stdio.h>
#include "MVertex_jmp.h"
#include "MSTK_private.h"

/* topoflag = 1 means respect model topology given by mesh entity
   classification and don't allow for dimensional reduction
   
   topoflag = 0 means the caller doesn't care about classifition and
   would like to merge the two vertices under all circumstances */

MVertex_ptr MVs_Merge(MVertex_ptr v1, MVertex_ptr v2, int topoflag) {
  RepType rtype = MESH_RepType(MV_Mesh(v1));

  return (*MVs_Merge_jmp[rtype])(v1,v2,topoflag);
}
