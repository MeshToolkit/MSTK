/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdlib.h>
#include <stdio.h>
#include "MEdge_jmp.h"
#include "MSTK_private.h"


MEdge_ptr MEs_Merge_RN(MEdge_ptr e1, MEdge_ptr e2, int topoflag) {

  /* Edges are not explicit entities in reduced representations -
     nothing to do */

  ME_Delete(e2,0);
  return e1;
}
