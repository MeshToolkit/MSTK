#include <stdlib.h>
#include <stdio.h>
#include "MEdge_jmp.h"
#include "MSTK_private.h"


MEdge_ptr MEs_Merge_RN(MEdge_ptr e1, MEdge_ptr e2) {

  /* Edges are not explicit entities in reduced representations -
     nothing to do */

  ME_Delete(e2,0);
  return e1;
}
