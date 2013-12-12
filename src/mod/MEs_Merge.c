#include <stdlib.h>
#include <stdio.h>
#include "MEdge_jmp.h"
#include "MSTK_private.h"


MEdge_ptr MEs_Merge(MEdge_ptr e1, MEdge_ptr e2) {
  RepType rtype = MESH_RepType(ME_Mesh(e1));

  return (*MEs_Merge_jmp[rtype])(e1,e2);
}
