#include <stdlib.h>
#include <stdio.h>
#include "MEdge_jmp.h"
#include "MSTK_private.h"

/* topoflag = 1 means respect model topology (as given by mesh entity
   classification) and do not allow dimensional reduction in the mesh

   topoflag = 0 means the caller does not care about classification
   and would like the two edges to be merged as requested */

MEdge_ptr MEs_Merge(MEdge_ptr e1, MEdge_ptr e2, int topoflag) {
  RepType rtype = MESH_RepType(ME_Mesh(e1));

  return (*MEs_Merge_jmp[rtype])(e1,e2,topoflag);
}
