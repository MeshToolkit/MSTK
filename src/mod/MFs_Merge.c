#include <stdlib.h>
#include <stdio.h>
#include "MFace_jmp.h"
#include "MSTK_private.h"

/* topoflag = 1 means respect model topology (as given by mesh entity
   classification) and do not allow dimensional reduction in the mesh

   topoflag = 0 means the caller does not care about classification
   and would like the two faces to be merged as requested */

MFace_ptr MFs_Merge(MFace_ptr f1, MFace_ptr f2, int topoflag) {
  RepType rtype = MESH_RepType(MF_Mesh(f1));

  return (*MFs_Merge_jmp[rtype])(f1,f2,topoflag);
}
