#include <stdlib.h>
#include <stdio.h>
#include "MVertex_jmp.h"
#include "MSTK_private.h"


MVertex_ptr MVs_Merge(MVertex_ptr v1, MVertex_ptr v2) {
  RepType rtype = MESH_RepType(MV_Mesh(v1));

  return (*MVs_Merge_jmp[rtype])(v1,v2);
}
