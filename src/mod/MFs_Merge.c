#include <stdlib.h>
#include <stdio.h>
#include "MFace_jmp.h"
#include "MSTK_private.h"

MFace_ptr MFs_Merge(MFace_ptr f1, MFace_ptr f2) {
  RepType rtype = MESH_RepType(MF_Mesh(f1));

  return (*MFs_Merge_jmp[rtype])(f1,f2);
}
