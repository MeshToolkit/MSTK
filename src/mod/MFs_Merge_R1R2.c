#include <stdlib.h>
#include <stdio.h>
#include "MFace_jmp.h"
#include "MSTK_private.h"


MFace_ptr MFs_Merge_R1R2(MFace_ptr f1, MFace_ptr f2) {

  /* Faces are not explicit entities in reduced representations -
     nothing to do */

  MF_Delete(f2,0);
  return f1;
}
