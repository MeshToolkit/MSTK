#define _H_MFace_Private

#include "MFace.h"
#include "MSTK_malloc.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  int MF_Num_Vertices_R1R2(MFace_ptr f) {
    MSTK_Report("MF_Num_Vertices_R1R2","Not implemented",ERROR);
    return 0;
  }

  Set_ptr MF_Vertices_R1R2(MFace_ptr f, int dir) {
    MSTK_Report("MF_Vertices_R1R2","Not implemented",ERROR);
    return 0;
  }
	

  int MF_UsesVertex_R1R2(MFace_ptr f, MVertex_ptr v) {
    MSTK_Report("MF_UsesVertex_R1R2","Not implemented",ERROR);
    return 0;
  }


#ifdef __cplusplus
}
#endif
