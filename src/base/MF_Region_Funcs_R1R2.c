#define _H_MFace_Private

#include "MFace.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  Set_ptr MF_Regions_R1R2(MFace_ptr f) {
    MSTK_Report("MF_Regions_R1R2","Not implemented",ERROR);
    return 0;    
  }

  MRegion_ptr MF_Region_R1R2(MFace_ptr f, int dir) {
    MSTK_Report("MF_Region_R1R2","Not implemented",ERROR);
    return 0;
  }


#ifdef __cplusplus
}
#endif
