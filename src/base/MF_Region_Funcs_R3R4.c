#define _H_MFace_Private

#include "MFace.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif


  Set_ptr MF_Regions_R3R4(MFace_ptr f) {
    Set_ptr fregs;
    MFace_UpAdj_R3R4 *upadj;

    upadj = (MFace_UpAdj_R3R4 *) f->upadj;

    fregs = Set_New(2);
    if (upadj->fregions[0])
      Set_Add(fregs,upadj->fregions[0]);
    if (upadj->fregions[1])
      Set_Add(fregs,upadj->fregions[1]);

    return fregs;
  }

  MRegion_ptr MF_Region_R3R4(MFace_ptr f, int i) {
    return ((MFace_UpAdj_R3R4 *) f->upadj)->fregions[i];
  }

  void MF_Add_Region_R3R4(MFace_ptr f, MRegion_ptr r, int side) {
    MFace_UpAdj_R3R4 *upadj;

    upadj = (MFace_UpAdj_R3R4 *) f->upadj;

    if (side == MSTK_UNKNOWN) {
      if (MR_FaceDir(r,f))
	upadj->fregions[0] = r;
      else
	upadj->fregions[1] = r;
    }
    else
      upadj->fregions[side] = r;
  }

  void MF_Rem_Region_R3R4(MFace_ptr f, MRegion_ptr r) {
    MFace_UpAdj_R3R4 *upadj;

    upadj = (MFace_UpAdj_R3R4 *) f->upadj;

    if (upadj->fregions[0] == r)
      upadj->fregions[0] = NULL;
    else if (upadj->fregions[1] == r)
      upadj->fregions[1] = NULL;
    else
      MSTK_Report("MF_Rem_Region_R3R4","Region not connected to face",ERROR);
  }


#ifdef __cplusplus
}
#endif
