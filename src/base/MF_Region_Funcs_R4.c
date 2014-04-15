#define _H_MFace_Private

#include "MFace.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif


  List_ptr MF_Regions_R4(MFace_ptr f) {
    List_ptr fregs;
    MFace_Adj_R3 *adj;

    adj = (MFace_Adj_R3 *) f->adj;

    fregs = List_New(2);
    if (adj->fregions[0])
      List_Add(fregs,adj->fregions[0]);
    if (adj->fregions[1])
      List_Add(fregs,adj->fregions[1]);

    return fregs;
  }

  void MF_RegionIDs_R4(MFace_ptr f, int *nfr, int *fregids) {
    MFace_Adj_F1F3 *adj = (MFace_Adj_F1F3 *) f->adj;

    *nfr = 0;
    if (adj->fregions[0])
      fregids[(*nfr)++] = MEnt_ID(adj->fregions[0]);
    if (adj->fregions[1])
      fregids[(*nfr)++] = MEnt_ID(adj->fregions[1]);

  }


  MRegion_ptr MF_Region_R4(MFace_ptr f, int i) {
    return ((MFace_Adj_R3 *) f->adj)->fregions[i];
  }

  int MF_RegionID_R4(MFace_ptr f, int i) {
    return MEnt_ID(((MFace_Adj_R4 *) f->adj)->fregions[i]);
  }

  void MF_Add_Region_R4(MFace_ptr f, MRegion_ptr r, int side) {
    MFace_Adj_R3 *adj;

    adj = (MFace_Adj_R3 *) f->adj;

    if (side == MSTK_UNKNOWN) {
      if (MR_FaceDir(r,f))
	adj->fregions[0] = r;
      else
	adj->fregions[1] = r;
    }
    else
      adj->fregions[side] = r;
  }

  void MF_Rem_Region_R4(MFace_ptr f, MRegion_ptr r) {
    MFace_Adj_R3 *adj;

    adj = (MFace_Adj_R3 *) f->adj;

    if (adj->fregions[0] == r)
      adj->fregions[0] = NULL;
    else if (adj->fregions[1] == r)
      adj->fregions[1] = NULL;
    else
      MSTK_Report("MF_Rem_Region_R4","Region not connected to face",MSTK_ERROR);
  }


#ifdef __cplusplus
}
#endif
