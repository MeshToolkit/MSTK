#define _H_MEdge_Private

#include "MSTK_malloc.h"
#include "MEdge.h"
#include "MEdge_jmp.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  int ME_Num_Regions_R3R4(MEdge_ptr e) {
    MSTK_Report("ME_Num_Regions_R3R4","Not implemented",ERROR);
    return 0;
  }

  List_ptr ME_Regions_R3R4(MEdge_ptr e) {
    int i, k, mkr, ner, nef;
    List_ptr eregions, efaces;
    MFace_ptr eface;
    MRegion_ptr freg;

    ner = 0;
    eregions = List_New(10);
    mkr = MSTK_GetMarker();

    efaces = ME_Faces(e);
    nef = List_Num_Entries(efaces);

    for (i = 0; i < nef; i++) {
      eface = List_Entry(efaces,i);
      for (k = 0; k < 2; k++) {
	freg = MF_Region(eface,k);
	if (freg && !MEnt_IsMarked(freg,mkr)) {
	  MEnt_Mark(freg,mkr);
	  List_Add(eregions,freg);
	  ner++;
	}
      }
    }
    List_Delete(efaces);
    List_Unmark(eregions,mkr);
    MSTK_FreeMarker(mkr);

    if (ner)
      return eregions;
    else {
      List_Delete(eregions);
      return 0;
    }    
  }

#ifdef __cplusplus
}
#endif

