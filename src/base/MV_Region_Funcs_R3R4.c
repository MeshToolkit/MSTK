#define _H_MVertex_Private

#include "MVertex.h"
#include "MVertex_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Functions */

  int MV_Num_Regions_R3R4(MVertex_ptr v) {
    Set_ptr vregions;
    int nr;

#ifdef DEBUG
    MSTK_Report("MV_Num_Regions_R3R4",
		"Inefficient to call this routine with this representation",
		MESG);
#endif
    vregions = MV_Regions_R3R4(v);
    nr = Set_Num_Entries(vregions);
    Set_Delete(vregions);

    return nr;
  }


  Set_ptr MV_Regions_R3R4(MVertex_ptr v) {
    int i, j, nr, nf, mkr;
    MFace_ptr face;
    MRegion_ptr reg;
    Set_ptr vregions, vfaces;
 
    vregions = Set_New(10);
    nr = 0;
    mkr = MSTK_GetMarker();

    vfaces = MV_Faces_R4(v);
    nf = Set_Num_Entries(vfaces);
    for (i = 0; i < nf; i++) {
      face = Set_Entry(vfaces,i);
      for (j = 0; j < 2; j++) {
	reg = MF_Region(face,j);
	if (reg && !MEnt_IsMarked(reg,mkr)) {
	  MEnt_Mark(reg,mkr);
	  Set_Add(vregions,reg);
	  nr++;
	}
      }
    }
    Set_Delete(vfaces);
    Set_Unmark(vregions,mkr);
    MSTK_FreeMarker(mkr);

    if (nr > 0)
      return vregions;
    else {
      Set_Delete(vregions);
      return 0;
    }
  }

#ifdef __cplusplus
}
#endif
