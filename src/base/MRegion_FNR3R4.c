#define _H_MRegion_Private

#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MR_Set_RepType_FNR3R4(MRegion_ptr r) {
    MRegion_DownAdj_FN *downadj;

    r->downadj = (MRegion_DownAdj_FN *) MSTK_malloc(sizeof(MRegion_DownAdj_FN));
    downadj = (MRegion_DownAdj_FN *) r->downadj;

    downadj->fdirs = NULL;
    downadj->rfaces = NULL;
  }

  void MR_Delete_F1F3R3R4(MRegion_ptr r, int keep) {
    MRegion_DownAdj_FN *downadj;
    MFace_ptr f;
    int idx;

    downadj = (MRegion_DownAdj_FN *) r->downadj;
      
    if (MEnt_Dim(r) != MDELETED) { /* if regn hasnt been temporarily deleted */
      if (downadj) {
	idx = 0;
	while ((f = List_Next_Entry(downadj->rfaces,&idx)))
	  MF_Rem_Region(f,r);
      }
    }

    if (!keep) {
      if (downadj) {	
	if (downadj->rfaces) {
	  MSTK_free(downadj->fdirs);
	  List_Delete(downadj->rfaces);
	}
	MSTK_free(downadj);
      }
    }
  }

  void MR_Restore_F1F3R3R4(MRegion_ptr r) {
    MRegion_DownAdj_FN *downadj;
    MFace_ptr f;
    int i, j, k, nf, side;

    downadj = (MRegion_DownAdj_FN *) r->downadj;

    nf = List_Num_Entries(downadj->rfaces);
    for (i = 0; i < nf; i++) {
      f = List_Entry(downadj->rfaces,i);
 
      j = (int) i/(8*sizeof(unsigned int));
      k = i%(8*sizeof(unsigned int));
      side = !((downadj->fdirs[j])>>k & 1);

      MF_Add_Region(f,r,side);
    }
  }

  void MR_Destroy_For_MESH_Delete_FNR3R4(MRegion_ptr r) {
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;
    if (downadj) {
      if (downadj->rfaces) {
	MSTK_free(downadj->fdirs);
	List_Delete(downadj->rfaces);
      }
      MSTK_free(downadj);
    }
  }

  int MR_Num_AdjRegions_FNR3R4(MRegion_ptr r) {
    List_ptr adjr;
    int nr;

#ifdef DEBUG
    MSTK_Report("MR_Num_AdjRegions",
		"Inefficient to call this routine with this representation",
		WARN);
#endif

    adjr = MR_AdjRegions(r);
    if (adjr) {
      nr = List_Num_Entries(adjr);
      List_Delete(adjr);
      return nr;
    }
    else
      return 0;
  }


  void MR_Add_AdjRegion_FNR3R4(MRegion_ptr r, int facenum, MRegion_ptr aregion) {
#ifdef DEBUG
    MSTK_Report("MR_Add_AdjRegion",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MR_Rem_AdjRegion_FNR3R4(MRegion_ptr r, MRegion_ptr aregion) {
#ifdef DEBUG
    MSTK_Report("MR_Rem_AdjRegion",
		"Function call not suitable for this representation",WARN);
#endif
  }


#ifdef __cplusplus
}
#endif
