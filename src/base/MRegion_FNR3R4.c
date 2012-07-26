#define _H_MRegion_Private

#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MR_Set_RepType_FNR3R4(MRegion_ptr r) {
    MRegion_Adj_FN *adj;

    r->adj = (MRegion_Adj_FN *) MSTK_malloc(sizeof(MRegion_Adj_FN));
    adj = (MRegion_Adj_FN *) r->adj;

    adj->fdirs = NULL;
    adj->rfaces = NULL;
  }

  void MR_Delete_F1F3R3R4(MRegion_ptr r, int keep) {
    MRegion_Adj_FN *adj;
    MFace_ptr f;
    int idx;

    adj = (MRegion_Adj_FN *) r->adj;
      
    if (MEnt_Dim((MEntity_ptr) r) != MDELETED) { /* if regn hasnt been temporarily deleted */
      if (adj) {
        if (adj->rfaces) {
          idx = 0;
          while ((f = List_Next_Entry(adj->rfaces,&idx)))
            MF_Rem_Region(f,r);
        }
      }
    }

    if (!keep) {
      if (adj) {	
	if (adj->rfaces) {
	  MSTK_free(adj->fdirs);
	  List_Delete(adj->rfaces);
	}
	MSTK_free(adj);
      }
    }
  }

  void MR_Restore_F1F3R3R4(MRegion_ptr r) {
    MRegion_Adj_FN *adj;
    MFace_ptr f;
    int i, j, k, nf, side;

    adj = (MRegion_Adj_FN *) r->adj;

    nf = List_Num_Entries(adj->rfaces);
    for (i = 0; i < nf; i++) {
      f = List_Entry(adj->rfaces,i);
 
      j = (int) i/(8*sizeof(unsigned int));
      k = i%(8*sizeof(unsigned int));
      side = !((adj->fdirs[j])>>k & 1);

      MF_Add_Region(f,r,side);
    }
  }

  void MR_Destroy_For_MESH_Delete_FNR3R4(MRegion_ptr r) {
    MRegion_Adj_FN *adj;

    adj = (MRegion_Adj_FN *) r->adj;
    if (adj) {
      if (adj->rfaces) {
	MSTK_free(adj->fdirs);
	List_Delete(adj->rfaces);
      }
      MSTK_free(adj);
    }
  }

  int MR_Num_AdjRegions_FNR3R4(MRegion_ptr r) {
    List_ptr adjr;
    int nr;

#ifdef DEBUG
    MSTK_Report("MR_Num_AdjRegions",
		"Inefficient to call this routine with this representation",
		MSTK_WARN);
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
		"Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void MR_Rem_AdjRegion_FNR3R4(MRegion_ptr r, MRegion_ptr aregion) {
#ifdef DEBUG
    MSTK_Report("MR_Rem_AdjRegion",
		"Function call not suitable for this representation",MSTK_WARN);
#endif
  }


#ifdef __cplusplus
}
#endif
