#define _H_MRegion_Private

#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MR_Set_RepType_R1(MRegion_ptr r) {
    MRegion_Adj_R1 *adj;

    r->adj = (MRegion_Adj_R1 *) MSTK_malloc(sizeof(MRegion_Adj_R1));
    adj = r->adj;
    adj->rvertices = NULL;
    adj->fvtemplate = NULL;
  }

  void MR_Delete_R1(MRegion_ptr r, int keep) {
    MRegion_Adj_R1 *adj;
    int idx;
    MVertex_ptr v;

    adj = (MRegion_Adj_R1 *) r->adj;

    if (MEnt_Dim((MEntity_ptr) r) != MDELETED) { /* if regn hasnt been temporarily deleted */
      if (adj) {
	idx = 0;
	while ((v = List_Next_Entry(adj->rvertices,&idx)))
	  MV_Rem_Region(v,r);
      }
    }
    
    if (!keep) {
      if (adj) {
	if (adj->rvertices)
	  List_Delete(adj->rvertices);
	if (adj->fvtemplate) {
	  int nf = adj->fvtemplate[0][0], i;
	  for (i = 0; i < nf; i++)
	    MSTK_free(adj->fvtemplate[i]);
	  MSTK_free(adj->fvtemplate);
	}
	MSTK_free(adj);
      }

      MSTK_free(r);
    }
  }

  void MR_Restore_R1(MRegion_ptr r) {
    MRegion_Adj_R1 *adj;
    int idx;
    MVertex_ptr v;

    adj = (MRegion_Adj_R1 *) r->adj;
    idx = 0;
    while ((v = List_Next_Entry(adj->rvertices,&idx)))
      MV_Add_Region(v,r);
  }

  void MR_Destroy_For_MESH_Delete_R1(MRegion_ptr r) {
    MRegion_Adj_R1 *adj;

    adj = (MRegion_Adj_R1 *) r->adj;
    if (adj) {
      if (adj->rvertices)
	List_Delete(adj->rvertices);
      if (adj->fvtemplate) {
	int nf = adj->fvtemplate[0][0], i;
	for (i = 0; i < nf; i++)
	  MSTK_free(adj->fvtemplate[i]);
	MSTK_free(adj->fvtemplate);
      }
      MSTK_free(adj);
    }
  }


  void MR_Set_Faces_R1(MRegion_ptr r, int nf, MFace_ptr *mfaces, int *dirs) {
#ifdef DEBUG
    MSTK_Report("MR_Set_Faces_R1",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MR_Replace_Face_R1(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf, int dir) {
#ifdef DEBUG
    MSTK_Report("MR_Replace_Face_R1",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MR_Replace_Face_i_R1(MRegion_ptr r, int i, MFace_ptr nuf, int dir) {
#ifdef DEBUG
    MSTK_Report("MR_Replace_Face_i_R1",
		"Function call not suitable for this representation",WARN);
#endif
  }


  int MR_Num_AdjRegions_R1(MRegion_ptr r) {
    List_ptr adjr;
    int nr;

#ifdef DEBUG
    MSTK_Report("MR_Num_AdjRegions",
		"Inefficient to call this routine with this representation",
		WARN);
#endif

    adjr = MR_AdjRegions_R1(r);
    if (adjr) {
      nr = List_Num_Entries(adjr);
      List_Delete(adjr);
      return nr;
    }
    else
      return 0;
  }


  List_ptr MR_AdjRegions_R1(MRegion_ptr r) {
    MSTK_Report("MR_AdjRegions",
		"Not yet implemented for this representation",WARN);
    return NULL;
  }


  void MR_Add_AdjRegion_R1(MRegion_ptr r, int facenum, MRegion_ptr ar) {
#ifdef DEBUG
    MSTK_Report("MR_Add_AdjRegion_R1",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MR_Rem_AdjRegion_R1(MRegion_ptr r, MRegion_ptr ar) {
#ifdef DEBUG
    MSTK_Report("MR_Rem_AdjRegion_R1",
		"Function call not suitable for this representation",WARN);
#endif
  }

#ifdef __cplusplus
}
#endif
