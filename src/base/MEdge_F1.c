#define _H_MEdge_Private

#include "MEdge.h"
#include "MEdge_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void ME_Set_RepType_F1(MEdge_ptr e) {
    MEdge_UpAdj_F1 *upadj;

    upadj = e->upadj = (MEdge_UpAdj_F1 *) MSTK_malloc(sizeof(MEdge_UpAdj_F1));
    upadj->nf = 0;
    upadj->efaces = Set_New(10);
  }

  void ME_Delete_F1(MEdge_ptr e) {
    MEdge_UpAdj_F1 *upadj;

    upadj = (MEdge_UpAdj_F1 *) e->upadj;
    Set_Delete(upadj->efaces);
    MSTK_free(upadj);

    MV_Rem_Edge(e->vertex[0],e);
    MV_Rem_Edge(e->vertex[1],e);

    MSTK_free(e);
  }

  int ME_Num_Faces_F1(MEdge_ptr e) {
    return ((MEdge_UpAdj_F1 *)e->upadj)->nf;
  }

  int ME_Num_Regions_F1(MEdge_ptr e) {
    Set_ptr eregs;
    int nr;

#ifdef DEBUG
    MSTK_Report("ME_Num_Regions",
		"Inefficient to call this routine with this representation",
		MESG);
#endif
    
    eregs = ME_Regions_F1(e);
    nr = eregs ? Set_Num_Entries(eregs) : 0;
    if (eregs) Set_Delete(eregs);

    return nr;
  }

  Set_ptr ME_Faces_F1(MEdge_ptr e) {
    return Set_Copy(((MEdge_UpAdj_F1 *)e->upadj)->efaces);
  }

  Set_ptr ME_Regions_F1(MEdge_ptr e) {
    MEdge_UpAdj_F1 *upadj; 
    int i, j, nr, mkr;
    Set_ptr eregs;
    MFace_ptr eface;
    MRegion_ptr freg;
 
    upadj = (MEdge_UpAdj_F1 *) e->upadj;
    nr = 0;
    eregs = Set_New(10);
    mkr = MSTK_GetMarker();

    for (i = 0; i < upadj->nf; i++) {
      eface = Set_Entry(upadj->efaces,i);
      for (j = 0; j < 2; j++) {
	freg = MF_Region(eface,j);
	if (freg && !MEnt_IsMarked(freg,mkr)) {
	  MEnt_Mark(freg,mkr);
	  Set_Add(eregs,freg);
	  nr++;
	}
      }
    }
    Set_Unmark(eregs,mkr);
    MSTK_FreeMarker(mkr);
    
    if (nr) 
      return eregs;
    else {
      Set_Delete(eregs);
      return 0;
    }
  }


  void ME_Add_Face_F1(MEdge_ptr e, MFace_ptr f) {
    MEdge_UpAdj_F1 *upadj;

    upadj = (MEdge_UpAdj_F1 *) e->upadj;

    if (upadj->efaces == NULL)
      upadj->efaces = Set_New(10);
    Set_Add(upadj->efaces,f);
    (upadj->nf)++;
  }

  void ME_Rem_Face_F1(MEdge_ptr e, MFace_ptr f) {
    MEdge_UpAdj_F1 *upadj; 
    int ok;

    upadj = (MEdge_UpAdj_F1 *) e->upadj;

    ok = Set_Rem(upadj->efaces,f);
    if (ok) (upadj->nf)--;
  }

  void ME_Add_Region_F1(MEdge_ptr e, MRegion_ptr r) {
#ifdef DEBUG
    MSTK_Report("ME_Add_Region","Function call not suitable for this representation",WARN);
#endif
  }

  void ME_Rem_Region_F1(MEdge_ptr e, MRegion_ptr r) {
#ifdef DEBUG
    MSTK_Report("ME_Rem_Region","Function call not suitable for this representation",WARN);
#endif
  }

#ifdef __cplusplus
}
#endif

