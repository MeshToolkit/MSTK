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
    upadj->efaces = List_New(2);
  }

  void ME_Delete_F1(MEdge_ptr e, int keep) {
    MEdge_UpAdj_F1 *upadj;

    if (MEnt_Dim(e) != MDELETED) { /* if edge hasnt been temporarily deleted */
      MV_Rem_Edge(e->vertex[0],e);
      MV_Rem_Edge(e->vertex[1],e);
    }

    if (!keep) {
      upadj = (MEdge_UpAdj_F1 *) e->upadj;
      if (upadj) {
	if (upadj->efaces)
	  List_Delete(upadj->efaces);
	MSTK_free(upadj);
      }
    }
  }

  void ME_Restore_F1(MEdge_ptr e) {

    MEnt_Set_Dim(e,MEDGE);
    MV_Add_Edge(e->vertex[0],e);
    MV_Add_Edge(e->vertex[1],e);
  }

  void ME_Destroy_For_MESH_Delete_F1(MEdge_ptr e) {
    MEdge_UpAdj_F1 *upadj;

    upadj = (MEdge_UpAdj_F1 *) e->upadj;
    if (upadj) {
      if (upadj->efaces)
	List_Delete(upadj->efaces);
      MSTK_free(upadj);
    }
  }

  int ME_Num_Faces_F1(MEdge_ptr e) {
    List_ptr efaces = ((MEdge_UpAdj_F1 *)e->upadj)->efaces;
    return List_Num_Entries(efaces);
  }

  int ME_Num_Regions_F1(MEdge_ptr e) {
    List_ptr eregs;
    int nr;

#ifdef DEBUG
    MSTK_Report("ME_Num_Regions",
		"Inefficient to call this routine with this representation",
		MESG);
#endif
    
    eregs = ME_Regions_F1(e);
    nr = eregs ? List_Num_Entries(eregs) : 0;
    if (eregs) List_Delete(eregs);

    return nr;
  }

  List_ptr ME_Faces_F1(MEdge_ptr e) {
    List_ptr efaces = ((MEdge_UpAdj_F1 *)e->upadj)->efaces; 
    if (efaces && List_Num_Entries(efaces))
      return List_Copy(efaces);
    else
      return NULL;
  }

  List_ptr ME_Regions_F1(MEdge_ptr e) {
    MEdge_UpAdj_F1 *upadj; 
    int i, j, nr, nf, mkr;
    List_ptr eregs;
    MFace_ptr eface;
    MRegion_ptr freg;
 
    upadj = (MEdge_UpAdj_F1 *) e->upadj;
    nf = List_Num_Entries(upadj->efaces);
    eregs = List_New(nf);
   
    mkr = MSTK_GetMarker();
    nr = 0;
    for (i = 0; i < nf; i++) {
      eface = List_Entry(upadj->efaces,i);

      for (j = 0; j < 2; j++) {
	freg = MF_Region(eface,j);
	if (freg && !MEnt_IsMarked(freg,mkr)) {
	  MEnt_Mark(freg,mkr);
	  List_Add(eregs,freg);
	  nr++;
	}
      }
    }
    List_Unmark(eregs,mkr);
    MSTK_FreeMarker(mkr);
    
    if (nr) 
      return eregs;
    else {
      List_Delete(eregs);
      return 0;
    }
  }


  void ME_Add_Face_F1(MEdge_ptr e, MFace_ptr f) {
    MEdge_UpAdj_F1 *upadj;

    upadj = (MEdge_UpAdj_F1 *) e->upadj;

    if (upadj->efaces == NULL)
      upadj->efaces = List_New(10);
    List_Add(upadj->efaces,f);
  }

  void ME_Rem_Face_F1(MEdge_ptr e, MFace_ptr f) {
    MEdge_UpAdj_F1 *upadj; 
    int ok;

    upadj = (MEdge_UpAdj_F1 *) e->upadj;

    List_Rem(upadj->efaces,f);
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

