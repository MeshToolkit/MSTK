#define _H_MEdge_Private

#include "MEdge.h"
#include "MEdge_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void ME_Set_RepType_F1(MEdge_ptr e) {
    MEdge_Adj_F1 *adj;

    adj = e->adj = (MEdge_Adj_F1 *) MSTK_malloc(sizeof(MEdge_Adj_F1));
    adj->efaces = List_New(2);
  }

  void ME_Delete_F1(MEdge_ptr e, int keep) {
    MEdge_Adj_F1 *adj;

    if (MEnt_Dim((MEntity_ptr) e) != MDELETED) { /* if edge hasnt been temporarily deleted */
      if (e->vertex[0]) MV_Rem_Edge(e->vertex[0],e);
      if (e->vertex[1]) MV_Rem_Edge(e->vertex[1],e);
    }

    if (!keep) {
      adj = (MEdge_Adj_F1 *) e->adj;
      if (adj) {
	if (adj->efaces)
	  List_Delete(adj->efaces);
	MSTK_free(adj);
      }
    }
  }

  void ME_Restore_F1(MEdge_ptr e) {

    MEnt_Set_Dim((MEntity_ptr) e,MEDGE);
    MV_Add_Edge(e->vertex[0],e);
    MV_Add_Edge(e->vertex[1],e);
  }

  void ME_Destroy_For_MESH_Delete_F1(MEdge_ptr e) {
    MEdge_Adj_F1 *adj;

    adj = (MEdge_Adj_F1 *) e->adj;
    if (adj) {
      if (adj->efaces)
	List_Delete(adj->efaces);
      MSTK_free(adj);
    }
  }

  int ME_Num_Faces_F1(MEdge_ptr e) {
    List_ptr efaces = ((MEdge_Adj_F1 *)e->adj)->efaces;
    return List_Num_Entries(efaces);
  }

  int ME_Num_Regions_F1(MEdge_ptr e) {
    List_ptr eregs;
    int nr;

#ifdef DEBUG
    MSTK_Report("ME_Num_Regions",
		"Inefficient to call this routine with this representation",
		MSTK_MESG);
#endif
    
    eregs = ME_Regions_F1(e);
    nr = eregs ? List_Num_Entries(eregs) : 0;
    if (eregs) List_Delete(eregs);

    return nr;
  }

  List_ptr ME_Faces_F1(MEdge_ptr e) {
    List_ptr efaces = ((MEdge_Adj_F1 *)e->adj)->efaces; 
    if (efaces && List_Num_Entries(efaces))
      return List_Copy(efaces);
    else
      return NULL;
  }

  void ME_FaceIDs_F1(MEdge_ptr e, int *nef, int *efaceids) {
    List_ptr efaces = ((MEdge_Adj_F1 *)e->adj)->efaces; 
    *nef = 0;
    if (efaces) {
      int i;
      *nef = List_Num_Entries(efaces);
      for (i = 0; i < *nef; i++) 
        efaceids[i] = MEnt_ID(List_Entry(efaces,i));
    }
  }

  List_ptr ME_Regions_F1(MEdge_ptr e) {
    MEdge_Adj_F1 *adj; 
    int i, j, nr, nf, mkr;
    List_ptr eregs;
    MFace_ptr eface;
    MRegion_ptr freg;
 
    adj = (MEdge_Adj_F1 *) e->adj;
    nf = List_Num_Entries(adj->efaces);
    eregs = List_New(nf);
   
    mkr = MSTK_GetMarker();
    nr = 0;
    for (i = 0; i < nf; i++) {
      eface = List_Entry(adj->efaces,i);

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

  void ME_RegionIDs_F1(MEdge_ptr e, int *ner, int *eregionids) {
    int i;
    List_ptr eregs = ME_Regions_F1(e);
    *ner = 0;
    if (eregs) {
      *ner = List_Num_Entries(eregs);
      for (i = 0; i < *ner; i++)
        eregionids[i] = MEnt_ID(List_Entry(eregs,i));
      List_Delete(eregs);
    }
  }

  void ME_Add_Face_F1(MEdge_ptr e, MFace_ptr f) {
    MEdge_Adj_F1 *adj;

    adj = (MEdge_Adj_F1 *) e->adj;

    if (adj->efaces == NULL)
      adj->efaces = List_New(10);
    List_ChknAdd(adj->efaces,f);
  }

  void ME_Rem_Face_F1(MEdge_ptr e, MFace_ptr f) {
    MEdge_Adj_F1 *adj; 

    adj = (MEdge_Adj_F1 *) e->adj;

    List_Rem(adj->efaces,f);
  }

  void ME_Add_Region_F1(MEdge_ptr e, MRegion_ptr r) {
#ifdef DEBUG
    MSTK_Report("ME_Add_Region","Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void ME_Rem_Region_F1(MEdge_ptr e, MRegion_ptr r) {
#ifdef DEBUG
    MSTK_Report("ME_Rem_Region","Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  MEdge_ptr ME_NextInHash_F1(MEdge_ptr e) {
#ifdef DEBUG
    MSTK_Report("ME_NextInHash", "Function call not suitable for this representation", MSTK_WARN);
#endif
    return NULL;
  }

  void ME_Set_NextInHash_F1(MEdge_ptr e, MEdge_ptr next) {
#ifdef DEBUG
    MSTK_Report("ME_NextInHash", "Function call not suitable for this representation", MSTK_WARN);
#endif
  }

  int ME_IsLocked_F1(MEdge_ptr e) {
    return 0;
  }

#ifdef __cplusplus
}
#endif

