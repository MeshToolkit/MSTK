#define _H_MRegion_Private

#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif


  List_ptr MR_Edges_FN(MRegion_ptr r) {
    int i, j, n, mkr, fdir, fecheck, nf;
    MFace_ptr face;
    MEdge_ptr edge;
    List_ptr redges, fedges;
    MRegion_Adj_FN *adj;

    adj = (MRegion_Adj_FN *) r->adj;
    nf = List_Num_Entries(adj->rfaces);

    switch (nf) {
    case 4: /* Tet */
      mkr = MSTK_GetMarker();

      face = List_Entry(adj->rfaces,0); /* first face */
      fdir = adj->fdirs[0] & 1;    /* Sense in which face is used in region */

      redges = MF_Edges(face,!fdir,0);
      List_Mark(redges,mkr);
      n = 3;

      face = List_Entry(adj->rfaces,1);
      fedges = MF_Edges(face,1,0);
      for (i = 0; i < 3 && n < 5; i++) { 
	edge = List_Entry(fedges,i);
	if (!MEnt_IsMarked(edge,mkr)) {
	  List_Add(redges,edge);
	  MEnt_Mark(edge,mkr);
	  n++;
	}
      }
      List_Delete(fedges);

      face = List_Entry(adj->rfaces,2);
      fedges = MF_Edges(face,1,0);
      for (i = 0; i < 3 && n < 6; i++) { 
	edge = List_Entry(fedges,i);
	if (!MEnt_IsMarked(edge,mkr)) {
	  List_Add(redges,edge);
	  MEnt_Mark(edge,mkr);
	  n++;
	}
      }
      List_Delete(fedges);

      List_Unmark(redges,mkr);
      MSTK_FreeMarker(mkr);
      
      return redges;
      break;
    case 6: /* Hex ? */

      /* All faces must have 4 edges each */
      fecheck = 1;
      for (i = 0; i < nf; i++) {
	face = List_Entry(adj->rfaces,i);
	if (MF_Num_Edges(face) != 4)
	  fecheck = 0;
      }

      if (!fecheck)
	break;

      n = 0;
      mkr = MSTK_GetMarker();

      /* Add edges of first face */
      face = List_Entry(adj->rfaces,0); /* first face */
      fdir = adj->fdirs[0] & 1;    /* Sense in which face is used in region */

      redges = MF_Edges(face,!fdir,0);
      List_Mark(redges,mkr);
      n = 4;

      for (i = 1; i < nf-1 && n < 12; i++) {
	face = List_Entry(adj->rfaces,i);
	fedges = MF_Edges(face,1,0);
	for (j = 0; j < 4 && n < 12; j++) {
	  edge = List_Entry(fedges,j);
	  if (!MEnt_IsMarked(edge,mkr)) {
	    List_Add(redges,edge);
	    MEnt_Mark(edge,mkr);
	    n++;
	  }
	}
	List_Delete(fedges);
      }
      List_Unmark(redges,mkr);
      MSTK_FreeMarker(mkr);

      return redges;
      break;
    default:
      break;
    }

    /* Pyramids, Prisms, General Polyhedra */
    /* We should do separate procedures for pyramids and prisms */
    mkr = MSTK_GetMarker();
    
    /* Add edges of first face */
    face = List_Entry(adj->rfaces,0); /* first face */
    fdir = adj->fdirs[0] & 1;   /* Sense in which face is used in region */
    
    redges = MF_Edges(face,!fdir,0);
    List_Mark(redges,mkr);
    
    for (i = 1; i < nf-1; i++) {
      face = List_Entry(adj->rfaces,i);
      fedges = MF_Edges(face,1,0);
      n = List_Num_Entries(fedges);
      for (j = 0; j < n; j++) {
	edge = List_Entry(fedges,j);
	if (!MEnt_IsMarked(edge,mkr)) {
	  List_Add(redges,edge);
	  MEnt_Mark(edge,mkr);
	}
      }
      List_Delete(fedges);
    }
    List_Unmark(redges,mkr);
    MSTK_FreeMarker(mkr);
    
    return redges;
  }


  int MR_UsesEdge_FN(MRegion_ptr r, MEdge_ptr e) {
    int i, nf;
    MFace_ptr face;
    MRegion_Adj_FN *adj;

    adj = (MRegion_Adj_FN *) r->adj;
    nf = List_Num_Entries(adj->rfaces);

    for (i = 0; i < nf; i++) {
      face = List_Entry(adj->rfaces,i);
      if (MF_UsesEntity(face,e,1))
	return 1;
    }
    return 0;
  }


#ifdef __cplusplus
}
#endif
