#define _H_MRegion_Private

#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif


  List_ptr MR_Edges_R3R4(MRegion_ptr r) {
    int i, j, idx, n, fdir, fecheck, nf, nfv, found;
    MFace_ptr face;
    MEdge_ptr edge;
    MVertex_ptr v0, v1, ev0, ev1;
    List_ptr redges, fverts;
    MRegion_DownAdj_FN *downadj;
    Mesh_ptr mesh = MEnt_Mesh(r);

    downadj = (MRegion_DownAdj_FN *) r->downadj;
    nf = List_Num_Entries(downadj->rfaces);

    switch (nf) {
    case 4: /* Tet */
      face = List_Entry(downadj->rfaces,0); /* first face */
      fdir = downadj->fdirs[0] & 1;    /* Sense in which face is used in region */

      redges = MF_Edges(face,!fdir,0);
      n = 3;

      face = List_Entry(downadj->rfaces,1);
      fverts = MF_Vertices(face,1,0);

      for (i = 0; i < 3 && n < 5; i++) {
	v0 = List_Entry(fverts,i);
	v1 = List_Entry(fverts,(i+1)%3);
            
	found = 0; idx = 0;
	while ((edge = List_Next_Entry(redges,&idx))) {
	  ev0 = ME_Vertex(edge,0);
	  ev1 = ME_Vertex(edge,1);
	  if ((v0 == ev0 && v1 == ev1) || (v0 == ev1 && v1 == ev0)) {
	    found = 1;
	    break;
	  }
	}
	
	if (!found) {
	  edge = ME_New(mesh);
	  MEnt_Set_Volatile(edge);
	  ME_Set_Vertex(edge,0,v0);
	  ME_Set_Vertex(edge,1,v1);
	  ME_Set_GInfo_Auto(edge);
	  List_Add(redges,edge);
	  n++;
	}
      }
      List_Delete(fverts);

      face = List_Entry(downadj->rfaces,2);
      fverts = MF_Vertices(face,1,0);

      for (i = 0; i < 3 && n < 6; i++) {
	v0 = List_Entry(fverts,i);
	v1 = List_Entry(fverts,(i+1)%3);
            
	found = 0; idx = 0;
	while ((edge = List_Next_Entry(redges,&idx))) {
	  ev0 = ME_Vertex(edge,0);
	  ev1 = ME_Vertex(edge,1);
	  if ((v0 == ev0 && v1 == ev1) || (v0 == ev1 && v1 == ev0)) {
	    found = 1;
	    break;
	  }
	}
	
	if (!found) {
	  edge = ME_New(mesh);
	  MEnt_Set_Volatile(edge);
	  ME_Set_Vertex(edge,0,v0);
	  ME_Set_Vertex(edge,1,v1);
	  ME_Set_GInfo_Auto(edge);
	  List_Add(redges,edge);
	  n++;
	}
      }
      List_Delete(fverts);
      
      return redges;
      break;
    case 6: /* Hex ? */

      /* All faces must have 4 edges each */
      fecheck = 1;
      for (i = 0; i < nf; i++) {
	face = List_Entry(downadj->rfaces,i);
	if (MF_Num_Edges(face) != 4)
	  fecheck = 0;
      }

      if (!fecheck)
	break;

      n = 0;

      /* Add edges of first face */
      face = List_Entry(downadj->rfaces,0); /* first face */
      fdir = downadj->fdirs[0] & 1;    /* Sense in which face is used in region */

      redges = MF_Edges(face,!fdir,0);
      n = 4;

      for (i = 1; i < nf-1 && n < 12; i++) {
	face = List_Entry(downadj->rfaces,i);

	fverts = MF_Edges(face,1,0);

	for (j = 0; j < 4 && n < 12; j++) {
	  v0 = List_Entry(fverts,j);
	  v1 = List_Entry(fverts,(j+1)%4);

	  found = 0; idx = 0;
	  while ((edge = List_Next_Entry(redges,&idx))) {
	    ev0 = ME_Vertex(edge,0);
	    ev1 = ME_Vertex(edge,1);
	    if ((v0 == ev0 && v1 == ev1) || (v0 == ev1 && v1 == ev0)) {
	      found = 1;
	      break;
	    }
	  }
	  
	  if (!found) {
	    edge = ME_New(mesh);
	    MEnt_Set_Volatile(edge);
	    ME_Set_Vertex(edge,0,v0);
	    ME_Set_Vertex(edge,1,v1);
	    ME_Set_GInfo_Auto(edge);
	    List_Add(redges,edge);
	    n++;
	  }
	}
	List_Delete(fverts);
      }

      return redges;
      break;
    default:
      break;
    }

    /* Pyramids, Prisms, General Polyhedra */
    /* We should do separate procedures for pyramids and prisms */
    
    /* Add edges of first face */
    face = List_Entry(downadj->rfaces,0); /* first face */
    fdir = downadj->fdirs[0] & 1;   /* Sense in which face is used in region */
    
    redges = MF_Edges(face,!fdir,0);
    
    for (i = 1; i < nf-1; i++) {
      face = List_Entry(downadj->rfaces,i);

      fverts = MF_Edges(face,1,0);
      nfv = List_Num_Entries(fverts);

      for (j = 0; j < nfv; j++) {
	v0 = List_Entry(fverts,j);
	v1 = List_Entry(fverts,(j+1)%nfv);

	found = 0; idx = 0;
	while ((edge = List_Next_Entry(redges,&idx))) {
	  ev0 = ME_Vertex(edge,0);
	  ev1 = ME_Vertex(edge,1);
	  if ((v0 == ev0 && v1 == ev1) || (v0 == ev1 && v1 == ev0)) {
	    found = 1;
	    break;
	  }
	}
	  
	if (!found) {
	  edge = ME_New(mesh);
	  MEnt_Set_Volatile(edge);
	  ME_Set_Vertex(edge,0,v0);
	  ME_Set_Vertex(edge,1,v1);
	  ME_Set_GInfo_Auto(edge);
	  List_Add(redges,edge);
	  n++;
	}
      }
      List_Delete(fverts);

    }
    
    return redges;
  }

  int MR_UsesEdge_R3R4(MRegion_ptr r, MEdge_ptr e) {
    int i, idx;
    MFace_ptr face;
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;

    idx = 0;
    while ((face = List_Next_Entry(downadj->rfaces,&idx))) {
      if (MF_UsesEntity(face,e,1))
	return 1;
    }
    return 0;
  }

#ifdef __cplusplus
}
#endif
