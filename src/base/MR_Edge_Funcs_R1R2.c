#define _H_MRegion_Private

#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  List_ptr MR_Edges_R1R2(MRegion_ptr r) {
    MRegion_DownAdj_R1R2 *downadj;
    int i, j, k, i2, j2, iv0, iv1, iv02, iv12, nv, ne, nf, nfv, found, nfv2;
    MRType regtype;
    MVertex_ptr v0, v1;
    MEdge_ptr edge;
    List_ptr redges;
    Mesh_ptr mesh = MEnt_Mesh(r);

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    
    ne = 0;
    if (downadj->fvtemplate) {
      redges = List_New(0);

      nf = downadj->fvtemplate[0][0];

      for (i = 0; i < nf; i++) {
	nfv = downadj->fvtemplate[i+1][0];

	for (j = 0; j < nfv; j++) {

	  /* vertex indices for edge to be constructed */

	  iv0 = downadj->fvtemplate[i+1][j+1];
	  iv1 = downadj->fvtemplate[i+1][(j+1)%nfv+1];

	  /* Check if this edge occurs in a previously processed face */
	  
	  found = 0;
	  for (i2 = 0; i2 < i; i2++) {
	    nfv2 = downadj->fvtemplate[k+1][0];
	    for (j2 = 0; j2 < nfv2; j2++) {
	      iv02 = downadj->fvtemplate[i2+1][j2+1];
	      iv12 = downadj->fvtemplate[i2+1][(j2+1)%nfv2+1];
	      
	      if ((iv02 == iv1 && iv12 == iv0) ||
		  (iv02 == iv0 && iv12 == iv1)) {
		found = 1;
		break;
	      }		
	    }
	    if (found) break;
	  }

	  if (found) continue;

	  edge = ME_New(MEnt_Mesh(r));
	  MEnt_Set_Volatile(edge);
	  v0 = List_Entry(downadj->rvertices,iv0);
	  ME_Set_Vertex(edge,0,v0);
	  v1 = List_Entry(downadj->rvertices,iv1);	  
	  ME_Set_Vertex(edge,0,v1);
	  ME_Set_GInfo_Auto(edge);
	  List_Add(redges,edge);
	  ne++;
	}
      }
    }
    else {
      nv = List_Num_Entries(downadj->rvertices);
      regtype = MSTK_nv2rtype[nv];
      if (regtype == RUNKNOWN)
	return NULL;

      ne = MSTK_nre_template[regtype];
      redges = List_New(ne);

      for (i = 0; i < ne; i++) {

	/* vertex indices for edge to be constructed */
	iv0 = MSTK_rev_template[regtype][i][0];
	iv1 = MSTK_rev_template[regtype][i][1];

	/* construct the edge */
	edge = ME_New(mesh);
	MEnt_Set_Volatile(edge);
	v0 = List_Entry(downadj->rvertices,iv0);
	ME_Set_Vertex(edge,0,v0);
	v1 = List_Entry(downadj->rvertices,iv1);	  
	ME_Set_Vertex(edge,0,v1);
	ME_Set_GInfo_Auto(edge);
	List_Add(redges,edge);
      }
    }

    return NULL;
  }

  int MR_UsesEdge_R1R2(MRegion_ptr r, MEdge_ptr e) {
    MRegion_DownAdj_R1R2 *downadj;
    int i, j, j0, j1, nv, ne, nf, nfv;
    MRType regtype;
    MVertex_ptr v0, v1;

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;

    v0 = ME_Vertex(e,0); j0 = List_Locate(downadj->rvertices,v0);
    v1 = ME_Vertex(e,1); j1 = List_Locate(downadj->rvertices,v1);

    if (downadj->fvtemplate) {
      nf = downadj->fvtemplate[0][0];
      for (i = 0; i < nf; i++) {
	nfv = downadj->fvtemplate[i+1][0];
	for (j = 0; j < nfv; j++) {
	  if ((downadj->fvtemplate[i+1][j+1] == j0 &&
	       downadj->fvtemplate[i+1][(j+1)%nfv+1] == j1) 
	      ||
	      (downadj->fvtemplate[i+1][j+1] == j1 &&
	       downadj->fvtemplate[i+1][(j+1)%nfv+1] == j0))
	    return 1;
	}
      }
    }
    else {
      nv = List_Num_Entries(downadj->rvertices);      
      regtype = MSTK_nv2rtype[nv];
      
      ne = MSTK_nre_template[regtype];
      for (i = 0; i < ne; i++) {
	if ((MSTK_rev_template[regtype][i][0] == j0 && 
	     MSTK_rev_template[regtype][i][1] == j1) 
	    ||
	    (MSTK_rev_template[regtype][i][0] == j1 && 
	     MSTK_rev_template[regtype][i][1] == j0))
	  return 1;
      }
    }

    return 0;
  }

#ifdef __cplusplus
}
#endif
