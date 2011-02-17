#define _H_MRegion_Private

#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  List_ptr MR_Edges_R2(MRegion_ptr r) {
    MRegion_Adj_R2 *adj;
    int i, j, i2, j2, iv0, iv1, iv02, iv12, nv, ne, nf, nfv, found, nfv2;
    MRType regtype;
    MVertex_ptr v[2], vtmp;
    MEdge_ptr edge;
    List_ptr redges;
    Mesh_ptr mesh = MEnt_Mesh((MEntity_ptr) r);

    adj = (MRegion_Adj_R2 *) r->adj;
    
    ne = 0;
    if (adj->fvtemplate) {
      redges = List_New(0);

      nf = adj->fvtemplate[0][0];

      for (i = 0; i < nf; i++) {
	nfv = adj->fvtemplate[i+1][0];

	for (j = 0; j < nfv; j++) {

	  /* vertex indices for edge to be constructed */

	  iv0 = adj->fvtemplate[i+1][j+1];
	  iv1 = adj->fvtemplate[i+1][(j+1)%nfv+1];

	  /* Check if this edge occurs in a previously processed face */
	  
	  found = 0;
	  for (i2 = 0; i2 < i; i2++) {
	    nfv2 = adj->fvtemplate[i2+1][0];
	    for (j2 = 0; j2 < nfv2; j2++) {
	      iv02 = adj->fvtemplate[i2+1][j2+1];
	      iv12 = adj->fvtemplate[i2+1][(j2+1)%nfv2+1];
	      
	      if ((iv02 == iv1 && iv12 == iv0) ||
		  (iv02 == iv0 && iv12 == iv1)) {
		found = 1;
		break;
	      }		
	    }
	    if (found) break;
	  }

	  if (found) continue;

	  v[0] = List_Entry(adj->rvertices,iv0);
	  v[1] = List_Entry(adj->rvertices,iv1);	  
#ifdef HASHTABLE
	  if (v[0]>v[1]) {
	    vtmp = v[0];
	    v[0] = v[1];
	    v[1] = vtmp;
	  }

	  edge = Hash_Entry(MESH_Hash_Edges(MEnt_Mesh((MEntity_ptr) r)), 2, v);
	  if (edge == NULL) {
	    edge = ME_New(MEnt_Mesh((MEntity_ptr) r));
	    MEnt_Set_Volatile(edge);
	    ME_Set_Vertex(edge,0,v[0]);
	    ME_Set_Vertex(edge,1,v[1]);
	    ME_Set_GInfo_Auto(edge);

	    Hash_Add(MESH_Hash_Edges(MEnt_Mesh((MEntity_ptr) r)), edge, 2, v);
	  }
#else
	  edge = ME_New(MEnt_Mesh((MEntity_ptr)r));
	  MEnt_Set_Volatile(edge);
	  ME_Set_Vertex(edge,0,v[0]);
	  ME_Set_Vertex(edge,1,v[1]);
	  ME_Set_GInfo_Auto(edge);
#endif
	  List_Add(redges,edge);
	  ME_Lock(edge);
	  ne++;
	}
      }
    }
    else {
      nv = List_Num_Entries(adj->rvertices);
      regtype = MSTK_nv2rtype[nv];
      if (regtype == RUNKNOWN)
	return NULL;

      ne = MSTK_nre_template[regtype];
      redges = List_New(ne);

      for (i = 0; i < ne; i++) {

	/* vertex indices for edge to be constructed */
	iv0 = MSTK_rev_template[regtype][i][0];
	iv1 = MSTK_rev_template[regtype][i][1];

	v[0] = List_Entry(adj->rvertices,iv0);
	v[1] = List_Entry(adj->rvertices,iv1);	  
#ifdef HASHTABLE
	  if (v[0]>v[1]) {
	    vtmp = v[0];
	    v[0] = v[1];
	    v[1] = vtmp;
	  }

	edge = Hash_Entry(MESH_Hash_Edges(mesh), 2, v);
	if (edge == NULL) {
	  /* construct the edge */
	  edge = ME_New(mesh);
	  MEnt_Set_Volatile(edge);
	  ME_Set_Vertex(edge,0,v[0]);
	  ME_Set_Vertex(edge,1,v[1]);
	  ME_Set_GInfo_Auto(edge);
	  Hash_Add(MESH_Hash_Edges(mesh), edge, 2, v);
	}
#else
	/* construct the edge */
	edge = ME_New(mesh);
	MEnt_Set_Volatile(edge);
	ME_Set_Vertex(edge,0,v[0]);
	ME_Set_Vertex(edge,1,v[1]);
	ME_Set_GInfo_Auto(edge);
#endif
	List_Add(redges,edge);
	ME_Lock(edge);
      }
    }

    if (!MESH_AutoLock(mesh)) {
       i = 0;
       while ((edge = List_Next_Entry(redges, &i))) {
	 ME_UnLock(edge);
       }
    }

    return redges;
  }

  int MR_UsesEdge_R2(MRegion_ptr r, MEdge_ptr e) {
    MRegion_Adj_R2 *adj;
    int i, j, j0, j1, nv, ne, nf, nfv;
    MRType regtype;
    MVertex_ptr v0, v1;

    adj = (MRegion_Adj_R2 *) r->adj;

    v0 = ME_Vertex(e,0); j0 = List_Locate(adj->rvertices,v0);
    v1 = ME_Vertex(e,1); j1 = List_Locate(adj->rvertices,v1);

    if (adj->fvtemplate) {
      nf = adj->fvtemplate[0][0];
      for (i = 0; i < nf; i++) {
	nfv = adj->fvtemplate[i+1][0];
	for (j = 0; j < nfv; j++) {
	  if ((adj->fvtemplate[i+1][j+1] == j0 &&
	       adj->fvtemplate[i+1][(j+1)%nfv+1] == j1) 
	      ||
	      (adj->fvtemplate[i+1][j+1] == j1 &&
	       adj->fvtemplate[i+1][(j+1)%nfv+1] == j0))
	    return 1;
	}
      }
    }
    else {
      nv = List_Num_Entries(adj->rvertices);      
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
