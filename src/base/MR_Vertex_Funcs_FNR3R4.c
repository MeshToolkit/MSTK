#define _H_MRegion_Private

#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif


  List_ptr MR_Vertices_FNR3R4(MRegion_ptr r) {
    int i, j, n, nf, mkr, found, diradj0=0, diropp=0, fdir, fdir0, fecheck;
    MFace_ptr face=NULL, face0=NULL, fadj0=NULL, fopp=NULL;
    MEdge_ptr edge, fedge00, upedge;
    MVertex_ptr vert, rv0, rvopp0=NULL;
    List_ptr rvertices, fverts, fedges0, adjfedges;
    MRegion_Adj_FN *adj;

    adj = (MRegion_Adj_FN *) r->adj;

    nf = List_Num_Entries(adj->rfaces);
    switch (nf) {
    case 4: /* Tet! */
      /* Add vertices of first face to list of region vertices */
      face0 = List_Entry(adj->rfaces,0); /* first face */
      fdir0 = adj->fdirs[0] & 1;    /* Sense in which face is used in region */
      rvertices = MF_Vertices(face0,!fdir0,0);
      n = 3;

      face = List_Entry(adj->rfaces,1);
      fverts = MF_Vertices(face,1,0);
      for (i = 0; i < 3 && n < 4; i++) { 
	vert = List_Entry(fverts,i);
	if (!List_Contains(rvertices,vert)) {
	  List_Add(rvertices,vert);
	  n++;
	}
      }
      List_Delete(fverts);

      return rvertices;
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

      face0 = List_Entry(adj->rfaces,0); /* face 0 */
      fdir0 = adj->fdirs[0] & 1; /* dir of face w.r.t. region */

      
      /* Add vertices of first face */
      rvertices = MF_Vertices(face0,!fdir0,0);

      /* vertex 0 of region and of face 0 */
      rv0 = List_Entry(rvertices,0);      
      
      /* Edges of face 0 */
      fedges0 = MF_Edges(face0,!fdir0,rv0);

      /* edge 0 of face 0 (wrt face dir pointing into the region) */
      fedge00 = List_Entry(fedges0,0); 

      /* Get the face adjacent to edge 0 of face 0 and also the
	 opposite face, a face that has no edge common face 0 */

      fopp = NULL; fadj0 = NULL;
      for (i = 1; i < nf; i++) {
	face = List_Entry(adj->rfaces,i);

	/* Check if face uses any edge of face 0 */

	for (j = 0, found = 0; j < 4; j++) {
	  edge = List_Entry(fedges0,j);

	  if (MF_UsesEntity(face,edge,1)) {
	    if (edge == fedge00) {
	      /* face uses edge 0 of face 0 (w.r.t. face dir pointing
		 into region) */
	      fadj0 = face;
	      diradj0 = (adj->fdirs[0])>>i & 1;  
	    }
	    found = 1;
	    break;
	  }
	}
	if (!found) {
	  fopp = face;
	  diropp = (adj->fdirs[0])>>i & 1;
	}
	if (fopp && fadj0)
	  break;
      }

      if (!fopp) {
	MSTK_Report("MR_Vertices_FNR3R4","Could not find opposite face",ERROR);
	List_Delete(fedges0);
	List_Delete(rvertices);
	return (void *) NULL;
      }

      /* edges of face adjacent to edge 0 of face 0 */
      adjfedges = MF_Edges(fadj0,1,rv0);

      /* edge connecting rv0 and vertex on opposite face */
      upedge = diradj0 ? List_Entry(adjfedges,3): List_Entry(adjfedges,0);

#ifdef DEBUG
      if (!ME_UsesEntity(upedge,rv0,0))
	MSTK_Report("MR_Vertices_FNR3R4","Cannot find correct vertical edge",
		    ERROR);
#endif
	
      if (ME_Vertex(upedge,0) == rv0)
	rvopp0 = ME_Vertex(upedge,1);
      else {
	rvopp0 = ME_Vertex(upedge,0);
      }

      List_Delete(fedges0);
      List_Delete(adjfedges);

      fverts = MF_Vertices(fopp,diropp,rvopp0);

      for (i = 0; i < 4; i++)
	List_Add(rvertices,List_Entry(fverts,i));
      List_Delete(fverts);

      return rvertices;
      break;
    default: 
      break;
    }


    /* Pyramids, Prisms, General Polyhedra */
    /* We should do separate procedures for pyramids and prisms */
    mkr = MSTK_GetMarker();
    
    /* Add vertices of first face */
    face = List_Entry(adj->rfaces,0); /* first face */
    fdir = adj->fdirs[0] & 1;    /* Sense in which face is used in region */
    
    rvertices = MF_Vertices(face,!fdir,0); 
    List_Mark(rvertices,mkr);
    
    for (i = 1; i < nf-1; i++) {
      face = List_Entry(adj->rfaces,i);
      fverts = MF_Vertices(face,1,0);
      n = List_Num_Entries(fverts);
      for (j = 0; j < n; j++) {
	vert = List_Entry(fverts,j);
	if (!MEnt_IsMarked(vert,mkr)) {
	  List_Add(rvertices,vert);
	  MEnt_Mark(vert,mkr);
	}
      }
      List_Delete(fverts);
    }
    List_Unmark(rvertices,mkr);
    MSTK_FreeMarker(mkr);
    
    return rvertices;
  }

  void MR_Replace_Vertex_FNR3R4(MRegion_ptr r, MVertex_ptr v, MVertex_ptr nuv){
#ifdef DEBUG
    MSTK_Report("MR_Replace_Vertex",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MR_Replace_Vertex_i_FNR3R4(MRegion_ptr r, int i, MVertex_ptr nuv) {
#ifdef DEBUG
    MSTK_Report("MF_Replace_Vertex_i",
		"Function call not suitable for this representation",WARN);
#endif
  }


  int MR_UsesVertex_FNR3R4(MRegion_ptr r, MVertex_ptr v) {
    int i, nf;
    MFace_ptr face;
    MRegion_Adj_FN *adj;

    adj = (MRegion_Adj_FN *) r->adj;
    nf = List_Num_Entries(adj->rfaces);

    for (i = 0; i < nf; i++) {
      face = List_Entry(adj->rfaces,i);
      if (MF_UsesEntity(face,v,MVERTEX))
	return 1;
    }
    return 0;
  }

#ifdef __cplusplus
}
#endif
