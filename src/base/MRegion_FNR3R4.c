#define _H_MRegion_Private

#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MR_Set_RepType_FNR3R4(MRegion_ptr r) {
    MRegion_DownAdj_FN *downadj;

    r->downadj = (MRegion_DownAdj_FN *) MSTK_malloc(sizeof(MRegion_DownAdj_FN));
    downadj = (MRegion_DownAdj_FN *) r->downadj;

    downadj->nf = (unsigned char) 0;
    downadj->fdirs = 0;
    downadj->rfaces = NULL;
  }

  void MR_Delete_F1F3R3R4(MRegion_ptr r) {
    MRegion_DownAdj_FN *downadj;
    MFace_ptr f;
    int i, nf;

    downadj = (MRegion_DownAdj_FN *) r->downadj;

    nf = Set_Num_Entries(downadj->rfaces);
    for (i = 0; i < nf; i++) {
      f = Set_Entry(downadj->rfaces,i);
      MF_Rem_Region(f,r);
    }
    Set_Delete(downadj->rfaces);
    MSTK_free(downadj);

    MSTK_free(r);
  }

  void MR_Set_Faces_F1F3R3R4(MRegion_ptr r, int nf, MFace_ptr *rfaces,int *dirs){
    int i;
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;
    downadj->nf = nf;
    downadj->fdirs = 0;
    downadj->rfaces = Set_New(nf);

    for (i = 0; i < nf; i++) {
      downadj->fdirs = downadj->fdirs | (dirs[i] << i);
      Set_Add(downadj->rfaces,rfaces[i]);
      MF_Add_Region(rfaces[i],r,!dirs[i]);
    }
  }

  void MR_Set_Vertices_FNR3R4(MRegion_ptr r, int nv, MFace_ptr *mvertices) {
#ifdef DEBUG
    MSTK_Report("MR_Set_Vertices_FNR3R4",
		"Function call not suitable for this representation",WARN);
#endif    
  }

  int MR_Num_Faces_FNR3R4(MRegion_ptr r) {
    return ((MRegion_DownAdj_FN *)r->downadj)->nf;
  }

  int MR_Num_AdjRegions_FNR3R4(MRegion_ptr r) {
    Set_ptr adjr;
    int nr;

#ifdef DEBUG
    MSTK_Report("MR_Num_AdjRegions",
		"Inefficient to call this routine with this representation",
		WARN);
#endif

    adjr = MR_AdjRegions(r);
    if (adjr) {
      nr = Set_Num_Entries(adjr);
      Set_Delete(adjr);
      return nr;
    }
    else
      return 0;
  }

  Set_ptr MR_Vertices_FNR3R4(MRegion_ptr r) {
    int i, j, n, mkr, found, diradj0, diropp, edir, fdir;
    MFace_ptr face, fadj0, fopp;
    MEdge_ptr edge;
    MVertex_ptr vert, rv0, rvopp0;
    Set_ptr rvertices, fverts, fedges;
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;

    switch (downadj->nf) {
    case 4: /* Tet */
      /* Add vertices of first face to list of region vertices */
      face = Set_Entry(downadj->rfaces,0); /* first face */
      fdir = downadj->fdirs & 1;    /* Sense in which face is used in region */
      rvertices = MF_Vertices(face,!fdir,0);
      n = 3;

      face = Set_Entry(downadj->rfaces,1);
      fverts = MF_Vertices(face,1,0);
      for (i = 0; i < 3 && n < 4; i++) { 
	vert = Set_Entry(fverts,i);
	if (!Set_Contains(rvertices,vert)) {
	  Set_Add(rvertices,vert);
	  n++;
	}
      }
      Set_Delete(fverts);

      return rvertices;
      break;
    case 6: /* Hex */
      /* Add vertices of first face */
      face = Set_Entry(downadj->rfaces,0); /* first face */
      fdir = downadj->fdirs & 1;    /* Sense in which face is used in region */
      rvertices = MF_Vertices(face,!fdir,0);
      rv0 = Set_Entry(rvertices,0);
      
      
      fedges = MF_Edges(face,1,rv0);
      
      /* Get the opposite face, a face that has no edge in common with
	 the first face */
      fopp = NULL; fadj0 = NULL;
      for (i = 1; i < downadj->nf; i++) {
	face = Set_Entry(downadj->rfaces,i);
	
	for (j = 0, found = 0; j < 4; j++) {
	  edge = Set_Entry(fedges,j);
	  if (MF_UsesEntity(face,edge,1)) {
	    if (j == 0) {
	      fadj0 = face;
	      diradj0 = (downadj->fdirs)>>i & 1;  
	    }
	    found = 1;
	    break;
	  }
	}
	if (!found) {
	  fopp = face;
	  diropp = (downadj->fdirs)>>i & 1;
	}
	if (fopp && fadj0)
	  break;
      }
      Set_Delete(fedges);

      if (!fopp) {
	MSTK_Report("MR_Vertices_FNR3R4","Could not find opposite face",ERROR);
	Set_Delete(rvertices);
	return (void *) NULL;
      }

      fedges = MF_Edges(fadj0,diradj0,0);
      for (i = 0; i < 4; i++) {
	edge = Set_Entry(fedges,i);
	edir = MF_EdgeDir_i(fadj0,i);
	if (ME_Vertex(edge,edir) == rv0) {
	  rvopp0 = ME_Vertex(edge,!edir);
	  break;
	}
      }
      Set_Delete(fedges);

      fverts = MF_Vertices(fopp,diropp,rvopp0);

      for (i = 0; i < 4; i++)
	Set_Add(rvertices,Set_Entry(fverts,i));
      Set_Delete(fverts);


      return rvertices;
      break;
    default: /* Pyramids, Prisms, General Polyhedra */
      /* We should do separate procedures for pyramids and prisms */
      mkr = MSTK_GetMarker();

      /* Add vertices of first face */
      face = Set_Entry(downadj->rfaces,0); /* first face */
      fdir = downadj->fdirs & 1;    /* Sense in which face is used in region */

      rvertices = MF_Vertices(face,!fdir,0); 
      Set_Mark(rvertices,mkr);

      for (i = 1; i < downadj->nf-1; i++) {
	face = Set_Entry(downadj->rfaces,i);
	fverts = MF_Vertices(face,1,0);
	n = Set_Num_Entries(fverts);
	for (j = 0; j < n; j++) {
	  vert = Set_Entry(fverts,j);
	  if (!MEnt_IsMarked(vert,mkr)) {
	    Set_Add(rvertices,vert);
	    MEnt_Mark(vert,mkr);
	  }
	}
	Set_Delete(fverts);
      }
      Set_Unmark(rvertices,mkr);
      MSTK_FreeMarker(mkr);

      return rvertices;
      break;
    }

    return 0;
  }

  Set_ptr MR_Edges_FNR3R4(MRegion_ptr r) {
    int i, j, n, mkr, fdir;
    MFace_ptr face;
    MEdge_ptr edge;
    Set_ptr redges, fedges;
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;

    switch (downadj->nf) {
    case 4: /* Tet */
      mkr = MSTK_GetMarker();

      face = Set_Entry(downadj->rfaces,0); /* first face */
      fdir = downadj->fdirs & 1;    /* Sense in which face is used in region */

      redges = MF_Edges(face,!fdir,0);
      Set_Mark(redges,mkr);
      n = 3;

      face = Set_Entry(downadj->rfaces,1);
      fedges = MF_Edges(face,1,0);
      for (i = 0; i < 3 && n < 5; i++) { 
	edge = Set_Entry(fedges,i);
	if (!MEnt_IsMarked(edge,mkr)) {
	  Set_Add(redges,edge);
	  MEnt_Mark(edge,mkr);
	  n++;
	}
      }
      Set_Delete(fedges);

      face = Set_Entry(downadj->rfaces,2);
      fedges = MF_Edges(face,1,0);
      for (i = 0; i < 3 && n < 6; i++) { 
	edge = Set_Entry(fedges,i);
	if (!MEnt_IsMarked(edge,mkr)) {
	  Set_Add(redges,edge);
	  MEnt_Mark(edge,mkr);
	  n++;
	}
      }
      Set_Delete(fedges);

      Set_Unmark(redges,mkr);
      MSTK_FreeMarker(mkr);
      
      return redges;
      break;
    case 6: /* Hex */
      n = 0;
      mkr = MSTK_GetMarker();

      /* Add edges of first face */
      face = Set_Entry(downadj->rfaces,0); /* first face */
      fdir = downadj->fdirs & 1;    /* Sense in which face is used in region */

      redges = MF_Edges(face,!fdir,0);
      Set_Mark(redges,mkr);
      n = 4;

      for (i = 1; i < (downadj->nf)-1 && n < 12; i++) {
	face = Set_Entry(downadj->rfaces,i);
	fedges = MF_Edges(face,1,0);
	for (j = 0; j < 4 && n < 12; j++) {
	  edge = Set_Entry(fedges,j);
	  if (!MEnt_IsMarked(edge,mkr)) {
	    Set_Add(redges,edge);
	    MEnt_Mark(edge,mkr);
	    n++;
	  }
	}
	Set_Delete(fedges);
      }
      Set_Unmark(redges,mkr);
      MSTK_FreeMarker(mkr);

      return redges;
      break;
    default: /* Pyramids, Prisms, General Polyhedra */
      /* We should do separate procedures for pyramids and prisms */
      mkr = MSTK_GetMarker();

      /* Add edges of first face */
      face = Set_Entry(downadj->rfaces,0); /* first face */
      fdir = downadj->fdirs & 1;    /* Sense in which face is used in region */

      redges = MF_Edges(face,!fdir,0);
      Set_Mark(redges,mkr);

      for (i = 1; i < downadj->nf-1; i++) {
	face = Set_Entry(downadj->rfaces,i);
	fedges = MF_Edges(face,1,0);
	n = Set_Num_Entries(fedges);
	for (j = 0; j < n; j++) {
	  edge = Set_Entry(fedges,j);
	  if (!MEnt_IsMarked(edge,mkr)) {
	    Set_Add(redges,edge);
	    MEnt_Mark(edge,mkr);
	  }
	}
	Set_Delete(fedges);
      }
      Set_Unmark(redges,mkr);
      MSTK_FreeMarker(mkr);

      return redges;
      break;
    }

    return 0;
  }

  Set_ptr MR_Faces_FNR3R4(MRegion_ptr r) {
    MRegion_DownAdj_FN *downadj;
    downadj = (MRegion_DownAdj_FN *) r->downadj;
    return Set_Copy(downadj->rfaces);
  }

  Set_ptr MR_AdjRegions_FNR3R4(MRegion_ptr r) {
    int i;
    MRegion_ptr freg;
    MFace_ptr face;
    Set_ptr adjr;
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;

    adjr = Set_New(4);
    for (i = 0; i < downadj->nf; i++) {
      face = Set_Entry(downadj->rfaces,i);
      freg = MF_Region(face,0);
      if (freg) {
	if (freg == r) {
	  freg = MF_Region(face,1);
	  Set_Add(adjr,freg);
	}
	else 
	  Set_Add(adjr,freg);
      }
    }
    return adjr;
  }

  int MR_FaceDir_FNR3R4(MRegion_ptr r, MFace_ptr f) {
    int i;
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;
    for (i = 0; i < downadj->nf; i++)
      if (f == (MFace_ptr) Set_Entry(downadj->rfaces,i))
	return ((downadj->fdirs)>>i & 1);
    return -1;
  }

  int MR_FaceDir_i_FNR3R4(MRegion_ptr r, int i) {
    MRegion_DownAdj_FN *downadj;
    downadj = (MRegion_DownAdj_FN *) r->downadj;
    return ((downadj->fdirs)>>i & 1);
  }

  void MR_Replace_Face_FNR3R4(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf, int nudir) {
    int i;
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;
    for (i = 0; i < downadj->nf; i++)
      if (f == (MFace_ptr) Set_Entry(downadj->rfaces,i)) {
	downadj->fdirs = (downadj->fdirs & ~(1<<i)); /* set bit i to 0 */
	downadj->fdirs = (downadj->fdirs | (nudir<<i)); /* set to nudir*/
	Set_Replacei(downadj->rfaces,i,nuf);

	MF_Rem_Region(f,r);
	MF_Add_Region(nuf,r,!nudir);
	return;
      }
  }

  void MR_Replace_Face_i_FNR3R4(MRegion_ptr r, int i, MFace_ptr nuf, int nudir) {
    MRegion_DownAdj_FN *downadj;
    MFace_ptr f;

    f = Set_Entry(downadj->rfaces,i);
    downadj = (MRegion_DownAdj_FN *) r->downadj;
    downadj->fdirs = (downadj->fdirs & ~(1<<i)); /* set bit i to 0 */
    downadj->fdirs = (downadj->fdirs | (nudir<<i)); /* set to nudir*/
    Set_Replacei(downadj->rfaces,i,nuf);

    MF_Rem_Region(f,r);
    MF_Add_Region(nuf,r,!nudir);
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

  void MR_Add_AdjRegion_FNR3R4(MRegion_ptr r, int facenum, MRegion_ptr aregion) {
#ifdef DEBUG
    MSTK_Report("MR_Add_AdjRegion",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MR_Rem_AdjRegion_FNR3R4(MRegion_ptr r, MRegion_ptr aregion) {
#ifdef DEBUG
    MSTK_Report("MR_Rem_AdjRegion",
		"Function call not suitable for this representation",WARN);
#endif
  }


  int MR_UsesFace_FNR3R4(MRegion_ptr r, MFace_ptr f) {
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;
    return Set_Contains(downadj->rfaces,f);
  }

  int MR_UsesEdge_FNR3R4(MRegion_ptr r, MEdge_ptr e) {
    int i;
    MFace_ptr face;
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;
    for (i = 0; i < downadj->nf; i++) {
      face = Set_Entry(downadj->rfaces,i);
      if (MF_UsesEntity(face,e,1))
	return 1;
    }
    return 0;
  }

  int MR_UsesVertex_FNR3R4(MRegion_ptr r, MVertex_ptr v) {
    int i;
    MFace_ptr face;
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;
    for (i = 0; i < downadj->nf; i++) {
      face = Set_Entry(downadj->rfaces,i);
      if (MF_UsesEntity(face,v,0))
	return 1;
    }
    return 0;
  }

#ifdef __cplusplus
}
#endif
