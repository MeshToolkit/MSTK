#define _H_MRegion_Private

#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif


  List_ptr MR_Edges_FN(MRegion_ptr r) {
    MFace_ptr botface;      /* bottom face */
    int botfdir;            /* dir in which bottom face used by region */
    MEdge_ptr botedge;      /* edge of bottom face */
    MVertex_ptr botvert;    /* vertex of bottom face */
    MEdge_ptr latedge;      /* lateral edge */    
    MFace_ptr topface;      /* top face for prisms and hexes */
    MEdge_ptr topedge;      /* edge of top face */
    MVertex_ptr topvert;    /* vertex of top face */
    int topfdir;            /* dir in which top face used by region */
    MVertex_ptr vapex;      /* Apex vertex for tet/pyramid */

    int i, j, n, fdir, nf;
    int idx, found, fedir, dir;
    MFace_ptr face;
    MEdge_ptr edge;
    MVertex_ptr vert;
    List_ptr redges, fedges, fverts;
    MRegion_Adj_FN *adj;

    adj = (MRegion_Adj_FN *) r->adj;
    nf = List_Num_Entries(adj->rfaces);

    switch (nf) {
    case 4: { /* Tet */
      
      int alltri;
      
      /* Make sure all faces have 3 edges each */
      alltri = 1;
      idx = 0;
      while ((face = List_Next_Entry(adj->rfaces,&idx))) {
	if (MF_Num_Vertices(face) != 3) {
	  alltri = 0;
	  break;
	}
      }
      
      if (alltri) { /* Tet */
	
	botface = List_Entry(adj->rfaces,0); /* first face */
	botfdir = adj->fdirs[0] & 1; /* Sense in which face is used in region */
	
	redges = MF_Edges(botface,!botfdir,0);
	
	face = List_Entry(adj->rfaces,1);	
	fverts = MF_Vertices(face,1,0);
	
	vapex = NULL;
	for (i = 0; i < 3; i++) {
	  vert = List_Entry(fverts,i);
	  if (!MF_UsesEntity(botface,vert,MVERTEX)) {
	    vapex = vert;
	    break;
	  }
	}					
	List_Delete(fverts);
	
	if (vapex == NULL) {
	  MSTK_Report("MR_Edges_FN","Could not find fourth vertex of tet",MSTK_ERROR);
	  return redges;
	}
	
	
	for (i = 0; i < 3; i++) {
	  botedge = List_Entry(redges,i);
	  fedir = MF_EdgeDir(botface,botedge);
	  dir = !(botfdir^fedir);
	  botvert = ME_Vertex(botedge,dir);
	  
	  latedge = MVs_CommonEdge(botvert,vapex);
	  if (latedge)
	    List_Add(redges,latedge);
	  else {
	    MSTK_Report("MR_Edges_FN","Could not find edge between bottom face and apex vertex of tet",MSTK_ERROR);
	    continue;
	  }
	}

	return redges;
      }
      break;
    }
    case 5:  {
      int ne;
      int nquads = 0, ntris = 0;
      int itri0 = -1, itri1 = -1, iquad0 = -1;
      MFace_ptr tri0 = NULL, tri1 = NULL, quad0 = NULL;

      for (i = 0; i < nf; i++) {
	face = List_Entry(adj->rfaces,i);
	ne = MF_Num_Edges(face);
	if (ne == 3) {
	  ntris++;
	  if (!tri0) {
	    tri0 = face;
	    itri0 = i;
	  }
	  else {
	    tri1 = face;
	    itri1 = i;
	  }
	}
	else if (ne == 4) {
	  nquads++;
	  if (!quad0) {
	    quad0 = face;
	    iquad0 = i;
	  }
	}
      }
      
      if (nquads == 3 && ntris == 2) {

	botfdir = MR_FaceDir_i(r,itri0); /* Sense in which face is used in region */
	redges = MF_Edges(tri0,!botfdir,0);

	/* find the vertical edge between vertex 0 of bottom triangle and
	   top triangle */

	botedge = List_Entry(redges,0);
	fedir = MF_EdgeDir(tri0,botedge);
	dir = !(botfdir^fedir);
	botvert = ME_Vertex(botedge,dir); 

	topfdir = MR_FaceDir_i(r,itri1);
	fverts = MF_Vertices(tri1,topfdir,0);

	found = 0;
	for (i = 0; i < 3; i++) {
	  topvert = List_Entry(fverts,i);

	  if ((latedge = MVs_CommonEdge(botvert,topvert))) {
	    found = 1;
	    break;
	  }
	}

	if (!found) {
	  MSTK_Report("MR_Edges_FN","Could not find edge connecting bottom and top faces of prism",MSTK_ERROR);
	  return redges;
	}

	fedges = MF_Edges(tri1,topfdir,topvert);
	List_Cat(redges,fedges);
	List_Delete(fedges);


	List_Add(redges,latedge);

	botedge = List_Entry(redges,1);
	fedir = MF_EdgeDir(tri0,botedge);
	dir = !(botfdir^fedir);
	botvert = ME_Vertex(botedge,dir);
	topvert = List_Entry(fverts,(i+1)%3);
	latedge = MVs_CommonEdge(botvert,topvert);
	if (latedge)
	  List_Add(redges,latedge);
	else {
	  MSTK_Report("MR_Edges_FN","Could not find edge connecting bottom and top faces of prism",MSTK_ERROR);
	  return redges;
	}

	botedge = List_Entry(redges,2);
	fedir = MF_EdgeDir(tri0,botedge);
	dir = !(botfdir^fedir);
	botvert = ME_Vertex(botedge,dir);
	topvert = List_Entry(fverts,(i+2)%3);
	latedge = MVs_CommonEdge(botvert,topvert);
	if (latedge)
	  List_Add(redges,latedge);
	else {
	  MSTK_Report("MR_Edges_FN","Could not find edge connecting bottom and top faces of prism",MSTK_ERROR);
	  return redges;
	}

	List_Delete(fverts);
	return redges;
      }
      else if (nquads == 1 && ntris == 4) {

	botfdir = MR_FaceDir_i(r,iquad0);  /* Sense in which face is used in region */
	redges = MF_Edges(quad0,!botfdir,0);

        vapex = NULL;
	fverts = MF_Vertices(tri0,1,0);
	for (i = 0; i < 3; i++) {
	  vert = List_Entry(fverts,i);
	  if (!MF_UsesEntity(quad0,vert,MVERTEX)) { /* found apex vtx of pyramid */
	    vapex = vert;
	    break;
	  }
	}
	List_Delete(fverts);


	for (i = 0; i < 4; i++) {	  
	  botedge = List_Entry(redges,i);
	  fedir = MF_EdgeDir(quad0,botedge);	  
	  dir = !(botfdir^fedir);
	  botvert = ME_Vertex(botedge,dir);
	  latedge = MVs_CommonEdge(botvert,vapex);
	  if (latedge)
	    List_Add(redges,latedge);
	  else {
	    MSTK_Report("MR_Edges_FN","Could not find edge between bottom face and apex vertex for pyramid",MSTK_ERROR);
	    return redges;
	  }
	}
	  
	return redges;
      }
      break;
    }
    case 6: { /* Hex ? */

      int allquad;

      /* All faces must have 4 edges each */
      allquad = 1;
      for (i = 0; i < nf; i++) {
	face = List_Entry(adj->rfaces,i);
	if (MF_Num_Edges(face) != 4) {
	  allquad = 0;
	  break;
	}
      }

      if (allquad) { /* Hex */

	n = 0;
	
	/* Add edges of first face */
	botface = List_Entry(adj->rfaces,0); /* first face */
	botfdir = adj->fdirs[0] & 1;    /* Sense in which face is used in region */
	
	redges = MF_Edges(botface,!botfdir,0);

	/* Find the opposite face */

	topface = NULL;
	for (i = 1; i < 6; i++) {
	  face = List_Entry(adj->rfaces,i);
	  fedges = MF_Edges(face,1,0);
	  int none_used = 1;
	  for (j = 0; j < 4; j++) {
	    MEdge_ptr edge = List_Entry(fedges,j);
	    if (MF_UsesEntity(botface,edge,MEDGE)) {
	      none_used = 0;
	      break;
	    }
	  }
	  List_Delete(fedges);
	  if (none_used) {
	    topface = face;
	    break;
	  }
	}

	if (!topface) {
	  MSTK_Report("MR_Edges_FN","Cannot find top face for hex",MSTK_ERROR);
	  return redges;
	}
	

	botedge = List_Entry(redges,0);
	fedir = MF_EdgeDir(botface,botedge);
	dir = !(botfdir^fedir);
	botvert = ME_Vertex(botedge,dir); 

	topfdir = MR_FaceDir_i(r,i);
	fverts = MF_Vertices(topface,topfdir,0);

	found = 0;
	for (i = 0; i < 4; i++) {
	  topvert = List_Entry(fverts,i);

	  if ((latedge = MVs_CommonEdge(botvert,topvert))) {
	    found = 1;
	    break;
	  }
	}

	if (!found) {
	  MSTK_Report("MR_Edges_FN","Could not find edge connecting bottom and top faces of hex",MSTK_ERROR);
	  return redges;
	}

	fedges = MF_Edges(topface,topfdir,topvert);
	List_Cat(redges,fedges);
	List_Delete(fedges);


	List_Add(redges,latedge);

	botedge = List_Entry(redges,1);
	fedir = MF_EdgeDir(botface,botedge);
	dir = !(botfdir^fedir);
	botvert = ME_Vertex(botedge,dir);
	topvert = List_Entry(fverts,(i+1)%4);
	latedge = MVs_CommonEdge(botvert,topvert);
	if (latedge)
	  List_Add(redges,latedge);
	else {
	  MSTK_Report("MR_Edges_FN","Could not find edge connecting bottom and top faces of hex",MSTK_ERROR);
	  return redges;
	}

	botedge = List_Entry(redges,2);
	fedir = MF_EdgeDir(botface,botedge);
	dir = !(botfdir^fedir);
	botvert = ME_Vertex(botedge,dir);
	topvert = List_Entry(fverts,(i+2)%4);
	latedge = MVs_CommonEdge(botvert,topvert);
	if (latedge)
	  List_Add(redges,latedge);
	else {
	  MSTK_Report("MR_Edges_FN","Could not find edge connecting bottom and top faces of prism",MSTK_ERROR);
	  return redges;
	}


	botedge = List_Entry(redges,3);
	fedir = MF_EdgeDir(botface,botedge);
	dir = !(botfdir^fedir);
	botvert = ME_Vertex(botedge,dir);
	topvert = List_Entry(fverts,(i+3)%4);
	latedge = MVs_CommonEdge(botvert,topvert);
	if (latedge)
	  List_Add(redges,latedge);
	else {
	  MSTK_Report("MR_Edges_FN","Could not find edge connecting bottom and top faces of prism",MSTK_ERROR);
	  return redges;
	}

	
	List_Delete(fverts);

	return redges;
      }
      break;
    }
    default:
      break;
    }

    /* General Polyhedra */
    
    /* Add edges of first face */
    face = List_Entry(adj->rfaces,0); /* first face */
    fdir = adj->fdirs[0] & 1;   /* Sense in which face is used in region */
    
    redges = MF_Edges(face,!fdir,0);
    
    for (i = 1; i < nf-1; i++) {
      face = List_Entry(adj->rfaces,i);
      fedges = MF_Edges(face,1,0);
      n = List_Num_Entries(fedges);
      for (j = 0; j < n; j++) {
	edge = List_Entry(fedges,j);
        int inlist;
        inlist = List_Contains(redges,edge);
        if (!inlist) {
	  List_Add(redges,edge);
	}
      }
      List_Delete(fedges);
    }
    
    return redges;
  }


  void MR_EdgeIDs_FN(MRegion_ptr r, int *nre, int *redgeids) {
    int i;
    List_ptr redges = MR_Edges(r);
    *nre = List_Num_Entries(redges);
    for (i = 0; i < *nre; i++)
      redgeids[i] = MEnt_ID(List_Entry(redges,i));
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
