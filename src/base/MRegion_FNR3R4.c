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
    downadj->fdirs = NULL;
    downadj->rfaces = NULL;
  }

  void MR_Delete_F1F3R3R4(MRegion_ptr r, int keep) {
    MRegion_DownAdj_FN *downadj;
    MFace_ptr f;
    int i, nf;

      downadj = (MRegion_DownAdj_FN *) r->downadj;
      
    if (r->dim != MDELREGION) { /* if region has not been temporarily deleted */
      if (downadj) {
	nf = List_Num_Entries(downadj->rfaces);
	for (i = 0; i < nf; i++) {
	  f = List_Entry(downadj->rfaces,i);
	  MF_Rem_Region(f,r);
	}
      }
    }

    if (keep) {
      MSTK_KEEP_DELETED = 1;
      r->dim = MDELREGION;
    }
    else {
#ifdef DEBUG
      r->dim = MDELREGION;
#endif

      if (downadj) {	
	if (downadj->rfaces) {
	  MSTK_free(downadj->fdirs);
	  List_Delete(downadj->rfaces);
	}
	MSTK_free(downadj);
      }

      MSTK_free(r);
    }
  }

  void MR_Restore_F1F3R3R4(MRegion_ptr r) {
    MRegion_DownAdj_FN *downadj;
    MFace_ptr f;
    int i, j, k, nf, side;

    if (r->dim != MDELREGION)
      return;

    r->dim = MREGION;

    downadj = (MRegion_DownAdj_FN *) r->downadj;

    nf = List_Num_Entries(downadj->rfaces);
    for (i = 0; i < nf; i++) {
      f = List_Entry(downadj->rfaces,i);
 
      j = (int) i/(8*sizeof(unsigned int));
      k = i%(8*sizeof(unsigned int));
      side = !((downadj->fdirs[j])>>k & 1);

      MF_Add_Region(f,r,side);
    }
  }

  void MR_Destroy_For_MESH_Delete_FNR3R4(MRegion_ptr r) {
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;
    if (downadj) {
      if (downadj->rfaces) {
	MSTK_free(downadj->fdirs);
	List_Delete(downadj->rfaces);
      }
      MSTK_free(downadj);
    }

    MSTK_free(r);
  }

  void MR_Set_Faces_FNR3R4(MRegion_ptr r, int nf, MFace_ptr *rfaces,int *dirs){
    int i, j, k;
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;
    downadj->nf = nf;
    k = 1 + ((int) nf/(8*sizeof(unsigned int)));
    downadj->fdirs = (unsigned int *) MSTK_calloc(k,sizeof(unsigned int));
    downadj->rfaces = List_New(nf);

    for (i = 0; i < nf; i++) {
      j = (int) i/(8*sizeof(unsigned int));
      k = i%(8*sizeof(unsigned int));
      downadj->fdirs[j] = downadj->fdirs[j] | (dirs[i] << k);

      List_Add(downadj->rfaces,rfaces[i]);

      /* could change to MESH_RepType(r->mesh) != F4 */
      if (r->repType != F4) 
	MF_Add_Region(rfaces[i],r,!dirs[i]);
    }
  }

  void MR_Set_Vertices_FNR3R4(MRegion_ptr r, int nv, MVertex_ptr *mvertices, int nf, int **rfvtemplate) {
    MVertex_ptr verts[MAXPV2], everts[2];
    MEdge_ptr redges[12], fedges[4]; /* will be used for std. elements only */
    MFace_ptr rfaces[MAXPF3];
    MRegion_ptr oregion;
    int rfdirs[MAXPF3], odir, i, j, fgdim, fgid, nfv, ind, nre;
    int vgdim[MAXPV3], vgid[MAXPV3], egdim[12], egid[12], sgn;
    int fedirs[4], redirs[12], evgdim[2], evgid[2], nfe;
    List_ptr fverts, fregs;
    Mesh_ptr mesh = r->mesh;
    MRType regtype;


    /* Templates for TETS, PYRAMIDS, PRISMS and HEXS (move to global file?) */

    /* Face-vertex templates for standard regions */
    static int rfvtmpl[4][6][5] = 
      {{{3,0,1,2,-1},{3,0,1,3,-1},{3,1,2,3,-1},{3,0,3,2,-1},{0,-1,-1,-1,-1},
	{0,-1,-1,-1,-1}},
       {{4,0,1,2,3},{3,0,1,4,-1},{3,1,2,4,-1},{3,2,3,4,-1},{3,3,0,4,-1},
	{0,-1,-1,-1,-1}},
       {{3,0,1,2,-1},{3,3,4,5,-1},{4,0,1,4,3},{4,1,2,5,4},{4,2,0,3,5},
	{0,-1,-1,-1,-1}},
       {{4,0,1,2,3},{4,4,5,6,7},{4,0,1,5,4},{4,1,2,6,5},{4,2,3,7,6},
	{4,3,0,4,7}}
      };

    /* Face direction templates for standard regions */
    static int rfdirtmpl[4][6] = 
      {{0,1,1,1,-1,-1},{0,1,1,1,1,-1},{0,1,1,1,1,-1},{0,1,1,1,1,1}};

    /* Face-edge templates for standard regions */
    static int rfetmpl[4][6][5] = 
      {{{3,0,1,2,-99},{3,0,4,-3,-99},{3,1,5,-4,-99},{3,2,3,-5,-99},
	{0,-99,-99,-99,-99},{0,-99,-99,-99,-99}},
       {{4,0,1,2,3},{3,0,5,-4,-99},{3,1,6,-5,-99},{3,2,7,-6,-99},
	{3,3,4,-7,-99},{0,-99,-99,-99,-99}},
       {{3,0,1,2,-99},{3,3,4,5,-99},{4,0,7,-3,-6},{4,1,8,-4,-7},{4,2,6,-5,-8},
	{0,-99,-99,-99,-99}},
       {{4,0,1,2,3},{4,4,5,6,7},{4,0,9,-4,-8},{4,1,10,-5,-9},{4,2,11,-6,-10},
	{4,3,8,-7,-11}}};

    /* Edge-vertex templates for standard regions */
    static int revtmpl[4][12][2] = 
      {{{0,1},{1,2},{2,0},{0,3},{1,3},{2,3},{-99,-99},{-99,-99},{-99,-99},
	{-99,-99},{-99,-99},{-99,-99}},
       {{0,1},{1,2},{2,3},{3,0},{0,4},{1,4},{2,4},{3,4},{-99,-99},{-99,-99},
	{-99,-99},{-99,-99}},
       {{0,1},{1,2},{2,0},{3,4},{4,5},{5,3},{0,3},{1,4},{2,5},{-99,-99},
	{-99,-99},{-99,-99}},
       {{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},
	{3,7}}};

    /* Number of edges for standard regions */
    static int nretmpl[4] = {6,8,9,12};

    for (i = 0; i < nv; i++) {
      vgdim[i] = MV_GEntDim(mvertices[i]);
      vgid[i] = MV_GEntID(mvertices[i]);
    }

    if (rfvtemplate) {
      for (i = 0; i < nf; i++) {
	fgdim = 4;
	fgid = 0;

	nfv = rfvtemplate[i][0];
	for (j = 0; j < nfv; j++) {
	  ind = rfvtemplate[i][j+1];
	  verts[j] = mvertices[ind];
	  if (vgdim[ind] > fgdim) {
	    fgdim = vgdim[ind];
	    fgid = vgid[ind];
	  }
	}

       
	rfaces[i] = MVs_CommonFace(nfv, verts);

	if (!rfaces[i]) {
	  rfaces[i] = MF_New(mesh);

	  if (fgdim >= 2) {
	    MF_Set_GEntDim(rfaces[i],fgdim);
	    MF_Set_GEntID(rfaces[i],fgid);
	  }
	  else {
	    /* Cannot determine unique classification for face */
	    MSTK_Report("MR_Set_Vertices",
			"Cannot determine unique classification for face",
			WARN);
	    MF_Set_GEntDim(rfaces[i],4);
	    MF_Set_GEntID(rfaces[i],0);
	  }


	  MF_Set_Vertices(rfaces[i], nfv, verts);
	  rfdirs[i] = 1;
	}
	else {
	  fregs = MF_Regions(rfaces[i]);
	  if (!fregs) {
	    fverts = MF_Vertices(rfaces[i],1,verts[0]);
	    rfdirs[i] = (verts[1] == List_Entry(fverts,1)) ? 1 : 0;
	    List_Delete(fverts);
	  }
	  else {
#ifdef DEBUG
	    if (List_Num_Entries(fregs) > 1)
	      MSTK_Report("MR_Set_Vertices",
			  "Face already connected on both sides",FATAL);
#endif
	    oregion = List_Entry(fregs,0);
	    odir = MR_FaceDir(oregion,rfaces[i]);
	    rfdirs[i] = !odir;

	    List_Delete(fregs);
	  }
	}

      }

      MR_Set_Faces(r, nf, rfaces, rfdirs);
    }
    else {
      switch (nv) {
      case 4:
	regtype = TET;
	nf = 4;
	break;
      case 5:
	regtype = PYRAMID;
	nf = 5;
	break;
      case 6:
	regtype = PRISM;
	nf = 5;
	break;
      case 8:
	regtype = HEX;
	nf = 6;
	break;
      default:
	MSTK_Report("MR_Set_Vertices",
		    "Polyhedron: Need number of faces and vertex template",
		    FATAL);
	regtype = RUNKNOWN;
	break;
      }

      /* First collect the edges */
      nre = nretmpl[regtype-1];
      for (i = 0; i < nre; i++) {

	/* Check if edge exists */

	ind = revtmpl[regtype-1][i][0];
	everts[0] = mvertices[ind];
	evgdim[0] = MV_GEntDim(everts[0]);
	evgid[0]  = MV_GEntID(everts[0]);
	
	ind = revtmpl[regtype-1][i][1];
	everts[1] = mvertices[ind];
	evgdim[1] = MV_GEntDim(everts[1]);
	evgid[1]  = MV_GEntID(everts[1]);

	redges[i] = MVs_CommonEdge(everts[0],everts[1]);

	/* Create the edge if it does not exist */

	if (!redges[i]) {
	  redges[i] = ME_New(mesh);

	  ME_Set_Vertex(redges[i],0,everts[0]);
	  ME_Set_Vertex(redges[i],1,everts[1]);
	  redirs[i] = 1;

	  egdim[i] = 4;
	  egid[i] = 0;

	  if (evgdim[0] > evgdim[1]) {
	    egdim[i] = evgdim[0];
	    egid[i] = evgid[0];
	  }
	  else if (evgdim[1] > evgdim[0]) {
	    egdim[i] = evgdim[1];
	    egid[i] = evgid[1];
	  }
	  else { /* evgdim[0] == evgdim[1] */
	    if (evgdim[0] == 0) {
	      MSTK_Report("MR_Set_Vertices_FNR3R4",
			  "Cannot determine edge classification. Guessing...",
			  WARN);
	      egdim[i] = 1;
	      egid[i] = 0;
	    }
	    else {
	      if (evgid[0] == evgid[1]) {
		/* Assume that edge is classified on same entity as
		   its two vertices. Could be wrong but will be right
		   more often than not */
		egdim[i] = evgdim[0];
		egid[i] = evgid[0];
	      }
	      else {
		/* Vertices are classified on two different edges,
		   faces or regions */
		MSTK_Report("MR_Set_Vertices_FNR3R4",
			    "Cannot determine edge classification. Guessing..",
			    WARN);
		egdim[i] = evgdim[0]+1;
		egid[i] = 0;
	      }
	    }
	  }

	  ME_Set_GEntDim(redges[i],egdim[i]);
	  ME_Set_GEntID(redges[i],egid[i]);
	}
	else {
	  redirs[i] = (ME_Vertex(redges[i],0) == everts[0]) ? 1 : 0;
	  egdim[i] = ME_GEntDim(redges[i]);
	  egid[i] = ME_GEntID(redges[i]);
	}
      }


      for (i = 0; i < nf; i++) {
	fgdim = 4;
	fgid = 0;

	nfe = rfetmpl[regtype-1][i][0];
	for (j = 0; j < nfe; j++) {
	  ind = rfetmpl[regtype-1][i][j+1];
	  sgn = (ind >= 0) ? 1 : 0;
	  ind = abs(ind);
	  fedges[j] = redges[ind];
	  fedirs[j] = !(redirs[ind]^sgn);

	  if (egdim[ind] > fgdim) {
	    fgdim = egdim[ind];
	    fgid = egid[ind];
	  }
	}

	rfaces[i] = MEs_CommonFace(nfe, fedges);

	if (!rfaces[i]) {
	  rfaces[i] = MF_New(mesh);
	  rfdirs[i] = rfdirtmpl[regtype-1][i];

	  MF_Set_Edges(rfaces[i], nfe, fedges, fedirs);

	  if (fgdim >= 2) {
	    MF_Set_GEntDim(rfaces[i],fgdim);
	    MF_Set_GEntID(rfaces[i],fgid);
	  }
	  else {
	    /* Cannot determine unique classification for face */
	    MSTK_Report("MR_Set_Vertices",
			"Cannot determine face classification. Guessing...",
			WARN);
	    MF_Set_GEntDim(rfaces[i],(fgdim+1));
	    MF_Set_GEntID(rfaces[i],0);
	  }


	}
	else {
	  fregs = MF_Regions(rfaces[i]);
	  if (!fregs) {
	    ind = rfvtmpl[regtype-1][i][1];
	    verts[0] = mvertices[ind];
	    ind = rfvtmpl[regtype-1][i][2];
	    verts[1] = mvertices[ind];

	    fverts = MF_Vertices(rfaces[i],1,verts[0]);
	    rfdirs[i] = ((verts[1] == List_Entry(fverts,1)) ? 
			 rfdirtmpl[regtype-1][i] : !rfdirtmpl[regtype-1][i]);
	    List_Delete(fverts);
	  }
	  else {
#ifdef DEBUG
	    if (List_Num_Entries(fregs) > 1)
	      MSTK_Report("MR_Set_Vertices",
			  "Face already connected on both sides",FATAL);
#endif
	    oregion = List_Entry(fregs,0);
	    odir = MR_FaceDir(oregion,rfaces[i]);
	    rfdirs[i] = !odir;

	    List_Delete(fregs);
	  }
	}

      }

      MR_Set_Faces(r, nf, rfaces, rfdirs);
    }
  }

  int MR_Num_Faces_FNR3R4(MRegion_ptr r) {
    return ((MRegion_DownAdj_FN *)r->downadj)->nf;
  }

  int MR_Num_AdjRegions_FNR3R4(MRegion_ptr r) {
    List_ptr adjr;
    int nr;

#ifdef DEBUG
    MSTK_Report("MR_Num_AdjRegions",
		"Inefficient to call this routine with this representation",
		WARN);
#endif

    adjr = MR_AdjRegions(r);
    if (adjr) {
      nr = List_Num_Entries(adjr);
      List_Delete(adjr);
      return nr;
    }
    else
      return 0;
  }

  List_ptr MR_Vertices_FNR3R4(MRegion_ptr r) {
    int i, j, n, mkr, found, diradj0=0, diropp=0, edir, fdir, fdir0, fecheck;
    MFace_ptr face=NULL, face0=NULL, fadj0=NULL, fopp=NULL;
    MEdge_ptr edge, fedge00, upedge;
    MVertex_ptr vert, rv0, rvopp0=NULL;
    List_ptr rvertices, fverts, fedges0, adjfedges;
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;


    switch (downadj->nf) {
    case 4: /* Tet! */
      /* Add vertices of first face to list of region vertices */
      face0 = List_Entry(downadj->rfaces,0); /* first face */
      fdir0 = downadj->fdirs[0] & 1;    /* Sense in which face is used in region */
      rvertices = MF_Vertices(face0,!fdir0,0);
      n = 3;

      face = List_Entry(downadj->rfaces,1);
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
      for (i = 0; i < downadj->nf; i++) {
	face = List_Entry(downadj->rfaces,i);
	if (MF_Num_Edges(face) != 4)
	  fecheck = 0;
      }

      if (!fecheck)
	break;

      face0 = List_Entry(downadj->rfaces,0); /* face 0 */
      fdir0 = downadj->fdirs[0] & 1; /* dir of face w.r.t. region */

      
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
      for (i = 1; i < downadj->nf; i++) {
	face = List_Entry(downadj->rfaces,i);

	/* Check if face uses any edge of face 0 */

	for (j = 0, found = 0; j < 4; j++) {
	  edge = List_Entry(fedges0,j);

	  if (MF_UsesEntity(face,edge,1)) {
	    if (edge == fedge00) {
	      /* face uses edge 0 of face 0 (w.r.t. face dir pointing
		 into region) */
	      fadj0 = face;
	      diradj0 = (downadj->fdirs[0])>>i & 1;  
	    }
	    found = 1;
	    break;
	  }
	}
	if (!found) {
	  fopp = face;
	  diropp = (downadj->fdirs[0])>>i & 1;
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
    face = List_Entry(downadj->rfaces,0); /* first face */
    fdir = downadj->fdirs[0] & 1;    /* Sense in which face is used in region */
    
    rvertices = MF_Vertices(face,!fdir,0); 
    List_Mark(rvertices,mkr);
    
    for (i = 1; i < downadj->nf-1; i++) {
      face = List_Entry(downadj->rfaces,i);
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

  List_ptr MR_Edges_FNR3R4(MRegion_ptr r) {
    int i, j, n, mkr, fdir, fecheck;
    MFace_ptr face;
    MEdge_ptr edge;
    List_ptr redges, fedges;
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;

    switch (downadj->nf) {
    case 4: /* Tet */
      mkr = MSTK_GetMarker();

      face = List_Entry(downadj->rfaces,0); /* first face */
      fdir = downadj->fdirs[0] & 1;    /* Sense in which face is used in region */

      redges = MF_Edges(face,!fdir,0);
      List_Mark(redges,mkr);
      n = 3;

      face = List_Entry(downadj->rfaces,1);
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

      face = List_Entry(downadj->rfaces,2);
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
      for (i = 0; i < downadj->nf; i++) {
	face = List_Entry(downadj->rfaces,i);
	if (MF_Num_Edges(face) != 4)
	  fecheck = 0;
      }

      if (!fecheck)
	break;

      n = 0;
      mkr = MSTK_GetMarker();

      /* Add edges of first face */
      face = List_Entry(downadj->rfaces,0); /* first face */
      fdir = downadj->fdirs[0] & 1;    /* Sense in which face is used in region */

      redges = MF_Edges(face,!fdir,0);
      List_Mark(redges,mkr);
      n = 4;

      for (i = 1; i < (downadj->nf)-1 && n < 12; i++) {
	face = List_Entry(downadj->rfaces,i);
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
    face = List_Entry(downadj->rfaces,0); /* first face */
    fdir = downadj->fdirs[0] & 1;   /* Sense in which face is used in region */
    
    redges = MF_Edges(face,!fdir,0);
    List_Mark(redges,mkr);
    
    for (i = 1; i < downadj->nf-1; i++) {
      face = List_Entry(downadj->rfaces,i);
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

  List_ptr MR_Faces_FNR3R4(MRegion_ptr r) {
    MRegion_DownAdj_FN *downadj;
    downadj = (MRegion_DownAdj_FN *) r->downadj;
    return List_Copy(downadj->rfaces);
  }

  List_ptr MR_AdjRegions_FNR3R4(MRegion_ptr r) {
    int i;
    MRegion_ptr freg;
    MFace_ptr face;
    List_ptr adjr;
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;

    adjr = List_New(4);
    for (i = 0; i < downadj->nf; i++) {
      face = List_Entry(downadj->rfaces,i);
      freg = MF_Region(face,0);
      if (freg) {
	if (freg == r) {
	  freg = MF_Region(face,1);
	  List_Add(adjr,freg);
	}
	else 
	  List_Add(adjr,freg);
      }
    }
    return adjr;
  }

  int MR_FaceDir_FNR3R4(MRegion_ptr r, MFace_ptr f) {
    int i,j,k;
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;

    for (i = 0; i < downadj->nf; i++) {

      if (f == (MFace_ptr) List_Entry(downadj->rfaces,i)) {

	j = (int) i/(8*sizeof(unsigned int));
	k = i%(8*sizeof(unsigned int));

	return ((downadj->fdirs[j])>>k & 1);
      }

    }

    return -1;
  }

  int MR_FaceDir_i_FNR3R4(MRegion_ptr r, int i) {
    int j, k;
    MRegion_DownAdj_FN *downadj;
    downadj = (MRegion_DownAdj_FN *) r->downadj;
    
    j = (int) i/(8*sizeof(unsigned int));
    k = i%(8*sizeof(unsigned int));

    return ((downadj->fdirs[j])>>k & 1);
  }

  /* If we allow replacing a face with multiple faces, we have to check if 
     we need to expand downadj->fdirs */

  void MR_Replace_Face_FNR3R4(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf, int nudir) {
    int i, j, k;
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;
    for (i = 0; i < downadj->nf; i++)
      if (f == (MFace_ptr) List_Entry(downadj->rfaces,i)) {

	j = (int) i/(8*sizeof(unsigned int));
	k = i%(8*sizeof(unsigned int));

	downadj->fdirs[j] = (downadj->fdirs[j] & ~(1<<k)); /* set bit i to 0 */
	downadj->fdirs[j] = (downadj->fdirs[j] | (nudir<<k)); /* set to nudir*/

	List_Replacei(downadj->rfaces,i,nuf);

	MF_Rem_Region(f,r);
	MF_Add_Region(nuf,r,!nudir);
	return;
      }
  }

  void MR_Replace_Face_i_FNR3R4(MRegion_ptr r, int i, MFace_ptr nuf, int nudir) {
    int j, k;
    MRegion_DownAdj_FN *downadj;
    MFace_ptr f;

    downadj = (MRegion_DownAdj_FN *) r->downadj;
    f = List_Entry(downadj->rfaces,i);

    j = (int) i/(8*sizeof(unsigned int));
    k = i%(8*sizeof(unsigned int));
    downadj->fdirs[j] = (downadj->fdirs[j] & ~(1<<k)); /* set bit i to 0 */
    downadj->fdirs[j] = (downadj->fdirs[j] | (nudir<<k)); /* set to nudir*/

    List_Replacei(downadj->rfaces,i,nuf);

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
    return List_Contains(downadj->rfaces,f);
  }

  int MR_UsesEdge_FNR3R4(MRegion_ptr r, MEdge_ptr e) {
    int i;
    MFace_ptr face;
    MRegion_DownAdj_FN *downadj;

    downadj = (MRegion_DownAdj_FN *) r->downadj;
    for (i = 0; i < downadj->nf; i++) {
      face = List_Entry(downadj->rfaces,i);
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
      face = List_Entry(downadj->rfaces,i);
      if (MF_UsesEntity(face,v,0))
	return 1;
    }
    return 0;
  }

#ifdef __cplusplus
}
#endif
