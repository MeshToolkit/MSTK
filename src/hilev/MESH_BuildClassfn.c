#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Routine to build classification for mesh entities given a mesh no
   information or model ID info only for the highest level entities
   (typically regions in volume meshes and faces in surface meshes) */


int MESH_BuildClassfn(Mesh_ptr mesh) {
  int i, k, idx, gid, gdim;
  int grid0, grid1, gfid0, gfid1, geid0, geid1;
  int max_greg_id, max_gface_id, max_gedge_id, max_gvertex_id, zeroid;
  int nfr, nef, nbf, nve, nbe, nbe2, gid2, gdim2;
  double PI=3.141592, ang, angr;
  MVertex_ptr vertex;
  MEdge_ptr edge;
  MFace_ptr face;
  MRegion_ptr region, freg0, freg1;
  List_ptr vedges, efaces, GFfaces, GEedges, GFedges, fregs;

  ang = 5*PI/6.0;  /* 75 degrees */

  /* Verify that mesh regions have classification information; if
     not, assign all regions to the same model region */

  zeroid = 0;
  max_greg_id = 0;
  idx = 0;
  while ((region = MESH_Next_Region(mesh,&idx))) {
    gid = MR_GEntID(region);
    if (gid) {
      if (gid > max_greg_id)
	max_greg_id = gid;
    }
    else
      zeroid = 1;
  }

  if (zeroid) {
    idx = 0; 
    while ((region = MESH_Next_Region(mesh,&idx))) {
      /* Make sure region is not classified on zero or -ve ID model region */

      if (MR_GEntID(region) <= 0)
	MR_Set_GEntID(region,(max_greg_id+1));
    }
  }

    
  /* Verify that mesh faces on the boundary have classification
     information; if not, assign all faces to the same model faces */
  
  max_gface_id = 0;
  idx = 0;
  while ((face = MESH_Next_Face(mesh,&idx))) {
    gdim = MF_GEntDim(face);
    if (gdim != 2)
      continue;

    gid = MF_GEntID(face);
    if (gid) {
      if (gid > max_gface_id)
	max_gface_id = gid;
    }
  }

  idx = 0;
  while ((face = MESH_Next_Face(mesh,&idx))) {
    
    gdim = MF_GEntDim(face);
    gid = MF_GEntID(face);

    /* Face has no classification? Assign classification info */
    /* Face classified as interior face? Verify */
    
    freg0 = freg1 = NULL;
    grid0 = grid1 = 0;
    fregs = MF_Regions(face);
    if (fregs) {
      nfr = List_Num_Entries(fregs);
      
      freg0 = List_Entry(fregs,0);
      grid0 = MR_GEntID(freg0);
      
      if (nfr == 2) {
	freg1 = List_Entry(fregs,1); 
	grid1 = MR_GEntID(freg1);
      }
      
      List_Delete(fregs);
    }
    else
      nfr = 0;
    
    if (nfr == 2 && (grid0 == grid1)) {
      /* Interior face */
      
      MF_Set_GEntDim(face,3);
      MF_Set_GEntID(face,grid0);
    }
    else {
      /* Boundary face */
      
      if (gdim > 2 || (gdim == 2 && gid <= 0)) { 
	MF_Set_GEntDim(face,2);
	MF_Set_GEntID(face,(max_gface_id+1)); /* for now assign 1 ID */
      }
    }
  }


  max_gedge_id = 0;
  idx = 0;
  while ((edge = MESH_Next_Edge(mesh,&idx))) {
    gid = ME_GEntID(edge);
    gdim = ME_GEntDim(edge);
    if (gdim == 1 && gid > max_gedge_id)
      max_gedge_id = gid;
  }

  idx = 0;
  while ((edge = MESH_Next_Edge(mesh,&idx))) {
    gdim = ME_GEntDim(edge);
    gid = ME_GEntID(edge);

    efaces = ME_Faces(edge);
    if (!efaces) {
      ME_Set_GEntDim(edge,1);
      ME_Set_GEntID(edge,(max_gedge_id+1));
      continue;
    }
    
    nef = List_Num_Entries(efaces);
    if (nef > 2) {
      GFfaces = List_New(nef);
      for (i = 0; i < nef; i++) {
	face = List_Entry(efaces,i);
	if (MF_GEntDim(face) == 2)
	  List_Add(GFfaces,face);
      }
    }
    else 
      GFfaces = List_Copy(efaces);
    
    nbf = List_Num_Entries(GFfaces);
    
    switch (nbf) {
    case 0: /* Interior edge */
      ME_Set_GEntDim(edge,3);
      gid2 = MF_GEntID(List_Entry(efaces,0));
      ME_Set_GEntID(edge,gid2);
      break;
    case 1: /* Edge of a non-manifold face */
      if (gdim > 1 || (gdim == 1 && gid <= 0)) {
	ME_Set_GEntDim(edge,1);
	ME_Set_GEntID(edge,(max_gedge_id+1));
      }
      break;
    case 2: 
      if (gdim > 2 || (gdim == 2 && gid <= 0)) {
	gfid0 = MF_GEntID(List_Entry(GFfaces,0));
	gfid1 = MF_GEntID(List_Entry(GFfaces,1));
	
	if (gfid0 == gfid1) { /* Edge on a model face */
	  ME_Set_GEntDim(edge,2);
	  ME_Set_GEntID(edge,gfid0);
	}
	else { /* Edge on a model edge b/w two model faces */
	  ME_Set_GEntDim(edge,1);
	  ME_Set_GEntID(edge,(max_gedge_id+1));
	}
	break;
      }
    default: /* Edge at the junction of several model faces */
      if (gdim > 1 || (gdim == 1 && gid <= 0)) {	 
	ME_Set_GEntDim(edge,1);
	ME_Set_GEntID(edge,(max_gedge_id+1));
      }
    }
    List_Delete(efaces);
    List_Delete(GFfaces);
  }

  /* Reclassify vertices in a surface mesh as being on an edge if they
     do not have a full cycle of faces around them */

  max_gvertex_id = 0;
  idx = 0;
  while ((vertex = MESH_Next_Vertex(mesh,&idx))) {
    gid = MV_GEntID(vertex);
    gdim = MV_GEntDim(vertex);
    if (gdim == 0 && gid > max_gvertex_id)
      max_gvertex_id = gid;
  }

  idx = 0;
  while ((vertex = MESH_Next_Vertex(mesh,&idx))) {

    gdim = MV_GEntDim(vertex);
    gid  = MV_GEntID(vertex);

    vedges = MV_Edges(vertex);
    if (!vedges) {
      MV_Set_GEntDim(vertex,0);
      MV_Set_GEntID(vertex,(max_gvertex_id+1));
      MSTK_Report("buildclass","Found an isolated mesh vertex",WARN);
      continue;
    }
    
    nve = List_Num_Entries(vedges);
    
    GEedges = List_New(nve);
    GFedges = List_New(nve);
    for (i = 0; i < nve; i++) {
      edge = List_Entry(vedges,i);
      gdim2 = ME_GEntDim(edge);
      if (gdim2 == 1)
	List_Add(GEedges,edge);
      if (gdim2 == 2)
	List_Add(GFedges,edge);
    }

    nbe = List_Num_Entries(GEedges);
    nbe2 = List_Num_Entries(GFedges);
    
    switch (nbe) {
    case 0:
      if (nbe2) { 
	
	/* No edges on model edges, only on model face; Classify on
	   model face */
	if (gdim > 2 || (gdim == 2 && gid <= 0)) {
	  MV_Set_GEntDim(vertex,2);
	  gid2 = ME_GEntID(List_Entry(GFedges,0));
	  MV_Set_GEntID(vertex,gid2);
	}
      }
      else {
	
	/* No boundary edges connected to vertex at all; must be interior */
	
	MV_Set_GEntDim(vertex,3);
	gid2 = ME_GEntID(List_Entry(vedges,0));
	MV_Set_GEntID(vertex,gid2);
      }
      break;

    case 1:
      /* Must be the end of a model edge; so its a model vertex */

      if (gdim > 0 || (gdim == 0 && gid <= 0)) {
	MV_Set_GEntDim(vertex,0);
	MV_Set_GEntID(vertex,(max_gvertex_id+1));
	max_gvertex_id++;
      }
      break;

    case 2:      
      geid0 = ME_GEntID(List_Entry(GEedges,0));
      geid1 = ME_GEntID(List_Entry(GEedges,1));
      
      if (geid0 == geid1) { /* Classified on the model edge */
	if (gdim > 1 || (gdim == 1 && gid <= 0)) {
	  MV_Set_GEntDim(vertex,1);
	  MV_Set_GEntID(vertex,geid0);
	}
      }
      else { /* Classified on a model vertex */
	if (gdim > 0 || (gdim == 0 && gid <= 0)) {
	  MV_Set_GEntDim(vertex,0);
	  MV_Set_GEntID(vertex,(max_gvertex_id+1));
	  max_gvertex_id++;
	}
      }
      break;

    default:
      if (gdim > 0 || (gdim == 0 && gid <= 0)) {
	MV_Set_GEntDim(vertex,0);
	MV_Set_GEntID(vertex,(max_gvertex_id+1));
	max_gvertex_id++;
      }
      break;
    }

    List_Delete(vedges);
    List_Delete(GEedges);
    List_Delete(GFedges);
  }

  return 1;
}


#ifdef __cplusplus
}
#endif
