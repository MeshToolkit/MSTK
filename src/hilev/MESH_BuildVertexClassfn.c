#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Routine to build classification for mesh vertexs given no
   classification or classification information for mesh regions only
   (if they exist)
 */

int MESH_BuildVertexClassfn(Mesh_ptr mesh) {
  int i, j, k, idx, idx2, fnd, fnd2, gvid, gvid2, gdim;
  int ngverts, ngvalloc, ngve, max_loc_geids, *loc_geids;
  int max_gvert_id, processedmk, submk;
  int nbe, nve, nsub, *gvids, **gvedgeids, gfid;
  MVertex_ptr vertex;
  MEdge_ptr vedge;
  List_ptr vedges, vbedges;


  /* Verify that mesh edges on the boundary have classification
     information; if not, assign all edges to the same model vertex */

  ngverts = 0; ngvalloc = 10;
  gvids = (int *) malloc(ngvalloc*sizeof(int)); /* model vertex ids */
  /* Number of model faces connected to vertex followed by model face ids */
  gvedgeids = (int **) malloc(ngvalloc*sizeof(int *));


  /* Take stock of existing model vertex information */

  max_gvert_id = 0;
  idx = 0;
  while ((vertex = MESH_Next_Vertex(mesh,&idx))) {
    gdim = MV_GEntDim(vertex);
    if (gdim != 0)
      continue;

    gvid = MV_GEntID(vertex);
    if (gvid) {
      if (gvid > max_gvert_id)
	max_gvert_id = gvid;
      
      if (ngvalloc == ngverts) {
	ngvalloc *= 2;
	gvids = (int *) realloc(gvids,ngvalloc*sizeof(int));
	gvedgeids = (int **)realloc(gvedgeids,ngvalloc*sizeof(int *));
      }

      gvids[ngverts++] = gvid;
      
      vedges = MV_Edges(vertex);
      nve = List_Num_Entries(vedges);

      gvedgeids[i] = (int *) malloc((1+nve)*sizeof(int));
      ngve = 0;
      for (i = 0; i < nve; i++) {
	vedge = List_Entry(vedges,i);
	if (ME_GEntDim(vedge) == 1) {
	  gvedgeids[i][1+ngve] = ME_GEntID(vedge);
	  ngve++;
	}
      }
      gvedgeids[i][0] = ngve;

      List_Delete(vedges);
    }
  }


  /* Build new model vertex information based on adjacent model edge info */
  
  /* NOTE: The following cases involve some repetition and can be
     "cleverly" folded into a shorter piece of code, but then the
     method by which the individual cases are handled gets a littled
     obscured. So leave as is */

  max_loc_geids = 10;
  loc_geids = (int *) malloc(max_loc_geids*sizeof(int));
  idx = 0;
  while ((vertex = MESH_Next_Vertex(mesh,&idx))) {
    
    gdim = MV_GEntDim(vertex);
    gvid = MV_GEntID(vertex);

    /* Vertex has no classification? Assign classification info */
    /* Vertex classified on edge, face or region? Verify */
    
    vedges = MV_Edges(vertex);
    nve = vedges ? List_Num_Entries(vedges) : 0;

    if (!nve) {
      /* Isolated vertex (must be a model vertex) */

      MV_Set_GEntDim(vertex,0);

      max_gvert_id++;
      MV_Set_GEntID(vertex,max_gvert_id);
	
      if (ngvalloc == ngverts) {
	ngvalloc *= 2;
	gvids = (int *) realloc(gvids,ngvalloc*sizeof(int));
	gvedgeids = (int **) realloc(gvedgeids,ngvalloc*sizeof(int *));
      }
	    
      gvids[ngverts] = max_gvert_id;
      gvedgeids[ngverts] = malloc(1*sizeof(int *));
      gvedgeids[ngverts][0] = 0;

      continue; /* nothing else to do */
    }

    if (nve > max_loc_geids) {      
      loc_geids = (int *) realloc(loc_geids,nve*sizeof(int));
      max_loc_geids = nve;
    }

    ngve = 0;
    for (i = 0; i < nve; i++) {
      vedge = List_Entry(vedges,i);
      if (ME_GEntDim(vedge) == 1) {
	loc_geids[ngve] = MV_GEntID(vedge);
	ngve++;
      }
    }
    
    switch (ngve) {
    case 0:
      /* Vertex on model face or in the interior - we took care of the
	 case of the isolated vertex b4 */

      fnd = 0; gfid = 0;
      for (i = 0; i < nve; i++) {
	vedge = List_Entry(vedges,i);
	if (ME_GEntDim(vedge) == 2) {
	  fnd = 1;
	  gfid = ME_GEntID(vedge);
	  break;
	}
      }

      if (fnd) {
	/* found at least one edge classified on model face - must be
	   classified on model face */
	MV_Set_GEntDim(vertex,2);
	MV_Set_GEntID(vertex,gfid);
      }
      else { /* no boundary edges, interior vertex */
	MV_Set_GEntDim(vertex,3);
	vedge = List_Entry(vedges,0);
	MV_Set_GEntID(vertex,ME_GEntID(vedge));
      }

      break;
    case 2: 
      if (loc_geids[0] == loc_geids[1]) {

	if (gdim == 0) {
	  /* Looks like during edge classification, this was tagged as
	     being a sharp corner. */

	  MV_Set_GEntDim(vertex,0);

	  max_gvert_id++;
	  MV_Set_GEntID(vertex,max_gvert_id);
	    
	  if (ngvalloc == ngverts) {
	    ngvalloc *= 2;
	    gvids = (int *) realloc(gvids,ngvalloc*sizeof(int));
	    gvedgeids = (int **) realloc(gvedgeids,ngvalloc*sizeof(int *));
	  }
	    
	  gvids[ngverts] = max_gvert_id;
	  gvedgeids[ngverts] = malloc((1+ngve)*sizeof(int *));
	  gvedgeids[ngverts][0] = ngve;
	  for (i = 0; i < ngve; i++)
	    gvedgeids[ngverts][1+i] = loc_geids[i];
	  ngverts++;
	}
	else {
	  /* vertex must be on model edge */

	  MV_Set_GEntDim(vertex,1);

	  MV_Set_GEntID(vertex,loc_geids[0]);
	}
      }
      else {
	/* vertex is on model vertex between two model edges */

	MV_Set_GEntDim(vertex,0);

	/* Check if such a vertex was encountered before (if so, we
	   have coincident mesh vertices) */
	
	for (i = 0, fnd = 0; i < ngverts; i++) {
	  if (gvedgeids[i][0] != 2)
	    continue;
	  
	  if (((loc_geids[0] == gvedgeids[i][1]) && 
	       (loc_geids[1] == gvedgeids[i][2])) ||
	      ((loc_geids[1] == gvedgeids[i][1]) &&
	       (loc_geids[0] == gvedgeids[i][2])))   {
	    fnd = 1;
	    break;
	  }
	}

	if (fnd) 
	  MV_Set_GEntID(vertex,gvids[i]);
	else {
	  max_gvert_id++;
	  MV_Set_GEntID(vertex,max_gvert_id);
	  
	  if (ngvalloc == ngverts) {
	    ngvalloc *= 2;
	    gvids = (int *) realloc(gvids,ngvalloc*sizeof(int));
	    gvedgeids = (int **) realloc(gvedgeids,ngvalloc*sizeof(int *));
	  }
	  
	  gvids[ngverts] = max_gvert_id;
	  gvedgeids[ngverts] = malloc((1+ngve)*sizeof(int *));
	  gvedgeids[ngverts][0] = ngve;
	  for (i = 0; i < ngve; i++)
	    gvedgeids[ngverts][1+i] = loc_geids[i];
	  ngverts++;
	}
      }
      break;

    default:
      /* if ngve is 1, vertex is on model vertex at the end of an edge
	 (edge surrounded by same model face on both sides) if ngve is
	 >= 3, vertex is on model vertex at junction of many model
	 edges
      */
      
      MV_Set_GEntDim(vertex,0);
      
      /* Check if a previously encountered model vertex has this
	 combination of connected edges (if so, then we have
	 coincident mesh vertices) */
      
      for (i = 0, fnd = 0; i < ngverts; i++) {
	if (ngve != gvedgeids[i][0]) 
	  continue; /* number of connected model edges is different */
	
	/* see if all the model edges of the current vertex can be
	   found in the i'th model vertex's edges */

	fnd2 = 0;
	for (j = 0; j < ngve; j++) {
	  
	  fnd2 = 0;
	  for (k = 0; k < ngve; k++) {
	    if (loc_geids[j] == gvedgeids[i][1+k]) {
	      fnd2 = 1;
	      break;
	    }
	  }

	  /* Did not find loc_gfid[j] in the list gvedgeids[i] */

	  if (!fnd2)
	    break;
	}

	/* if a model edge connected to this vertex was not found in
	   the model edge list of the previously processed vertex,
	   then the two model vertices are obviously different */

	if (!fnd2)
	  continue;
	else {
	  fnd = 1;
	  break;
	}
      }

      if (fnd)
	MV_Set_GEntID(vertex,gvids[i]);
      else {
	max_gvert_id++;
	MV_Set_GEntID(vertex,max_gvert_id);
	
	if (ngvalloc == ngverts) {
	  ngvalloc *= 2;
	  gvids = (int *) realloc(gvids,ngvalloc*sizeof(int));
	  gvedgeids = (int **) realloc(gvedgeids,ngvalloc*sizeof(int *));
	}
	
	gvids[ngverts] = max_gvert_id;
	gvedgeids[ngverts] = malloc((1+ngve)*sizeof(int *));
	gvedgeids[ngverts][0] = ngve;
	for (i = 0; i < ngve; i++)
	  gvedgeids[ngverts][1+i] = loc_geids[i];
	ngverts++;
      }
      break;
    }

    List_Delete(vedges);    /* needed vedges in case 0 */
  }
  free(loc_geids);


  free(gvids);
  for (i = 0; i < ngverts; i++)
    free(gvedgeids[i]);
  free(gvedgeids);
       
  return 1;
}

#ifdef __cplusplus
}
#endif
