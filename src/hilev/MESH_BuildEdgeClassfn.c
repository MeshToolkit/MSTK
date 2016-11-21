#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Routine to build classification for mesh edges given no
   classification or classification information for mesh regions only
   (if they exist)

   if use_geometry is 1, then the code tries to refine the
   classification of mesh edges on model edges based on geometric
   measures

 */

  int MESH_BuildEdgeClassfn(Mesh_ptr mesh, int use_geometry) {
  int i, j, k, idx, idx2, fnd, fnd2, geid, geid2, gdim;
  int ngedges, ngealloc, ngef, max_loc_gfids, *loc_gfids;
  int max_gedge_id;
  int nve, nbe, nef, nsub, *geids, **gefaceids;
  double PI=3.141592, cosang, COSSHARPANG;
  MVertex_ptr ev[2];
  MEdge_ptr edge, subedge, adjedge;
  MFace_ptr eface;
  List_ptr efaces, vedges, vbedges, geedges, subedges;
  double rval;
  void *pval;

  COSSHARPANG = cos(5*PI/6);  /* 165 degrees */


  /* Verify that mesh edges on the boundary have classification
     information; if not, assign all edges to the same model edge */

  ngedges = 0; ngealloc = 10;
  geids = (int *) malloc(ngealloc*sizeof(int)); /* model edge ids */
  /* Number of model faces connected to edge followed by model face ids */
  gefaceids = (int **) malloc(ngealloc*sizeof(int *));


  /* Take stock of existing model edge information */

  max_gedge_id = 0;
  idx = 0;
  while ((edge = MESH_Next_Edge(mesh,&idx))) {
    gdim = ME_GEntDim(edge);
    if (gdim != 1)
      continue;

    geid = ME_GEntID(edge);
    if (geid) {

      /* Has this model edge been encountered? If not, add it to list
	 of model edges */

      for (i = 0, fnd = 0; i < ngedges; i++)
	if (geids[i] == geid) {
	  fnd = 1;
	  break;
	}

      if (!fnd) {
	if (geid > max_gedge_id)
	  max_gedge_id = geid;

	if (ngealloc == ngedges) {
	  ngealloc *= 2;
	  geids = (int *) realloc(geids,ngealloc*sizeof(int));
	  gefaceids = (int **)realloc(gefaceids,ngealloc*sizeof(int *));
	}

	geids[ngedges] = geid;

	efaces = ME_Faces(edge);
	nef = List_Num_Entries(efaces);

	gefaceids[ngedges] = (int *) malloc((1+nef)*sizeof(int));
        ngef = 0;
	for (i = 0; i < nef; i++) {
	  eface = List_Entry(efaces,i);
	  if (MF_GEntDim(eface) == 2) {
	    gefaceids[ngedges][1+ngef] = MF_GEntID(eface);
	    ngef++;
	  }
	}
	gefaceids[ngedges][0] = ngef;

	ngedges++;

	List_Delete(efaces);
      }

    }
  }


  /* Build new model edge information based on adjacent model region info */
  
  /* NOTE: The following cases involve some repetition and can be
     "cleverly" folded into a shorter piece of code, but then the
     method by which the individual cases are handled gets a littled
     obscured. So leave as is */

  max_loc_gfids = 10;
  loc_gfids = (int *) malloc(max_loc_gfids*sizeof(int));
  idx = 0;
  while ((edge = MESH_Next_Edge(mesh,&idx))) {
    
    gdim = ME_GEntDim(edge);
    geid = ME_GEntID(edge);

    /* Edge has no classification? Assign classification info */
    /* Edge classified as face or interior edge? Verify */
    
    efaces = ME_Faces(edge);
    nef = efaces ? List_Num_Entries(efaces) : 0;

    if (!nef) {
      /* Isolated edge (must be on model edge). Did we encounter such
	 an edge before? */

      ME_Set_GEntDim(edge,1);

      for (i = 0, fnd = 0; i < ngedges; i++) {
	if (gefaceids[i][0] == 0) {
	  fnd = 1;
	  break;
	}
      }

      if (fnd) 
	ME_Set_GEntID(edge,geids[i]);
      else {
	max_gedge_id++;
	ME_Set_GEntID(edge,max_gedge_id);
	    
	if (ngealloc == ngedges) {
	  ngealloc *= 2;
	  geids = (int *) realloc(geids,ngealloc*sizeof(int));
	  gefaceids = (int **) realloc(gefaceids,ngealloc*sizeof(int *));
	}
	    
	geids[ngedges] = max_gedge_id;
	gefaceids[ngedges] = malloc(1*sizeof(int *));
	gefaceids[ngedges][0] = 0;
      }

      continue; /* nothing else to do */
    }

    if (nef > max_loc_gfids) {      
      loc_gfids = (int *) realloc(loc_gfids,nef*sizeof(int));
      max_loc_gfids = nef;
    }

    ngef = 0;
    for (i = 0; i < nef; i++) {
      eface = List_Entry(efaces,i);
      if (MF_GEntDim(eface) == 2) {
	loc_gfids[ngef] = MF_GEntID(eface);
	ngef++;
      }
    }
    
    switch (ngef) {
    case 0:
      /* Interior edge - we took care of the case of the isolated edge b4 */

      ME_Set_GEntDim(edge,3);
      eface = List_Entry(efaces,0);
      ME_Set_GEntID(edge,MF_GEntID(eface));

      break;
    case 2: 
      if (loc_gfids[0] == loc_gfids[1]) {

	if (gdim == 1) {
	  /* Looks like during face classification, this was tagged as
	     being a sharp edge. This means that it is a sharp edge on
	     the interior of a model face (same model face on both
	     sides) */

	  ME_Set_GEntDim(edge,1);

	  /* Check if such an edge was encountered before */

	  for (i = 0, fnd = 0; i < ngedges; i++) {
	    if (gefaceids[i][0] != 2)
	      continue;

	    if ((loc_gfids[0] == gefaceids[i][0]) && 
		(gefaceids[i][0] == gefaceids[i][1])) {
	      fnd = 1;
	      break;
	    }
	  }

	  if (fnd) 
	    ME_Set_GEntID(edge,geids[i]);
	  else {
	    max_gedge_id++;
	    ME_Set_GEntID(edge,max_gedge_id);
	    
	    if (ngealloc == ngedges) {
	      ngealloc *= 2;
	      geids = (int *) realloc(geids,ngealloc*sizeof(int));
	      gefaceids = (int **) realloc(gefaceids,ngealloc*sizeof(int *));
	    }
	    
	    geids[ngedges] = max_gedge_id;
	    gefaceids[ngedges] = malloc((1+ngef)*sizeof(int *));
	    gefaceids[ngedges][0] = ngef;
	    for (i = 0; i < ngef; i++)
	      gefaceids[ngedges][1+i] = loc_gfids[i];
	    ngedges++;
	  }
	}
	else {
	  /* edge must be on model face */

	  ME_Set_GEntDim(edge,2);

	  ME_Set_GEntID(edge,loc_gfids[0]);
	}
      }
      else {
	/* edge is on model edge between two faces */

	ME_Set_GEntDim(edge,1);

	/* Check if such an edge was encountered before */
	
	for (i = 0, fnd = 0; i < ngedges; i++) {
	  if (gefaceids[i][0] != 2)
	    continue;
	  
	  if (((loc_gfids[0] == gefaceids[i][1]) && 
	       (loc_gfids[1] == gefaceids[i][2])) ||
	      ((loc_gfids[1] == gefaceids[i][1]) &&
	       (loc_gfids[0] == gefaceids[i][2])))   {
	    fnd = 1;
	    break;
	  }
	}

	if (fnd) 
	  ME_Set_GEntID(edge,geids[i]);
	else {
	  max_gedge_id++;
	  ME_Set_GEntID(edge,max_gedge_id);
	  
	  if (ngealloc == ngedges) {
	    ngealloc *= 2;
	    geids = (int *) realloc(geids,ngealloc*sizeof(int));
	    gefaceids = (int **) realloc(gefaceids,ngealloc*sizeof(int *));
	  }
	  
	  geids[ngedges] = max_gedge_id;
	  gefaceids[ngedges] = malloc((1+ngef)*sizeof(int *));
	  gefaceids[ngedges][0] = ngef;
	  for (i = 0; i < ngef; i++)
	    gefaceids[ngedges][1+i] = loc_gfids[i];
	  ngedges++;
	}
      }
      break;

    default:
      /* if ngef is 1, edge is on model edge of non-manifold model
	 face 
	 if ngef is >= 3, edge is on model edge at junction of
	 many model faces
      */
      
      ME_Set_GEntDim(edge,1);
      
      /* Check if a previously encountered model edge has this
	 combination of connected faces */
      
      for (i = 0, fnd = 0; i < ngedges; i++) {
	if (ngef != gefaceids[i][0]) 
	  continue; /* number of connected model faces is different */
	
	/* see if all the model faces of the current edge can be found
	   in the i'th model edge's faces */

	fnd2 = 0;
	for (j = 0; j < ngef; j++) {
	  
	  fnd2 = 0;
	  for (k = 0; k < ngef; k++) {
	    if (loc_gfids[j] == gefaceids[i][1+k]) {
	      fnd2 = 1;
	      break;
	    }
	  }

	  /* Did not find loc_gfid[j] in the list gefaceids[i] */

	  if (!fnd2)
	    break;
	}

	/* if a model face connected to this edge was not found in the
	   model face list of the previously processed, then the two model
	   edges are obviously different */

	if (!fnd2)
	  continue;
	else {
	  fnd = 1;
	  break;
	}
      }

      if (fnd)
	ME_Set_GEntID(edge,geids[i]);
      else {
	max_gedge_id++;
	ME_Set_GEntID(edge,max_gedge_id);
	
	if (ngealloc == ngedges) {
	  ngealloc *= 2;
	  geids = (int *) realloc(geids,ngealloc*sizeof(int));
	  gefaceids = (int **) realloc(gefaceids,ngealloc*sizeof(int *));
	}
	
	geids[ngedges] = max_gedge_id;
	gefaceids[ngedges] = malloc((1+ngef)*sizeof(int *));
	gefaceids[ngedges][0] = ngef;
	for (i = 0; i < ngef; i++)
	  gefaceids[ngedges][1+i] = loc_gfids[i];
	ngedges++;
      }
      break;
    }

    List_Delete(efaces);    /* needed efaces in case 0 */
  }
  free(loc_gfids);


  if (use_geometry == 1) {


#ifdef MSTK_USE_MARKERS
    int processedmk = MSTK_GetMarker();
    int submk = MSTK_GetMarker();
#else
    MAttrib_ptr processedatt = MAttrib_New(mesh, "processed", INT, MEDGE);
    MAttrib_ptr sublistatt = MAttrib_New(mesh, "sublist", INT, MEDGE);
#endif

    /* Now assign model edge IDs based on whether a sharp set of edges
       enclose a set of edges */

    for (i = 0; i < ngedges; i++) {

      /* Find all mesh edges with this model edge id */

      geedges = List_New(10);
      idx = 0; 
      while ((edge = MESH_Next_Edge(mesh,&idx))) {
        if (ME_GEntDim(edge) == 1 && ME_GEntID(edge) == geids[i])
          List_Add(geedges,edge);
      }

      /* Process edges of this list and subdivide them into subedges */
      /* The way we do that is 
       
         1) we put an unprocessed edge from the original list in a subedge
         list

         2) we then add its neighboring edges to subedge list if they are
         of the same color (same model edge id) and do not have a sharp
         edge separating them from the current edge

         3) we then process the next edge in the subedge list 

         4) we are done if we cannot find any more neighbors of edges in
         the subedge list to add to the subedge list

         5) we then repeat steps 1 through 4 until we are left with no
         more edges to process from the original list
       
      */
      nsub = 0;
      idx = 0;
      while ((edge = List_Next_Entry(geedges,&idx))) {
        int emarked;
#ifdef MSTK_USE_MARKERS
        emarked = MEnt_IsMarked(edge,processedmk);
#else
        MEnt_Get_AttVal(edge, processedatt, &emarked, &rval, &pval);
#endif
        if (emarked)
          continue;
        
        /* Found a edge in geedges that has not been processed */
#ifdef MSTK_USE_MARKERS
        MEnt_Mark(edge,processedmk);
#else
        MEnt_Set_AttVal(edge, processedatt, 1, 0.0, NULL);
#endif

        subedges = List_New(10);
        List_Add(subedges,edge);
#ifdef MSTK_USE_MARKERS
        MEnt_Mark(edge,submk);
#else
        MEnt_Set_AttVal(edge, sublistatt, 1, 0.0, NULL);
#endif

        idx2 = 0;
        while ((subedge = List_Next_Entry(subedges,&idx2))) {
          geid = ME_GEntID(subedge);

          ev[0] = ME_Vertex(subedge,0);
          ev[1] = ME_Vertex(subedge,1);

          for (j = 0; j < 2; j++) {
            vedges = MV_Edges(ev[j]);
            nve = List_Num_Entries(vedges);

            vbedges = List_New(nve); /* list of boundary edges cnctd 2 vert */
            for (k = 0; k < nve; k++) {
              adjedge = List_Entry(vedges,k);
              if (ME_GEntDim(adjedge) == 1)
                List_Add(vbedges,adjedge);
            }

            nbe = List_Num_Entries(vbedges);
            if (nbe == 2) {
              /* we might be on a model vertex or on a model edge */

              adjedge = List_Entry(vbedges,0);
              if (adjedge == subedge)
                adjedge = List_Entry(vbedges,1);
              geid2 = ME_GEntID(adjedge);

              if (geid == geid2) {
                /* The two edges are of the same ID. If the angle
                   between them is not sharp they can be classified as
                   being on the same subedge */
	      
                cosang = MEs_Angle(subedge,adjedge);

                if (cosang <= COSSHARPANG) {
                  /* Add edge2 to subedge list unless its already there */
                  int adjemarked;
#ifdef MSTK_USE_MARKERS
                  adjemarked = MEnt_IsMarked(adjedge,submk);
#else
                  MEnt_Get_AttVal(adjedge, sublistatt, &adjemarked, &rval, &pval);
#endif
                  if (!adjemarked) {
                    List_Add(subedges,adjedge);
#ifdef MSTK_USE_MARKERS
                    MEnt_Mark(adjedge,submk);
#else
                    MEnt_Set_AttVal(adjedge, sublistatt, 1, 0.0, NULL);
#endif
                  }
                }
                else {
                  /* The two edges make a very sharp angle. We will
                     consider the edge b/w them to be a model vertex */
                  /* Tag the edge as being on a model vertex (we don't
                     know the model vertex ID as yet) and continue */

                  MV_Set_GEntDim(ev[j],0);
                  MV_Set_GEntID(ev[j],0);
                }
              }
              else {
                /* we reached a model vertex */
                /* Tag the edge as being on a model vertex (we don't know
                   the model vertex ID as yet) and continue */

                MV_Set_GEntDim(ev[j],0);
                MV_Set_GEntID(ev[j],0);
              }
            }
            else {
              /* we reached a a model vertex */
              /* Tag the edge as being on a model vertex (we don't know
                 the model vertex ID as yet) and continue */

              ME_Set_GEntDim(ev[j],0);
              ME_Set_GEntID(ev[j],0);
            }
            List_Delete(vedges);
            List_Delete(vbedges);	  
          }

          /* Finished processing all neighbors of the edge */
        }

        /* Now we have a list of edges which we believe constitutes a
           model edge by itself. If this is the first subedge (which
           means it could also be the entire model edge originally
           considered), leave the model edge tag as it is. If not,
           assign the edges in the subedge a new model edge ID */

        if (nsub != 0) {
          max_gedge_id++;
          idx2 = 0;
          while ((subedge = List_Next_Entry(subedges,&idx2))) 
            ME_Set_GEntID(subedge,max_gedge_id);
        }
        nsub++;

        /* Done with this subedge */

#ifdef MSTK_USE_MARKERS
        idx2 = 0;
        while ((subedge = List_Next_Entry(subedges,&idx2))) {
          MEnt_Mark(subedge,processedmk);
          MEnt_Unmark(subedge,submk);
        }
#else
        idx2 = 0;
        while ((subedge = List_Next_Entry(subedges,&idx2))) {
          MEnt_Set_AttVal(subedge, processedatt, 1, 0.0, NULL);
          MEnt_Set_AttVal(subedge, sublistatt, 0, 0.0, NULL);
        }
#endif
        List_Delete(subedges);
      }

#ifdef MSTK_USE_MARKERS
      List_Unmark(geedges,processedmk);
#else
      idx2 = 0;
      while ((edge = List_Next_Entry(geedges, &idx2)))
        MEnt_Set_AttVal(edge, processedatt, 1, 0.0, NULL);
#endif
      List_Delete(geedges);
    }

#ifdef MSTK_USE_MARKERS
      MSTK_FreeMarker(processedmk);
      MSTK_FreeMarker(submk);
#else
      MAttrib_Delete(processedatt);
      MAttrib_Delete(sublistatt);
#endif
  } /* if use_geometry == 1 */

  free(geids);
  for (i = 0; i < ngedges; i++)
    free(gefaceids[i]);
  free(gefaceids);
       
  return 1;
}


#ifdef __cplusplus
}
#endif
