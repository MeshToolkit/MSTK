/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Routine to build classification for mesh faces given no
   classification or classification information for mesh regions only
   (if they exist)

   if use_geometry is 1, the code will try to use geometric measures
   to refine the geometric entity that boundary mesh faces are
   classified on
 */

  int MESH_BuildFaceClassfn(Mesh_ptr mesh, int use_geometry) {
  int i, j, k, idx, idx2, fnd, gfid, gfid2, gdim;
  int ngfaces, ngfalloc, grid0, grid1;
  int max_gface_id;
  int nfe, nef, nbf, nfr, nsub, *gfids, (*gfregids)[2];
  double PI=3.141592, ang, COSSHARPANG;
  MEdge_ptr edge;
  MFace_ptr face, subface, adjface;
  MRegion_ptr freg0, freg1;
  List_ptr fregs, fedges, efaces, ebfaces, gffaces, subfaces;
#ifdef MSTK_USE_MARKERS
  int processedmk, submk;
#else
  MAttrib_ptr processedatt, sublistatt;
  double rval;
  void *pval;
#endif

  COSSHARPANG = cos(5*PI/6);  /* 150 degrees */


  /* Verify that mesh faces on the boundary have classification
     information; if not, assign all faces to the same model faces */

  ngfaces = 0; ngfalloc = 10;
  gfids = (int *) malloc(ngfalloc*sizeof(int));
  gfregids = (int (*)[2]) malloc(ngfalloc*sizeof(int [2]));


  /* Take stock of existing model face information */

  max_gface_id = 0;
  idx = 0;
  while ((face = MESH_Next_Face(mesh,&idx))) {
    gdim = MF_GEntDim(face);
    if (gdim != 2)
      continue;

    gfid = MF_GEntID(face);
    if (gfid) {

      /* Has this model face been encountered? If not, add it to list
	 of model faces */

      for (i = 0, fnd = 0; i < ngfaces; i++)
	if (gfids[i] == gfid) {
	  fnd = 1;
	  break;
	}

      if (!fnd) {
	if (gfid > max_gface_id)
	  max_gface_id = gfid;

	if (ngfalloc == ngfaces) {
	  ngfalloc *= 2;
	  gfids = (int *) realloc(gfids,ngfalloc*sizeof(int));
	  gfregids = (int (*)[2])realloc(gfregids,ngfalloc*sizeof(int [2]));
	}

	gfids[ngfaces] = gfid;
	gfregids[ngfaces][0] = gfregids[ngfaces][1] = 0;
	fregs = MF_Regions(face);
	if (fregs) {
	  nfr = List_Num_Entries(fregs);
      
	  freg0 = List_Entry(fregs,0);
	  gfregids[ngfaces][0] = MR_GEntID(freg0); /* NOTE 1 (see EOF) */
      
	  if (nfr == 2) {
	    freg1 = List_Entry(fregs,1); 
	    gfregids[ngfaces][1] = MR_GEntID(freg1);
	  }	  
	  
	  List_Delete(fregs);
	}

	ngfaces++;
      }

    }
  }


  /* Build new model face information based on adjacent model region info */

  idx = 0;
  while ((face = MESH_Next_Face(mesh,&idx))) {
    
    gdim = MF_GEntDim(face);
    gfid = MF_GEntID(face);

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
      
      if (gdim > 2 || (gdim == 2 && gfid <= 0)) { 
	MF_Set_GEntDim(face,2);

	/* Check if this type of face with these adjacent regions has
	   been encountered before; if it has, we just use the
	   existing model face id */

	for (i = 0, fnd = 0; i < ngfaces; i++)
	  if ((gfregids[i][0] == grid0 && gfregids[i][1] == grid1) ||
	      (gfregids[i][0] == grid1 && gfregids[i][1] == grid0)) {
	    fnd = 1;
	    break;
	  }

	if (fnd)
	  MF_Set_GEntID(face,gfids[i]);
	else {
	  max_gface_id++;
	  MF_Set_GEntID(face,max_gface_id);

	  if (ngfalloc == ngfaces) {
	    ngfalloc *= 2;
	    gfids = (int *) realloc(gfids,ngfalloc*sizeof(int));
	    gfregids = (int (*)[2]) realloc(gfregids,ngfalloc*sizeof(int [2]));
	  }

	  gfids[ngfaces] = max_gface_id;
	  gfregids[ngfaces][0] = grid0;
	  gfregids[ngfaces][1] = grid1;
	  ngfaces++;
	}
      }
    }
  }

  if (use_geometry == 1) {

#ifdef MSTK_USE_MARKERS
    processedmk = MSTK_GetMarker();
    submk = MSTK_GetMarker();
#else
    processedatt = MAttrib_New(mesh, "processedatt", INT, MALLTYPE);
    sublistatt = MAttrib_New(mesh, "sublist", INT, MFACE);
#endif

    /* Now assign model face IDs based on whether a sharp set of edges
       enclose a set of faces */

    for (i = 0; i < ngfaces; i++) {

      /* Find all mesh faces with this model face id */

      gffaces = List_New(10);
      idx = 0; 
      while ((face = MESH_Next_Face(mesh,&idx))) {
        if (MF_GEntDim(face) == 2 && MF_GEntID(face) == gfids[i])
          List_Add(gffaces,face);
      }

      /* Process faces of this list and subdivide them into subfaces */
      /* The way we do that is 
       
         1) we put an unprocessed face from the original list in a subface
         list

         2) we then add its neighboring faces to subface list if they are
         of the same color (same model face id) and do not have a sharp
         edge separating them from the current face

         3) we then process the next face in the subface list 

         4) we are done if we cannot find any more neighbors of faces in
         the subface list to add to the subface list

         5) we then repeat steps 1 through 4 until we are left with no
         more faces to process from the original list
       
      */

      nsub = 0;
      idx = 0;
      while ((face = List_Next_Entry(gffaces,&idx))) {
        int processed;
#ifdef MSTK_USE_MARKERS
        processed = MEnt_IsMarked(face,processedmk);
#else
        MEnt_Get_AttVal(face, processedatt, &processed, &rval, &pval);
#endif
        if (processed) continue;

        /* Found a face in gffaces that has not been processed */
#ifdef MSTK_USE_MARKERS
        MEnt_Mark(face,processedmk);
#else
        MEnt_Set_AttVal(face, processedatt, 1, 0.0, NULL);
#endif

        subfaces = List_New(10);
        List_Add(subfaces,face);
#ifdef MSTK_USE_MARKERS
        MEnt_Mark(face,submk);
#else
        MEnt_Set_AttVal(face, sublistatt, 1, 0.0, NULL);
#endif

        idx2 = 0;
        while ((subface = List_Next_Entry(subfaces,&idx2))) {
          gfid = MF_GEntID(subface);

          fedges = MF_Edges(subface,1,0);
          nfe = List_Num_Entries(fedges);

          for (j = 0; j < nfe; j++) {
            edge = List_Entry(fedges,j);
            efaces = ME_Faces(edge);
            nef = List_Num_Entries(efaces);

            ebfaces = List_New(nef); /* list of boundary faces cnctd 2 edge */
            for (k = 0; k < nef; k++) {
              adjface = List_Entry(efaces,k);
              if (MF_GEntDim(adjface) == 2)
                List_Add(ebfaces,adjface);
            }
            List_Delete(efaces);

            nbf = List_Num_Entries(ebfaces);
            if (nbf == 2) {
              /* we might be on a model face or on a model edge */
            
              adjface = List_Entry(ebfaces,0);
              if (adjface == subface)
                adjface = List_Entry(ebfaces,1);
              gfid2 = MF_GEntID(adjface);

            
              if (gfid == gfid2) {
                /* The two faces are of the same ID. If the angle
                   between them is not sharp they can be classified as
                   being on the same subface */
	      
                ang = MFs_DihedralAngle(subface,adjface,edge);
              
                if (ang <= COSSHARPANG) {
                  /* Add face2 to subface list unless its already there */
                  int inlist;
#ifdef MSTK_USE_MARKERS
                  inlist = MEnt_IsMarked(adjface,submk);
#else
                  MEnt_Get_AttVal(adjface, sublistatt, &inlist, &rval, &pval);
#endif                  
                  if (!inlist) {
                    List_Add(subfaces,adjface);
#ifdef MSTK_USE_MARKERS
                    MEnt_Mark(adjface,submk);
#else
                    MEnt_Set_AttVal(adjface, sublistatt, 1, 0.0, NULL);
#endif
                  }
                }
                else {
                  /* The two faces make a very sharp angle. We will
                     consider the edge b/w them to be a model edge */
                  /* Tag the edge as being on a model edge (we don't
                     know the model edge ID as yet) and continue */
                
                  ME_Set_GEntDim(edge,1);
                  ME_Set_GEntID(edge,0);
                }
              }
              else {
                /* we reached a model edge */
                /* Tag the edge as being on a model edge (we don't know
                   the model edge ID as yet) and continue */

                ME_Set_GEntDim(edge,1);
                ME_Set_GEntID(edge,0);
              }
            }
            else {
              /* we reached a a model edge */
              /* Tag the edge as being on a model edge (we don't know
                 the model edge ID as yet) and continue */

              ME_Set_GEntDim(edge,1);
              ME_Set_GEntID(edge,0);
            }
            List_Delete(ebfaces);
          }

          /* Finished processing all neighbors of the face */
          List_Delete(fedges);
        }

        /* Now we have a list of faces which we believe constitutes a
           model face by itself. If this is the first subface (which
           means it could also be the entire model face originally
           considered), leave the model face tag as it is. If not,
           assign the faces in the subface a new model face ID */

        if (nsub != 0) {
          max_gface_id++;
          idx2 = 0;
          while ((subface = List_Next_Entry(subfaces,&idx2))) 
            MF_Set_GEntID(subface,max_gface_id);
        }
        nsub++;

        /* Done with this subface */

        idx2 = 0;
        while ((subface = List_Next_Entry(subfaces,&idx2))) {
#ifdef MSTK_USE_MARKERS
          MEnt_Mark(subface,processedmk);
          MEnt_Unmark(subface,submk);
#else
          MEnt_Set_AttVal(subface, processedatt, 1, 0.0, NULL);
          MEnt_Set_AttVal(subface, sublistatt, 0, 0.0, NULL);
#endif
        }
        List_Delete(subfaces);
      }

#ifdef MSTK_USE_MARKERS
      List_Unmark(gffaces,processedmk);
#else
      idx = 0;
      while ((face = List_Next_Entry(gffaces,&idx)))
        MEnt_Set_AttVal(face, processedatt, 0, 0.0, NULL);
#endif
      List_Delete(gffaces);
    }
  
#ifdef MSTK_USE_MARKERS
    MSTK_FreeMarker(processedmk);
    MSTK_FreeMarker(submk);
#else
    MAttrib_Delete(processedatt);
    MAttrib_Delete(sublistatt);
#endif
  } /* if use_geometry == 1 */

  free(gfids);
  free(gfregids);
       
  return 1;
}

  /* NOTE 1: Did not take into account which region is on which side
     of face in this code. I think it is unimportant */			       


#ifdef __cplusplus
}
#endif
