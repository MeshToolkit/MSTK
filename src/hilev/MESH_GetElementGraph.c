/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Mesh.h"
#include "MSTK.h"
#include "MSTK_private.h"

/* Get the element adjacency graph for a mesh */
/* adjbeg - pointer to array containing the offsets of the adjelems list at */
/*          which the adjacent elements are enumerated for each element     */

int MESH_GetElementGraph(Mesh_ptr mesh, int *nelems, int **adjbeg,
			 int **adjelems) {

  *nelems = MESH_Num_Regions(mesh);
  if (*nelems) {  /* 3D mesh */

    /* Now build the graph */
    *adjbeg = (int *) malloc((*nelems+1)*sizeof(int));
    int nalloc = 4*(*nelems);  /* start with assuming 4 nbrs per entity */
    *adjelems = (int *) malloc(nalloc*sizeof(int));

    int pos = 0;
    int i = 0;
    int idx = 0;
    MRegion_ptr mr;
    while ((mr = MESH_Next_Region(mesh, &idx))) {
      (*adjbeg)[i] = pos;

      List_ptr rfaces = MR_Faces(mr);
      int idx2 = 0;
      MFace_ptr rf;
      while ((rf = List_Next_Entry(rfaces, &idx2))) {
	List_ptr fregs = MF_Regions(rf);
	if (List_Num_Entries(fregs) == 1) continue;

	MRegion_ptr oppr;
	if (List_Entry(fregs,0) == mr)
	  oppr = List_Entry(fregs,1);
	else
	  oppr = List_Entry(fregs,0);
	if (pos == nalloc) {
	  nalloc *= 2;
	  *adjelems = (int *) realloc(adjelems, nalloc*sizeof(int));
	}
	(*adjelems)[pos] = MR_ID(oppr)-1;
	pos++;
	
	List_Delete(fregs);
      }
      List_Delete(rfaces);
      i++;
    }
    (*adjbeg)[*nelems] = pos;

  } else {

    *nelems = MESH_Num_Faces(mesh);

    *adjbeg = (int *) malloc((*nelems+1)*sizeof(int));
    int nalloc = 4*(*nelems);  /* start with assuming 4 nbrs per entity */
    *adjelems = (int *) malloc(nalloc*sizeof(int));

    int pos = 0;
    int i = 0;
    int idx = 0;
    MFace_ptr mf;
    while ((mf = MESH_Next_Face(mesh, &idx))) {
      (*adjbeg)[i] = pos;

      List_ptr fedges = MF_Edges(mf, 1, 0);
      int idx2 = 0;
      MEdge_ptr fe;
      while ((fe = List_Next_Entry(fedges, &idx2))) {
	List_ptr efaces = ME_Faces(fe);
	if (List_Num_Entries(efaces) == 1) continue;

	/* have to account for non-manifold configurations */
	int idx3 = 0;
	MFace_ptr ef;
	while ((ef = List_Next_Entry(efaces, &idx3))) {
	  if (ef != mf) {
	    if (pos == nalloc) {
	      nalloc *= 2;
	      *adjelems = (int *) realloc(*adjelems, nalloc*sizeof(int));
	    }
	    (*adjelems)[pos] = MF_ID(ef)-1;
	    pos++;
	  }
	}
	List_Delete(efaces);
      }
      List_Delete(fedges);
      i++;
    }
    (*adjbeg)[*nelems] = pos;
  }

  return 1;
}
