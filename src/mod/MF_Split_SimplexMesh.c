/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>

#include "MSTK.h"

/* Split a face in a simplex mesh and also subdivide elements
   (triangles or tetrahedra connected to this edge. If the point is
   inside the face three smaller faces are created in its place.
   There are no checks to ensure that the point is actually inside the
   face and not on its boundary or outside the face. This is a
   strictly topological routine. Geometric checks for validity must be
   done outside this routine */


MVertex_ptr MF_Split_SimplexMesh(MFace_ptr fsplit, double *splitxyz) {
  int i, j, k, rfdir=1, ntets=0, ntris=0, *rid=NULL, fgdim, fgid, found;
  MVertex_ptr vsplit, ev[2], (*tetverts)[4]=NULL, triverts[3], fv;
  MVertex_ptr fvarr[3], rvarr[4];
  MFace_ptr f;
  MRegion_ptr r;
  List_ptr fedges, ftets, rfaces, fverts;
  Mesh_ptr mesh = MF_Mesh(fsplit);

  /* point is not on the boundary of the face */
    
  ftets = MF_Regions(fsplit);
  if (ftets) {
    ntets = List_Num_Entries(ftets);   
    tetverts = (MVertex_ptr (*)[4]) malloc(ntets*sizeof(MVertex_ptr [4]));
    rid = (int *) malloc(ntets*sizeof(int));
  }

  for (i = 0; i < ntets; i++) {

    r = List_Entry(ftets,i);

    rfaces = MR_Faces(r);

    /* Find the face to be split and get the first three vertices
       in a suitable order from it. Also, find another face and 
       get a vertex that is not in the face to be split. This vertex
       forms the fourth vertex of the tet */

    for (j = 0; j < 4; j++) {
      f = List_Entry(rfaces,j);

      if (f == fsplit) {
        rfdir = MR_FaceDir_i(r,j);        
        fverts = MF_Vertices(f,!rfdir,0);

        for (k = 0; k < 3; k++) {
          tetverts[i][0] = List_Entry(fverts,0);
          tetverts[i][1] = List_Entry(fverts,1);
          tetverts[i][2] = List_Entry(fverts,2);
        }
        List_Delete(fverts);
        break;
      }
    }

    found = 0;
    for (j = 0; j < 4; j++) {
      f = List_Entry(rfaces,j);
      if (f != fsplit) {
        fverts = MF_Vertices(f,!rfdir,0);
        for (k = 0; k < 3; k++) {
          fv = List_Entry(fverts,k);
          if (!MF_UsesEntity(fsplit,fv,MVERTEX)) {
            tetverts[i][3] = fv;
            found = 1;
            break;
          }
        }
        List_Delete(fverts);
      }
      if (found) break;
    }

    List_Delete(rfaces);
  }

  /* Now that we finished collecting info about the connected tets we
     can delete them */

  if (ftets) {
    for (i = 0; i < ntets; i++)
      MR_Delete(List_Entry(ftets,i),0);

    List_Delete(ftets);
  }

  /* Delete the face itself */

  fverts = MF_Vertices(fsplit,1,0);
  for (i = 0; i < 3; i++)
    triverts[i] = List_Entry(fverts,i);
  List_Delete(fverts);

  fgdim = MF_GEntDim(fsplit);
  fgid  = MF_GEntID(fsplit);


  /* Split the face */

  vsplit = MF_Split(fsplit, splitxyz);

  
  /* Create three tets for each tet that was deleted */

  for (i = 0; i < ntets; i++) {
    for (j = 0; j < 3; j++) {
      r = MR_New(mesh);
      rvarr[0] = vsplit;
      rvarr[1] = tetverts[i][j];
      rvarr[2] = tetverts[i][(j+1)%3];
      rvarr[3] = tetverts[i][3];
      MR_Set_Vertices(r, 4, rvarr, 0, NULL);
      MR_Set_GEntID(r,rid[i]);
    }
  }

  if (ntets) {
    free(tetverts);
    free(rid);
  }

  return vsplit;
}
