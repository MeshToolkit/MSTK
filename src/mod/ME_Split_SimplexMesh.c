#include <stdio.h>
#include <stdlib.h>

#include "MSTK.h"

/* Split an edge in a simplex mesh and also subdivide elements (triangles or tetrahedra connected to this edge */

MVertex_ptr ME_Split_SimplexMesh(MEdge_ptr esplit, double *splitxyz) {
  int i, j, k, rfdir, ntets=0, ntris=0, *fdim, *fid, *rid, found;
  MVertex_ptr vsplit, ev[2], (*tetverts)[4], (*triverts)[3], fv;
  MVertex_ptr fvarr[3], rvarr[4];
  MFace_ptr f;
  MRegion_ptr r;
  List_ptr etets, rfaces, etris, fverts;
  Mesh_ptr mesh = ME_Mesh(esplit);

  ev[0] = ME_Vertex(esplit,0);
  ev[1] = ME_Vertex(esplit,1);

  etets = ME_Regions(esplit);
  if (etets) {
    ntets = List_Num_Entries(etets);   
    tetverts = (MVertex_ptr (*)[4]) malloc(ntets*sizeof(MVertex_ptr [4]));
    rid = (int *) malloc(ntets*sizeof(int));
  }

  for (i = 0; i < ntets; i++) {

    r = List_Entry(etets,i);

    rfaces = MR_Faces(r);

    /* Find a tet face that uses ev[0] but not ev[1] */

    found = 0;
    for (j = 0; !found && j < 4; j++) {
      f = List_Entry(rfaces,j);

      fverts = MF_Vertices(f,1,0);
      if (List_Contains(fverts,ev[0]) &&
          !List_Contains(fverts,ev[1])) {

        found = 1;
        
        /* Get the two vertices (a,b) of this face excluding ev[0] in
           such an order that ev[0],a,b,ev[1] will form a valid
           tet. This requires checking whether the face points into or
           out of this tet (look at rfdir) */

        rfdir = MR_FaceDir_i(r,j);
        for (k = 0; k < 3; k++) {
          fv = List_Entry(fverts,k);
          if (fv == ev[0]) {
            tetverts[i][0] = ev[0];
            tetverts[i][1] = rfdir ? List_Entry(fverts,(k+2)%3) : List_Entry(fverts,(k+1)%3);
            tetverts[i][2] = rfdir ? List_Entry(fverts,(k+1)%3) : List_Entry(fverts,(k+2)%3);
            tetverts[i][3] = ev[1];
          }
        }
      }
      List_Delete(fverts);

      if (found) break;
    }

    List_Delete(rfaces);
  }

  /* Now that we finished collecting info about the connected tets we
     can delete them */

  if (etets) {
    for (i = 0; i < ntets; i++)
      MR_Delete(List_Entry(etets,i),0);

    List_Delete(etets);
  }

  /* Now get the triangular face connected to the edge. For each
     triangular face, record the vertex opposite to edge esplit and
     delete the triangular face */

  etris = ME_Faces(esplit);
  if (etris) {
    ntris = List_Num_Entries(etris);
    triverts = (MVertex_ptr (*)[3]) malloc(ntris*sizeof(MVertex_ptr[3]));
    fdim = (int *) malloc(ntris*sizeof(int));
    fid = (int *) malloc(ntris*sizeof(int));
  }
  
  for (i = 0; i < ntris; i++) {
    f = List_Entry(etris,i);

    fverts = MF_Vertices(f,1,0);
    for (j = 0; j < 3; j++) {
      fv = List_Entry(fverts,j);
      if (fv != ev[0] && fv != ev[1]) {        
        triverts[i][0] = fv;
        triverts[i][1] = List_Entry(fverts,(j+1)%3);
        triverts[i][2] = List_Entry(fverts,(j+2)%3);
        fdim[i] = MF_GEntDim(f);
        fid[i] = MF_GEntID(f);
        break;
      }
    }
    List_Delete(fverts);

    MF_Delete(f,0);
  }

  if (etris) List_Delete(etris);

  /* Now split the edge itself */

  vsplit = ME_Split(esplit, splitxyz);

  /* Now for each tri face that we deleted, create two tri faces that
     incorporate the split vertex, one of the split edge vertices and
     opposite vertex */

  for (i = 0; i < ntris; i++) {

    /* First triangle */

    fvarr[0] = triverts[i][0]; 
    fvarr[1] = triverts[i][1];
    fvarr[2] = vsplit;

    f = MF_New(mesh);
    MF_Set_Vertices(f,3,fvarr);
    MF_Set_GEntDim(f,fdim[i]);
    MF_Set_GEntID(f,fid[i]);
    
    /* Second triangle */

    fvarr[0] = triverts[i][0];
    fvarr[1] = vsplit;
    fvarr[2] = triverts[i][2];

    f = MF_New(mesh);
    MF_Set_Vertices(f,3,fvarr);
    MF_Set_GEntDim(f,fdim[i]);
    MF_Set_GEntID(f,fid[i]);
  }
  if (ntris) {
    free(triverts);
    free(fdim);
    free(fid);
  }

  /* Now for each tet that we deleted, create two tets (these will use
     the split faces that are already created */

  for (i = 0; i < ntets; i++) {
    rvarr[0] = vsplit;
    rvarr[1] = tetverts[i][2];
    rvarr[2] = tetverts[i][1];
    rvarr[3] = tetverts[i][0];

    r = MR_New(mesh);
    MR_Set_Vertices(r,4,rvarr,0,NULL);
    MR_Set_GEntID(r,rid[i]);

    rvarr[0] = vsplit;
    rvarr[1] = tetverts[i][1];
    rvarr[2] = tetverts[i][2];
    rvarr[3] = tetverts[i][3];

    r = MR_New(mesh);
    MR_Set_Vertices(r,4,rvarr,0,NULL);
    MR_Set_GEntID(r,rid[i]);
  }

  if (ntets) {
    free(tetverts);
    free(rid);
  }

  return vsplit;
}
