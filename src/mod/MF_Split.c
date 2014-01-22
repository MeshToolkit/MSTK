#include <stdio.h>
#include <stdlib.h>

#include "MSTK.h"

/* Routine to split a face at a point in the face. New triangular
   faces are created using the split vertex and each edge of the
   face. The new faces are incorporated into connected regions. There
   is no check to ensure that the point is actually in the face and
   not on its boundary or outside the face. This is a strictly
   topological routine and geometric checks for validity must be done
   elsewhere */

#ifdef __cplusplus
extern "C" {
#endif

MVertex_ptr MF_Split(MFace_ptr fsplit, double *splitxyz) {
  Mesh_ptr mesh;
  MVertex_ptr vsplit, fv[3];
  MEdge_ptr enew, *fedges0, *fedges1, fe;
  MFace_ptr fnew[MAXPV2];
  MRegion_ptr fr;
  int gid, gdim, i, j, idx, nnew;
  int nfv, fnewdir[MAXPV2], rfdir;
  List_ptr fregs, fverts;

  gdim = MF_GEntDim(fsplit);
  gid = MF_GEntID(fsplit);
 
  /* Collect information */

  mesh = MF_Mesh(fsplit);
  
  /* Regions connected to face */
  
  fregs = MF_Regions(fsplit);

  /* Vertices of face */

  fverts = MF_Vertices(fsplit,1,0);
  nfv = List_Num_Entries(fverts);

  /* Create the splitting vertex */

  vsplit = MV_New(mesh);
  MV_Set_Coords(vsplit,splitxyz);
  MV_Set_GEntDim(vsplit,gdim);
  MV_Set_GEntID(vsplit,gid);

  /* Create the 'nfe' faces */

  for (i = 0; i < nfv; i++) {
    fv[0] = vsplit;
    fv[1] = List_Entry(fverts,i);
    fv[2] = List_Entry(fverts,(i+1)%nfv);
    
    fnew[i] = MF_New(mesh);
    MF_Set_GEntDim(fnew[i],gdim);
    MF_Set_GEntID(fnew[i],gid);
    MF_Set_Vertices(fnew[i],3,fv);
  }
  List_Delete(fverts);
  nnew = nfv;


  for (i = 0; i < List_Num_Entries(fregs); i++) {
    fr = List_Entry(fregs,i);
    rfdir = MR_FaceDir(fr,fsplit);
    for (j = 0; j < nnew; j++)
      fnewdir[j] = rfdir;
    MR_Replace_Faces(fr,1,&fsplit,nnew,fnew,fnewdir);
  }
  List_Delete(fregs);

  MF_Delete(fsplit,0);

  return vsplit;
}

#ifdef __cplusplus
}
#endif
