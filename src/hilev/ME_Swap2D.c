#include <stdio.h>
#include <stdlib.h>
#include "MSTK.h"


int ME_Swap2D(MEdge_ptr eswp, MEdge_ptr *enew, MFace_ptr fnew[2]) {
  Mesh_ptr mesh;
  MEdge_ptr fe[4], fenew[6];
  MVertex_ptr v[4];
  MFace_ptr f[2];
  int dir[3], gid, k, nf;
  Set_ptr fedges, efaces;

  if (ME_GEntDim(eswp) != 2) {
    MSTK_Report("ME_Swap2D","Can swap edges on model faces only",ERROR);
    return 0;
  } 
  /*  gid = ME_GEntID(eswp); */
  gid = 0;
 
  /* Collect information */

  mesh = ME_Mesh(eswp);
  v[0] = ME_Vertex(eswp,0);
  v[1] = ME_Vertex(eswp,1);

  efaces = ME_Faces(eswp);
  nf = Set_Num_Entries(efaces);

  if (nf != 2) {
    MSTK_Report("ME_Swap2D","Edge must be connected to exactly two faces for 2D swap",ERROR);
    Set_Delete(efaces);
    return 0;
  }


  k = 0;
  f[0] = Set_Entry(efaces,0);
  if (!MF_EdgeDir(f[0],eswp)) {
    k = 1;
    f[0] = Set_Entry(efaces,1);
  }

  fedges = MF_Edges(f[0],1,v[0]);
  if (Set_Num_Entries(fedges) != 3) {
    MSTK_Report("ME_Swap2D","Cannot swap edge for non-triangular meshes",ERROR);
    Set_Delete(efaces);
    Set_Delete(fedges);
    return 0;
  }
  fe[0] = Set_Entry(fedges,1);
  fe[1] = Set_Entry(fedges,2);
  Set_Delete(fedges);

  v[3] = ME_Vertex(fe[0],1);
  if (v[3] == v[1])
    v[3] = ME_Vertex(fe[0],0);

  f[1] = Set_Entry(efaces,!k);
  fedges = MF_Edges(f[1],1,v[1]);
  if (Set_Num_Entries(fedges) != 3) {
    MSTK_Report("ME_Swap2D","Cannot swap edge for non-triangular meshes",ERROR);
    Set_Delete(efaces);
    Set_Delete(fedges);
    return 0;
  }
  fe[2] = Set_Entry(fedges,1);
  fe[3] = Set_Entry(fedges,2);
  Set_Delete(fedges);

  v[2] = ME_Vertex(fe[3],0);
  if (v[2] == v[1])
    v[2] = ME_Vertex(fe[3],1);


  Set_Delete(efaces);


  /* Delete old configuration */

  MF_Delete(f[0]);
  MF_Delete(f[1]);
  ME_Delete(eswp);

  /* Create new configuration */
  *enew = ME_New(mesh);
  ME_Set_Vertex(*enew,0,v[2]);
  ME_Set_Vertex(*enew,1,v[3]);
  ME_Set_GEntDim(*enew,2);
  ME_Set_GEntID(*enew,gid);

  fnew[0] = MF_New(mesh);
  MF_Set_GEntDim(fnew[0],2);
  MF_Set_GEntID(fnew[0],gid);
  fenew[0] = *enew; dir[0] = 1;
  fenew[1] = fe[1];  dir[1] = (ME_Vertex(fe[1],0) == v[3]) ? 1 : 0;
  fenew[2] = fe[2];  dir[2] = (ME_Vertex(fe[2],0) == v[0]) ? 1 : 0;
  MF_Set_Edges(fnew[0],3,fenew,dir);

  fnew[1] = MF_New(mesh);
  MF_Set_GEntDim(fnew[1],2);
  MF_Set_GEntID(fnew[1],gid);
  fenew[0] = *enew; dir[0] = 0;
  fenew[1] = fe[3];  dir[1] = (ME_Vertex(fe[3],0) == v[2]) ? 1 : 0;
  fenew[2] = fe[0];  dir[2] = (ME_Vertex(fe[0],0) == v[1]) ? 1 : 0;
  MF_Set_Edges(fnew[1],3,fenew,dir);

  return 1;
}
