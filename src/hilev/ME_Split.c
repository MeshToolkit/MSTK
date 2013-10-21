#include <stdio.h>
#include <stdlib.h>
#include "MSTK.h"

/* Routine to split an edge by inserting a vertex - two new edges are
   created in place of this edge and the faces connected to the old
   edge are updated to reflect the new topology. The split edge is
   deleted */

int ME_Split(MEdge_ptr esplit, MVertex_ptr vnew, MEdge_ptr enew[2]) {
  Mesh_ptr mesh;
  MEdge_ptr fe[4], fenew[6], enewrev[2];
  MVertex_ptr v[4];
  MFace_ptr f[2], ef;
  int dir[3], gid, gdim, k, nf, idx, fedir;
  List_ptr fedges, efaces;

  gdim = ME_GEntDim(esplit);
  gid = ME_GEntID(esplit);
 
  /* Collect information */

  mesh = ME_Mesh(esplit);
  v[0] = ME_Vertex(esplit,0);
  v[1] = ME_Vertex(esplit,1);

  enew[0] = ME_New(mesh);
  ME_Set_Vertex(enew[0],0,v[0]);
  ME_Set_Vertex(enew[0],1,vnew);
  ME_Set_GEntDim(enew[0],gdim);
  ME_Set_GEntID(enew[0],gid);

  enew[1] = ME_New(mesh);
  ME_Set_Vertex(enew[1],0,vnew);
  ME_Set_Vertex(enew[1],1,v[1]);
  ME_Set_GEntDim(enew[1],gdim);
  ME_Set_GEntID(enew[1],gid);

  enewrev[0] = enew[1];
  enewrev[1] = enew[0];

  efaces = ME_Faces(esplit);
  idx = 0;
  while ((ef = List_Next_Entry(efaces,&idx))) {
    fedir = MF_EdgeDir(ef,esplit);
    if (fedir) 
      MF_Replace_Edges(ef,1,&esplit,2,enew);
    else
      MF_Replace_Edges(ef,1,&esplit,2,enewrev);
  }

  return 1;
}
