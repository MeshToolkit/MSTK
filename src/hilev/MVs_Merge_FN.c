#include <stdlib.h>
#include <stdio.h>
#include "MVertex_jmp.h"
#include "MSTK_private.h"


MVertex_ptr MVs_Merge_FN(MVertex_ptr v1, MVertex_ptr v2) {
  int idx, gdim, gid;
  MEdge_ptr   edge;
  Mesh_ptr    mesh;
  List_ptr    vedges2;

  mesh = MV_Mesh(v1);
  gid = MV_GEntID(v1);
  gdim = MV_GEntDim(v1);

  if (mesh != MF_Mesh(v2)) {
    MSTK_Report("MVs_Join","Vertices not from same mesh",ERROR);
    return 0;
  }
  else if (gid != MV_GEntID(v2) || gdim != MV_GEntDim(v2)) {
    MSTK_Report("MFs_Join","Faces not from same geometric entity",ERROR);
    return 0;
  }

  vedges2 = MV_Edges(v2);

  if (vedges2) {
    idx = 0;
    while ((edge = List_Next_Entry(vedges2,&idx)))
      ME_Replace_Vertex(edge,v2,v1);

    List_Delete(vedges2);
  }

  MV_Delete(v2,0);

  return v1;
}
