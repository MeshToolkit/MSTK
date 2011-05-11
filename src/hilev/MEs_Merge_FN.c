#include <stdlib.h>
#include <stdio.h>
#include "MEdge_jmp.h"
#include "MSTK_private.h"

/* Replace e2 with e1 in all the faces using e2 and delete e2 */

MEdge_ptr MEs_Merge_FN(MEdge_ptr e1, MEdge_ptr e2) {
  int i, idx, gdim1, gid1, gdim2, gid2;
  MVertex_ptr v11, v12, v21, v22;
  MFace_ptr   face;
  Mesh_ptr    mesh;
  List_ptr    efaces2;

  mesh = ME_Mesh(e1);
  gid1 = ME_GEntID(e1);
  gdim1 = ME_GEntDim(e1);
  gid2 = ME_GEntID(e2);
  gdim2 = ME_GEntDim(e2);

  if (mesh != MF_Mesh(e2)) {
    MSTK_Report("MEs_Merge","Edges not from same mesh - Cannot merge",ERROR);
    return 0;
  }
  else {
    if (gdim1 == gdim2) {
      if (gid1 != gid2) {
	MSTK_Report("MEs_Merge","Edges are on different geometric entities of the same dimension - Cannot merge",ERROR);
	return 0;
      }
    }
    else if (gdim2 < gdim1) {
      MSTK_Report("MEs_Merge","Cannot merge edge on lower dimensional geometric entity to edge on higher dimensional geometric entity",ERROR);
      return 0;
    }
  }
  
  v11 = ME_Vertex(e1,0); v12 = ME_Vertex(e1,1);
  v21 = ME_Vertex(e2,0); v22 = ME_Vertex(e2,1);

  if ((v11 != v21 && v11 != v22) ||
      (v12 != v21 && v12 != v22)) {
    MSTK_Report("MEs_Merge","Vertices of edges must be merged before attempting to merge edges",ERROR);
    return 0;
  }


  efaces2 = ME_Faces(e2);

  if (efaces2) {
    idx = 0;
    while ((face = List_Next_Entry(efaces2,&idx)))
      MF_Replace_Edges(face,1,&e2,1,&e1);

    List_Delete(efaces2);
  }

  ME_Delete(e2,0);

  return e1;
}
