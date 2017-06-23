#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MSTK.h"
#include "MSTK_private.h"

/* Routine to split an edge at 'n' given locations - n+1 edges are
   created in place of this edge and the faces connected to the old
   edge are updated to reflect the new topology. The split edge is
   deleted. The list returns the new edges */

#ifdef __cplusplus
extern "C" {
#endif

List_ptr ME_MultiSplit(MEdge_ptr esplit, int n, double (*xyz)[3]) {
  Mesh_ptr mesh;
  MEdge_ptr *enew, *enewrev;
  MVertex_ptr *verts;
  int i, j;
  int gid, gdim;
  List_ptr efaces, eregions, enewlist;
  double (*vxyz)[3];

  gdim = ME_GEntDim(esplit);
  gid = ME_GEntID(esplit);

  efaces = ME_Faces(esplit);
  eregions = ME_Regions(esplit);

  verts = (MVertex_ptr *) malloc((n+2)*sizeof(MVertex_ptr));
  vxyz = (double (*)[3]) malloc((n+2)*sizeof(double [3]));

  /* Collect information */

  mesh = ME_Mesh(esplit);
  verts[0] = ME_Vertex(esplit,0);
  verts[n+1] = ME_Vertex(esplit,1);
  
  MV_Coords(verts[0],vxyz[0]);
  MV_Coords(verts[n+1],vxyz[n+1]);

  if (xyz) { /* user specified split coordinates */
    for (i = 1; i <= n; i++)
      memcpy(vxyz[i],xyz[i-1],sizeof(double[3]));
  }
  else { /* equispaced split coordinates */
    for (i = 1; i <= n; i++)
      for (j = 0; j < 3; j++) 
        vxyz[i][j] = vxyz[0][j] + i*(vxyz[n+1][j]-vxyz[0][j])/(n+1);
  }

  for (i = 0; i < n; i++) {
    verts[i+1] = MV_New(mesh);    
    MV_Set_Coords(verts[i+1],vxyz[i+1]);
    MV_Set_GEntDim(verts[i+1],gdim);
    MV_Set_GEntID(verts[i+1],gid);
  }

  enew = (MEdge_ptr *) malloc((n+1)*sizeof(MEdge_ptr));
  enewrev = (MEdge_ptr *) malloc((n+1)*sizeof(MEdge_ptr));

  for (i = 0; i <= n; i++) {
    enew[i] = ME_New(mesh);
    ME_Set_Vertex(enew[i],0,verts[i]);
    ME_Set_Vertex(enew[i],1,verts[i+1]);
    ME_Set_GEntDim(enew[i],gdim);
    ME_Set_GEntID(enew[i],gid);

    enewrev[n-i] = enew[i];
  }

  if (efaces) {
    MFace_ptr ef;
    int idx = 0;
    while ((ef = List_Next_Entry(efaces,&idx))) {
      int fedir = MF_EdgeDir(ef,esplit);
      if (fedir) 
        MF_Replace_Edges(ef,1,&esplit,n+1,enew);
      else
        MF_Replace_Edges(ef,1,&esplit,n+1,enewrev);
    }
    List_Delete(efaces);
  }

  if (eregions) {
    MRegion_ptr er;
    int idx = 0;
    while ((er = List_Next_Entry(eregions,&idx))) 
      MR_Update_ElementType(er);
    List_Delete(eregions);
  }

  ME_Delete(esplit,0);

  enewlist = List_New(n+1);
  for (i = 0; i < n+1; i++) 
    List_Add(enewlist,enew[i]);

  free(verts);
  free(vxyz);
  free(enew);
  free(enewrev);

  return enewlist;
}

#ifdef __cplusplus
}
#endif
