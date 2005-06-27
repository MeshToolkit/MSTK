#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MSTK.h"
#include "MSTK_private.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Routine to build classification for mesh entities given a mesh no
   information or model ID info only for the highest level entities
   (typically regions in volume meshes and faces in surface meshes) */


int MESH_BuildClassfn(Mesh_ptr mesh) {
  int i, k, ok, idx, gid, gdim;
  int grid0, grid1, gfid0, gfid1, geid0, geid1;
  int max_greg_id, max_gface_id, max_gedge_id, max_gvertex_id, zeroid;
  int nfr, nef, nbf, nve, nbe, nbe2, gid2, gdim2, *gfids, (*gfregids)[2];
  double PI=3.141592, ang, angr;
  MVertex_ptr vertex;
  MEdge_ptr edge;
  MFace_ptr face;
  MRegion_ptr region, freg0, freg1;
  List_ptr vedges, efaces, GFfaces, GEedges, GFedges, fregs;

  ang = 5*PI/6.0;  /* 75 degrees */

  /* MESH REGIONS >>>>>>>>>>>>>>> */

  /* Verify that mesh regions have classification information; if
     not, assign all regions to the same model region */

  zeroid = 0;
  max_greg_id = 0;
  idx = 0;
  while ((region = MESH_Next_Region(mesh,&idx))) {
    gid = MR_GEntID(region);
    if (gid) {
      if (gid > max_greg_id)
	max_greg_id = gid;
    }
    else
      zeroid = 1;
  }

  if (zeroid) {
    idx = 0; 
    while ((region = MESH_Next_Region(mesh,&idx))) {
      /* Make sure region is not classified on zero or -ve ID model region */

      if (MR_GEntID(region) <= 0)
	MR_Set_GEntID(region,(max_greg_id+1));
    }
  }

  /* <<<<<<<<<<<<<<<< end MESH REGIONS */



  /* MESH FACES >>>>>>>>>>>>>>> */
    
  ok = MESH_BuildFaceClassfn(mesh);


  /* <<<<<<<<<<<<<<< MESH FACES */


  /* MESH EDGES >>>>>>>>>>>>>>> */
    
  ok = MESH_BuildEdgeClassfn(mesh);


  /* <<<<<<<<<<<<<<< MESH EDGES */


  /* MESH EDGES >>>>>>>>>>>>>>> */
    
  ok = MESH_BuildVertexClassfn(mesh);


  /* <<<<<<<<<<<<<<< MESH EDGES */


  return 1;
}


#ifdef __cplusplus
}
#endif
