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
  int ok, idx, gid;
  int max_greg_id,  zeroid;
  double PI=3.141592, ang;
  MRegion_ptr region;

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
