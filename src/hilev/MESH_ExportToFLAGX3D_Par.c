#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif



  /* Function to export a partitioned MSTK mesh to the Parallel FLAG
     X3D format (LA-UR-04-9033) */

  /* if natt = 0, all attributes are written out if natt = -1, no
     attributes are written out if natt > 0, only attributes specified
     in attnames are written out Except for cell material IDs,
     attributes are ignored by FLAG; still we will write them out for
     completeness */

  /* opts is an array of flags that controls how mesh is exported to FLAG 
     Currently, it is a dummy argument
  */


int MESH_ExportToFLAGX3D_Par(Mesh_ptr mesh, const char *filename,  
			     const int nparts, const int natt, 
			     const char **attnames, int *opts, int *procids) {
  int nv, ne, nf, nr;

  if (nparts < 2) {
    MESH_ExportToFLAGX3D(mesh,filename,natt,attnames,opts);
    return 1;
  }

  nv = MESH_Num_Vertices(mesh);
  ne = MESH_Num_Edges(mesh);
  nf = MESH_Num_Faces(mesh);
  nr = MESH_Num_Regions(mesh);
  if (!nr && !nf) {
    fprintf(stderr,
	    "Format does not support meshes with only edges or nodes\n");
    return 0;
  }


  if (nr)
    MESH_Vol_ExportToFLAGX3D_Par(mesh,filename,nparts,natt,attnames,opts,procids);
  else
    MESH_Surf_ExportToFLAGX3D_Par(mesh,filename,nparts,natt,attnames,opts,procids);

  return 1;
}

#ifdef __cplusplus
  }
#endif

