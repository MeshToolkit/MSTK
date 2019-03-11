/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Update attributes on partitions 

     Authors: Duo Wang
              Rao Garimella
  */

  int MESH_UpdateAttributes(Mesh_ptr mesh, MSTK_Comm comm) {
    int i, natt, attr_type;
    char attr_name[256];
    MAttrib_ptr attrib;

    natt = MESH_Num_Attribs(mesh);
    for(i = 0; i < natt; i++) {
      attrib = MESH_Attrib(mesh,i);
      attr_type = MAttrib_Get_Type(attrib);
      if (attr_type == POINTER) continue; 

      MAttrib_Get_Name(attrib,attr_name);
      MPI_Barrier(comm);
      MESH_Update1Attribute(mesh,attrib,comm);
    }

    return 1;
  }

  /* Keep this for backward compatibility */

  int MSTK_UpdateAttr(Mesh_ptr mesh, MSTK_Comm comm) {

#ifdef DEBUG
     MSTK_Report("MSTK_UpdateAttr","MSTK_UpdateAttr is deprecated. Use MESH_UpdateAttributes",MSTK_WARN);
#endif

     return MESH_UpdateAttributes(mesh,comm);     
  }


#ifdef __cplusplus
}
#endif

