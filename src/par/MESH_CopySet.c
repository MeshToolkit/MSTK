/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "MSTK.h"

#ifdef __cplusplus
extern "C" {
#endif

 /* 
    this function copy set information from mesh to submesh
    mset_global - Set in the global_mesh
    Create the set in the submesh if it is not there

 */

  int MESH_CopySet(Mesh_ptr mesh, int num, Mesh_ptr *submeshes, MSet_ptr gmset) {
  int i, idx, idx2, found;
  MType mtype;
  MEntity_ptr lment, gment;
  MSet_ptr lmset;
  char msetname[256];
  void *loc;
  MAttrib_ptr g2latt;
  List_ptr lmentlist;
  MSet_ptr *lmset_array;
  Mesh_ptr submesh;

  lmset_array = (MSet_ptr *) calloc(num,sizeof(MSet_ptr));

  g2latt = MESH_AttribByName(mesh,"Global2Local");

  MSet_Name(gmset,msetname);
  mtype = MSet_EntDim(gmset);

  for (i = 0; i < num; ++i)
    lmset_array[i] = MESH_MSetByName(submeshes[i],msetname);
 
  

  idx = 0;
  while ((gment = MSet_Next_Entry(gmset,&idx))) {
    MEnt_Get_AttVal(gment,g2latt,0,0,&lmentlist);
    
    idx2 = 0;
    while ((lment = List_Next_Entry(lmentlist,&idx2))) {
      submesh = MEnt_Mesh(lment);

      for (i = 0; i < num; ++i) {
	if (submesh == submeshes[i]) {
	  lmset = lmset_array[i];
	  if (!lmset)
	    lmset = lmset_array[i] = MSet_New(submesh,msetname,mtype);
      
	  MSet_Add(lmset,lment);
	  break;
	}
      }

    }
  }

  free(lmset_array);
  return 1;
}
  
#ifdef __cplusplus
}
#endif

