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

int MESH_CopySet(Mesh_ptr mesh, Mesh_ptr submesh, MSet_ptr gmset) {
  int idx, idx2, found;
  MType mtype;
  MEntity_ptr lment, gment;
  MSet_ptr lmset;
  char msetname[256];
  void *loc;
  MAttrib_ptr g2latt;
  List_ptr lmentlist;

  g2latt = MESH_AttribByName(mesh,"Global2Local");

  MSet_Name(gmset,msetname);
  lmset = MESH_MSetByName(submesh,msetname);

  if (!lmset) {
    /* Create the mesh set */
    mtype = MSet_EntDim(gmset);
    lmset = MSet_New(submesh,msetname,mtype);
  }

  idx = 0;
  while ((gment = MSet_Next_Entry(gmset,&idx))) {
    MEnt_Get_AttVal(gment,g2latt,0,0,&lmentlist);
    
    idx2 = 0;
    found = 0;
    while ((lment = List_Next_Entry(lmentlist,&idx2))) {
      if (MEnt_Mesh(lment) == submesh) {
	found = 1;
	MSet_Add(lmset,lment);
        break;
      }
    }
  }

  return 1;
}
  
#ifdef __cplusplus
}
#endif

