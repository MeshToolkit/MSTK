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
  int *gids, gid, j, num, idx, nent;
  MType mtype;
  MEntity_ptr lment, gment;
  MSet_ptr lmset;
  char msetname[256];

  MSet_Name(gmset,msetname);
  lmset = MESH_MSetByName(submesh,msetname);

  if (!lmset) {
    /* Create the mesh set */
    mtype = MSet_EntDim(gmset);
    lmset = MSet_New(submesh,msetname,mtype);
  }

  nent = MSet_Num_Entries(gmset);
  gids = (int *) MSTK_malloc(nent*sizeof(int));
  j = 0;
  idx = 0;
  while ((gment = MSet_Next_Entry(gmset,&idx)))
    gids[j++] = MEnt_GlobalID(gment);

  if (mtype == MALLTYPE) {

    MSTK_Report("MESH_CopySet","Not implemented for MALLTYPE",ERROR);

  }
  else {
    
    switch (mtype) {
    case MVERTEX:
      idx = 0;
      while ((lment = MESH_Next_Vertex(submesh,&idx))) {
	gid = MV_GlobalID(lment);
	for (j = 0; j < nent; j++)
	  if (gid == gids[j]) {
	    MSet_Add(lmset,lment);
	    break;
	  }
      }
      
      break;
    case MEDGE:
      idx = 0;
      while ((lment = MESH_Next_Edge(submesh,&idx))) {
	gid = MV_GlobalID(lment);
	for (j = 0; j < nent; j++)
	  if (gid == gids[j]) {
	    MSet_Add(lmset,lment);
	    break;
	  }
      }
      break;
    case MFACE:
      idx = 0;
      while ((lment = MESH_Next_Face(submesh,&idx))) {
	gid = MV_GlobalID(lment);
	for (j = 0; j < nent; j++)
	  if (gid == gids[j]) {
	    MSet_Add(lmset,lment);
	    break;
	  }
      }
      break;
    case MREGION:
      idx = 0;
      while ((lment = MESH_Next_Region(submesh,&idx))) {
	gid = MV_GlobalID(lment);
	for (j = 0; j < nent; j++)
	  if (gid == gids[j]) {
	    MSet_Add(lmset,lment);
	    break;
	  }
      }
      break;
    }
  }

  MSTK_free(gids);
  
  return 1;
}
  
#ifdef __cplusplus
}
#endif

