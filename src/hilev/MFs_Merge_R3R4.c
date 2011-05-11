#include <stdlib.h>
#include <stdio.h>
#include "MFace_jmp.h"
#include "MSTK_private.h"

/* Replace f2 with f1 in all the regions using f2 and delete f2 */

MFace_ptr MFs_Merge_R3R4(MFace_ptr f1, MFace_ptr f2) {
  int i, idx, gdim, gid, dir;
  MRegion_ptr reg;
  Mesh_ptr    mesh;
  List_ptr    fregs2, fverts1, fverts2;


  /* Check that we are not violating topological conformity with geometric model */

  mesh = MF_Mesh(f1);
  gid = MF_GEntID(f1);
  gdim = MF_GEntDim(f1);

  if (mesh != MF_Mesh(f2)) {
    MSTK_Report("MFs_Merge","Faces not from same mesh - Cannot merge",ERROR);
    return 0;
  }
  else if (gid != MF_GEntID(f2) || gdim != MF_GEntDim(f2)) {
    MSTK_Report("MFs_Merge","Faces not from same geometric entity - Cannot merge",ERROR);
    return 0;
  }
  

  /* Debug check to make sure edges of faces are identical */

#ifdef DEBUG
  {
    List_ptr fverts1, fverts2;
    MVertex_ptr fv1;

    fverts1 = MF_Vertices(f1,1,0);
    fverts2 = MF_Vertices(f2,1,0);
    
    if (List_Num_Entries(f1) != List_Num_Entries(f2)) {
      MSTK_Report("MFs_Merge",
		  "To-be-merged faces have different number of edges",
		  ERROR);
      return 0;
    }

    idx = 0;
    while ((fv1 = List_Next_Entry(fverts1,&idx))) {
      if (!List_Contains(fverts2,fv2)) {
	MSTK_Report("MFs_Merge",
		    "Vertices of faces must be merged before merging faces",
		    ERROR);
	return 0;
      }
    }

    List_Delete(fverts1);
    List_Delete(fverts2);
  }
#endif


  /* Replace f2 with f1 in all regions connected to f2 */

  fregs2 = MF_Regions(f2);

  if (fregs2) {
    idx = 0;
    while ((reg = List_Next_Entry(fregs2,&idx))) {
      dir = MR_FaceDir(reg,f2);
      MR_Replace_Face(reg,f2,f1,dir);
    }

    List_Delete(fregs2);
  }


  MF_Delete(f2,0);

  return f1;
}
