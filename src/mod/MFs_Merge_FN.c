#include <stdlib.h>
#include <stdio.h>
#include "MFace_jmp.h"
#include "MSTK_private.h"

/* Replace f2 with f1 in all the regions using f2 and delete f2 */

MFace_ptr MFs_Merge_FN(MFace_ptr f1, MFace_ptr f2) {
  int i, idx, gdim1, gid1, gdim2, gid2, dir, nfe1, nfe2;
  MRegion_ptr reg;
  MVertex_ptr fv0;
  Mesh_ptr    mesh;
  List_ptr    fregs2, fedges1, fedges2;




  /* Check that we are not violating topological conformity with geometric model */

  mesh = MF_Mesh(f1);
  gid1 = MF_GEntID(f1);
  gdim1 = MF_GEntDim(f1);
  gid2 = MF_GEntID(f2);
  gdim2 = MF_GEntDim(f2);

  if (mesh != MF_Mesh(f2)) {
    MSTK_Report("MFs_Merge","Faces not from same mesh - Cannot merge",MSTK_ERROR);
    return 0;
  }
  else {
    if (gdim1 == gdim2) {
      if (gid1 != gid2) {
	MSTK_Report("MFs_Merge","Faces are on different geometric entities of the same dimension - Cannot merge",MSTK_ERROR);
	return 0;
      }
    }
    else if (gdim2 < gdim1) {
      MSTK_Report("MFs_Merge","Cannot merge face on lower dimensional geometric entity to face on higher dimensional geometric entity",MSTK_ERROR);
      return 0;
    }
  }
  

  /* check to make sure edges of faces are identical */

  fedges1 = MF_Edges(f1,1,0);
  nfe1 = List_Num_Entries(fedges1);
  dir = MF_EdgeDir_i(f1,0);
  fv0 = ME_Vertex(List_Entry(fedges1,0),!dir);
  
  fedges2 = MF_Edges(f2,1,fv0);
  nfe2 = List_Num_Entries(fedges2);
  
  if (nfe1 != nfe2) {
    MSTK_Report("MFs_Merge",
		"To-be-merged faces have different number of edges",
		MSTK_ERROR);
    return 0;
  }
  
  
  if (List_Entry(fedges1,0) == List_Entry(fedges2,0)) {
    dir = 1;

#ifdef DEBUG  

    /* Full edge-by-edge check */

    for (i = 0; i < nfe1; i++) {
      if (List_Entry(fedges1,i) != List_Entry(fedges2,i)) {
	MSTK_Report("MFs_Merge",
		    "Edges of faces must be merged before merging faces",
		    MSTK_ERROR);
	return 0;
      }
    }
#endif 

  }
  else if (List_Entry(fedges1,0) == List_Entry(fedges2,nfe1-1)) {
    dir = 0;
    
#ifdef DEBUG  

    /* Full edge-by-edge check */

    for (i = 0; i < nfe1; i++) {
      if (List_Entry(fedges1,i) != List_Entry(fedges2,nfe1-i-1)) {
	MSTK_Report("MFs_Merge",
		    "Edges of faces must be merged before merging faces",
		    MSTK_ERROR);
	return 0;
      }
    }

#endif 

  }
  else {
    MSTK_Report("MFs_Merge",
		"Edges of faces must be merged before merging faces",
		MSTK_ERROR);
    return 0;
  }

  List_Delete(fedges1);
  List_Delete(fedges2);



  /* Replace f2 with f1 in all regions connected to f2 */

  fregs2 = MF_Regions(f2);

  if (fregs2) {
    idx = 0;
    while ((reg = List_Next_Entry(fregs2,&idx))) {
      int f2dir = MR_FaceDir(reg,f2);
      if (dir)                      /* f1 and f2 oriented the same way */
	MR_Replace_Faces(reg,1,&f2,1,&f1,&f2dir);
      else {                         /* f1 and f2 oriented in opposite way */
        f2dir = !f2dir;
	MR_Replace_Faces(reg,1,&f2,1,&f1,&f2dir);
      }
    }

    List_Delete(fregs2);
  }


  MF_Delete(f2,0);

  return f1;
}
