#define _H_MFace_Private

#include "MFace.h"
#include "MSTK_malloc.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  int MFs_AreSame_R1R2(MFace_ptr f1, MFace_ptr f2) {
    
    if (f1 == f2)
      return 1;
    else {
      int nfv1, nfv2, i, j, k, dir=-1;
      MVertex_ptr fv;
      List_ptr fverts1, fverts2;
      
      fverts1 = MF_Vertices(f1,1,0); nfv1 = List_Num_Entries(fverts1);
      fverts2 = MF_Vertices(f2,1,0); nfv2 = List_Num_Entries(fverts2);
      
      if (nfv1 != nfv2) {
	List_Delete(fverts1);
	List_Delete(fverts2);
	return 0;
      }
      
      fv = List_Entry(fverts1,0);
      k = List_Locate(fverts2,fv);
      if (k == -1) { /* Could not be found */
	List_Delete(fverts1);
	List_Delete(fverts2);
	return 0;
      }
      
      fv = List_Entry(fverts1,1);
      if (fv == List_Entry(fverts2,(k+1)%nfv1))
	dir = 1;
      else if (fv == List_Entry(fverts2,(k-1+nfv1)%nfv1))
	dir = 0;
      else { /* Could not be found */
	List_Delete(fverts1);
	List_Delete(fverts2);
	return 0;
      }
      
      if (dir) {
	for (j = 2; j < nfv1; j++) {
	  fv = List_Entry(fverts1,j);
	  if (fv != List_Entry(fverts2,(k+j)%nfv1)) {
	    List_Delete(fverts1); 
	    List_Delete(fverts2);
	    return 0;
	  }
	}
      }
      else {
	for (j = 2; j < nfv1; j++) {
	  fv = List_Entry(fverts1,j);
	  if (fv != List_Entry(fverts2,(k-j+nfv1)%nfv1)) {
	    List_Delete(fverts1); 
	    List_Delete(fverts2);
	    return 0;
	  }
	}
      }
      List_Delete(fverts1);
      List_Delete(fverts2);    

      return 1;
    }
  }


#ifdef __cplusplus
}
#endif
