#define _H_MVertex_Private

#include "MVertex.h"
#include "MVertex_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Functions */

  int MV_Num_Faces_R1R2(MVertex_ptr v) {
    int i, nf;
    List_ptr vfaces;

#ifdef DEBUG
    MSTK_Report("MV_Num_Faces_R1R2",
		"Inefficient to call this routine with this representation",
		MESG);
#endif
    vfaces = MV_Faces_R1R2(v);
    nf = List_Num_Entries(vfaces);
    for (i = 0; i < nf; i++) {
      /* Must destroy the temporary faces */
    }
    List_Delete(vfaces);

    return nf;
  }


  List_ptr MV_Faces_R1R2(MVertex_ptr v) {
    MSTK_Report("MV_Faces","Not yet implemented for this representation",
		MESG);
    return NULL;
  }

#ifdef __cplusplus
}
#endif
