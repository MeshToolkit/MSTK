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
    MFace_ptr vface;

#ifdef DEBUG
    MSTK_Report("MV_Num_Faces_R1R2",
		"Inefficient to call this routine with this representation",
		MESG);
#endif
    vfaces = MV_Faces_R1R2(v);
    nf = List_Num_Entries(vfaces);
    for (i = 0; i < nf; i++) {
      vface = List_Entry(vfaces,i);
      if (MEnt_IsVolatile(vface))
	MF_Delete(vface,0);
    }
    List_Delete(vfaces);

    return nf;
  }


  List_ptr MV_Faces_R1R2(MVertex_ptr v) {
    MVertex_UpAdj_R1R2 *upadj;
    int idx, idx1, idx2, found;
    MEntity_ptr ent;
    MFace_ptr rface, lstface;
    List_ptr rfaces, vfaces;

    vfaces = List_New(0);

    upadj = (MVertex_UpAdj_R1R2 *) v->upadj;
    idx = 0;
    while ((ent = (MEntity_ptr) List_Next_Entry(upadj->velements,&idx))) {

      if (MEnt_Dim(ent) == MREGION) {

	rfaces = MR_Faces(ent);

	idx1 = 0;
	while ((rface = List_Next_Entry(rfaces,&idx1))) {
	  if (MF_UsesEntity(rface,v,MVERTEX)) {
	    
	    idx2 = 0; found = 0;
	    while ((lstface = List_Next_Entry(vfaces,&idx2))) {
	      if (MFs_AreSame(rface,lstface)) {
		found = 1;
		break;
	      }
	    }

	    if (!found)
	      List_Add(vfaces,rface);
	  }
	}
	
	List_Delete(rfaces);
      }
      else { /* Must be a face */
	List_Add(vfaces,ent);
      }
    }

    if (List_Num_Entries(vfaces))
      return vfaces;
    else {
      List_Delete(vfaces);
      return NULL;
    }
  }

  void MV_Add_Face_R1R2(MVertex_ptr v, MFace_ptr mface) {
    MVertex_UpAdj_R1R2 *upadj;

    if (MF_Region(mface,0) || MF_Region(mface,1)) {
      MSTK_Report("MV_Add_Face_R1",
		 "Can only add faces with no regions in this representation",
		 ERROR);
      return;
    }

    upadj = (MVertex_UpAdj_R1R2 *)v->upadj;
    List_Add(upadj->velements,mface);
  }

  void MV_Rem_Face_R1R2(MVertex_ptr v, MFace_ptr mface) {
    MVertex_UpAdj_R1R2 *upadj;

    if (MF_Region(mface,0) || MF_Region(mface,1)) {
      MSTK_Report("MV_Rem_Face_R1",
      "Set should contain only faces with no regions in this representation",
		 ERROR);
      return;
    }

    upadj = (MVertex_UpAdj_R1R2 *) v->upadj;
    List_Rem(upadj->velements,mface);
   }

#ifdef __cplusplus
}
#endif
