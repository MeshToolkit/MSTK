#define _H_MEdge_Private

#include "MEdge.h"
#include "MEdge_jmp.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  int ME_Num_Faces_R1R2(MEdge_ptr e) {
    int nf;
    List_ptr efaces = ME_Faces(e);

    if (efaces) {
      nf = List_Num_Entries(efaces);
      List_Delete(efaces);
    }
    else
      nf = 0;

    return nf;
  }

  List_ptr ME_Faces_R1R2(MEdge_ptr e) {
    int idx, idx1, idx2, found;
    MFace_ptr face, lstface;
    MRegion_ptr reg;
    List_ptr rfaces, vfaces1, efaces, vregs0, vregs1, cmnregs=NULL;

    efaces = List_New(0);

    vregs0 = MV_Regions(e->vertex[0]);
    vregs1 = MV_Regions(e->vertex[1]);
    if (vregs0 || vregs1) {
      if (vregs0 && vregs1) {
	cmnregs = List_New(0);
	
	idx = 0;
	while ((reg = List_Next_Entry(vregs0,&idx))) {
	  if (List_Contains(vregs1,reg))
	    List_Add(cmnregs,reg);
	}
      }
      if (vregs0)
	List_Delete(vregs0);
      if (vregs1)
	List_Delete(vregs1);
      
      if (!List_Num_Entries(cmnregs)) {
	List_Delete(cmnregs);
	return NULL;
      }

      idx = 0;
      while ((reg = List_Next_Entry(cmnregs,&idx))) {
	rfaces = MR_Faces(reg);

	idx1 = 0;
	while ((face = List_Next_Entry(rfaces,&idx1))) {
	  if (MF_UsesEntity(face,(MEntity_ptr) e,MEDGE)) {
	    
	    idx2 = 0; found = 0;
	    while ((lstface = List_Next_Entry(efaces,&idx2))) {
	      if (MFs_AreSame(face,lstface)) {
		found = 1;
		break;
	      }
	    }

	    if (!found)
	      List_Add(efaces,face);
	  }
	}
	
	List_Delete(rfaces);
      }
      List_Delete(cmnregs);
    }
    else { /* Must be only faces are connected to edge */

      vfaces1 = MV_Faces(e->vertex[0]);

      if (vfaces1) {
	idx = 0;
	while ((face = List_Next_Entry(vfaces1,&idx))) {
	  if (MF_UsesEntity(face,e->vertex[1],MVERTEX))
	    List_Add(efaces,face);
	}
      }
      
    }

    if (List_Num_Entries(efaces))
      return efaces;
    else {
      List_Delete(efaces);
      return NULL;
    }
  }

  void ME_FaceIDs_R1R2(MEdge_ptr e, int *nef, int *efaceids) {
    List_ptr efaces = ME_Faces_R1R2(e);
    if (efaces) {
      int i;
      *nef = List_Num_Entries(efaces);
      for (i = 0; i < *nef; i++) 
        efaceids[i] = MEnt_ID(List_Entry(efaces,i));
    }
    else
      *nef = 0;
  }

#ifdef __cplusplus
}
#endif

