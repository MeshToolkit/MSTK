#define _H_MVertex_Private

#include "MVertex.h"
#include "MVertex_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  int MV_Num_Regions_R1R2(MVertex_ptr v) {
    MVertex_UpAdj_R1R2 *upadj;
    int idx, nr = 0;
    MEntity_ptr ent;

    upadj = (MVertex_UpAdj_R1R2 *) v->upadj;
    idx = 0;
    while ((ent = (MEntity_ptr) List_Next_Entry(upadj->velements,&idx))) {
      if (MEnt_Dim(ent) == MREGION)
	nr++;
    }
    return nr;
  }

  List_ptr MV_Regions_R1R2(MVertex_ptr v) {
    MVertex_UpAdj_R1R2 *upadj;
    int idx, nel, nr = 0, dim;
    MEntity_ptr ent;
    List_ptr vregions;

    upadj = (MVertex_UpAdj_R1R2*) v->upadj;
    nel = List_Num_Entries(upadj->velements);
    vregions = List_New(nel);

    idx = 0;
    while ((ent = (MEntity_ptr) List_Next_Entry(upadj->velements,&idx))) {
      if (MEnt_Dim(ent) == MREGION) {
	List_Add(vregions,ent);
	nr++;
      }
    }
    if (nr)
      return vregions;
    else {
      List_Delete(vregions);
      return 0;
    }      
  }

  void MV_Add_Region_R1R2(MVertex_ptr v, MRegion_ptr mregion) {
    MVertex_UpAdj_R1R2 *upadj;

    upadj = (MVertex_UpAdj_R1R2 *) v->upadj;
    List_Add(upadj->velements,mregion);
  }

  void MV_Rem_Region_R1R2(MVertex_ptr v, MRegion_ptr mregion) {
   MVertex_UpAdj_R1R2 *upadj;

    upadj = (MVertex_UpAdj_R1R2 *) v->upadj;
    List_Rem(upadj->velements,mregion);
  }


#ifdef __cplusplus
}
#endif
