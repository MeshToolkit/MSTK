#define _H_MFace_Private

#include "MFace.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  List_ptr MF_Regions_R1(MFace_ptr f) {
    MFace_Adj_R1 *adj = (MFace_Adj_R1 *) f->adj;
    MVertex_ptr v, rv;
    MRegion_ptr r;
    List_ptr vregs, rverts, fregs;
    int nfv, nr, idx, i, fnd=0, idx1;
    
    v = List_Entry(adj->fvertices,0);
    vregs = MV_Regions(v);
    if (!vregs)
      return NULL;

    nfv = List_Num_Entries(adj->fvertices);
    nr = 0;
    fregs = List_New(2);

    idx = 0;
    while (nr < 2 && (r = List_Next_Entry(vregs,&idx))) {
      rverts = MR_Vertices(r);

      for (i = 1; i < nfv; i++) {
	v = List_Entry(adj->fvertices,i);
	idx1 = 0;
	fnd = 0;
	while ((rv = List_Next_Entry(rverts,&idx1))) {
	  if (rv == v) {
	    fnd = 1;
	    break;
	  }
	}
	if (!fnd)
	  break;
      }
      
      List_Delete(rverts);

      if (fnd) {
	List_Add(fregs,r);
	nr++;
      }
    }
    List_Delete(vregs);
    if (nr)
      return fregs;
    else {
      List_Delete(fregs);
      return 0;
    }      
  }

  void MF_RegionIDs_R1(MFace_ptr f, int *nfr, int *fregids) {
    int i;
    List_ptr fregs = MF_Regions_R1(f);
    *nfr = List_Num_Entries(fregs);
    for (i = 0; i < *nfr; i++)
      fregids[i] = MEnt_ID(List_Entry(fregs,i));
    List_Delete(fregs);
  }


  MRegion_ptr MF_Region_R1(MFace_ptr f, int dir) {
    List_ptr fregs;
    int nr;
    MRegion_ptr reg, ret = NULL;

#ifdef DEBUG
    MSTK_Report("MF_Region_R1",
		"May be more efficient to call MF_Regions",MSTK_WARN);
#endif

    fregs = MF_Regions_R1(f);
    if (!fregs)
      return NULL;

    nr = List_Num_Entries(fregs);

    reg = List_Entry(fregs,0);
    if (MR_FaceDir(reg,f) == !dir) {
      ret = reg;
    }
    else if (nr == 2) {
      reg = List_Entry(fregs,1);
      if (MR_FaceDir(reg,f) == !dir) {
	ret = reg;
      }
    }

    List_Delete(fregs);

    return ret;
  }

  int MF_RegionID_R1(MFace_ptr f, int ir) {
    MRegion_ptr r = MF_Region_R1(f,ir);
    return r ? MEnt_ID(r) : 0;
  }


#ifdef __cplusplus
}
#endif
