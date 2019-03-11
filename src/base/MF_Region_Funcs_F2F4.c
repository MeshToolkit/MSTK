/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#define _H_MFace_Private

#include "MFace.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif


  List_ptr MF_Regions_F2F4(MFace_ptr f) {
    List_ptr fregs, eregs;
    MRegion_ptr r;
    MEdge_ptr e;
    MFace_Adj_F2F4 *adj;
    int i, k=0, nr;
    
    adj = (MFace_Adj_F2F4 *) f->adj;
    
    fregs = List_New(2);
    e = List_Entry(adj->fedges,0);
    
    eregs = ME_Regions(e);
    if (eregs) {
      nr = List_Num_Entries(eregs);
      
      for (i = 0; i < nr; i++) {
	r = List_Entry(eregs,i);
	if (MR_UsesEntity(r,(MEntity_ptr) f,2)) {
	  List_Add(fregs,r);
	  k++;
	}
	if (k == 2)
	  break;
      }
      
      List_Delete(eregs);
    }
			
    
    if (k) 
      return fregs;
    else {
      List_Delete(fregs);
      return 0;
    }
  }

  void MF_RegionIDs_F2F4(MFace_ptr f, int *nfr, int *fregids) {
    int i;
    List_ptr fregs = MF_Regions_F2F4(f);
    *nfr = List_Num_Entries(fregs);
    for (i = 0; i < *nfr; i++)
      fregids[i] = MEnt_ID(List_Entry(fregs,i));
    List_Delete(fregs);
  }

  MRegion_ptr MF_Region_F2F4(MFace_ptr f, int ir) {
    List_ptr eregs;
    MRegion_ptr r, r1;
    MEdge_ptr e;
    int nr, i, dir, fdir;
    MFace_Adj_F2F4 *adj;
    
#ifdef DEBUG
    MSTK_Report("MF_Region_F4","More efficient to use MF_Regions",MSTK_MESG);
#endif

    dir = (ir) ? 0 : 1;
    
    adj = (MFace_Adj_F2F4 *) f->adj;
    e = List_Entry(adj->fedges,0);
    
    eregs = ME_Regions(e);
    if (eregs) {
      nr = List_Num_Entries(eregs);
      
      r1 = 0;
      for (i = 0; i < nr; i++) {
	r = List_Entry(eregs,i);
	fdir = MR_FaceDir(r,f);
	if (fdir == !dir) {
	  r1 = r;
	  break;
	}
      }
      List_Delete(eregs);
	
      return r1;
    }
    
    return 0;
  }

  int MF_RegionID_F2F4(MFace_ptr f, int ir) {
    MRegion_ptr r = MF_Region_F2F4(f,ir);
    return r ? MEnt_ID(r) : 0;
  }

  void MF_Add_Region_F2F4(MFace_ptr f, MRegion_ptr r, int side) {
#ifdef DEBUG
    MSTK_Report("MF_Add_Region",
		"Function call unsuitable for this representation",
		MSTK_WARN);
#endif
  }

  void MF_Rem_Region_F2F4(MFace_ptr f, MRegion_ptr r) {
#ifdef DEBUG
    MSTK_Report("MF_Rem_Region",
		"Function call unsuitable for this representation",
		MSTK_WARN);
#endif
  }


#ifdef __cplusplus
}
#endif
