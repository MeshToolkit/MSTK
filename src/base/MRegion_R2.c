#define _H_MRegion_Private

#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MR_Set_RepType_R2(MRegion_ptr r) {
    MRegion_DownAdj_R1R2 *downadj;
    MRegion_SameAdj_R2 *sameadj;

    r->downadj = (MRegion_DownAdj_R1R2 *) MSTK_malloc(sizeof(MRegion_DownAdj_R1R2));
    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    downadj->nv = (unsigned char) 0;
    downadj->rvertices = NULL;

    r->sameadj = (MRegion_SameAdj_R2 *) MSTK_malloc(sizeof(MRegion_SameAdj_R2));
    sameadj = (MRegion_SameAdj_R2 *) r->sameadj;
    sameadj->nradj = (unsigned char) 0;
    sameadj->aregions = NULL;
  }

  void MR_Delete_R2(MRegion_ptr r, int keep) {
    MRegion_DownAdj_R1R2 *downadj;
    MRegion_SameAdj_R2 *sameadj;
    int i, nv, nr;
    MVertex_ptr v;
    MRegion_ptr r2;

    if (r->dim != MDELREGION) {
      downadj = (MRegion_DownAdj_R1R2 *) r->downadj;    
      nv = downadj->nv;
      for (i = 0; i < nv; i++) {
	v = List_Entry(downadj->rvertices,i);
	MV_Rem_Region(v,r);
      }
      
      sameadj = (MRegion_SameAdj_R2 *) r->sameadj;
      nr = sameadj->nradj;
      for (i = 0; i < nr; i++) {
	r2 = List_Entry(sameadj->aregions,i);
	MR_Rem_AdjRegion_R2(r,r2);
      }
    }

    if (keep) {
      MSTK_KEEP_DELETED = 1;
      r->dim = MDELREGION;
    }
    else {
#ifdef DEBUG
      r->dim = MDELREGION;
#endif

      List_Delete(downadj->rvertices);
      MSTK_free(downadj);
      
      List_Delete(sameadj->aregions);
      MSTK_free(sameadj);

      MSTK_free(r);
    }
  }

  void MR_Restore_R2(MRegion_ptr r) {
    MRegion_DownAdj_R1R2 *downadj;
    MRegion_SameAdj_R2 *sameadj;
    int i, nv, nr;
    MVertex_ptr v;
    MRegion_ptr r2;

    if (r->dim != MDELREGION)
      return;

    r->dim = MREGION;

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    nv = downadj->nv;
    for (i = 0; i < nv; i++) {
      v = List_Entry(downadj->rvertices,i);
      MV_Add_Region(v,r);
    }

    sameadj = (MRegion_SameAdj_R2 *) r->sameadj;
    nr = sameadj->nradj;
    for (i = 0; i < nr; i++) {
      r2 = List_Entry(sameadj->aregions,i);
      MR_Add_AdjRegion_R2(r,i,r2);
    }
  }


  void MR_Set_Vertices_R2(MRegion_ptr r, int nv, MVertex_ptr *rvertices) {
    int i;
    MRegion_DownAdj_R1R2 *downadj;

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    downadj->nv = nv;
    downadj->rvertices = List_New(nv);
    for (i = 0; i < nv; i++)
      List_Add(downadj->rvertices,rvertices[i]);
  }

  void MR_Add_AdjRegion_R2(MRegion_ptr r, int facenum, MRegion_ptr aregion) {
    MRegion_SameAdj_R2 *sameadj;

    /* Is r->sameadj allocated ? */
    /* We should make use of the facenum info */
    sameadj = (MRegion_SameAdj_R2 *) r->sameadj;
    List_Add(sameadj->aregions,aregion);
    sameadj->nradj++;
  }

  void MR_Rem_AdjRegion_R2(MRegion_ptr r, MRegion_ptr aregion) {
    MRegion_SameAdj_R2 *sameadj;

    sameadj = (MRegion_SameAdj_R2 *) r->sameadj;
    if (List_Rem(sameadj->aregions,aregion))
      sameadj->nradj--;
  }

  void MR_Set_Faces_R2(MRegion_ptr r, int nf, MFace_ptr *mfaces, int *dirs) {
#ifdef DEBUG
    MSTK_Report("MR_Set_Faces_R2","Function call not suitable for this representation",WARN);
#endif
  }

  void MR_Replace_Face_R2(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf, int dir) {
#ifdef DEBUG
    MSTK_Report("MR_Replace_Face_R2","Function call not suitable for this representation",WARN);
#endif
  }

  void MR_Replace_Face_i_R2(MRegion_ptr r, int i, MFace_ptr nuf, int dir) {
#ifdef DEBUG
    MSTK_Report("MR_Replace_Face_i_R2","Function call not suitable for this representation",WARN);
#endif
  }

  int MR_Num_Faces_R2(MRegion_ptr r) {
    return mrtype_nf[MR_ElementType(r)];
  }

  int MR_Num_AdjRegions_R2(MRegion_ptr r) {
    return ((MRegion_SameAdj_R2 *) r->sameadj)->nradj;
  }

  List_ptr MR_Vertices_R2(MRegion_ptr r) {
    MRegion_DownAdj_R1R2 *downadj;

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    return List_Copy(downadj->rvertices);
  }

  List_ptr MR_Edges_R2(MRegion_ptr r) {
    MSTK_Report("MR_Edges","Not yet implemented for this representation",WARN);
    return NULL;
  }

  List_ptr MR_Faces_R2(MRegion_ptr r) {
    MSTK_Report("MR_Faces","Not yet implemented for this representation",WARN);
    return NULL;
  }

  List_ptr MR_AdjRegions_R2(MRegion_ptr r) {
    MRegion_SameAdj_R2 *sameadj;

    sameadj = (MRegion_SameAdj_R2 *) r->sameadj;
    return List_Copy(sameadj->aregions);
  }

  int MR_FaceDir_R2(MRegion_ptr r, MFace_ptr f) {
    /* Must check if vertices of face are in same order as dictated by
       template */
    return 1; 
  }

  int MR_FaceDir_i_R2(MRegion_ptr r, int i) {
    return 1; /* Is this the right thing to do? */
  }

  int MR_UsesFace_R2(MRegion_ptr r, MFace_ptr f) {
    MSTK_Report("MR_UsesFace_R1","Not yet implemented for this representation",WARN);
    return 0;
  }

  int MR_UsesEdge_R2(MRegion_ptr r, MEdge_ptr e) {
    MSTK_Report("MR_UsesEdge_R1","Not yet implemented for this representation",WARN);
    return 0;
  }

  int MR_UsesVertex_R2(MRegion_ptr r, MVertex_ptr v) {
    MSTK_Report("MR_UsesVertex_R1","Not yet implemented for this representation",WARN);
    return 0;
  }

  void MR_Replace_Vertex_R2(MRegion_ptr r, MVertex_ptr v, MVertex_ptr nuv) {
    int i;
    MRegion_DownAdj_R1R2 *downadj;

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    for (i = 0; i < downadj->nv; i++)
      if (v == (MVertex_ptr) List_Entry(downadj->rvertices,i)) {
	List_Replacei(downadj->rvertices,i,nuv);
	return;
      }
  }

  void MR_Replace_Vertex_i_R2(MRegion_ptr r, int i, MVertex_ptr nuv) {
    MRegion_DownAdj_R1R2 *downadj;

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    List_Replacei(downadj->rvertices,i,nuv);
  }

#ifdef __cplusplus
}
#endif
