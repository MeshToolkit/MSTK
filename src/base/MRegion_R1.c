#define _H_MRegion_Private

#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MR_Set_RepType_R1(MRegion_ptr r) {
    MRegion_DownAdj_R1R2 *downadj;

    r->downadj = (MRegion_DownAdj_R1R2 *) MSTK_malloc(sizeof(MRegion_DownAdj_R1R2));
    downadj = r->downadj;
    downadj->nv = (unsigned char) 0;
    downadj->rvertices = NULL;
  }

  void MR_Delete_R1(MRegion_ptr r) {
    MRegion_DownAdj_R1R2 *downadj;

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    Set_Delete(downadj->rvertices);
    MSTK_free(downadj);
  }

  void MR_Set_Vertices_R1(MRegion_ptr r, int nv, MVertex_ptr *rvertices) {
    int i;
    MRegion_DownAdj_R1R2 *downadj;

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    downadj->nv = nv;
    downadj->rvertices = Set_New(nv);
    for (i = 0; i < nv; i++)
      Set_Add(downadj->rvertices,rvertices[i]);
  }

  void MR_Set_Faces_R1(MRegion_ptr r, int nf, MFace_ptr *mfaces, int *dirs) {
#ifdef DEBUG
    MSTK_Report("MR_Set_Faces_R1","Function call not suitable for this representation",WARN);
#endif
  }

  void MR_Replace_Face_R1(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf, int dir) {
#ifdef DEBUG
    MSTK_Report("MR_Replace_Face_R1","Function call not suitable for this representation",WARN);
#endif
  }

  void MR_Replace_Face_i_R1(MRegion_ptr r, int i, MFace_ptr nuf, int dir) {
#ifdef DEBUG
    MSTK_Report("MR_Replace_Face_i_R1","Function call not suitable for this representation",WARN);
#endif
  }

  int MR_Num_Faces_R1(MRegion_ptr r) {
    return mrtype_nf[MR_ElementType(r)];
  }

  int MR_Num_AdjRegions_R1(MRegion_ptr r) {
    Set_ptr adjr;
    int nr;

#ifdef DEBUG
    MSTK_Report("MR_Num_AdjRegions",
		"Inefficient to call this routine with this representation",
		WARN);
#endif

    adjr = MR_AdjRegions(r);
    if (adjr) {
      nr = Set_Num_Entries(adjr);
      Set_Delete(adjr);
      return nr;
    }
    else
      return 0;
  }

  Set_ptr MR_Vertices_R1(MRegion_ptr r) {
    MRegion_DownAdj_R1R2 *downadj;
    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    return Set_Copy(downadj->rvertices);
  }

  Set_ptr MR_Edges_R1(MRegion_ptr r) {
    MSTK_Report("MR_Edges","Not yet implemented for this representation",WARN);
    return NULL;
  }

  Set_ptr MR_Faces_R1(MRegion_ptr r) {
    MSTK_Report("MR_Faces","Not yet implemented for this representation",WARN);
    return NULL;
  }

  Set_ptr MR_AdjRegions_R1(MRegion_ptr r) {
    MSTK_Report("MR_Faces","Not yet implemented for this representation",WARN);
    return NULL;
  }

  int MR_FaceDir_R1(MRegion_ptr r, MFace_ptr f) {
    /* Must check if vertices of face are in same order as dictated by
       template */
    MSTK_Report("MR_FaceDir_R1","Not yet implemented for this representation",WARN);
    return 1; 
  }

  int MR_FaceDir_i_R1(MRegion_ptr r, int i) {
    return 1; /* Is this the right thing to do? */
  }


  int MR_UsesFace_R1(MRegion_ptr r, MFace_ptr f) {
    MSTK_Report("MR_UsesFace_R1","Not yet implemented for this representation",WARN);
    return 0;
  }

  int MR_UsesEdge_R1(MRegion_ptr r, MEdge_ptr e) {
    MSTK_Report("MR_UsesEdge_R1","Not yet implemented for this representation",WARN);
    return 0;
  }

  int MR_UsesVertex_R1(MRegion_ptr r, MVertex_ptr v) {
    MSTK_Report("MR_UsesVertex_R1","Not yet implemented for this representation",WARN);
    return 0;
  }

  void MR_Replace_Vertex_R1(MRegion_ptr r, MVertex_ptr v, MVertex_ptr nuv) {
    int i;
    MRegion_DownAdj_R1R2 *downadj;

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    for (i = 0; i < downadj->nv; i++)
      if (v == (MVertex_ptr) Set_Entry(downadj->rvertices,i)) {
	Set_Replacei(downadj->rvertices,i,nuv);
	return;
      }
  }

  void MR_Replace_Vertex_i_R1(MRegion_ptr r, int i, MVertex_ptr nuv) {
    MRegion_DownAdj_R1R2 *downadj;

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    Set_Replacei(downadj->rvertices,i,nuv);
  }


  void MR_Add_AdjRegion_R1(MRegion_ptr r, int facenum, MRegion_ptr ar) {
#ifdef DEBUG
    MSTK_Report("MR_Add_AdjRegion_R1","Function call not suitable for this representation",WARN);
#endif
  }

  void MR_Rem_AdjRegion_R1(MRegion_ptr r, MRegion_ptr ar) {
#ifdef DEBUG
    MSTK_Report("MR_Rem_AdjRegion_R1","Function call not suitable for this representation",WARN);
#endif
  }

#ifdef __cplusplus
}
#endif
