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
    downadj->rvertices = NULL;

    r->sameadj = (MRegion_SameAdj_R2 *) MSTK_malloc(sizeof(MRegion_SameAdj_R2));
    sameadj = (MRegion_SameAdj_R2 *) r->sameadj;
    sameadj->aregions = NULL;
  }

  void MR_Delete_R2(MRegion_ptr r, int keep) {
    MRegion_DownAdj_R1R2 *downadj;
    MRegion_SameAdj_R2 *sameadj;
    int idx, nr;
    MVertex_ptr v;
    MRegion_ptr r2;

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;    
    sameadj = (MRegion_SameAdj_R2 *) r->sameadj;

    if (MEnt_Dim(r) != MDELETED) {
      if (downadj) {
	idx = 0;
	while ((v = List_Next_Entry(downadj->rvertices,&idx)))
	  MV_Rem_Region(v,r);
      }
      
      if (sameadj) {
	idx = 0;
	while ((r2 = List_Next_Entry(sameadj->aregions,&idx)))
	  MR_Rem_AdjRegion_R2(r,r2);
      }
    }

    if (!keep) {
      if (downadj) {
	if (downadj->rvertices)
	  List_Delete(downadj->rvertices);
	if (downadj->fvtemplate) {
	  int nf = downadj->fvtemplate[0][0], i;
	  for (i = 0; i < nf; i++)
	    MSTK_free(downadj->fvtemplate[i]);
	  MSTK_free(downadj->fvtemplate);
	}
	MSTK_free(downadj);
      }
      
      if (sameadj) {
	if (sameadj->aregions)
	  List_Delete(sameadj->aregions);
	MSTK_free(sameadj);
      }
    }
  }

  void MR_Restore_R2(MRegion_ptr r) {
    MRegion_DownAdj_R1R2 *downadj;
    MRegion_SameAdj_R2 *sameadj;
    int idx, nr,i;
    MVertex_ptr v;
    MRegion_ptr r2;

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    idx = 0;
    while ((v = List_Next_Entry(downadj->rvertices,&idx)))
      MV_Add_Region(v,r);

    sameadj = (MRegion_SameAdj_R2 *) r->sameadj;
    idx = 0; i = 0;
    while ((r2 = List_Next_Entry(sameadj->aregions,&idx))) {
      MR_Add_AdjRegion_R2(r,i,r2);
      i++;
    }
  }


  void MR_Destroy_For_MESH_Delete_R2(MRegion_ptr r) {
    MRegion_DownAdj_R1R2 *downadj;
    MRegion_SameAdj_R2 *sameadj;

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    if (downadj) {
      if (downadj->rvertices)
	List_Delete(downadj->rvertices);
      if (downadj->fvtemplate) {
	int nf = downadj->fvtemplate[0][0], i;
	for (i = 0; i < nf; i++)
	  MSTK_free(downadj->fvtemplate[i]);
	MSTK_free(downadj->fvtemplate);
      }
      MSTK_free(downadj);
    }
    
    sameadj = (MRegion_SameAdj_R2 *) r->sameadj;
    if (sameadj) {
      if (sameadj->aregions)
	List_Delete(sameadj->aregions);
      MSTK_free(sameadj);
    }
  }

  int MR_Set_GInfo_Auto_R2(MRegion_ptr r) {
    return MR_Set_GInfo_Auto_R1R2(r);
  }

  void MR_Set_Vertices_R2(MRegion_ptr r, int nv, MVertex_ptr *rvertices, 
			  int nf, int **rfvtemplate) {
    MR_Set_Vertices_R1R2(r,nv,rvertices,nf,rfvtemplate);
  }

  void MR_Add_AdjRegion_R2(MRegion_ptr r, int facenum, MRegion_ptr aregion) {
    MRegion_SameAdj_R2 *sameadj;

    /* Is r->sameadj allocated ? */
    /* We should make use of the facenum info */
    sameadj = (MRegion_SameAdj_R2 *) r->sameadj;
    List_Add(sameadj->aregions,aregion);
  }

  void MR_Rem_AdjRegion_R2(MRegion_ptr r, MRegion_ptr aregion) {
    MRegion_SameAdj_R2 *sameadj;

    sameadj = (MRegion_SameAdj_R2 *) r->sameadj;
    List_Rem(sameadj->aregions,aregion);
  }

  void MR_Set_Faces_R2(MRegion_ptr r, int nf, MFace_ptr *mfaces, int *dirs) {
#ifdef DEBUG
    MSTK_Report("MR_Set_Faces_R2",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MR_Replace_Face_R2(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf, int dir) {
#ifdef DEBUG
    MSTK_Report("MR_Replace_Face_R2",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MR_Replace_Face_i_R2(MRegion_ptr r, int i, MFace_ptr nuf, int dir) {
#ifdef DEBUG
    MSTK_Report("MR_Replace_Face_i_R2",
		"Function call not suitable for this representation",WARN);
#endif
  }

  int MR_Num_Faces_R2(MRegion_ptr r) {
    return MR_Num_Faces_R1R2(r);
  }

  int MR_Num_AdjRegions_R2(MRegion_ptr r) {
    List_ptr adjregs = ((MRegion_SameAdj_R2 *) r->sameadj)->aregions;
    return List_Num_Entries(adjregs);
  }

  List_ptr MR_Vertices_R2(MRegion_ptr r) {
    return MR_Vertices_R1R2(r);
  }

  List_ptr MR_Edges_R2(MRegion_ptr r) {
    return MR_Edges_R1R2(r);
  }

  List_ptr MR_Faces_R2(MRegion_ptr r) {
    return MR_Faces_R1R2(r);
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
    return MR_FaceDir_i_R1R2(r,i);
  }

  int MR_UsesFace_R2(MRegion_ptr r, MFace_ptr f) {
    return MR_UsesFace_R1R2(r,f);
  }

  int MR_UsesEdge_R2(MRegion_ptr r, MEdge_ptr e) {
    return MR_UsesEdge_R1R2(r,e);
  }

  int MR_UsesVertex_R2(MRegion_ptr r, MVertex_ptr v) {
    return MR_UsesVertex_R1R2(r,v);
  }

  void MR_Replace_Vertex_R2(MRegion_ptr r, MVertex_ptr v, MVertex_ptr nuv) {
    MR_Replace_Vertex_R1R2(r,v,nuv);
  }

  void MR_Replace_Vertex_i_R2(MRegion_ptr r, int i, MVertex_ptr nuv) {
    MR_Replace_Vertex_i_R1R2(r,i,nuv);
  }

#ifdef __cplusplus
}
#endif
