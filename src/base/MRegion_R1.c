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
    downadj->rvertices = NULL;
  }

  void MR_Delete_R1(MRegion_ptr r, int keep) {
    MRegion_DownAdj_R1R2 *downadj;
    int idx;
    MVertex_ptr v;

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;

    if (MEnt_Dim(r) != MDELETED) { /* if regn hasnt been temporarily deleted */
      if (downadj) {
	idx = 0;
	while ((v = List_Next_Entry(downadj->rvertices,&idx)))
	  MV_Rem_Region(v,r);
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

      MSTK_free(r);
    }
  }

  void MR_Restore_R1(MRegion_ptr r) {
    MRegion_DownAdj_R1R2 *downadj;
    int idx;
    MVertex_ptr v;

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    idx = 0;
    while ((v = List_Next_Entry(downadj->rvertices,&idx)))
      MV_Add_Region(v,r);
  }

  void MR_Destroy_For_MESH_Delete_R1(MRegion_ptr r) {
    MRegion_DownAdj_R1R2 *downadj;

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
  }


  int MR_Set_GInfo_Auto_R1(MRegion_ptr r) {
    return MR_Set_GInfo_Auto_R1R2(r);
  }

  void MR_Set_Vertices_R1(MRegion_ptr r, int nv, MVertex_ptr *rvertices, 
			  int nf, int **template) {
    MR_Set_Vertices_R1R2(r,nv,rvertices,nf,template);
  }

  void MR_Set_Faces_R1(MRegion_ptr r, int nf, MFace_ptr *mfaces, int *dirs) {
#ifdef DEBUG
    MSTK_Report("MR_Set_Faces_R1",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MR_Replace_Face_R1(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf, int dir) {
#ifdef DEBUG
    MSTK_Report("MR_Replace_Face_R1",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MR_Replace_Face_i_R1(MRegion_ptr r, int i, MFace_ptr nuf, int dir) {
#ifdef DEBUG
    MSTK_Report("MR_Replace_Face_i_R1",
		"Function call not suitable for this representation",WARN);
#endif
  }

  int MR_Num_Faces_R1(MRegion_ptr r) {
    MRegion_DownAdj_R1R2 *downadj;
    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    
    if (downadj->fvtemplate) {
      return downadj->fvtemplate[0][0];
    }
    else {
      switch (List_Num_Entries(downadj->rvertices)) {
      case 4:
	return 4;  /* Tet */
      case 5:
	return 5;  /* Pyramid */
      case 6:
	return 5;  /* Prism */
      case 8:
	return 6;  /* Hex */
      default:
	return 0;
      }
    }
  }

  int MR_Num_AdjRegions_R1(MRegion_ptr r) {
    List_ptr adjr;
    int nr;

#ifdef DEBUG
    MSTK_Report("MR_Num_AdjRegions",
		"Inefficient to call this routine with this representation",
		WARN);
#endif

    adjr = MR_AdjRegions_R1(r);
    if (adjr) {
      nr = List_Num_Entries(adjr);
      List_Delete(adjr);
      return nr;
    }
    else
      return 0;
  }

  List_ptr MR_Vertices_R1(MRegion_ptr r) {
    return MR_Vertices_R1R2(r);
  }

  List_ptr MR_Edges_R1(MRegion_ptr r) {
    return MR_Edges_R1R2(r);
  }

  List_ptr MR_Faces_R1(MRegion_ptr r) {
    return MR_Faces_R1R2(r);
  }

  List_ptr MR_AdjRegions_R1(MRegion_ptr r) {
    MSTK_Report("MR_AdjRegions",
		"Not yet implemented for this representation",WARN);
    return NULL;
  }

  int MR_FaceDir_R1(MRegion_ptr r, MFace_ptr f) {
    return MR_FaceDir_R1R2(r,f);
  }

  int MR_FaceDir_i_R1(MRegion_ptr r, int i) {
    return 1; /* Is this the right thing to do? Yes! */
  }


  int MR_UsesFace_R1(MRegion_ptr r, MFace_ptr f) {
    return MR_UsesFace_R1R2(r,f);
  }

  int MR_UsesEdge_R1(MRegion_ptr r, MEdge_ptr e) {
    return MR_UsesEdge_R1R2(r,e);
  }

  int MR_UsesVertex_R1(MRegion_ptr r, MVertex_ptr v) {
    return MR_UsesVertex_R1R2(r,v);
  }

  void MR_Replace_Vertex_R1(MRegion_ptr r, MVertex_ptr v, MVertex_ptr nuv) {
    MR_Replace_Vertex_R1R2(r,v,nuv);
  }

  void MR_Replace_Vertex_i_R1(MRegion_ptr r, int i, MVertex_ptr nuv) {
    MR_Replace_Vertex_i_R1R2(r,i,nuv);
  }


  void MR_Add_AdjRegion_R1(MRegion_ptr r, int facenum, MRegion_ptr ar) {
#ifdef DEBUG
    MSTK_Report("MR_Add_AdjRegion_R1",
		"Function call not suitable for this representation",WARN);
#endif
  }

  void MR_Rem_AdjRegion_R1(MRegion_ptr r, MRegion_ptr ar) {
#ifdef DEBUG
    MSTK_Report("MR_Rem_AdjRegion_R1",
		"Function call not suitable for this representation",WARN);
#endif
  }

#ifdef __cplusplus
}
#endif
