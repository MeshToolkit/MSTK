#define _H_MRegion_Private

#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MR_Set_RepType_F4(MRegion_ptr r) {
    MR_Set_RepType_FNR3R4(r);
  }

  void MR_Delete_F4(MRegion_ptr r) {
    MRegion_DownAdj_FN *downadj;
    MFace_ptr f;
    MEdge_ptr e;
    Set_ptr fedges;
    int i, j, nf, ne;

    downadj = (MRegion_DownAdj_FN *) r->downadj;

    nf = Set_Num_Entries(downadj->rfaces);
    for (i = 0; i < nf; i++) {
      f = Set_Entry(downadj->rfaces,i);
      
      fedges = MF_Edges(f,1,0);
      ne = Set_Num_Entries(fedges);
      for (j = 0; j < ne; j++) {
	e = Set_Entry(fedges,j);
	ME_Rem_Region(e,r);
      }
      Set_Delete(fedges);
    }
    Set_Delete(downadj->rfaces);
    MSTK_free(downadj);

    MSTK_free(r);
  }

  void MR_Set_Faces_F4(MRegion_ptr r, int nf, MFace_ptr *rfaces, int *dirs) {
    int i, j, ne;
    MRegion_DownAdj_FN *downadj;
    MEdge_ptr e;
    Set_ptr fedges;

    downadj = (MRegion_DownAdj_FN *) r->downadj;
    downadj->nf = nf;
    downadj->fdirs = 0;
    downadj->rfaces = Set_New(nf);

    for (i = 0; i < nf; i++) {
      downadj->fdirs = downadj->fdirs | (dirs[i] << i);
      Set_Add(downadj->rfaces,rfaces[i]);
      
      fedges = MF_Edges(rfaces[i],1,0);
      ne = Set_Num_Entries(fedges);
      for (j = 0; j < ne; j++) {
	e = Set_Entry(fedges,j);
	ME_Add_Region(e,r);
	ME_Rem_Face(e,rfaces[i]);
      }
      Set_Delete(fedges);
    }
  }

  void MR_Set_Vertices_F4(MRegion_ptr r, int nv, MVertex_ptr *mvertices) {
    MR_Set_Vertices_FNR3R4(r,nv,mvertices);
  }

  int MR_Num_Faces_F4(MRegion_ptr r) {
    return MR_Num_Faces_FNR3R4(r);
  }

  /*
  int MR_Num_AdjRegions_F4(MRegion_ptr r) {
    return MR_Num_AdjRegions_FNR3R4(r);
  }
  */

  Set_ptr MR_Vertices_F4(MRegion_ptr r) {
    return MR_Vertices_FNR3R4(r);
  }

  Set_ptr MR_Edges_F4(MRegion_ptr r) {
    return MR_Edges_FNR3R4(r);
  }

  Set_ptr MR_Faces_F4(MRegion_ptr r) {
    return MR_Faces_FNR3R4(r);
  }

  Set_ptr MR_AdjRegions_F4(MRegion_ptr r) {
    return MR_AdjRegions_FNR3R4(r);
  }

  int MR_FaceDir_F4(MRegion_ptr r, MFace_ptr f) {
    return MR_FaceDir_FNR3R4(r,f);
  }

  int MR_FaceDir_i_F4(MRegion_ptr r, int i) {
    return MR_FaceDir_i_FNR3R4(r,i);
  }

  void MR_Replace_Face_F4(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf,int nudir){
    MR_Replace_Face_FNR3R4(r,f,nuf,nudir);
  }

  void MR_Replace_Face_i_F4(MRegion_ptr r, int i, MFace_ptr nuf, int nudir) {
    MR_Replace_Face_i_FNR3R4(r,i,nuf,nudir);
  }

  void MR_Replace_Vertex_F4(MRegion_ptr r, MVertex_ptr v, MVertex_ptr nuv){
    MR_Replace_Vertex_FNR3R4(r,v,nuv);
  }

  void MR_Replace_Vertex_i_F4(MRegion_ptr r, int i, MVertex_ptr nuv) {
    MR_Replace_Vertex_i_FNR3R4(r,i,nuv);
  }

  void MR_Add_AdjRegion_F4(MRegion_ptr r, int facenum, MRegion_ptr aregion) {
    MR_Add_AdjRegion_FNR3R4(r,facenum,aregion);
  }

  void MR_Rem_AdjRegion_F4(MRegion_ptr r, MRegion_ptr aregion) {
    MR_Rem_AdjRegion_FNR3R4(r,aregion);
  }

  int MR_UsesFace_F4(MRegion_ptr r, MFace_ptr f) {
    return MR_UsesFace_FNR3R4(r,f);
  }

  int MR_UsesEdge_F4(MRegion_ptr r, MEdge_ptr e) {
    return MR_UsesEdge_FNR3R4(r,e);
  }

  int MR_UsesVertex_F4(MRegion_ptr r, MVertex_ptr v) {
    return MR_UsesVertex_FNR3R4(r,v);
  }


#ifdef __cplusplus
}
#endif
