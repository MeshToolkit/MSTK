#define _H_MRegion_Private

#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MR_Set_RepType_F1(MRegion_ptr r) {
    MR_Set_RepType_FNR3R4(r);
  }

  void MR_Delete_F1(MRegion_ptr r) {
    MR_Delete_FNR3R4(r);
  }

  void MR_Set_Faces_F1(MRegion_ptr r, int nf, MFace_ptr *rfaces, int *dirs) {
    MR_Set_Faces_FNR3R4(r,nf,rfaces,dirs);
  }

  void MR_Set_Vertices_F1(MRegion_ptr r, int nv, MFace_ptr *mvertices) {
    MR_Set_Vertices_FNR3R4(r,nv,mvertices);
  }

  int MR_Num_Faces_F1(MRegion_ptr r) {
    return MR_Num_Faces_FNR3R4(r);
  }

  /*
  int MR_Num_AdjRegions_F1(MRegion_ptr r) {
    return MR_Num_AdjRegions_FNR3R4(r);
  }
  */

  List_ptr MR_Vertices_F1(MRegion_ptr r) {
    return MR_Vertices_FNR3R4(r);
  }

  List_ptr MR_Edges_F1(MRegion_ptr r) {
    return MR_Edges_FNR3R4(r);
  }

  List_ptr MR_Faces_F1(MRegion_ptr r) {
    return MR_Faces_FNR3R4(r);
  }

  List_ptr MR_AdjRegions_F1(MRegion_ptr r) {
    return MR_AdjRegions_FNR3R4(r);
  }

  int MR_FaceDir_F1(MRegion_ptr r, MFace_ptr f) {
    return MR_FaceDir_FNR3R4(r,f);
  }

  int MR_FaceDir_i_F1(MRegion_ptr r, int i) {
    return MR_FaceDir_i_FNR3R4(r,i);
  }

  void MR_Replace_Face_F1(MRegion_ptr r, MFace_ptr f, MFace_ptr nuf,int nudir){
    MR_Replace_Face_FNR3R4(r,f,nuf,nudir);
  }

  void MR_Replace_Face_i_F1(MRegion_ptr r, int i, MFace_ptr nuf, int nudir) {
    MR_Replace_Face_i_FNR3R4(r,i,nuf,nudir);
  }

  void MR_Replace_Vertex_F1(MRegion_ptr r, MVertex_ptr v, MVertex_ptr nuv){
    MR_Replace_Vertex_FNR3R4(r,v,nuv);
  }

  void MR_Replace_Vertex_i_F1(MRegion_ptr r, int i, MVertex_ptr nuv) {
    MR_Replace_Vertex_i_FNR3R4(r,i,nuv);
  }

  void MR_Add_AdjRegion_F1(MRegion_ptr r, int facenum, MRegion_ptr aregion) {
    MR_Add_AdjRegion_FNR3R4(r,facenum,aregion);
  }

  void MR_Rem_AdjRegion_F1(MRegion_ptr r, MRegion_ptr aregion) {
    MR_Rem_AdjRegion_FNR3R4(r,aregion);
  }

  int MR_UsesFace_F1(MRegion_ptr r, MFace_ptr f) {
    return MR_UsesFace_FNR3R4(r,f);
  }

  int MR_UsesEdge_F1(MRegion_ptr r, MEdge_ptr e) {
    return MR_UsesEdge_FNR3R4(r,e);
  }

  int MR_UsesVertex_F1(MRegion_ptr r, MVertex_ptr v) {
    return MR_UsesVertex_FNR3R4(r,v);
  }

#ifdef __cplusplus
}
#endif
