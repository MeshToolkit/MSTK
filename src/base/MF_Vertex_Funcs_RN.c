#define _H_MFace_Private

#include "MFace.h"
#include "MFace_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif


  int MF_Set_GInfo_Auto_RN(MFace_ptr f) {
    int i, same, nv, fgdim, fgid, vgdim, vgid, vgdim0, vgid0;
    MVertex_ptr v;
    MFace_DownAdj_RN *downadj;

    downadj = (MFace_DownAdj_RN *) f->downadj;
    nv = List_Num_Entries(downadj->fvertices);

    same = 1;
    fgdim = -1;
    fgid = -1;

    v = List_Entry(downadj->fvertices,0);    
    vgid0 = MV_GEntID(v);
    vgdim0 = MV_GEntDim(v);

    for (i = 1; i < nv; i++) {
      v = List_Entry(downadj->fvertices,i);
      vgid = MV_GEntID(v);
      vgdim = MV_GEntDim(v);
      if (vgdim == vgdim0 && vgid == vgid0)
	continue; /* all vertices have same classification so far */
      else {
	same = 0;
	if (vgdim > fgdim) {
	  fgdim = vgdim;
	  fgid = vgid;
	}
      }
    }
    if (same) {
      fgdim = vgdim0;
      fgid = vgid;
    }
     
    if (fgdim == -1)
      fgdim = 4;
    MEnt_Set_GEntDim(f,fgdim);
    MEnt_Set_GEntID(f,fgid);

    if (fgdim == 4)
      return 0;
    else
      return 1;
  }

  void MF_Set_Vertices_RN(MFace_ptr f, int n, MVertex_ptr *v) {
    MFace_DownAdj_RN *downadj;
    int i;

    downadj = (MFace_DownAdj_RN *) f->downadj;
    downadj->fvertices = List_New(n);

    for (i = 0; i < n; i++)
      List_Add(downadj->fvertices,v[i]);
  }

  void MF_Replace_Vertex_i_RN(MFace_ptr f, int i, MVertex_ptr v) {
    MFace_DownAdj_RN *downadj;

    downadj = (MFace_DownAdj_RN *) f->downadj;
    if (downadj->fvertices == NULL)
      MSTK_Report("MF_Replace_Vertex_RN",
		  "No initial set of vertices for face",ERROR);

    List_Replacei(downadj->fvertices,i,v);
  }

  void MF_Replace_Vertex_RN(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv) {
    MFace_DownAdj_RN *downadj;

    downadj = (MFace_DownAdj_RN *) f->downadj;
    if (downadj->fvertices == NULL)
      MSTK_Report("MF_Replace_Vertex_RN",
		  "No initial set of vertices for face",ERROR);

    List_Replace(downadj->fvertices,v,nuv);
  }

  void MF_Insert_Vertex_RN(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v) {
    MFace_DownAdj_RN *downadj;

    downadj = (MFace_DownAdj_RN *) f->downadj;
    if (downadj->fvertices == NULL)
      downadj->fvertices = List_New(4);

    List_Insert(downadj->fvertices,nuv,b4v);
  }

  void MF_Insert_Vertex_i_RN(MFace_ptr f, MVertex_ptr nuv, int i) {
    MFace_DownAdj_RN *downadj;

    downadj = (MFace_DownAdj_RN *) f->downadj;
    if (downadj->fvertices == NULL)
      downadj->fvertices = List_New(4);

    List_Inserti(downadj->fvertices,nuv,i);
  }

  int MF_Num_Vertices_RN(MFace_ptr f) {
    MFace_DownAdj_RN *downadj;
    downadj = (MFace_DownAdj_RN *) f->downadj;
    return List_Num_Entries(downadj->fvertices);
  }

  List_ptr MF_Vertices_RN(MFace_ptr f, int dir, MVertex_ptr v0) {
    MFace_DownAdj_RN *downadj;
    List_ptr fverts;
    int i, k=0, nv, fnd=0;

    downadj = (MFace_DownAdj_RN *) f->downadj;
    nv = List_Num_Entries(downadj->fvertices);

    if (!v0) {
      if (dir) 
	fverts = List_Copy(downadj->fvertices);
      else {
	fverts = List_New(nv);

	for (i = nv; i >= 0; i--)
	  List_Add(fverts,List_Entry(downadj->fvertices,i));
      }
    }
    else {
      for (i = 0; i < nv; i++) {
	if (List_Entry(downadj->fvertices,i) == v0) {
	  fnd = 1;
	  k = i;
	  break;
	}
      }

      if (!fnd) {
	MSTK_Report("MF_Num_Vertices_RN","Cannot find starting vertex",ERROR);
	return 0;
      }

      fverts = List_New(nv);

      for (i = 0; i < nv; i++) {
	if (dir)
	  List_Add(fverts,List_Entry(downadj->fvertices,(k+i)%nv));
	else
	  List_Add(fverts,List_Entry(downadj->fvertices,(k+nv-i)%nv));
      }
    }

    return fverts;      
  }
	
  int MF_UsesVertex_RN(MFace_ptr f, MVertex_ptr v) {
    MFace_DownAdj_RN *downadj;
    downadj = (MFace_DownAdj_RN *) f->downadj;
    return List_Contains(downadj->fvertices,v);
  }

#ifdef __cplusplus
}
#endif
