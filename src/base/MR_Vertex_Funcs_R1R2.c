#define _H_MRegion_Private

#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  int MR_Set_GInfo_Auto_R1R2(MRegion_ptr r) {
    int i, same, nv, rgdim, rgid, vgdim, vgid, vgdim0, vgid0;
    MVertex_ptr v;
    MRegion_DownAdj_R1R2 *downadj;

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    nv = List_Num_Entries(downadj->rvertices);

    same = 1;
    rgdim = -1;
    rgid = -1;
    vgid = -1;

    v = List_Entry(downadj->rvertices,0);    
    vgid0 = MV_GEntID(v);
    vgdim0 = MV_GEntDim(v);

    for (i = 1; i < nv; i++) {
      v = List_Entry(downadj->rvertices,i);
      vgid = MV_GEntID(v);
      vgdim = MV_GEntDim(v);
      if (vgdim == vgdim0 && vgid == vgid0)
	continue; /* all vertices have same classification so far */
      else {
	same = 0;
	if (vgdim > rgdim) {
	  rgdim = vgdim;
	  rgid = vgid;
	}
      }
    }
    if (same) {
      rgdim = vgdim0;
      rgid = vgid;
    }
     
    if (rgdim == -1)
      rgdim = 4;
    MEnt_Set_GEntDim(r,rgdim);
    MEnt_Set_GEntID(r,rgid);

    if (rgdim == 4)
      return 0;
    else
      return 1;
  }

  void MR_Set_Vertices_R1R2(MRegion_ptr r, int nv, MVertex_ptr *rvertices, 
			    int nf, int **rfvtemplate) {
    int i;
    MRegion_DownAdj_R1R2 *downadj;

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    downadj->rvertices = List_New(nv);
    for (i = 0; i < nv; i++)

#ifdef DEBUG
      if (MR_Mesh(r) != MV_Mesh(rvertices[i]))
	MSTK_Report("MR_Set_Vertices_R1R2",
		    "Region and vertex belong to different meshes",FATAL);
#endif

      List_Add(downadj->rvertices,rvertices[i]);

    /* If this is a non-standard element and the face vertex template
       has been specified store this information for later
       retrieval */
    if (nf && rfvtemplate) {
      downadj->fvtemplate = (int **) MSTK_malloc((nf+1)*sizeof(int *));
      downadj->fvtemplate[0] = (int *) MSTK_malloc(1*sizeof(int));
      for (i = 0; i < nf; i++) {
	int j, nfv;
	nfv = rfvtemplate[i][0];
	downadj->fvtemplate[i] = (int *) MSTK_malloc((nfv+1)*sizeof(int));
	downadj->fvtemplate[i][0] = nfv;
	for (j = 0; j < nfv; j++)
	  downadj->fvtemplate[i][j+1] = rfvtemplate[i][j+1];
      }
    }
  }


  List_ptr MR_Vertices_R1R2(MRegion_ptr r) {
    MRegion_DownAdj_R1R2 *downadj;
    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    return List_Copy(downadj->rvertices);
  }


  int MR_UsesVertex_R1R2(MRegion_ptr r, MVertex_ptr v) {
    MRegion_DownAdj_R1R2 *downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    return List_Contains(downadj->rvertices,v);
  }

  void MR_Replace_Vertex_R1R2(MRegion_ptr r, MVertex_ptr v, MVertex_ptr nuv) {
    int i, nv;
    MRegion_DownAdj_R1R2 *downadj;

#ifdef DEBUG
    if (MR_Mesh(r) != MV_Mesh(nuv))
      MSTK_Report("MR_Replace_Vertex_R1R2",
		  "Region and vertex belong to different meshes",
		  FATAL);
#endif

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    nv = List_Num_Entries(downadj->rvertices);
    for (i = 0; i < nv; i++)
      if (v == (MVertex_ptr) List_Entry(downadj->rvertices,i)) {
	List_Replacei(downadj->rvertices,i,nuv);
	return;
      }
  }

  void MR_Replace_Vertex_i_R1R2(MRegion_ptr r, int i, MVertex_ptr nuv) {
    MRegion_DownAdj_R1R2 *downadj;

#ifdef DEBUG
    if (MR_Mesh(r) != MV_Mesh(nuv))
      MSTK_Report("MR_Replace_Vertex_R1R2",
		  "Region and vertex belong to different meshes",
		  FATAL);
#endif

    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    List_Replacei(downadj->rvertices,i,nuv);
  }


#ifdef __cplusplus
}
#endif
