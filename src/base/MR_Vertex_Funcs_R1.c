#define _H_MRegion_Private

#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  int MR_Set_GInfo_Auto_R1(MRegion_ptr r) {
    int i, same, nv, rgdim, rgid, vgdim, vgid, vgdim0, vgid0;
    MVertex_ptr v;
    MRegion_Adj_R1 *adj;

    adj = (MRegion_Adj_R1 *) r->adj;
    nv = List_Num_Entries(adj->rvertices);

    same = 1;
    rgdim = -1;
    rgid = -1;
    vgid = -1;

    v = List_Entry(adj->rvertices,0);    
    vgid0 = MV_GEntID(v);
    vgdim0 = MV_GEntDim(v);

    for (i = 1; i < nv; i++) {
      v = List_Entry(adj->rvertices,i);
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
    MEnt_Set_GEntDim((MEntity_ptr) r,rgdim);
    MEnt_Set_GEntID((MEntity_ptr) r,rgid);

    if (rgdim == 4)
      return 0;
    else
      return 1;
  }

  void MR_Set_Vertices_R1(MRegion_ptr r, int nv, MVertex_ptr *rvertices, 
			    int nf, int **rfvtemplate) {
    int i;
    MRegion_Adj_R1 *adj;

    adj = (MRegion_Adj_R1 *) r->adj;
    adj->rvertices = List_New(nv);
    for (i = 0; i < nv; i++) {

#ifdef DEBUG
      if (MR_Mesh(r) != MV_Mesh(rvertices[i]))
	MSTK_Report("MR_Set_Vertices_R1",
		    "Region and vertex belong to different meshes",MSTK_FATAL);
#endif

      List_Add(adj->rvertices,rvertices[i]);
      MV_Add_Region(rvertices[i], r);
    }

    /* If this is a non-standard element and the face vertex template
       has been specified store this information for later
       retrieval */
    if (nf && rfvtemplate) {
      adj->fvtemplate = (int **) MSTK_malloc((nf+1)*sizeof(int *));
      adj->fvtemplate[0] = (int *) MSTK_malloc(1*sizeof(int));
      for (i = 0; i < nf; i++) {
	int j, nfv;
	nfv = rfvtemplate[i][0];
	adj->fvtemplate[i] = (int *) MSTK_malloc((nfv+1)*sizeof(int));
	adj->fvtemplate[i][0] = nfv;
	for (j = 0; j < nfv; j++)
	  adj->fvtemplate[i][j+1] = rfvtemplate[i][j+1];
      }

      switch (nf) {
      case 4: {
        int alltri = 1;
        for (i = 0; i < nf; i++) 
          if (rfvtemplate[i][0] != 3) {
            alltri = 0;
            break;
          }
        if (alltri)
          r->mrtype = TET;
        else
          r->mrtype = POLYHED;
        break;
      }
      case 5: {
        int nquads = 0; 
        int ntris = 0;
        for (i = 0; i < nf; i++)
          if (rfvtemplate[i][0] == 3)
            ntris++;
          else if (rfvtemplate[i][1] == 4)
            nquads++;

        if (nquads == 3 && ntris == 2)
          r->mrtype = PRISM;
        else if (nquads == 1 && ntris == 4)
          r->mrtype = PYRAMID;
        else
          r->mrtype = POLYHED;
        break;
      }
      case 6: {
        int allquad = 1;
        for (i = 0; i < nf; i++)
          if (rfvtemplate[i][0] != 4) {
            allquad = 0;
            break;
          }
        if (allquad)
          r->mrtype = HEX;
        else
          r->mrtype = POLYHED;
        break;
      }
      default:
        r->mrtype = POLYHED;
      }
    }
    else {
      switch (nv) {
      case 4:
	r->mrtype = TET;
	nf = 4;
	break;
      case 5:
	r->mrtype = PYRAMID;
	nf = 5;
	break;
      case 6:
	r->mrtype = PRISM;
	nf = 5;
	break;
      case 8:
	r->mrtype = HEX;
	nf = 6;
	break;
      default:
	MSTK_Report("MR_Set_Vertices",
		    "Polyhedron: Need number of faces and vertex template",
		    MSTK_FATAL);
	r->mrtype = RUNKNOWN;
	break;
      }
    }
  }


  List_ptr MR_Vertices_R1(MRegion_ptr r) {
    MRegion_Adj_R1 *adj = (MRegion_Adj_R1 *) r->adj;
    return List_Copy(adj->rvertices);
  }


  void MR_VertexIDs_R1(MRegion_ptr r, int *nrv, int *rvertids) {
    int i;
    MRegion_Adj_R1 *adj = (MRegion_Adj_R1 *) r->adj;
    *nrv = List_Num_Entries(adj->rvertices);
    for (i = 0; i < *nrv; i++)
      rvertids[i] = MEnt_ID(List_Entry(adj->rvertices,i));
  }

  int MR_UsesVertex_R1(MRegion_ptr r, MVertex_ptr v) {
    MRegion_Adj_R1 *adj = (MRegion_Adj_R1 *) r->adj;
    return List_Contains(adj->rvertices,v);
  }

  void MR_Replace_Vertex_R1(MRegion_ptr r, MVertex_ptr v, MVertex_ptr nuv) {
    int i, nv;
    MRegion_Adj_R1 *adj;

#ifdef DEBUG
    if (MR_Mesh(r) != MV_Mesh(nuv))
      MSTK_Report("MR_Replace_Vertex_R1",
		  "Region and vertex belong to different meshes",
		  MSTK_FATAL);
#endif

    adj = (MRegion_Adj_R1 *) r->adj;
    nv = List_Num_Entries(adj->rvertices);
    for (i = 0; i < nv; i++)
      if (v == (MVertex_ptr) List_Entry(adj->rvertices,i)) {
	List_Replacei(adj->rvertices,i,nuv);
	return;
      }
  }

  void MR_Replace_Vertex_i_R1(MRegion_ptr r, int i, MVertex_ptr nuv) {
    MRegion_Adj_R1 *adj;

#ifdef DEBUG
    if (MR_Mesh(r) != MV_Mesh(nuv))
      MSTK_Report("MR_Replace_Vertex_R1",
		  "Region and vertex belong to different meshes",
		  MSTK_FATAL);
#endif

    adj = (MRegion_Adj_R1 *) r->adj;
    List_Replacei(adj->rvertices,i,nuv);
  }


#ifdef __cplusplus
}
#endif
