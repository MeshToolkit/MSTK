#define _H_MFace_Private

#include "MFace.h"
#include "MFace_jmp.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif


  int MF_Set_GInfo_Auto_R4(MFace_ptr f) {
    int i, same, nv, fgdim, fgid, vgdim, vgid, vgdim0, vgid0;
    MVertex_ptr v;
    MFace_Adj_R4 *adj;

    adj = (MFace_Adj_R4 *) f->adj;
    nv = List_Num_Entries(adj->fvertices);

    same = 1;
    fgdim = -1;
    fgid = -1;
    vgid = -1;

    v = List_Entry(adj->fvertices,0);    
    vgid0 = MV_GEntID(v);
    vgdim0 = MV_GEntDim(v);

    for (i = 1; i < nv; i++) {
      v = List_Entry(adj->fvertices,i);
      vgid = MV_GEntID(v);
      vgdim = MV_GEntDim(v);
      if (vgdim == vgdim0 && vgid == vgid0)
	continue; /* all vertices have same classification so far */
      else {
	same = 0;
	break;
      }
    }
    if (same) {
      fgdim = vgdim0;
      fgid = vgid;
    }
     
    if (fgdim == -1 || fgdim < 2) {
      List_ptr fregions;

      /* We are unable to find proper classification info from the
	 vertices. Lets look at the number of regions connected to the
	 face and their classification */

      fregions = MF_Regions(f);
      
      if (fregions == NULL || List_Num_Entries(fregions) == 1) {
	
	/* In a complete mesh, this face must be on a model face */

	fgdim = 2;

      }
      else {
	MRegion_ptr fregion0, fregion1;
	int rgid0, rgid1;

	/* Internal face. Check if it is a mesh region or on an
	   interior interface */

	fregion0 = List_Entry(fregions,0); rgid0 = MEnt_GEntID(fregion0);
	fregion1 = List_Entry(fregions,1); rgid1 = MEnt_GEntID(fregion1);

	if (rgid0 == -1 || rgid1 == -1) {
	  
	  /* One of the regions is not classified properly. Just
	     assume this is an internal face */

	  fgdim = 3;	  

	}
	else {

	  fgdim = (rgid0 == rgid1) ? 3 : 2;

	}
      }

      if (fregions) List_Delete(fregions);
    }

    MEnt_Set_GEntDim((MEntity_ptr) f,fgdim);
    MEnt_Set_GEntID((MEntity_ptr) f,fgid);

    if (fgdim == 4)
      return 0;
    else
      return 1;
  }

  void MF_Set_Vertices_R4(MFace_ptr f, int n, MVertex_ptr *v) {
    MFace_Adj_R4 *adj;
    int i;

    adj = (MFace_Adj_R4 *) f->adj;
    if (adj->fvertices)
      List_Delete(adj->fvertices);
    adj->fvertices = List_New(n);

    for (i = 0; i < n; i++) {

#ifdef DEBUG
      if (MF_Mesh(f) != MV_Mesh(v[i]))
	MSTK_Report("MF_Set_Vertices_R4",
		    "Face and Vertex are not from the same mesh",
		    MSTK_FATAL);
#endif

      List_Add(adj->fvertices,v[i]);
    }
  }


  void MF_Rem_Vertex_R4(MFace_ptr f, MVertex_ptr v) {
    MFace_Adj_R4 *adj;

    adj = (MFace_Adj_R4 *) f->adj;
    if (adj->fvertices == NULL)
      MSTK_Report("MF_Rem_Vertex_R4",
		  "No initial set of vertices for face",MSTK_ERROR);

#ifdef DEBUG
    if (MF_Mesh(f) != MV_Mesh(v))
      MSTK_Report("MF_Rem_Vertex_R4",
		  "Face and Vertex are not from the same mesh",
		  MSTK_FATAL);
#endif

    List_Rem(adj->fvertices,v);
  }


  void MF_Replace_Vertex_i_R4(MFace_ptr f, int i, MVertex_ptr v) {
    MFace_Adj_R4 *adj;

    adj = (MFace_Adj_R4 *) f->adj;
    if (adj->fvertices == NULL)
      MSTK_Report("MF_Replace_Vertex_R4",
		  "No initial set of vertices for face",MSTK_ERROR);

#ifdef DEBUG
    if (MF_Mesh(f) != MV_Mesh(v))
      MSTK_Report("MF_Set_Vertices_R4",
		  "Face and Vertex are not from the same mesh",
		  MSTK_FATAL);
#endif

    List_Replacei(adj->fvertices,i,v);
  }

  void MF_Replace_Vertex_R4(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv) {
    MFace_Adj_R4 *adj;

    adj = (MFace_Adj_R4 *) f->adj;
    if (adj->fvertices == NULL)
      MSTK_Report("MF_Replace_Vertex_R4",
		  "No initial set of vertices for face",MSTK_ERROR);

#ifdef DEBUG
    if (MF_Mesh(f) != MV_Mesh(v))
      MSTK_Report("MF_Set_Vertices_R4",
		  "Face and Vertex are not from the same mesh",
		  MSTK_FATAL);
#endif

    List_Replace(adj->fvertices,v,nuv);
  }

  void MF_Insert_Vertex_R4(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v) {
    MFace_Adj_R4 *adj;

    adj = (MFace_Adj_R4 *) f->adj;
    if (adj->fvertices == NULL)
      adj->fvertices = List_New(4);

#ifdef DEBUG
    if (MF_Mesh(f) != MV_Mesh(nuv))
      MSTK_Report("MF_Set_Vertices_R4",
		  "Face and Vertex are not from the same mesh",
		  MSTK_FATAL);
#endif

    List_Insert(adj->fvertices,nuv,b4v);
  }

  void MF_Insert_Vertex_i_R4(MFace_ptr f, MVertex_ptr nuv, int i) {
    MFace_Adj_R4 *adj;

    adj = (MFace_Adj_R4 *) f->adj;
    if (adj->fvertices == NULL)
      adj->fvertices = List_New(4);

#ifdef DEBUG
    if (MF_Mesh(f) != MV_Mesh(nuv))
      MSTK_Report("MF_Set_Vertices_R4",
		  "Face and Vertex are not from the same mesh",
		  MSTK_FATAL);
#endif

    List_Inserti(adj->fvertices,nuv,i);
  }

  int MF_Num_Vertices_R4(MFace_ptr f) {
    MFace_Adj_R4 *adj;
    adj = (MFace_Adj_R4 *) f->adj;
    return List_Num_Entries(adj->fvertices);
  }

  List_ptr MF_Vertices_R4(MFace_ptr f, int dir, MVertex_ptr v0) {
    MFace_Adj_R4 *adj;
    List_ptr fverts;
    int i, k=0, nv, fnd=0;

    adj = (MFace_Adj_R4 *) f->adj;
    nv = List_Num_Entries(adj->fvertices);

    if (!v0) {
      if (dir) 
	fverts = List_Copy(adj->fvertices);
      else {
	fverts = List_New(nv);

	for (i = nv-1; i >= 0; i--)
	  List_Add(fverts,List_Entry(adj->fvertices,i));
      }
    }
    else {
      for (i = 0; i < nv; i++) {
	if (List_Entry(adj->fvertices,i) == v0) {
	  fnd = 1;
	  k = i;
	  break;
	}
      }

      if (!fnd) {
	MSTK_Report("MF_Num_Vertices_R4","Cannot find starting vertex",MSTK_ERROR);
	return 0;
      }

      fverts = List_New(nv);

      for (i = 0; i < nv; i++) {
	if (dir)
	  List_Add(fverts,List_Entry(adj->fvertices,(k+i)%nv));
	else
	  List_Add(fverts,List_Entry(adj->fvertices,(k+nv-i)%nv));
      }
    }

    return fverts;      
  }
	
  void MF_VertexIDs_R4(MFace_ptr f, int dir, int startvertid, int *nfv,
                       int *fvertids) {
    MFace_Adj_R2 *adj;
    List_ptr fverts;
    int i, k=0, nv, fnd=0;

    adj = (MFace_Adj_R2 *) f->adj;
    nv = List_Num_Entries(adj->fvertices);

    if (!startvertid) {
      if (dir) {
        for (i = 0; i < nv; i++)
          fvertids[i] = MEnt_ID(List_Entry(adj->fvertices,i));
      }
      else {
	for (i = nv-1; i >= 0; i--)
          fvertids[i] = MEnt_ID(List_Entry(adj->fvertices,i));
      }
    }
    else {
      for (i = 0; i < nv; i++) {
	if (MEnt_ID(List_Entry(adj->fvertices,i)) == startvertid) {
	  fnd = 1;
	  k = i;
	  break;
	}
      }

      if (!fnd) {
	MSTK_Report("MF_VertexIDs_R4","Cannot find starting vertex",MSTK_FATAL);
      }

      if (dir) {
        for (i = 0; i < nv; i++)
          fvertids[i] = MEnt_ID(List_Entry(adj->fvertices,(k+i)%nv));
      }
      else {
        for (i = 0; i < nv; i++)
          fvertids[i] = MEnt_ID(List_Entry(adj->fvertices,(k+nv-i)%nv));
      }
    }
  }	
	
  int MF_UsesVertex_R4(MFace_ptr f, MVertex_ptr v) {
    MFace_Adj_R4 *adj;
    adj = (MFace_Adj_R4 *) f->adj;
    return List_Contains(adj->fvertices,v);
  }

#ifdef __cplusplus
}
#endif
