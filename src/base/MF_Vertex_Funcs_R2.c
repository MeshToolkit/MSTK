#define _H_MFace_Private

#include "MFace.h"
#include "MSTK_malloc.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  int MFs_AreSame_R2(MFace_ptr f1, MFace_ptr f2) {
    
    if (f1 == f2)
      return 1;
    else {
      int nfv1, nfv2, i, j, k, dir=-1;
      MVertex_ptr fv;
      List_ptr fverts1, fverts2;
      
      fverts1 = MF_Vertices(f1,1,0); nfv1 = List_Num_Entries(fverts1);
      fverts2 = MF_Vertices(f2,1,0); nfv2 = List_Num_Entries(fverts2);
      
      if (nfv1 != nfv2) {
	List_Delete(fverts1);
	List_Delete(fverts2);
	return 0;
      }
      
      fv = List_Entry(fverts1,0);
      k = List_Locate(fverts2,fv);
      if (k == -1) { /* Could not be found */
	List_Delete(fverts1);
	List_Delete(fverts2);
	return 0;
      }
      
      fv = List_Entry(fverts1,1);
      if (fv == List_Entry(fverts2,(k+1)%nfv1))
	dir = 1;
      else if (fv == List_Entry(fverts2,(k-1+nfv1)%nfv1))
	dir = 0;
      else { /* Could not be found */
	List_Delete(fverts1);
	List_Delete(fverts2);
	return 0;
      }
      
      if (dir) {
	for (j = 2; j < nfv1; j++) {
	  fv = List_Entry(fverts1,j);
	  if (fv != List_Entry(fverts2,(k+j)%nfv1)) {
	    List_Delete(fverts1); 
	    List_Delete(fverts2);
	    return 0;
	  }
	}
      }
      else {
	for (j = 2; j < nfv1; j++) {
	  fv = List_Entry(fverts1,j);
	  if (fv != List_Entry(fverts2,(k-j+nfv1)%nfv1)) {
	    List_Delete(fverts1); 
	    List_Delete(fverts2);
	    return 0;
	  }
	}
      }
      List_Delete(fverts1);
      List_Delete(fverts2);    

      return 1;
    }
  }


  int MF_Set_GInfo_Auto_R2(MFace_ptr f) {
    int i, same, nv, fgdim, fgid, vgdim, vgid, vgdim0, vgid0;
    MVertex_ptr v;
    MFace_Adj_R2 *adj;

    adj = (MFace_Adj_R2 *) f->adj;
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

  void MF_Set_Vertices_R2(MFace_ptr f, int n, MVertex_ptr *v) {
    MFace_Adj_R2 *adj;
    int i;

    adj = (MFace_Adj_R2 *) f->adj;
    if (adj->fvertices)
      List_Delete(adj->fvertices);
    adj->fvertices = List_New(n);

    for (i = 0; i < n; i++) {

#ifdef DEBUG
      if (MF_Mesh(f) != MV_Mesh(v[i]))
	MSTK_Report("MF_Set_Vertices_R2",
		    "Face and Vertex are not from the same mesh",
		    FATAL);
#endif

      List_Add(adj->fvertices,v[i]);
    }
  }

  void MF_Replace_Vertex_i_R2(MFace_ptr f, int i, MVertex_ptr v) {
    MFace_Adj_R2 *adj;

    adj = (MFace_Adj_R2 *) f->adj;
    if (adj->fvertices == NULL)
      MSTK_Report("MF_Replace_Vertex_R2",
		  "No initial set of vertices for face",ERROR);

#ifdef DEBUG
    if (MF_Mesh(f) != MV_Mesh(v))
      MSTK_Report("MF_Set_Vertices_R2",
		  "Face and Vertex are not from the same mesh",
		  FATAL);
#endif

    List_Replacei(adj->fvertices,i,v);
  }

  void MF_Replace_Vertex_R2(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv) {
    MFace_Adj_R2 *adj;

    adj = (MFace_Adj_R2 *) f->adj;
    if (adj->fvertices == NULL)
      MSTK_Report("MF_Replace_Vertex_R2",
		  "No initial set of vertices for face",ERROR);

#ifdef DEBUG
    if (MF_Mesh(f) != MV_Mesh(v))
      MSTK_Report("MF_Set_Vertices_R2",
		  "Face and Vertex are not from the same mesh",
		  FATAL);
#endif

    List_Replace(adj->fvertices,v,nuv);
  }

  void MF_Insert_Vertex_R2(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v) {
    MFace_Adj_R2 *adj;

    adj = (MFace_Adj_R2 *) f->adj;
    if (adj->fvertices == NULL)
      adj->fvertices = List_New(4);

#ifdef DEBUG
    if (MF_Mesh(f) != MV_Mesh(nuv))
      MSTK_Report("MF_Set_Vertices_R2",
		  "Face and Vertex are not from the same mesh",
		  FATAL);
#endif

    List_Insert(adj->fvertices,nuv,b4v);
  }

  void MF_Insert_Vertex_i_R2(MFace_ptr f, MVertex_ptr nuv, int i) {
    MFace_Adj_R2 *adj;

    adj = (MFace_Adj_R2 *) f->adj;
    if (adj->fvertices == NULL)
      adj->fvertices = List_New(4);

#ifdef DEBUG
    if (MF_Mesh(f) != MV_Mesh(nuv))
      MSTK_Report("MF_Set_Vertices_R2",
		  "Face and Vertex are not from the same mesh",
		  FATAL);
#endif

    List_Inserti(adj->fvertices,nuv,i);
  }

  int MF_Num_Vertices_R2(MFace_ptr f) {
    MFace_Adj_R2 *adj;
    adj = (MFace_Adj_R2 *) f->adj;
    return List_Num_Entries(adj->fvertices);
  }

  List_ptr MF_Vertices_R2(MFace_ptr f, int dir, MVertex_ptr v0) {
    MFace_Adj_R2 *adj;
    List_ptr fverts;
    int i, k=0, nv, fnd=0;

    adj = (MFace_Adj_R2 *) f->adj;
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
	MSTK_Report("MF_Num_Vertices_R2","Cannot find starting vertex",ERROR);
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
	
  int MF_UsesVertex_R2(MFace_ptr f, MVertex_ptr v) {
    MFace_Adj_R2 *adj;
    adj = (MFace_Adj_R2 *) f->adj;
    return List_Contains(adj->fvertices,v);
  }

#ifdef __cplusplus
}
#endif
