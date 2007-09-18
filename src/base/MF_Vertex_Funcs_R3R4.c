#define _H_MFace_Private

#include "MFace.h"
#include "MFace_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MF_Set_Vertices_R3R4(MFace_ptr f, int n, MVertex_ptr *v) {
    MFace_Adj_R1 *adj;
    int i;

    adj = f->adj;
    adj->fvertices = List_New(n);

    for (i = 0; i < n; i++)
      List_Add(adj->fvertices,v[i]);
  }

  void MF_Replace_Vertex_i_R3R4(MFace_ptr f, int i, MVertex_ptr v) {
    MFace_Adj_R1 *adj;

    adj = f->adj;
    if (!adj->fvertices)
      MSTK_Report("MF_Replace_Vertex_R3R4","No initial set of vertices for face",ERROR);

    List_Replacei(adj->fvertices,i,v);
  }

  void MF_Replace_Vertex_R3R4(MFace_ptr f, MVertex_ptr v, MVertex_ptr nuv) {
    MFace_Adj_R1 *adj;

    adj = f->adj;
    if (!adj->fvertices)
      MSTK_Report("MF_Replace_Vertex_R3R4","No initial set of vertices for face",ERROR);

    List_Replace(adj->fvertices,v,nuv);
  }

  void MF_Insert_Vertex_R3R4(MFace_ptr f, MVertex_ptr nuv, MVertex_ptr b4v) {
    MFace_Adj_R1 *adj;

    adj = f->adj;
    if (!adj->fvertices)
      adj->fvertices = List_New(4);

    List_Insert(adj->fvertices,nuv,b4v);
  }

  void MF_Insert_Vertex_i_R3R4(MFace_ptr f, MVertex_ptr nuv, int i) {
    MFace_Adj_R1 *adj;

    adj = f->adj;
    if (!adj->fvertices)
      adj->fvertices = List_New(4);

    List_Inserti(adj->fvertices,nuv,i);
  }

  int MF_Num_Vertices_R3R4(MFace_ptr f) {
    MFace_Adj_R1 *adj;
    adj = (MFace_Adj_R1 *) f->adj;
    if (adj->fvertices)
      return List_Num_Entries(adj->fvertices);
    else
      return 0;
  }

  List_ptr MF_Vertices_R3R4(MFace_ptr f, int dir, MVertex_ptr v0) {
    MFace_Adj_R1 *adj;
    List_ptr fverts;
    int i, k=0, nv, fnd=0;

    adj = (MFace_Adj_R1 *) f->adj;
    if (!adj->fvertices)
      return 0;

    nv = List_Num_Entries(adj->fvertices);

    if (!v0) {
      if (dir) 
	fverts = List_Copy(adj->fvertices);
      else {
	fverts = List_New(nv);

	for (i = 0; i < nv; i++)
	  List_Add(fverts,List_Entry(adj->fvertices,i));
      }
    }
    else {
      fverts = List_New(nv);

      for (i = 0; i < nv; i++) {
	if (List_Entry(adj->fvertices,i) == v0) {
	  fnd = 1;
	  k = i;
	  break;
	}
      }

      if (!fnd) {
	MSTK_Report("MF_Num_Vertices_R3R4","Cannot find starting vertex",ERROR);
	return 0;
      }

      for (i = 0; i < nv; i++) {
	if (dir)
	  List_Add(fverts,List_Entry(adj->fvertices,(k+i)%nv));
	else
	  List_Add(fverts,List_Entry(adj->fvertices,(k+nv-i)%nv));
      }
    }

    return fverts;      
  }
	
  int MF_UsesVertex_R3R4(MFace_ptr f, MVertex_ptr v) {
    MFace_Adj_R1 *adj;
    adj = (MFace_Adj_R1 *) f->adj;
    return List_Contains(adj->fvertices,v);
  }

#ifdef __cplusplus
}
#endif
