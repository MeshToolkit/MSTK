#define _H_MFace_Private

#include "MFace.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MF_Set_Edges_FN(MFace_ptr f, int n, MEdge_ptr *e, int *dir) {
    int i;
    MFace_DownAdj_FN *downadj;

    downadj = (MFace_DownAdj_FN *) f->downadj;

    downadj->ne = n;
    downadj->edirs = 0;
    downadj->fedges = Set_New(n);
    
    for (i = 0; i < n; i++) {
      downadj->edirs = downadj->edirs | (dir[i] << i);
      Set_Add(downadj->fedges,e[i]);
      ME_Add_Face(e[i],f);
    }
  }

  void MF_Replace_Edge_i_FN(MFace_ptr f, int i, MEdge_ptr e, int dir) {
    MFace_DownAdj_FN *downadj;
    MEdge_ptr olde;

    downadj = (MFace_DownAdj_FN *) f->downadj;

    if (downadj->ne == 0)
      MSTK_Report("MF_Replace_Edge_i","No initial set of edges for face",
		  ERROR);

    olde = Set_Entry(downadj->fedges,i);
      
    downadj->edirs = (downadj->edirs & ~(1<<i)); /* set bit i to 0 */
    downadj->edirs = (downadj->edirs | (dir<<i)); /* set to dir */
    Set_Replacei(downadj->fedges,i,e);

    ME_Rem_Face(olde,f);
    ME_Add_Face(e,f);
  }

  void MF_Replace_Edge_FN(MFace_ptr f, MEdge_ptr e, MEdge_ptr nue, int dir) {
    int i, found = 0;
    MFace_DownAdj_FN *downadj;
    MEdge_ptr olde;

    downadj = (MFace_DownAdj_FN *) f->downadj;


    if (downadj->ne == 0)
      MSTK_Report("MF_Replace_Edge_F1","No initial set of edges for face",
		  ERROR);
      
    for (i = 0; i < downadj->ne; i++)
      if ((olde = Set_Entry(downadj->fedges,i)) == e) {
	found = 1;
	break;
      }
    if (!found) {
      MSTK_Report("MF_Replace_Edge","Edge not found in face",ERROR);
      return;
    }

    downadj->edirs = (downadj->edirs & ~(1<<i)); /* set bit i to 0 */
    downadj->edirs = (downadj->edirs | (dir<<i)); /* set to dir */
    Set_Replacei(downadj->fedges,i,e);

    ME_Rem_Face(olde,f);
    ME_Add_Face(e,f);
  }

  int MF_Num_Edges_FN(MFace_ptr f) {
    return ((MFace_DownAdj_FN *) f->downadj)->ne;
  }	

  Set_ptr MF_Edges_FN(MFace_ptr f, int dir, MVertex_ptr v0) {
    int i, k, n, fnd, edir;
    Set_ptr fedges;
    MEdge_ptr e;
    MFace_DownAdj_FN *downadj;

    downadj = (MFace_DownAdj_FN *) f->downadj;


    if (v0 == NULL) {
      if (dir)
	return Set_Copy(downadj->fedges);
      else {
	n = downadj->ne;
	fedges = Set_New(n);
	for (i = n-1; i >= 0; i--)
	  Set_Add(fedges,Set_Entry(downadj->fedges,i));
	return fedges;
      }
    }
    else {
      n = downadj->ne;
      fnd = 0;
      for (i = 0; i < n; i++) {
	e = Set_Entry(downadj->fedges,i);
	edir = ((downadj->edirs)>>i) & 1;
	if (ME_Vertex(e,edir^dir) == v0) {
	  fnd = 1;
	  k = i;
	  break;
	}
      }

      if (!fnd)
	MSTK_Report("MF_Edges_F1","Cannot find vertex in face!!",FATAL);
	
      fedges = Set_New(n);
      for (i = 0; i < n; i++) {
	e = dir ? Set_Entry(downadj->fedges,(k+i)%n) :
	  Set_Entry(downadj->fedges,(k+n-i)%n);
	Set_Add(fedges,e);
      }	
    }

    return fedges;
  }

  int MF_EdgeDir_FN(MFace_ptr f, MEdge_ptr e) {
    int n, i;
    MFace_DownAdj_FN *downadj;

    downadj = (MFace_DownAdj_FN *) f->downadj;
    
    n = downadj->ne;
    for (i = 0; i < n; i++) {
      if (Set_Entry(downadj->fedges,i) == e)
	return ((downadj->edirs)>>i) & 1;
    }

    MSTK_Report("MF_Edges_FN","Cannot find edge in face!!",FATAL);

    return 0;
  }

  int MF_EdgeDir_i_FN(MFace_ptr f, int i) {
    MFace_DownAdj_FN *downadj;

    downadj = (MFace_DownAdj_FN *) f->downadj;

    return ((downadj->edirs)>>i) & 1;
  }
			
  int MF_UsesEdge_FN(MFace_ptr f, MEdge_ptr e) {
    MFace_DownAdj_FN *downadj;

    downadj = (MFace_DownAdj_FN *) f->downadj;

    return Set_Contains(downadj->fedges,e);
  }

#ifdef __cplusplus
}
#endif
