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
    downadj->fedges = List_New(n);
    
    for (i = 0; i < n; i++) {
      downadj->edirs = downadj->edirs | (dir[i] << i);
      List_Add(downadj->fedges,e[i]);
      ME_Add_Face(e[i],f);
    }
  }

  void MF_Replace_Edge_i_FN(MFace_ptr f, int i, int nnu, MEdge_ptr *nuedges, int *nudirs) {
    MFace_DownAdj_FN *downadj;
    MEdge_ptr olde;
    int j, k, dir;

    downadj = (MFace_DownAdj_FN *) f->downadj;

    if (downadj->ne == 0)
      MSTK_Report("MF_Replace_Edge_i","No initial set of edges for face",
		  ERROR);

    olde = List_Entry(downadj->fedges,i);
    downadj->edirs = (downadj->edirs & ~(1<<i)); /* set bit i to 0 */

    /* First one is easy - Just replace the old edge with the new one */
    List_Replacei(downadj->fedges,i,nuedges[0]);
    downadj->edirs = (downadj->edirs | (nudirs[0]<<i)); /* set to dir */

    /* Insert the rest of the nuedges at the right place */
    for (j = 1; j < nnu; j++)
      List_Inserti(downadj->fedges,nuedges[j],i+j);

    /* Move the direction bits after the i'th bit to the left by
       (nnu-1) spaces and insert the direction bits for the rest of
       the inserted edges */

    for (k = downadj->ne-1; k > i; k--) {

      /* set bit (k+nnu-1) to 0 */
      downadj->edirs = downadj->edirs & ~(1<<(k+nnu-1));

      /* get bit k */
      dir = (downadj->edirs>>k) & 1;

      /* move bit k to (k+nnu-1)'th position */
      downadj->edirs = downadj->edirs | (dir<<(k+nnu-1));

      /* set bit k to 0 */
      downadj->edirs = downadj->edirs & ~(1<<k);
    }

    /* set bit j according to input */
    for (j = 1; j < nnu; j++)
      downadj->edirs = downadj->edirs | (nudirs[j]<<(i+j));

    /* One edge replaced an existing edge, others were added */
    downadj->ne += nnu-1; 

    /* Update upward adjacencies */

    ME_Rem_Face(olde,f);
    for (i = 0; i < nnu; i++)
      ME_Add_Face(nuedges[i],f);
  }

  void MF_Replace_Edge_FN(MFace_ptr f, MEdge_ptr e, int nnu, MEdge_ptr *nuedges, int *nudirs) {
    int i, found = 0;
    MFace_DownAdj_FN *downadj;
    MEdge_ptr olde;

    downadj = (MFace_DownAdj_FN *) f->downadj;


    if (downadj->ne == 0)
      MSTK_Report("MF_Replace_Edge_F1","No initial set of edges for face",
		  ERROR);
      
    for (i = 0; i < downadj->ne; i++)
      if ((olde = List_Entry(downadj->fedges,i)) == e) {
	found = 1;
	break;
      }
    if (!found) {
      MSTK_Report("MF_Replace_Edge","Edge not found in face",ERROR);
      return;
    }

    MF_Replace_Edge_i_FN(f, i, nnu, nuedges, nudirs);
  }

  int MF_Num_Edges_FN(MFace_ptr f) {
    return ((MFace_DownAdj_FN *) f->downadj)->ne;
  }	

  List_ptr MF_Edges_FN(MFace_ptr f, int dir, MVertex_ptr v0) {
    int i, k, n, fnd, edir;
    List_ptr fedges;
    MEdge_ptr e;
    MFace_DownAdj_FN *downadj;

    downadj = (MFace_DownAdj_FN *) f->downadj;


    if (v0 == NULL) {
      if (dir)
	return List_Copy(downadj->fedges);
      else {
	n = downadj->ne;
	fedges = List_New(n);
	for (i = n-1; i >= 0; i--)
	  List_Add(fedges,List_Entry(downadj->fedges,i));
	return fedges;
      }
    }
    else {
      n = downadj->ne;
      fnd = 0;
      for (i = 0; i < n; i++) {
	e = List_Entry(downadj->fedges,i);
	edir = ((downadj->edirs)>>i) & 1;
	if (ME_Vertex(e,edir^dir) == v0) {
	  fnd = 1;
	  k = i;
	  break;
	}
      }

      if (!fnd)
	MSTK_Report("MF_Edges_F1","Cannot find vertex in face!!",FATAL);
	
      fedges = List_New(n);
      for (i = 0; i < n; i++) {
	e = dir ? List_Entry(downadj->fedges,(k+i)%n) :
	  List_Entry(downadj->fedges,(k+n-i)%n);
	List_Add(fedges,e);
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
      if (List_Entry(downadj->fedges,i) == e)
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

    return List_Contains(downadj->fedges,e);
  }

#ifdef __cplusplus
}
#endif
