#define _H_MVertex_Private

#include "MVertex.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Functions */


  List_ptr MV_Edges_R2R4(MVertex_ptr v) {
    int i;
    MVertex_ptr adjv;
    List_ptr vedges;
    MVertex_SameAdj_R2R4 *sameadj;

    /* Have to create volatile edges */

    vedges = List_New(sameadj->nvadj);
    for (i = 0; i < sameadj->nvadj; i++) {
      /* create a volatile edge from v and the adjacent vertex */
      adjv = List_Entry(sameadj->adjverts,i);
    }

    return vedges;
  }

#ifdef __cplusplus
}
#endif
