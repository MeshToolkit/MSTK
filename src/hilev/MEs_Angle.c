#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MSTK.h"
#include "MSTK_VecFuncs.h"


#ifdef __cplusplus
extern "C" {
#endif

  /* ***COSINE*** of Angle between two mesh edges if they share a
     vertex. The angle will be the smaller of the two possible
     angles. */

  double MEs_Angle(MEdge_ptr edge1, MEdge_ptr edge2) {
    int i, j, fnd;
    double exyz[2][3], vec1[3], vec2[3], dp;
    MVertex_ptr ev1[2], ev2[2];

    ev1[0] = ME_Vertex(edge1,0);
    ev1[1] = ME_Vertex(edge1,1);
    ev2[0] = ME_Vertex(edge2,0);
    ev2[1] = ME_Vertex(edge2,1);
    
    fnd = 1;
    if (ev1[0] == ev2[0]) {
      i = 0; j = 0;
    }
    else if (ev1[0] == ev2[1]) {
      i = 0; j = 1;
    }
    else if (ev1[1] == ev2[0]) {
      i = 1; j = 0;
    }
    else if (ev1[1] == ev2[1]) {
      i = 1; j = 1;
    }
    else {
      MSTK_Report("MEs_Angle","Edges do not share common vertex",ERROR);
      return 0.0;
    }

    /* Edge vector 1 going away from common vertex */

    MV_Coords(ev1[i],exyz[0]);
    MV_Coords(ev1[1-i],exyz[1]);

    MSTK_VDiff3(exyz[1],exyz[0],vec1);
    MSTK_VNormalize3(vec1);

    /* Edge vector 2 going away from common vertex */

    MV_Coords(ev2[j],exyz[0]);
    MV_Coords(ev2[1-j],exyz[1]);
    MSTK_VDiff3(exyz[1],exyz[0],vec2);
    MSTK_VNormalize3(vec2);

    /* Angle between vectors */

    dp = MSTK_VDot3(vec1,vec2);
    
    return dp;
  }

#ifdef __cplusplus
}
#endif
