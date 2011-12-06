#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MSTK.h"
#include "MSTK_VecFuncs.h"


#ifdef __cplusplus
extern "C" {
#endif

  /* ***COSINE*** of Dihedral Angle between two mesh faces if they
     share an edge. The dihedral angle will be the smaller of the two
     possible angles. Clearly, the result will only be an estimate if
     either of the faces is not planar */

  double MFs_DihedralAngle(MFace_ptr face1, MFace_ptr face2, MEdge_ptr edge) {
    int i, fnd, nfe1, fedir, nfv1, nfv2;
    double fxyz[MAXPV2][3], vec1[3], vec2[3];
    double normal1[3], normal2[3], dp, mid[3];
    MVertex_ptr fv, ev0, ev1;
    List_ptr fedges1, fedges2, fverts1, fverts2;
    

    if (!edge) {
      fedges1 = MF_Edges(face1,1,0);      
      nfe1 = List_Num_Entries(fedges1);

      fedges2 = MF_Edges(face2,1,0);

      for (i = 0, fnd = 0; i < nfe1 && !fnd; i++) {
	edge = List_Entry(fedges1,i);
	if (List_Contains(fedges2,edge)) 
	  fnd = 1;	
      }

      List_Delete(fedges1);
      List_Delete(fedges2);

      if (!fnd) {
	MSTK_Report("MFs_DihedralAngle","Faces do not share common edge",MSTK_ERROR);
	return 0.0;
      }
    }


    /* For non-convex faces, picking the three points from which to
       calculate the normal is an issue. We will pick the two points
       of the common edge and the geometric center of the face.

       THIS MAY NOT WORK IF THE GEOMETRIC CENTER IS OUTSIDE THE FACE

       Eventually we will have to find a point in the face such that
       all triangles formed by connecting the point and each of the
       edges have consistent normals */


    ev0 = ME_Vertex(edge,0); 
    ev1 = ME_Vertex(edge,1); 

    /* Normal of face 1 */

    fedir = MF_EdgeDir(face1,edge);

    fverts1 = MF_Vertices(face1,fedir,ev0);
    nfv1 = List_Num_Entries(fverts1);

    for (i = 0; i < nfv1; i++) {
      fv = List_Entry(fverts1,i); 
      MV_Coords(fv,fxyz[i]);
    }
    List_Delete(fverts1);

    MSTK_VDiff3(fxyz[1],fxyz[0],vec1);

    if (nfv1 == 3) { /* Triangles - always convex */
      MSTK_VDiff3(fxyz[2],fxyz[0],vec2);
    }
    else { /* Others - can be non-convex */

      /* use geometric center as third point */
      
      mid[0] = mid[1] = mid[2] = 0.0;
      for (i = 0; i < nfv1; i++) {
	mid[0] += fxyz[i][0];
	mid[1] += fxyz[i][1];
	mid[2] += fxyz[i][2];
      }
      mid[0] /= nfv1; mid[1] /= nfv1; mid[2] /= nfv1;

      MSTK_VDiff3(mid,fxyz[0],vec2);
    }
    
    MSTK_VCross3(vec1,vec2,normal1);
    MSTK_VNormalize3(normal1);
    

    /* Normal of face 2 */

    fedir = MF_EdgeDir(face2,edge);

    fverts2 = MF_Vertices(face2,!fedir,ev1);
    nfv2 = List_Num_Entries(fverts2);

    for (i = 0; i < nfv2; i++) {
      fv = List_Entry(fverts2,i); 
      MV_Coords(fv,fxyz[i]);
    }
    List_Delete(fverts2);

    MSTK_VDiff3(fxyz[1],fxyz[0],vec1);

    if (nfv1 == 3) { /* Triangles - always convex */
      MSTK_VDiff3(fxyz[2],fxyz[0],vec2);
    }
    else { /* Others - can be non-convex */

      /* use geometric center as third point */
      
      mid[0] = mid[1] = mid[2] = 0.0;
      for (i = 0; i < nfv2; i++) {
	mid[0] += fxyz[i][0];
	mid[1] += fxyz[i][1];
	mid[2] += fxyz[i][2];
      }
      mid[0] /= nfv2; mid[1] /= nfv2; mid[2] /= nfv2;

      MSTK_VDiff3(mid,fxyz[0],vec2);
    }
    
    MSTK_VCross3(vec1,vec2,normal2);
    MSTK_VNormalize3(normal2);

    /* We have to negate the second normal. Otherwise we will get a
       dihedral angle of 0 for two faces in a plane when it should
       be 180 */

    MSTK_VNeg3(normal2);
    

    /* Angle between normals */
    dp = MSTK_VDot3(normal1,normal2);
    
    return dp;
  }

#ifdef __cplusplus
}
#endif
