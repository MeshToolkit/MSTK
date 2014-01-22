#include <UnitTest++.h>

#include "../../include/MSTK.h"

// Test if we can split a face correctly in a 3D mesh

TEST(MF_Split_with_Edge) 
{
  int idx, i, n, gid, gdim, ok;
  unsigned int dir;
  Mesh_ptr mesh;
  MVertex_ptr v0, v1, vnew0, vnew1;
  MEdge_ptr enew, e0, e1;
  MFace_ptr mf, efaces[2];
  MRegion_ptr fr;
  List_ptr fedges, eflist, fregions, rfaces;
  double xyz0[3], xyz1[3], xyznew[3];

  MSTK_Init();

  mesh = MESH_New(UNKNOWN_REP);
  ok = MESH_InitFromFile(mesh,"serial/reghex3D.mstk",NULL);
  CHECK_EQUAL(ok,1);

  CHECK(MESH_Num_Vertices(mesh) > 0);
  
  idx = 0; n = 0;
  while ((mf = MESH_Next_Face(mesh,&idx))) {
    if (MF_GEntDim(mf) == 3) {

      /* Store regions of face before it gets deleted */
      fregions = MF_Regions(mf);
      
      /* Get edges of face */

      fedges = MF_Edges(mf,1,0);

      if (n == 0) {
        /* split two opposite edges at their midpoints and connect 
           the two vertices to split the face */

        e0 = List_Entry(fedges,0);
        v0 = ME_Vertex(e0,0);
        MV_Coords(v0,xyz0);
        v1 = ME_Vertex(e0,1);
        MV_Coords(v1,xyz1);
        for (i = 0; i < 3; i++) xyznew[i] = (xyz0[i]+xyz1[i])/2.0;

        vnew0 = ME_Split(e0,xyznew);

        e1 = List_Entry(fedges,2);
        v0 = ME_Vertex(e1,0);
        MV_Coords(v0,xyz0);
        v1 = ME_Vertex(e1,1);
        MV_Coords(v1,xyz1);
        for (i = 0; i < 3; i++) xyznew[i] = (xyz0[i]+xyz1[i])/2.0;

        vnew1 = ME_Split(e1,xyznew);
      }
      else if (n == 1) {
        /* split one edge and connect it to an opposite vertex to get 
           a splitting edge for the face */

        e0 = List_Entry(fedges,0);
        v0 = ME_Vertex(e0,0);
        MV_Coords(v0,xyz0);
        v1 = ME_Vertex(e0,1);
        MV_Coords(v1,xyz1);
        for (i = 0; i < 3; i++) xyznew[i] = (xyz0[i]+xyz1[i])/2.0;

        vnew0 = ME_Split(e0,xyznew);

        e1 = List_Entry(fedges,2);
        dir = MF_EdgeDir_i(mf,2);
        vnew1 = ME_Vertex(e1,!dir);
      }
      else if (n == 2) {
        /* Create a splitting edge by connecting two non-adjacent vertices */
        e0 = List_Entry(fedges,0);
        dir = MF_EdgeDir_i(mf,0);
        vnew0 = ME_Vertex(e0,!dir);

        e1 = List_Entry(fedges,2);
        dir = MF_EdgeDir_i(mf,2);
        vnew1 = ME_Vertex(e1,!dir);
      }
      else
        continue;
      List_Delete(fedges);

      enew = MF_Split_with_Edge(mf,vnew0,vnew1);
      CHECK(enew);
       
      eflist = ME_Faces(enew);
      CHECK_EQUAL(2,List_Num_Entries(eflist));
      efaces[0] = List_Entry(eflist,0);
      efaces[1] = List_Entry(eflist,1);
      List_Delete(eflist);

      idx = 0;
      while ((fr = List_Next_Entry(fregions,&idx))) {
        rfaces = MR_Faces(fr);
        CHECK(List_Contains(rfaces,efaces[0]));
        CHECK(List_Contains(rfaces,efaces[1]));
        CHECK(!List_Contains(rfaces,mf));
        List_Delete(rfaces);
      }
      List_Delete(fregions);

      n++;
    }
  }
}
