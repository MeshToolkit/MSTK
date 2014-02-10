#include <UnitTest++.h>

#include "../../include/MSTK.h"

// Test if we can correctly read a polygonal exodus mesh 
// Then test if we can correctly write it and read it back

TEST(Read_Write_ExodusII_Poly2) 
{
  int idx, ok, nfe;
  Mesh_ptr mesh, mesh2;
  MFace_ptr mf;
  int one_tri, one_quad, one_penta;

  MSTK_Init();

  mesh = MESH_New(UNKNOWN_REP);
  ok = MESH_ImportFromFile(mesh,"serial/poly2.exo","exo",NULL,NULL);
  CHECK_EQUAL(ok,1);

  CHECK(MESH_Num_Vertices(mesh) > 0);
  
  idx = 0; one_tri = 0; one_quad = 0; one_penta = 0;
  while ((mf = MESH_Next_Face(mesh,&idx))) {
    int ne = MF_Num_Edges(mf);
    if (ne == 3)
      one_tri++;
    else if (ne == 4)
      one_quad++;
    else if (ne == 5)
      one_penta++;
  }

  CHECK_EQUAL(one_tri,1);
  CHECK_EQUAL(one_quad,1);
  CHECK_EQUAL(one_penta,1);

  ok = MESH_ExportToFile(mesh,"./poly2-tmp.exo","exo",0,NULL,NULL,NULL);

  CHECK_EQUAL(ok,1);

  mesh2 = MESH_New(UNKNOWN_REP);
  ok = MESH_ImportFromFile(mesh2,"./poly2-tmp.exo","exo",NULL,NULL);

  idx = 0; one_tri = 0; one_quad = 0; one_penta = 0;
  while ((mf = MESH_Next_Face(mesh2,&idx))) {
    int ne = MF_Num_Edges(mf);
    if (ne == 3)
      one_tri++;
    else if (ne == 4)
      one_quad++;
    else if (ne == 5)
      one_penta++;
  }

  CHECK_EQUAL(one_tri,1);
  CHECK_EQUAL(one_quad,1);
  CHECK_EQUAL(one_penta,1);
}
