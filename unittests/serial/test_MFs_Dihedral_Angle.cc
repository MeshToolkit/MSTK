#include <UnitTest++.h>

#include "../../include/MSTK.h"

TEST(MFs_DihedralAngle) {


  MSTK_Init();

  Mesh_ptr mesh = MESH_New(UNKNOWN_REP);
  
  int ok = MESH_InitFromFile(mesh,"serial/twohex.mstk",NULL);
  CHECK_EQUAL(1,ok);


  // Get the top center edge of the mesh

  MVertex_ptr mv0 = MESH_VertexFromID(mesh,8);
  MVertex_ptr mv1 = MESH_VertexFromID(mesh,11);
  MEdge_ptr me = MVs_CommonEdge(mv0,mv1);
  CHECK(me);

  List_ptr efaces = ME_Faces(me);
  CHECK_EQUAL(3,List_Num_Entries(efaces));

  MFace_ptr f0=NULL, f1=NULL;
  for (int i = 0; i < 3; i++) {
    MFace_ptr ef = List_Entry(efaces,i);
    if (MF_GEntDim(ef) == 3) continue;
    if (f0)
      f1 = ef;
    else
      f0 = ef;
  }
  List_Delete(efaces);

  // The cosine of the dihedral angle between these two faces (being
  // the top horizontal faces) should be -1.0 (180 degrees angle)

  CHECK_CLOSE(-1.0,MFs_DihedralAngle(f0,f1,me),1.e-02);
  CHECK_CLOSE(-1.0,MFs_DihedralAngle(f1,f0,me),1.e-02);


  // Now get an edge bottom left corner

  mv0 = MESH_VertexFromID(mesh,1);
  mv1 = MESH_VertexFromID(mesh,4);
  me = MVs_CommonEdge(mv0,mv1);
  CHECK(me);

  efaces = ME_Faces(me);
  CHECK_EQUAL(2,List_Num_Entries(efaces));
  f0 = List_Entry(efaces,0);
  f1 = List_Entry(efaces,1);
  List_Delete(efaces);

  // The cosine of the dihedral angle between these two faces (being at
  // a corner) should be 0.0 (90 or 270 degrees angle)

  CHECK_CLOSE(0.0,MFs_DihedralAngle(f0,f1,me),1.0e-02);
  CHECK_CLOSE(0.0,MFs_DihedralAngle(f1,f0,me),1.0e-02);

  MESH_Delete(mesh);

}
