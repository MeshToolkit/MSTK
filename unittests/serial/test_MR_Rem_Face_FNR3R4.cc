#include <UnitTest++.h>

#include "../include/MSTK.h"

// Test if we can remove a face correctly from a region (check list of
// remaining faces and their directions)

TEST(MR_Rem_Face_FNR3R4)
{
  int ok, i, j, nrf;
  Mesh_ptr mesh;
  MRegion_ptr mr;
  List_ptr rflist;
  MFace_ptr rfaces[6];
  int rfdirs[6];


  MSTK_Init();

  mesh = MESH_New(UNKNOWN_REP);

  ok = MESH_InitFromFile(mesh,"serial/onehex.mstk");
  CHECK_EQUAL(1,ok);

  mr = MESH_Region(mesh,0);
  rflist = MR_Faces(mr);

  for (i = 0; i < 6; i++) {
    rfaces[i] = List_Entry(rflist,i);
    rfdirs[i] = MR_FaceDir_i(mr,i);
  }

  List_Delete(rflist);

  // Remove faces from region 1 by 1 and see if 
  // the leftover faces have the right IDs and dirs

  for (i = 0; i < 6; i++) {
    MR_Rem_Face(mr,rfaces[i]);
    
    rflist = MR_Faces(mr);
    nrf = List_Num_Entries(rflist);

    CHECK_EQUAL(6-(i+1),nrf);  // Right number of faces?

    for (j = 0; j < nrf; j++) {
      CHECK_EQUAL(rfaces[i+1+j],List_Entry(rflist,j));
      CHECK_EQUAL(rfdirs[i+1+j],MR_FaceDir_i(mr,j));
    }                  

    List_Delete(rflist);
  }

}
