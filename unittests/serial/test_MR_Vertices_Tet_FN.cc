#include <UnitTest++.h>

#include "../../include/MSTK.h"

TEST(MR_Vertices_Tet_FN)
{
  int ok, i, nv, vids[4];
  Mesh_ptr mesh;
  MRegion_ptr mr;
  List_ptr rverts;

  int exp_nv = 4;
  int exp_vids[4] = {2,4,3,1};

  MSTK_Init();

  mesh = MESH_New(UNKNOWN_REP);

  ok = MESH_InitFromFile(mesh,"serial/onetet.mstk",NULL);
  CHECK_EQUAL(1,ok);

  CHECK_EQUAL(F1,MESH_RepType(mesh));

  CHECK_EQUAL(1,MESH_Num_Regions(mesh));

  CHECK_EQUAL(exp_nv,MESH_Num_Vertices(mesh));

  mr = MESH_Region(mesh,0);

  /* Check if we can get the correct vertices of a tet */

  rverts = MR_Vertices(mr);
  
  for (i = 0; i < exp_nv; i++) {
    MVertex_ptr v = List_Entry(rverts,i);
    int vid = MV_ID(v);
    CHECK_EQUAL(exp_vids[i],vid);
  }
  
  List_Delete(rverts);

  /* Check if we get the correct result when we get the IDs directly */

  MR_VertexIDs(mr,&nv,vids);
  
  CHECK_EQUAL(nv,exp_nv);
  CHECK_ARRAY_EQUAL(exp_vids,vids,exp_nv);
}
