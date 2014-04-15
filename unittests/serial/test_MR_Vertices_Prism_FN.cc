#include <UnitTest++.h>

#include "../../include/MSTK.h"

TEST(MR_Vertices_Prism_FN)
{
  int ok, i, nv, vids[6];
  Mesh_ptr mesh;
  MRegion_ptr mr;
  List_ptr rverts;

  int exp_nv = 6;
  int exp_vids[6] = {2,6,3,5,1,4};

  MSTK_Init();

  mesh = MESH_New(UNKNOWN_REP);

  ok = MESH_InitFromFile(mesh,"serial/oneprism.mstk",NULL);
  CHECK_EQUAL(1,ok);

  CHECK_EQUAL(F1,MESH_RepType(mesh));

  CHECK_EQUAL(1,MESH_Num_Regions(mesh));

  CHECK_EQUAL(exp_nv,MESH_Num_Vertices(mesh));

  mr = MESH_Region(mesh,0);

  /* Check if we can get the correct vertices of a prism */

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
