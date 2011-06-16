#include <UnitTest++.h>

#include "../include/MSTK.h"

TEST(MR_Edges_Tet_FN)
{
  int ok, i;
  Mesh_ptr mesh;
  MRegion_ptr mr;
  List_ptr redges;

  int exp_ne = 6;
  int exp_eids[6] = {5,3,6,1,2,4};

  MSTK_Init();

  mesh = MESH_New(UNKNOWN_REP);

  ok = MESH_InitFromFile(mesh,"serial/onetet.mstk");
  CHECK_EQUAL(1,ok);

  CHECK_EQUAL(F1,MESH_RepType(mesh));

  CHECK_EQUAL(1,MESH_Num_Regions(mesh));

  CHECK_EQUAL(exp_ne,MESH_Num_Edges(mesh));

  mr = MESH_Region(mesh,0);
  redges = MR_Edges(mr);

  for (i = 0; i < exp_ne; i++) {
    MEdge_ptr e = List_Entry(redges,i);
    int eid = ME_ID(e);
    CHECK_EQUAL(exp_eids[i],eid);
  }
  
  List_Delete(redges);
}
