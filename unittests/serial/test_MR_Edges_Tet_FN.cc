#include <UnitTest++.h>

#include "MSTK.h"

TEST(MR_Edges_Tet_FN)
{
  int ok, i, ne, eids[6];
  Mesh_ptr mesh;
  MRegion_ptr mr;
  List_ptr redges;

  int exp_ne = 6;
  int exp_eids[6] = {5,3,6,1,2,4};

  MSTK_Init();

  mesh = MESH_New(UNKNOWN_REP);

  ok = MESH_InitFromFile(mesh,"serial/onetet.mstk",NULL);
  CHECK_EQUAL(1,ok);

  CHECK_EQUAL(F1,MESH_RepType(mesh));

  CHECK_EQUAL(1,MESH_Num_Regions(mesh));

  CHECK_EQUAL(exp_ne,MESH_Num_Edges(mesh));

  mr = MESH_Region(mesh,0);

  /* Check if we can retrieve edges of a tet */

  redges = MR_Edges(mr);

  for (i = 0; i < exp_ne; i++) {
    MEdge_ptr e = List_Entry(redges,i);
    int eid = ME_ID(e);
    CHECK_EQUAL(exp_eids[i],eid);
  }
  
  List_Delete(redges);

  /* Check if we can directly retrieve the IDs of edges of a tet */

  MR_EdgeIDs(mr,&ne,eids);
  
  CHECK_EQUAL(ne,exp_ne);
  CHECK_ARRAY_EQUAL(exp_eids,eids,exp_ne);
}
