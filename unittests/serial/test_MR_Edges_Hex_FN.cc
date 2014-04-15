#include <UnitTest++.h>

#include "../../include/MSTK.h"

TEST(MR_Edges_Hex_FN)
{
  int ok, i, ne, eids[12];
  Mesh_ptr mesh;
  MRegion_ptr mr;
  List_ptr redges;

  int exp_ne = 12;
  int exp_eids[12] = {8,7,6,3,
		      5,1,2,12,
		      4,10,9,11};

  MSTK_Init();

  mesh = MESH_New(UNKNOWN_REP);

  ok = MESH_InitFromFile(mesh,"serial/onehex.mstk",NULL);
  CHECK_EQUAL(1,ok);

  CHECK_EQUAL(F1,MESH_RepType(mesh));

  CHECK_EQUAL(1,MESH_Num_Regions(mesh));

  CHECK_EQUAL(exp_ne,MESH_Num_Edges(mesh));

  mr = MESH_Region(mesh,0);

  /* Check if we can retrieve edges of a hex */

  redges = MR_Edges(mr);

  for (i = 0; i < exp_ne; i++) {
    MEdge_ptr e = List_Entry(redges,i);
    int eid = ME_ID(e);
    CHECK_EQUAL(exp_eids[i],eid);
  }
  
  List_Delete(redges);

  /* Check if we can directly retrieve the IDs of edges of a hex */

  MR_EdgeIDs(mr,&ne,eids);
  
  CHECK_EQUAL(ne,exp_ne);
  CHECK_ARRAY_EQUAL(exp_eids,eids,exp_ne);
}
