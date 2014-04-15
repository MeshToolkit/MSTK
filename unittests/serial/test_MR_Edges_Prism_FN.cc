#include <UnitTest++.h>

#include "../../include/MSTK.h"

TEST(MR_Edges_Prism_FN)
{
  int ok, i, ne, eids[9];
  Mesh_ptr mesh;
  MRegion_ptr mr;
  List_ptr redges;

  int exp_ne = 9;
  int exp_eids[9] = {8, 4, 3,
		     1, 5, 7,
		     2, 6, 9};

  MSTK_Init();

  mesh = MESH_New(UNKNOWN_REP);

  ok = MESH_InitFromFile(mesh,"serial/oneprism.mstk",NULL);
  CHECK_EQUAL(1,ok);

  CHECK_EQUAL(F1,MESH_RepType(mesh));

  CHECK_EQUAL(1,MESH_Num_Regions(mesh));

  CHECK_EQUAL(exp_ne,MESH_Num_Edges(mesh));

  mr = MESH_Region(mesh,0);

  /* Check if acn retrieve edges of a prism */

  redges = MR_Edges(mr);

  for (i = 0; i < exp_ne; i++) {
    MEdge_ptr e = List_Entry(redges,i);
    int eid = ME_ID(e);
    CHECK_EQUAL(exp_eids[i],eid);
  }
  
  List_Delete(redges);

  /* Check if we can directly retrieve the IDs of edges of a prism */

  MR_EdgeIDs(mr,&ne,eids);
  
  CHECK_EQUAL(ne,exp_ne);
  CHECK_ARRAY_EQUAL(exp_eids,eids,exp_ne);
  
}
