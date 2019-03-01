#include <UnitTest++.h>

#include "MSTK.h"

TEST(MR_Edges_polyhed_FN)
{
  int ok, i, ne, eids[15];
  Mesh_ptr mesh;
  MRegion_ptr mr;
  List_ptr redges;

  int exp_ne = 15;
  int exp_eids[15] = {6,3,11,5,2,14,15,12,4,7,8,10,13,1,9};

  MSTK_Init();

  mesh = MESH_New(UNKNOWN_REP);

  ok = MESH_InitFromFile(mesh,"serial/onepolyhed.mstk",NULL);
  CHECK_EQUAL(1,ok);

  CHECK_EQUAL(F1,MESH_RepType(mesh));

  CHECK_EQUAL(1,MESH_Num_Regions(mesh));

  CHECK_EQUAL(exp_ne,MESH_Num_Edges(mesh));

  mr = MESH_Region(mesh,0);

  /* Check if we can retrieve the edges of a polyhedron */

  redges = MR_Edges(mr);

  for (i = 0; i < exp_ne; i++) {
    MEdge_ptr e = List_Entry(redges,i);
    int eid = ME_ID(e);
    CHECK_EQUAL(exp_eids[i],eid);
  }
  
  List_Delete(redges);

  /* Check if we can directly retrieve the IDs of edges of a polyhedron */

  MR_EdgeIDs(mr,&ne,eids);
  
  CHECK_EQUAL(ne,exp_ne);
  CHECK_ARRAY_EQUAL(exp_eids,eids,exp_ne);

}
