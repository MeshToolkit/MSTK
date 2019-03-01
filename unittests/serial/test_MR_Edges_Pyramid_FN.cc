/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <UnitTest++.h>

#include "MSTK.h"

TEST(MR_Edges_Pyramid_FN)
{
  int ok, i, ne, eids[8];
  Mesh_ptr mesh;
  MRegion_ptr mr;
  List_ptr redges;

  int exp_ne = 8;
  int exp_eids[8] = {3,4,8,7,5,2,1,6};

  MSTK_Init();

  mesh = MESH_New(UNKNOWN_REP);

  ok = MESH_InitFromFile(mesh,"serial/onepyramid.mstk",NULL);
  CHECK_EQUAL(1,ok);

  CHECK_EQUAL(F1,MESH_RepType(mesh));

  CHECK_EQUAL(1,MESH_Num_Regions(mesh));

  CHECK_EQUAL(exp_ne,MESH_Num_Edges(mesh));

  mr = MESH_Region(mesh,0);

  /* Check if we can retrieve edges of a pyramid */

  redges = MR_Edges(mr);

  for (i = 0; i < exp_ne; i++) {
    MEdge_ptr e = List_Entry(redges,i);
    int eid = ME_ID(e);
    CHECK_EQUAL(exp_eids[i],eid);
  }
  
  List_Delete(redges);

  /* Check if we can directly retrieve the IDs of edges of a pyramid */

  MR_EdgeIDs(mr,&ne,eids);
  
  CHECK_EQUAL(ne,exp_ne);
  CHECK_ARRAY_EQUAL(exp_eids,eids,exp_ne);

}
