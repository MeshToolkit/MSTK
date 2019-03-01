/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <UnitTest++.h>

#include "MSTK.h"

TEST(MR_Vertices_Hex_FN)
{
  int ok, i, nv, vids[8];
  Mesh_ptr mesh;
  MRegion_ptr mr;
  List_ptr rverts;

  int exp_nv = 8;
  int exp_vids[8] = {6,4,7,3,5,1,8,2};

  MSTK_Init();

  mesh = MESH_New(UNKNOWN_REP);

  ok = MESH_InitFromFile(mesh,"serial/onehex.mstk",NULL);
  CHECK_EQUAL(1,ok);

  CHECK_EQUAL(F1,MESH_RepType(mesh));

  CHECK_EQUAL(1,MESH_Num_Regions(mesh));

  CHECK_EQUAL(exp_nv,MESH_Num_Vertices(mesh));

  mr = MESH_Region(mesh,0);

  /* Check if we can get the correct vertices of a hex */

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
