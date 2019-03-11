/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <UnitTest++.h>

#include "MSTK.h"

// Test if we can collapse degenerate edges correctly in a 3D mesh
// Here we are not worrying about correct classification of entities with 
// respect to a model, only if the mesh topology is correct at the end

TEST(ME_MultiCollapse) {
  Mesh_ptr mesh;

  MSTK_Init();

  mesh = MESH_New(F1);
  int ok = MESH_ImportFromExodusII(mesh, "serial/mixed_pinchout.exo", NULL,
                                   NULL);
  CHECK_EQUAL(1, ok);

  CHECK_EQUAL(1, MESH_CheckTopo(mesh));

  int idx = 0;
  MEdge_ptr me;
  while ((me = MESH_Next_Edge(mesh, &idx))) {
    if (ME_LenSqr(me) < 1.0e-12) {
      MVertex_ptr ev0 = ME_Vertex(me, 0);
      MVertex_ptr ev1 = ME_Vertex(me, 1);
      int gdim0 = MV_GEntDim(ev0);
      int gdim1 = MV_GEntDim(ev1);
      int gid0 = MV_GlobalID(ev0);
      int gid1 = MV_GlobalID(ev1);

      MVertex_ptr vkeep, vdel;
      if (gid0 < gid1) {
        vkeep = ev0;
        vdel = ev1;
      } else {
        vkeep = ev1;
        vdel = ev0;
      }
    
      int topoflag = 0;  // don't worry about violation of consistency with some
                         // some notion of a geometric model
      List_ptr deleted_ents = NULL, merged_entity_pairs = NULL;
      vkeep = ME_Collapse(me, vkeep, topoflag, &deleted_ents,
                          &merged_entity_pairs);

      if (!vkeep) {
        vkeep = vdel;
        vdel = (vkeep == ev0) ? ev1 : ev0;
        vkeep = ME_Collapse(me, vkeep, topoflag, &deleted_ents,
                            &merged_entity_pairs);
      }
      CHECK(vkeep != NULL);
      
      if (deleted_ents)
        List_Delete(deleted_ents);
      if (merged_entity_pairs)
        List_Delete(merged_entity_pairs);
    }
  }

  CHECK_EQUAL(1, MESH_CheckTopo(mesh));

  MESH_Delete(mesh);
}
