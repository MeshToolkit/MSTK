/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include "MSTK.h"

#include <vector>
#include <UnitTest++.h>


// Test if we can collapse an edge correctly in a 2D and 3D mesh

TEST(ME_MultiCollapse_Par) 
{
  Mesh_ptr mesh0, mesh;
  List_ptr deleted_ents_all = List_New(10);
  List_ptr merged_entity_pairs_all = List_New(10);

  MSTK_Init();

  mesh = MESH_New(F1);
  int opts[6] = {1, 1, 0, 0, 0, 0};  // weave meshes together with 1 ghost layer
  int ok = MESH_ImportFromNemesisI(mesh, "parallel/4proc/mixed_pinchout.par",
                                   opts, MPI_COMM_WORLD);
  CHECK_EQUAL(1, ok);

  int idx = 0;
  MEdge_ptr me;
  int entities_deleted;
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
    
      int topoflag = 0;  // Don't worry about messing up consistency of the
                         // mesh classification w.r.t. a model
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
      
      if (deleted_ents) {
        List_Cat(deleted_ents_all, deleted_ents);
        List_Delete(deleted_ents);

        entities_deleted = 1;
      }
      if (merged_entity_pairs) {
        List_Cat(merged_entity_pairs_all, merged_entity_pairs);
        List_Delete(merged_entity_pairs);
      }
    }
  }

  CHECK_EQUAL(1, MESH_CheckTopo(mesh));

  int numprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  if (numprocs > 1) {
    // Now we have to update master entity info across processors

    std::vector<int> merged_ents_info;

    int nmerged = List_Num_Entries(merged_entity_pairs_all)/2;
    for (int j = 0; j < nmerged; j++) {
      MEntity_ptr delent = List_Entry(merged_entity_pairs_all, 2*j);
      MEntity_ptr keepent = List_Entry(merged_entity_pairs_all, 2*j+1);
      merged_ents_info.push_back(static_cast<int>(MEnt_Dim(keepent)));
      merged_ents_info.push_back(MEnt_GlobalID(delent));
      merged_ents_info.push_back(MEnt_GlobalID(keepent));
    }
    
    int *nmerged_proc = new int[numprocs];
    int *nmerged_proc_x3 = new int[numprocs];
    MPI_Allgather(&nmerged, 1, MPI_INT, nmerged_proc, 1, MPI_INT,
                  MPI_COMM_WORLD);
    
    int *offset = new int[numprocs];
    int nmerged_global = 0;
    for (int p = 0; p < numprocs; p++) {
      offset[p] = 3*nmerged_global;
      nmerged_global += nmerged_proc[p];
      nmerged_proc_x3[p] = 3*nmerged_proc[p];
    }
    
    // We probably can make this more efficient by using point-to-point
    // communication
    
    int *merged_ents_info_global = new int[3*nmerged_global];
    MPI_Allgatherv(&(merged_ents_info[0]), 3*nmerged, MPI_INT,
                   merged_ents_info_global, nmerged_proc_x3, offset,
                   MPI_INT, MPI_COMM_WORLD);
    

    idx = 0;
    MVertex_ptr vertex;
    while ((vertex = MESH_Next_Vertex(mesh, &idx))) {
      if (MV_PType(vertex) == PGHOST) {
        int vgid = MV_GlobalID(vertex);
        for (int i = 0; i < nmerged_global; i++) {
          if (merged_ents_info_global[3*i] == MVERTEX &&
              merged_ents_info_global[3*i+1] == vgid) {
            // Found vertex that got deleted and replaced by another vtx
            // on a different proc
            MV_Set_GlobalID(vertex, merged_ents_info_global[3*i+2]);
            break;
          }
        }
      }
    }
    
    
    idx = 0;
    MEdge_ptr edge;
    while ((edge = MESH_Next_Edge(mesh, &idx))) {
      if (ME_PType(edge) == PGHOST) {
        int egid = ME_GlobalID(edge);
        for (int i = 0; i < nmerged_global; i++) {
          if (merged_ents_info_global[3*i] == MEDGE &&
              merged_ents_info_global[3*i+1] == egid) {
            // Found edge that got deleted and replaced by another edge
            // on a different proc
            ME_Set_GlobalID(edge, merged_ents_info_global[3*i+2]);
            break;
          }
        }
      }
    }
    
    
    idx = 0;
    MFace_ptr face;
    while ((face = MESH_Next_Face(mesh, &idx))) {
      if (MF_PType(face) == PGHOST) {
        int fgid = MF_GlobalID(face);
        for (int i = 0; i < nmerged_global; i++) {
          if (merged_ents_info_global[3*i] == MFACE &&
              merged_ents_info_global[3*i+1] == fgid) {
            // Found face that got deleted and replaced by another face
            // on a different proc
            MF_Set_GlobalID(face, merged_ents_info_global[3*i+2]);
            break;
          }
        }
      }
    }
    
    delete [] nmerged_proc;
    delete [] nmerged_proc_x3;
    delete [] merged_ents_info_global;
    delete [] offset;

    // ME_Collapse only marked these entities as DELETED but now
    // delete them for good
    idx = 0;
    MEntity_ptr delent;
    while ((delent = List_Next_Entry(deleted_ents_all, &idx)))
      MEnt_Delete(delent, 0);
    
    List_Delete(deleted_ents_all);
    List_Delete(merged_entity_pairs_all);
    
    // Now renumber global IDs to make them contiguous
    
    if (entities_deleted)
      MESH_Renumber_GlobalIDs(mesh, MALLTYPE, 0, NULL, MPI_COMM_WORLD);

    CHECK_EQUAL(1, MESH_Parallel_Check(mesh, MPI_COMM_WORLD));
  }

  MESH_Delete(mesh);
}
