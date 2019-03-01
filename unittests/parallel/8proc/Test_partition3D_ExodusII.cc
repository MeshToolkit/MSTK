#include <UnitTest++.h>

#include "MSTK.h"

SUITE(Parallel) {

TEST(Partition3D_sym_1ring_exo) {
  
  int opts[10], status, nproc, rank;
  Mesh_ptr mesh;
  
  MSTK_Init();
  MSTK_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_size(comm,&nproc);
  MPI_Comm_rank(comm,&rank);

  mesh = MESH_New(F1);

  opts[0] = 1;  // Partition the input mesh
  opts[1] = 0;  // use the default method for distributing mesh
  opts[2] = 1;  // build 1 layer of ghosts
  opts[3] = 0;  // use Metis as the partitioner if available

  status = MESH_ImportFromFile(mesh,"parallel/8proc/hex_3x3x3_ss.exo",NULL,opts,comm);

  CHECK(status);

  return;
}

TEST(ParRead_3D_exo) {
  
  int opts[10], status, nproc, rank;
  Mesh_ptr mesh1, mesh2;
  
  MSTK_Init();
  MSTK_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_size(comm,&nproc);
  MPI_Comm_rank(comm,&rank);

  mesh1 = MESH_New(F1);

  opts[0] = 1;  // Partition the input mesh
  opts[2] = 1;  // build 1 layer of ghosts
  opts[3] = 0;  // use metis as the partitioner if available

  opts[1] = 0;  // Read on processor 0, partition and distribute

  status = MESH_ImportFromFile(mesh1,"parallel/8proc/hex_3x3x3_ss.exo",NULL,
			       opts,comm);

  CHECK(status);

  mesh2 = MESH_New(F1);

  opts[1] = 1;  // Partition the graph and read portions of mesh

  status = MESH_ImportFromFile(mesh2,"parallel/8proc/hex_3x3x3_ss.exo",NULL,
			       opts,comm);

  CHECK(status);


  // Since we are partitioning the same graph either from a full mesh
  // or from sparse elemental information read from the file, we
  // should get the same counts (we cannot guarantee the numbering though)

  CHECK(MESH_Num_Vertices(mesh1) == MESH_Num_Vertices(mesh2));
  CHECK(MESH_Num_Edges(mesh1) == MESH_Num_Edges(mesh2));
  CHECK(MESH_Num_Faces(mesh1) == MESH_Num_Faces(mesh2));
  CHECK(MESH_Num_Regions(mesh1) == MESH_Num_Regions(mesh2));

  CHECK(MESH_Num_GhostVertices(mesh1) == MESH_Num_GhostVertices(mesh2));
  CHECK(MESH_Num_GhostEdges(mesh1) == MESH_Num_GhostEdges(mesh2));
  CHECK(MESH_Num_GhostFaces(mesh1) == MESH_Num_GhostFaces(mesh2));
  CHECK(MESH_Num_GhostRegions(mesh1) == MESH_Num_GhostRegions(mesh2));

  CHECK(MESH_Num_OverlapVertices(mesh1) == MESH_Num_OverlapVertices(mesh2));
  CHECK(MESH_Num_OverlapEdges(mesh1) == MESH_Num_OverlapEdges(mesh2));
  CHECK(MESH_Num_OverlapFaces(mesh1) == MESH_Num_OverlapFaces(mesh2));
  CHECK(MESH_Num_OverlapRegions(mesh1) == MESH_Num_OverlapRegions(mesh2));

  return;
}
}
