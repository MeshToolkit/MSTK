#include <UnitTest++.h>

#include "../../../include/MSTK.h"

SUITE(Parallel) {

TEST(Partition2D_1ring_exo) {
  
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

  status = MESH_ImportFromFile(mesh,"parallel/4proc/mesh5x5-skewed.exo",NULL,opts,comm);

  CHECK(status);

  return;
}

TEST(ParRead_2D_exo) {
  
  int opts[10], status, nproc, rank;
  Mesh_ptr mesh;
  
  MSTK_Init();
  MSTK_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_size(comm,&nproc);
  MPI_Comm_rank(comm,&rank);

  mesh = MESH_New(F1);

  opts[0] = 1;  // Partition the input mesh
  opts[1] = 1;  // Partition the graph and read portions of mesh
  opts[2] = 1;  // build 1 layer of ghosts
  opts[3] = 0;  // use metis as the partitioner if available

  status = MESH_ImportFromFile(mesh,"parallel/4proc/mesh5x5-skewed.exo",NULL,opts,comm);

  CHECK(status);

  return;
}
}
