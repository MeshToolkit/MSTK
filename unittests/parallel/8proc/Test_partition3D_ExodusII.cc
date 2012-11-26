#include <UnitTest++.h>

#include "../../../include/MSTK.h"

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
  opts[3] = 1;  // use zoltan as the partitioner if available

  status = MESH_ImportFromFile(mesh,"parallel/8proc/hex_3x3x3_ss.exo","exo",opts,comm);

  CHECK(status);

  return;
}
}
