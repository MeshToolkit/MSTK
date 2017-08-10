#include <UnitTest++.h>

#include "../../../include/MSTK.h"
#include "../../../include/MSTK_private.h"
#include <vector>
#include <algorithm>

SUITE(Parallel) {
TEST(RenumberGIDs_Dist) {
  Mesh_ptr mesh;
  char filename[256];
  int nproc, rank;
  MVertex_ptr mv;
  MEdge_ptr me;
  MFace_ptr mf;
  int dim, idx;

  MSTK_Init();

  MSTK_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm,&nproc);
  MPI_Comm_rank(comm,&rank);

  Mesh_ptr mesh0;
  if (rank == 0) {
    mesh0 = MESH_New(UNKNOWN_REP);

    sprintf(filename,"parallel/4proc/quad6x5.mstk");
    MESH_InitFromFile(mesh0, filename, comm);
  }

//  DebugWait=1;
//  while (DebugWait);

  int ring = 1; /* One ring ghosts */
  int with_attr = 1; /* Do allow exchange of attributes */
  int del_inmesh = 1; /* delete input mesh after partitioning */
  PartitionMethod method;

#if defined (_MSTK_HAVE_METIS)
  method = METIS;
#elif defined (_MSTK_HAVE_ZOLTAN)
  method = ZOLTAN_RCB;
#else
  fprintf(stderr,"Cannot find partitioner\n");
  status = 0;
  CHECK(status);
#endif

  dim = 2;
  mesh = NULL;
  MSTK_Mesh_Distribute(mesh0, &mesh, &dim, ring, with_attr, method, 
		       del_inmesh, comm);


  // Assign some random global IDs to the entities
  srand(time(NULL));

  idx = 0;
  while ((mv = MESH_Next_Vertex(mesh, &idx)))
    MV_Set_GlobalID(mv, 25*MV_GlobalID(mv));

  idx = 0;
  while ((me = MESH_Next_Edge(mesh, &idx)))
    ME_Set_GlobalID(me, 25*ME_GlobalID(me));

  idx = 0;
  while ((mf = MESH_Next_Face(mesh, &idx)))
    MF_Set_GlobalID(mf, 25*MF_GlobalID(mf));


  // Verify that the global IDs are discontinuous

  std::vector<int> vgids, egids, fgids;
  idx = 0;
  while ((mv = MESH_Next_Vertex(mesh, &idx)))
    if (MV_PType(mv) != PGHOST)
      vgids.push_back(MV_GlobalID(mv));
  std::sort(vgids.begin(), vgids.end());
  idx = 0;
  while ((me = MESH_Next_Edge(mesh, &idx)))
    if (ME_PType(me) != PGHOST)
      egids.push_back(ME_GlobalID(me));
  std::sort(egids.begin(), egids.end());
  idx = 0;
  while ((mf = MESH_Next_Face(mesh, &idx)))
    if (MF_PType(mf) != PGHOST)
      fgids.push_back(MF_GlobalID(mf));
  std::sort(fgids.begin(), fgids.end());

  bool vjump_found = false;
  for (int i = 0; i < vgids.size()-1; i++)
    if (vgids[i] != vgids[i+1]-1) {
      vjump_found = true;
      break;
    }
  bool ejump_found = false;
  for (int i = 0; i < egids.size()-1; i++)
    if (egids[i] != egids[i+1]-1) {
      ejump_found = true;
      break;
    }
  bool fjump_found = false;
  for (int i = 0; i < fgids.size()-1; i++)
    if (fgids[i] != fgids[i+1]-1) {
      fjump_found = true;
      break;
    }

  CHECK(vjump_found);
  CHECK(ejump_found);
  CHECK(fjump_found);

  // Call the renumbering routine
  
  MESH_Renumber_GlobalIDs(mesh, MALLTYPE, 0, NULL, comm);


  // Now verify that there are NO jumps in global IDs
  vgids.clear();
  idx = 0;
  while ((mv = MESH_Next_Vertex(mesh, &idx)))
    if (MV_PType(mv) != PGHOST)
      vgids.push_back(MV_GlobalID(mv));
  std::sort(vgids.begin(), vgids.end());

  vjump_found = false;
  for (int i = 0; i < vgids.size()-1; i++)
    if (vgids[i] != vgids[i+1]-1) {
      vjump_found = true;
      std::cout << "Vjump - (" << i << ", " << vgids[i] << ")  (" << i+1 << ", "
                << vgids[i+1] << ")\n";
      break;
    }
  CHECK(!vjump_found);

  egids.clear();
  idx = 0;
  while ((me = MESH_Next_Edge(mesh, &idx)))
    if (ME_PType(me) != PGHOST)
      egids.push_back(ME_GlobalID(me));
  std::sort(egids.begin(), egids.end());

  ejump_found = false;
  for (int i = 0; i < egids.size()-1; i++)
    if (egids[i] != egids[i+1]-1) {
      ejump_found = true;
      std::cout << "Ejump - (" << i << ", " << egids[i] << ")  (" << i+1 << ", "
                << egids[i+1] << ")\n";
      break;
    }
  CHECK(!ejump_found);

  fgids.clear();
  idx = 0;
  while ((mf = MESH_Next_Face(mesh, &idx)))
    if (MF_PType(mf) != PGHOST)
      fgids.push_back(MF_GlobalID(mf));
  std::sort(fgids.begin(), fgids.end());

  fjump_found = false;
  for (int i = 0; i < fgids.size()-1; i++)
    if (fgids[i] != fgids[i+1]-1) {
      fjump_found = true;
      std::cout << "Fjump - (" << i << ", " << fgids[i] << ")  (" << i+1 << ", "
                << fgids[i+1] << ")\n";
      break;
    }
  CHECK(!fjump_found);

  CHECK(MESH_Parallel_Check(mesh, comm) == 1);

  return;
}
}
