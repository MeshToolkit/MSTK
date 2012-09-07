#include <UnitTest++.h>

#include "../../../include/MSTK.h"

SUITE(Parallel) {
TEST(Partition3D_sym) {

  int i, j, idx;
  int nr, nf, ne, nv, ngr, nor, ngf, nof, nge, noe, ngv, nov;
  int nproc, rank, status, dim;
  int *regids, *gregids, *oregids;
  int *faceids, *gfaceids, *ofaceids;
  int *edgeids, *gedgeids, *oedgeids;
  int *vertexids, *gvertexids, *overtexids;
  Mesh_ptr mesh, mymesh;
  MRegion_ptr mr;
  MFace_ptr mf;
  MEdge_ptr me;
  MVertex_ptr mv;

  
  /* Expected centroids - to verify that the eight cells got assigned
     to the processors we expected */
  double exp_centroid[9][3] = {{0.75,0.75,0.75},
                               {0.75,0.75,0.25},
                               {0.25,0.75,0.75},
                               {0.25,0.75,0.25},
                               {0.75,0.25,0.25},
                               {0.75,0.25,0.75},
                               {0.25,0.25,0.75},
                               {0.25,0.25,0.25}};

  /* Total number of entities - ghost + owned */
  int expnr[8]={8,8,8,8,8,8,8,8};
  int expnf[8]={36,36,36,36,36,36,36,36};
  int expne[8]={54,54,54,54,54,54,54,54};
  int expnv[8]={27,27,27,27,27,27,27,27};

  /* Number of ghost entities */
  int expngr[8]={7,7,7,7,7,7,7,7};
  int expngf[8]={30,31,31,32,31,32,32,33};
  int expnge[8]={42,46,46,49,46,49,49,51};
  int expngv[8]={19,23,23,25,23,25,25,26};

  /* Number of overlap entities */
  int expnor[8]={1,1,1,1,1,1,1,1};
  int expnof[8]={6,5,5,4,5,4,4,3};
  int expnoe[8]={12,8,8,5,8,5,5,3};
  int expnov[8]={8,4,4,2,4,2,2,1};


  MSTK_Init();
  MSTK_Comm comm = MPI_COMM_WORLD;
  int debugwait=0;
  while (debugwait);


  MPI_Comm_size(comm,&nproc);
  MPI_Comm_rank(comm,&rank);

  Mesh_ptr globalmesh=NULL;
  if (rank == 0) {

    globalmesh = MESH_New(UNKNOWN_REP);
    status = MESH_InitFromFile(globalmesh,"parallel/8proc/hex2x2x2.mstk",comm);

    CHECK(status);

    CHECK(MESH_Num_Regions(globalmesh));

    dim = 3;
  }
    
  int method;

#if defined (_MSTK_HAVE_METIS)
  method = 0;
#elif defined (_MSTK_HAVE_ZOLTAN)
  method = 1;
#else
  fprintf(stderr,"No partitioner found");
  status = 0;
  CHECK(status);
#endif

  mymesh = NULL;
  status = MSTK_Mesh_Distribute(globalmesh, &mymesh, &dim, 1, 1, method, comm);

  CHECK(status);

  if (rank == 0) MESH_Delete(globalmesh);

  double centroid[3] = {0.0,0.0,0.0};
  double rxyz[8][3];
  int nrv;

  idx = 0; nr = 0;
  while (mr = MESH_Next_Region(mymesh,&idx)) {
    if (MR_PType(mr) != PGHOST) {
      MR_Coords(mr,&nrv,rxyz);
      for (i = 0; i < nrv; i++) {
        for (j = 0; j < 3; j++)
          centroid[j] += rxyz[i][j];
      }
      for (j = 0; j < 3; j++)
        centroid[j] /= nrv;
      nr++;
    }
  }
  CHECK_EQUAL(1,nr);
  CHECK_ARRAY_EQUAL(exp_centroid[rank],centroid,3);

  ngr = MESH_Num_GhostRegions(mymesh);
  CHECK_EQUAL(expngr[rank],ngr);

  nor = MESH_Num_OverlapRegions(mymesh);
  CHECK_EQUAL(expnor[rank],nor);

  nr = MESH_Num_Regions(mymesh);
  CHECK_EQUAL(expnr[rank],nr);

  ngf = MESH_Num_GhostFaces(mymesh);
  CHECK_EQUAL(expngf[rank],ngf);

  nof = MESH_Num_OverlapFaces(mymesh);
  CHECK_EQUAL(expnof[rank],nof);

  nf = MESH_Num_Faces(mymesh);
  CHECK_EQUAL(expnf[rank],nf);

  nge = MESH_Num_GhostEdges(mymesh);
  CHECK_EQUAL(expnge[rank],nge);

  noe = MESH_Num_OverlapEdges(mymesh);
  CHECK_EQUAL(expnoe[rank],noe);

  ne = MESH_Num_Edges(mymesh);
  CHECK_EQUAL(expne[rank],ne);

  ngv = MESH_Num_GhostVertices(mymesh);
  CHECK_EQUAL(expngv[rank],ngv);

  nov = MESH_Num_OverlapVertices(mymesh);
  CHECK_EQUAL(expnov[rank],nov);

  nv = MESH_Num_Vertices(mymesh);
  CHECK_EQUAL(expnv[rank],nv);

  return;
}
}
