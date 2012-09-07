#include <UnitTest++.h>

#include "../../../include/MSTK.h"

SUITE(Parallel) {
TEST(Partition2D) {

  int i, idx;
  int nf, nv, ngf, nof, ngv, nov;
  int nproc, rank, status, dim;
  int *faceids, *gfaceids, *ofaceids;
  int *vertexids, *gvertexids, *overtexids;
  Mesh_ptr mymesh;
  MFace_ptr mf;
  MVertex_ptr mv;

#if defined (_MSTK_HAVE_METIS)

  int expnf[4]={7,8,7,8};
  int expngf[4]={8,7,8,7};
  int expnof[4]={5,5,5,5};

  int expfaceids[4][8]={
    {1,2,3,4,5,6,7,0},
    {8,9,10,11,12,13,14,15},
    {16,17,18,19,20,21,22,0},
    {23,24,25,26,27,28,29,30}};

  int expgfaceids[4][8]={
    {9,11,12,15,20,21,22,30},
    {1,2,4,6,20,29,30,0},
    {1,2,3,9,25,27,28,30},
    {1,8,9,16,18,20,21,0}};

  int expofaceids[4][5]={
    {1,2,3,4,6},
    {8,9,11,12,15},
    {16,18,20,21,22},
    {25,27,28,29,30}};


  int expnv[4] = {14,10,10,8};

  int expvertexids[4][14]={
    {1,2,3,4,5,6,7,8,9,10,11,12,13,14},
    {15,16,17,18,19,20,21,22,23,24,0,0,0,0},
    {25,26,27,28,29,30,31,32,33,34,0,0,0,0},
    {35,36,37,38,39,40,41,42,10,0,0,0,0,0}};

#elif defined (_MSTK_HAVE_ZOLTAN)


  int expnf[4]={7,8,7,8};
  int expngf[4]={8,7,8,7};
  int expnof[4]={5,5,5,5};

  int expfaceids[4][8]={
    {1,2,3,4,5,6,7,0},
    {8,9,10,11,12,13,14,15},
    {16,17,18,19,20,21,22,0},
    {23,24,25,26,27,28,29,30}};

  int expgfaceids[4][8]={
    {8,10,11,13,20,21,29,30},
    {2,3,5,7,20,21,22,0},
    {2,3,8,9,25,27,28,30},
    {1,2,3,16,18,20,21,0}};

  int expofaceids[4][5]={
    {1,2,3,5,7},
    {8,9,10,11,13},
    {16,18,20,21,22},
    {25,27,28,29,30}};


  int expnv[4] = {14,10,10,8};

  int expvertexids[4][14]={
    {1,2,3,4,5,6,7,8,9,10,11,12,13,14},
    {15,16,17,18,19,20,21,22,23,24,0,0,0,0},
    {25,26,27,28,29,30,31,32,33,34,0,0,0,0},
    {35,36,37,38,39,40,41,42,10,0,0,0,0,0}};  

#endif

  MSTK_Init();
  MSTK_Comm comm = MPI_COMM_WORLD;

  int debugwait=0;
  while (debugwait);


  MPI_Comm_size(comm,&nproc);
  MPI_Comm_rank(comm,&rank);

  Mesh_ptr globalmesh=NULL;

  if (rank == 0) {

    globalmesh = MESH_New(UNKNOWN_REP);
    status = MESH_InitFromFile(globalmesh,"parallel/4proc/quad6x5.mstk", comm);

    CHECK(status);

    CHECK(!MESH_Num_Regions(globalmesh));
    CHECK(MESH_Num_Faces(globalmesh));

    dim = 2;
  }

  int method;

#if defined (_MSTK_HAVE_METIS)    
  method = 0;
#elif defined (_MSTK_HAVE_ZOLTAN)
  method = 1;
#else
  fprintf(stderr,"Cannot find partitioner");
  status = 0;
  CHECK(status);
#endif

  mymesh = NULL;
  status = MSTK_Mesh_Distribute(globalmesh,&mymesh, &dim, 1, 1, method, comm);

  CHECK(status);


  if (rank == 0)
    MESH_Delete(globalmesh);


  status = MESH_Parallel_Check(mymesh, comm);

  CHECK(status);


  nf = MESH_Num_Faces(mymesh)-MESH_Num_GhostFaces(mymesh);

  CHECK_EQUAL(expnf[rank],nf);

  faceids = (int *) malloc(nf*sizeof(int));

  idx = 0; i = 0;
  while ((mf = MESH_Next_Face(mymesh,&idx))) {
    if (MF_PType(mf) != PGHOST)
      faceids[i++] = MEnt_GlobalID(mf);
  }

  CHECK_ARRAY_EQUAL(expfaceids[rank],faceids,nf);

  free(faceids);


  nv = MESH_Num_Vertices(mymesh)-MESH_Num_GhostVertices(mymesh);

  CHECK_EQUAL(expnv[rank],nv);

  vertexids = (int *) malloc(nv*sizeof(int));

  idx = 0; i = 0;
  while ((mv = MESH_Next_Vertex(mymesh,&idx)))
    if (MV_PType(mv) != PGHOST)
      vertexids[i++] = MEnt_GlobalID(mv);

  CHECK_ARRAY_EQUAL(expvertexids[rank],vertexids,nv);

  free(vertexids);



  ngf = MESH_Num_GhostFaces(mymesh);
  CHECK_EQUAL(expngf[rank],ngf);

  gfaceids = (int *) malloc(ngf*sizeof(int));

  idx = 0; i = 0;
  while ((mf = MESH_Next_GhostFace(mymesh,&idx)))
    gfaceids[i++] = MEnt_GlobalID(mf);

  CHECK_ARRAY_EQUAL(expgfaceids[rank],gfaceids,ngf);

  free(gfaceids);

  nof = MESH_Num_OverlapFaces(mymesh);
  CHECK_EQUAL(expnof[rank],nof);


  ofaceids = (int *) malloc(nof*sizeof(int));

  idx = 0; i = 0;
  while ((mf = MESH_Next_OverlapFace(mymesh,&idx)))
    ofaceids[i++] = MEnt_GlobalID(mf);

  CHECK_ARRAY_EQUAL(expofaceids[rank],ofaceids,nof);

  free(ofaceids);


  /*
  ngv = MESH_Num_GhostVertices(mymesh);
  CHECK_EQUAL(expngv[rank],ngv);

  gvertexids = (int *) malloc(nv*sizeof(int));

  idx = 0; i = 0;
  while ((mv = MESH_Next_GhostVertex(mymesh,&idx)))
    gvertexids[i++] = MEnt_GlobalID(mv);

  free(gvertexids);

  CHECK_ARRAY_EQUAL(expgvertexids[rank],gvertexids);


  nov = MESH_NUM_OverlapVertices(mymesh);
  CHECK_EQUAL(expnov[rank],nov);

  overtexids = (int *) malloc(nov*sizeof(int));

  idx = 0; i = 0;
  while ((mv = MESH_Next_OverlapVertex(mymesh,&idx)))
    overtexids[i++] = MEnt_GlobalID(mv);

  free(overtexids);

  CHECK_ARRAY_EQUAL(expovertexids[rank],overtexids);
  */

  return;
}
}
