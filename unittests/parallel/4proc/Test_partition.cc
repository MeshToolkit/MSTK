#include <UnitTest++.h>

#include "../../../include/MSTK.h"

SUITE(Parallel) {
TEST(Partition2D) {

  int i, idx;
  int nf, nv, ngf, nof, ngv, nov;
  int nproc, rank, status, dim;
  int *faceids, *gfaceids, *ofaceids;
  int *vertexids, *gvertexids, *overtexids;
  Mesh_ptr mesh, mymesh;
  MFace_ptr mf;
  MVertex_ptr mv;


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


  MSTK_Init();

  int debugwait=0;
  while (debugwait);


  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if (rank == 0) {

    mesh = MESH_New(UNKNOWN_REP);
    status = MESH_InitFromFile(mesh,"parallel/4proc/quad6x5.mstk");

    CHECK(status);

    CHECK(!MESH_Num_Regions(mesh));
    CHECK(MESH_Num_Faces(mesh));

    dim = 2;
    mymesh = mesh;
    
  }
    
  status = MSTK_Mesh_Distribute(&mymesh, &dim, 1, 1, rank, nproc, MPI_COMM_WORLD);

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
