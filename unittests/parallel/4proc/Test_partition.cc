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
    {18,19,20,24,25,29,30,0},
    {16,17,21,22,23,26,27,28},
    {4,5,9,10,13,14,15,0},
    {1,2,3,6,7,8,11,12}};

  int expgfaceids[4][8]={
    {12,13,14,15,17,22,23,28},
    {11,12,13,18,19,24,29,0},
    {3,7,8,12,17,18,19,20},
    {4,9,13,14,16,17,18}};

  int expofaceids[4][5]={
    {18,19,20,24,29},
    {16,17,22,23,28},
    {4,9,13,14,15},
    {3,7,8,11,12}};


  int expnv[4] = {14,10,10,8};

  int expvertexids[4][14]={
    {18,19,26,25,33,32,40,39,27,34,41,28,35,42},
    {4,5,12,11,6,13,20,7,14,21,0,0,0,0},
    {22,23,30,29,37,36,24,31,38,17,0,0,0,0},
    {1,2,9,8,16,15,3,10,0,0,0,0,0,0}};


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
    
  int method = 0; /* with Metis */
  status = MSTK_Mesh_Distribute(&mymesh, &dim, 1, 1, method, rank, nproc, MPI_COMM_WORLD);

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
