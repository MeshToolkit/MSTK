#include <UnitTest++.h>

#include "../../../include/MSTK.h"

SUITE(Parallel) {
TEST(Weave2D) {

  int i, idx;
  int nf, nv, ngf, nof, ngv, nov;
  int nproc, rank, status, dim;
  int *faceids, *gfaceids, *ofaceids;
  int *vertexids, *gvertexids, *overtexids;
  Mesh_ptr mesh;
  MFace_ptr mf;
  MVertex_ptr mv;
  char filename[256];


  int expnf[4]={6,4,6,4};
  int expngf[4]={4,3,4,3};
  int expnof[4]={2,1,2,1};

  int expfaceids[4][6]={
    {1,2,3,4,5,6},
    {3,2,5,6,0,0},
    {4,5,1,2,3,6},
    {6,2,3,5,0,0}};

  int expgfaceids[4][4]={
    {3,4,5,6},
    {2,5,6,0},
    {1,2,3,6},
    {2,3,5,0}};

  int expofaceids[4][2]={
    {1,2},
    {3,0},
    {4,5},
    {6,0}};


  int expnv[4] = {12,9,12,9};

  int expvertexids[4][12]={
    {1,2,3,4,5,6,7,8,10,9,11,12},
    {3,7,6,8,2,5,11,10,12,0,0,0},
    {4,5,6,9,10,11,1,2,3,7,8,12},
    {6,8,11,12,2,3,5,7,10,0,0,0}};


  MSTK_Init();
  MSTK_Set_Comm(MPI_COMM_WORLD);

  int debugwait=0;
  while (debugwait);


  nproc = MSTK_Comm_size();
  rank = MSTK_Comm_rank();


  mesh = MESH_New(UNKNOWN_REP);

  sprintf(filename,"parallel/4proc/quad3x2.mstk.%-1d",rank);
  status = MESH_InitFromFile(mesh,filename);

  CHECK(status);

  CHECK(!MESH_Num_Regions(mesh));
  CHECK(MESH_Num_Faces(mesh));


  int input_type = 0;  /* no parallel info present in meshes */
  int num_ghost_layers = 1; /* always */
  status = MSTK_Weave_DistributedMeshes(mesh, num_ghost_layers, input_type);

  CHECK(status);


  nf = MESH_Num_Faces(mesh); /* includes ghost faces */

  CHECK_EQUAL(expnf[rank],nf);

  faceids = (int *) malloc(nf*sizeof(int));

  idx = 0; i = 0;
  while ((mf = MESH_Next_Face(mesh,&idx)))
    faceids[i++] = MEnt_GlobalID(mf);

  CHECK_ARRAY_EQUAL(expfaceids[rank],faceids,nf);

  free(faceids);


  nv = MESH_Num_Vertices(mesh); /* includes ghost vertices */

  CHECK_EQUAL(expnv[rank],nv);

  vertexids = (int *) malloc(nv*sizeof(int));

  idx = 0; i = 0;
  while ((mv = MESH_Next_Vertex(mesh,&idx)))
    vertexids[i++] = MEnt_GlobalID(mv);

  CHECK_ARRAY_EQUAL(expvertexids[rank],vertexids,nv);

  free(vertexids);



  ngf = MESH_Num_GhostFaces(mesh);
  CHECK_EQUAL(expngf[rank],ngf);

  gfaceids = (int *) malloc(ngf*sizeof(int));

  idx = 0; i = 0;
  while ((mf = MESH_Next_GhostFace(mesh,&idx)))
    gfaceids[i++] = MEnt_GlobalID(mf);

  CHECK_ARRAY_EQUAL(expgfaceids[rank],gfaceids,ngf);

  free(gfaceids);

  nof = MESH_Num_OverlapFaces(mesh);
  CHECK_EQUAL(expnof[rank],nof);


  ofaceids = (int *) malloc(nof*sizeof(int));

  idx = 0; i = 0;
  while ((mf = MESH_Next_OverlapFace(mesh,&idx)))
    ofaceids[i++] = MEnt_GlobalID(mf);

  CHECK_ARRAY_EQUAL(expofaceids[rank],ofaceids,nof);

  free(ofaceids);


  /*
  ngv = MESH_Num_GhostVertices(mesh);
  CHECK_EQUAL(expngv[rank],ngv);

  gvertexids = (int *) malloc(nv*sizeof(int));

  idx = 0; i = 0;
  while ((mv = MESH_Next_GhostVertex(mesh,&idx)))
    gvertexids[i++] = MEnt_GlobalID(mv);

  free(gvertexids);

  CHECK_ARRAY_EQUAL(expgvertexids[rank],gvertexids);


  nov = MESH_NUM_OverlapVertices(mesh);
  CHECK_EQUAL(expnov[rank],nov);

  overtexids = (int *) malloc(nov*sizeof(int));

  idx = 0; i = 0;
  while ((mv = MESH_Next_OverlapVertex(mesh,&idx)))
    overtexids[i++] = MEnt_GlobalID(mv);

  free(overtexids);

  CHECK_ARRAY_EQUAL(expovertexids[rank],overtexids);
  */

  return;
}
}
