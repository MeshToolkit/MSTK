/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <UnitTest++.h>
#include <cstdio>
#include <vector>
#include <algorithm>

#include "MSTK.h"

SUITE(Parallel) {
TEST(Weave2D) {

  int i, idx;
  int nf, nv, ngf, nof, ngv, nov;
  int nproc, rank, status, dim;
  Mesh_ptr mesh;
  MFace_ptr mf;
  MVertex_ptr mv;
  char filename[256];


  int expnf[4]={9,8,6,9};
  int expngf[4]={6,5,3,6};
  int expnof[4]={3,3,3,3};

  int expfaceids[4][9]={
    {1,2,3,4,5,6,10,11,12},
    {1,2,3,4,5,6,10,11,0},
    {7,8,9,10,11,12,0,0,0},
    {2,3,6,7,8,9,10,11,12}};

  int expgfaceids[4][6]={
    {4,5,6,10,11,12},
    {1,2,3,10,11,0},
    {10,11,12,0,0,0},
    {2,3,6,7,8,9}};

  int expofaceids[4][3]={
    {1,2,3},
    {4,5,6},
    {7,8,9},
    {10,11,12}};


  int expnv[4] = {16,15,12,16};

  int expvertexids[4][16]={
    {1,2,3,4,5,6,7,8,9,10,11,12,13,16,18,20},
    {1,2,3,4,5,6,7,8,9,10,11,12,13,16,18,0},
    {5,7,8,12,13,14,15,16,17,18,19,20,0,0,0,0},
    {1,3,5,6,7,8,10,12,13,14,15,16,17,18,19,20}};

  int expngv[4] = {8,11,4,16};
  int expgvertexids[4][16]={
    {9,10,11,12,13,16,18,20,0,0,0,0,0,0,0,0},
    {1,2,3,4,5,6,7,8,13,16,18,0,0,0,0,0},
    {5,7,8,12,0,0,0,0,0,0,0,0,0,0,0,0},
    {1,3,5,6,7,8,10,12,13,14,15,16,17,18,19,20}};

  int expnov[4] = {8,4,8,0};
  int expovertexids[4][8]={
    {1,2,3,4,5,6,7,8},
    {9,10,11,12,0,0,0,0},
    {13,14,15,16,17,18,19,20},
    {0,0,0,0,0,0,0,0}};
  
  MSTK_Init();
  MSTK_Comm comm = MPI_COMM_WORLD;

  int debugwait=0;
  while (debugwait);


  MPI_Comm_size(comm,&nproc);
  MPI_Comm_rank(comm,&rank);


  mesh = MESH_New(UNKNOWN_REP);

  sprintf(filename,"parallel/4proc/quad4x3.mstk.%-1d",rank);
  status = MESH_InitFromFile(mesh,filename,comm);

  CHECK(status);

  CHECK(!MESH_Num_Regions(mesh));
  CHECK(MESH_Num_Faces(mesh));

  int input_type = 0;  /* no parallel info present in meshes */
  int num_ghost_layers = 1; /* always */
  int topodim = 2;
  status = MSTK_Weave_DistributedMeshes(mesh, topodim, num_ghost_layers, input_type, comm);

  CHECK(status);

  CHECK(MESH_CheckTopo(mesh));


  nf = MESH_Num_Faces(mesh); /* includes ghost faces */

  CHECK_EQUAL(expnf[rank],nf);

  std::vector<int> faceids(nf);

  idx = 0; i = 0;
  while ((mf = MESH_Next_Face(mesh,&idx)))
    faceids[i++] = MEnt_GlobalID(mf);
  std::sort(faceids.begin(), faceids.end());

  CHECK_ARRAY_EQUAL(expfaceids[rank],faceids.data(),nf);



  nv = MESH_Num_Vertices(mesh); /* includes ghost vertices */

  CHECK_EQUAL(expnv[rank],nv);

  std::vector<int> vertexids(nv);

  idx = 0; i = 0;
  while ((mv = MESH_Next_Vertex(mesh,&idx)))
    vertexids[i++] = MEnt_GlobalID(mv);
  std::sort(vertexids.begin(), vertexids.end());

  CHECK_ARRAY_EQUAL(expvertexids[rank],vertexids.data(),nv);




  ngf = MESH_Num_GhostFaces(mesh);
  CHECK_EQUAL(expngf[rank],ngf);

  std::vector<int> gfaceids(ngf);

  idx = 0; i = 0;
  while ((mf = MESH_Next_GhostFace(mesh,&idx)))
    gfaceids[i++] = MEnt_GlobalID(mf);
  std::sort(gfaceids.begin(), gfaceids.end());

  CHECK_ARRAY_EQUAL(expgfaceids[rank],gfaceids,ngf);


  nof = MESH_Num_OverlapFaces(mesh);
  CHECK_EQUAL(expnof[rank],nof);


  std::vector<int> ofaceids(nof);

  idx = 0; i = 0;
  while ((mf = MESH_Next_OverlapFace(mesh,&idx)))
    ofaceids[i++] = MEnt_GlobalID(mf);
  std::sort(ofaceids.begin(), ofaceids.end());

  CHECK_ARRAY_EQUAL(expofaceids[rank],ofaceids,nof);



  ngv = MESH_Num_GhostVertices(mesh);
  CHECK_EQUAL(expngv[rank],ngv);

  std::vector<int> gvertexids(ngv);

  idx = 0; i = 0;
  while ((mv = MESH_Next_GhostVertex(mesh,&idx)))
    gvertexids[i++] = MEnt_GlobalID(mv);
  std::sort(gvertexids.begin(), gvertexids.end());

  CHECK_ARRAY_EQUAL(expgvertexids[rank],gvertexids.data(),ngv);


  nov = MESH_Num_OverlapVertices(mesh);
  CHECK_EQUAL(expnov[rank],nov);

  std::vector<int> overtexids(nov);

  idx = 0; i = 0;
  while ((mv = MESH_Next_OverlapVertex(mesh,&idx)))
    overtexids[i++] = MEnt_GlobalID(mv);
  std::sort(overtexids.begin(), overtexids.end());

  CHECK_ARRAY_EQUAL(expovertexids[rank],overtexids.data(),nov);

  return;
}
}
