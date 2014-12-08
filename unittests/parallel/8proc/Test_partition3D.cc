#include <UnitTest++.h>

#include "../../../include/MSTK.h"

SUITE(Parallel) {
#if defined (_MSTK_HAVE_METIS)
TEST(Partition3D_sym_0ring_METIS) {

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

  
  /* Total number of entities - ghost + owned */
  int expnr[8]={1,1,1,1,1,1,1,1};
  int expnf[8]={6,6,6,6,6,6,6,6};
  int expne[8]={12,12,12,12,12,12,12,12};
  int expnv[8]={8,8,8,8,8,8,8,8};

  /* Number of ghost entities */
  int expngr[8]={0,0,0,0};
  int expngf[8]={0,1,1,2,1,2,2,3};
  int expnge[8]={0,4,4,7,4,7,7,9};
  int expngv[8]={0,4,4,6,4,6,6,7};

  /* Number of overlap entities */
  int expnor[8]={0,0,0,0,0,0,0,0};
  int expnof[8]={3,2,2,1,2,1,1,0};
  int expnoe[8]={9,5,5,2,5,2,2,0};
  int expnov[8]={7,3,3,1,3,1,1,0};


  MSTK_Init();
  MSTK_Comm comm = MPI_COMM_WORLD;
  int debugwait=0;
  while (debugwait);


  MPI_Comm_size(comm,&nproc);
  MPI_Comm_rank(comm,&rank);

  int method=0;

  Mesh_ptr globalmesh=NULL;
  if (rank == 0) {
      
    globalmesh = MESH_New(UNKNOWN_REP);
    status = MESH_InitFromFile(globalmesh,"parallel/8proc/hex2x2x2.mstk",comm);
      
    CHECK(status);

    status = MESH_CheckTopo(globalmesh);
    CHECK(status);
      
    CHECK(MESH_Num_Regions(globalmesh));
      
    dim = 3;
  }
    
  mymesh = NULL;
  int num_ghost_layers = 0;
  int with_attr = 0;
  int del_inmesh = 1;
  status = MSTK_Mesh_Distribute(globalmesh, &mymesh, &dim, num_ghost_layers, 
                                with_attr, method, del_inmesh, comm);

  CHECK(status);

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

  MESH_Delete(mymesh);

  /*
  if (rank == 0) MESH_Delete(globalmesh);    
  */

  return;
}

TEST(Partition3D_sym_1ring_METIS) {

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
    
  int method=0;

  Mesh_ptr globalmesh=NULL;
  if (rank == 0) {
    
    globalmesh = MESH_New(UNKNOWN_REP);
    status = MESH_InitFromFile(globalmesh,"parallel/8proc/hex2x2x2.mstk",comm);
    
    CHECK(status);

    status = MESH_CheckTopo(globalmesh);
    CHECK(status);
    
    CHECK(MESH_Num_Regions(globalmesh));
    
    dim = 3;
  }
  
  mymesh = NULL;
  int num_ghost_layers = 1;
  int with_attr = 0;
  int del_inmesh = 1;
  status = MSTK_Mesh_Distribute(globalmesh, &mymesh, &dim, num_ghost_layers, 
				with_attr, method, del_inmesh, comm);
  
  CHECK(status);

  status = MESH_CheckTopo(mymesh);
  CHECK(status);
  
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
  
  MESH_Delete(mymesh);

  /*
  if (rank == 0) MESH_Delete(globalmesh);
  */

  return;
}
#endif


#if defined (_MSTK_HAVE_ZOLTAN)
TEST(Partition3D_sym_0ring_ZOLTAN_GRAPH) {

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

  
  /* Total number of entities - ghost + owned */
  int expnr[8]={1,1,1,1,1,1,1,1};
  int expnf[8]={6,6,6,6,6,6,6,6};
  int expne[8]={12,12,12,12,12,12,12,12};
  int expnv[8]={8,8,8,8,8,8,8,8};

  /* Number of ghost entities */
  int expngr[8]={0,0,0,0};
  int expngf[8]={0,1,1,2,1,2,2,3};
  int expnge[8]={0,4,4,7,4,7,7,9};
  int expngv[8]={0,4,4,6,4,6,6,7};

  /* Number of overlap entities */
  int expnor[8]={0,0,0,0,0,0,0,0};
  int expnof[8]={3,2,2,1,2,1,1,0};
  int expnoe[8]={9,5,5,2,5,2,2,0};
  int expnov[8]={7,3,3,1,3,1,1,0};


  MSTK_Init();
  MSTK_Comm comm = MPI_COMM_WORLD;
  int debugwait=0;
  while (debugwait);


  MPI_Comm_size(comm,&nproc);
  MPI_Comm_rank(comm,&rank);

  int method=1; // Zoltan with a graph partitioner

  Mesh_ptr globalmesh=NULL;
  if (rank == 0) {
      
    globalmesh = MESH_New(UNKNOWN_REP);
    status = MESH_InitFromFile(globalmesh,"parallel/8proc/hex2x2x2.mstk",comm);
      
    CHECK(status);

    status = MESH_CheckTopo(globalmesh);
    CHECK(status);
      
    CHECK(MESH_Num_Regions(globalmesh));
      
    dim = 3;
  }
    
  mymesh = NULL;
  int num_ghost_layers = 0;
  int with_attr = 0;
  int del_inmesh = 1;
  status = MSTK_Mesh_Distribute(globalmesh, &mymesh, &dim, num_ghost_layers, 
                                with_attr, method, del_inmesh, comm);

  CHECK(status);

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

  MESH_Delete(mymesh);

  /*
  if (rank == 0) MESH_Delete(globalmesh);    
  */

  return;
}

TEST(Partition3D_sym_1ring_ZOLTAN_GRAPH) {

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
    
  int method=1;   // Zoltan with a graph partitioner

  Mesh_ptr globalmesh=NULL;
  if (rank == 0) {
    
    globalmesh = MESH_New(UNKNOWN_REP);
    status = MESH_InitFromFile(globalmesh,"parallel/8proc/hex2x2x2.mstk",comm);
    
    CHECK(status);

    status = MESH_CheckTopo(globalmesh);
    CHECK(status);
    
    CHECK(MESH_Num_Regions(globalmesh));
    
    dim = 3;
  }
  
  mymesh = NULL;
  int num_ghost_layers = 1;
  int with_attr = 0;
  int del_inmesh = 1;
  status = MSTK_Mesh_Distribute(globalmesh, &mymesh, &dim, num_ghost_layers, 
				with_attr, method, del_inmesh, comm);
  
  CHECK(status);

  status = MESH_CheckTopo(mymesh);
  CHECK(status);
  
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
  
  MESH_Delete(mymesh);
  
  /* if (rank == 0) MESH_Delete(globalmesh); */

  return;
}
#endif


// This test is less of a check on the functioning of the mesh distribution
// and more of a check to see if we were able to generate a partitioning only
// in the XY plane using ZOLTAN RCB partitioner

#if defined (_MSTK_HAVE_ZOLTAN)
TEST(Partition3D_sym_0ring_ZOLTAN_RCB) {

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

  
  /* Total number of entities - ghost + owned */
  int expnr[8]={4,4,4,4,4,4,4,4};
  int expnf[8]={20,20,20,20,20,20,20,20};
  int expne[8]={33,33,33,33,33,33,33,33};
  int expnv[8]={18,18,18,18,18,18,18,18};

  /* Number of ghost entities */
  int expngr[8]={0,0,0,0};
  int expngf[8]={0,4,2,6,4,4,6,6};
  int expnge[8]={0,12,7,17,12,12,17,17};
  int expngv[8]={0,9,6,12,9,9,12,12};

  /* Number of overlap entities */
  int expnor[8]={0,0,0,0,0,0,0,0};
  int expnof[8]={6,6,4,4,6,2,4,0};
  int expnoe[8]={17,15,10,10,15,5,10,0};
  int expnov[8]={12,9,6,6,9,3,6,0};


  MSTK_Init();
  MSTK_Comm comm = MPI_COMM_WORLD;
  int debugwait=0;
  while (debugwait);


  MPI_Comm_size(comm,&nproc);
  MPI_Comm_rank(comm,&rank);

  int method=2;   // Zoltan with Recursive Coordinate Bisection

  Mesh_ptr globalmesh=NULL;
  if (rank == 0) {
      
    globalmesh = MESH_New(UNKNOWN_REP);
    status = MESH_InitFromFile(globalmesh,"parallel/8proc/hex4x4x2.mstk",comm);
      
    CHECK(status);

    status = MESH_CheckTopo(globalmesh);
    CHECK(status);
      
    CHECK(MESH_Num_Regions(globalmesh));
      
    dim = 3;
  }
    
  mymesh = NULL;
  int num_ghost_layers = 0;
  int with_attr = 0;
  int del_inmesh = 1;
  status = MSTK_Mesh_Distribute(globalmesh, &mymesh, &dim, num_ghost_layers, 
                                with_attr, method, del_inmesh, comm);

  CHECK(status);

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

  MESH_Delete(mymesh);

  /* if (rank == 0) MESH_Delete(globalmesh); */

  return;
}


#endif


}
