#include <UnitTest++.h>

#include "../../../include/MSTK.h"

SUITE(Parallel) {
TEST(UpdAtt2D_Dist) {

  int i, idx, ival, status;
  int vcol, fcol, mpid;
  int nproc, rank, dim;
  double rval;
  void *pval;
  Mesh_ptr mesh;
  MVertex_ptr mv;
  MFace_ptr mf;
  char filename[256];


  MSTK_Init();

  MSTK_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm,&nproc);
  MPI_Comm_rank(comm,&rank);

  Mesh_ptr mesh0;
  if (rank == 0) {
    mesh0 = MESH_New(UNKNOWN_REP);

    sprintf(filename,"parallel/4proc/quad10x10.mstk");
    MESH_InitFromFile(mesh0, filename, comm);
    
    if (MESH_Num_Regions(mesh0) > 0) {
      fprintf(stderr,"Code is for surface meshes only. Exiting...\n");
    }
    else if (MESH_Num_Faces(mesh0) > 0)
      dim = 2;
    else {
      fprintf(stderr,"Mesh is neither solid nor surface mesh. Exiting...\n");
      exit(-1);
    }

    status = MESH_CheckTopo(mesh0);
    CHECK(status);
  }

//  DebugWait=1;
//  while (DebugWait);

  int ring = 1; /* One ring ghosts */
  int with_attr = 1; /* Do allow exchange of attributes */
  int del_inmesh = 1; /* delete input mesh after partitioning */
  int method;

#if defined (_MSTK_HAVE_METIS)
  method = 0;
#elif defined (_MSTK_HAVE_ZOLTAN)
  method = 1;
#else
  fprintf(stderr,"Cannot find partitioner\n");
  status = 0;
  CHECK(status);
#endif

  mesh = NULL;
  MSTK_Mesh_Distribute(mesh0, &mesh, &dim, ring, with_attr, method, 
		       del_inmesh, comm);

  status = MESH_CheckTopo(mesh);
  CHECK(status);

  /*  if (rank == 0) MESH_Delete(mesh0); */

  MAttrib_ptr vcolatt, fcolatt;
  vcolatt = MAttrib_New(mesh,"vcolatt",INT,MVERTEX);

  idx = 0;
  while (mv = MESH_Next_OverlapVertex(mesh,&idx))
    MEnt_Set_AttVal(mv,vcolatt,rank+1,0.0,NULL);

  fcolatt = MAttrib_New(mesh,"fcolatt",INT,MFACE);

  idx = 0;
  while (mf = MESH_Next_OverlapFace(mesh,&idx))
    MEnt_Set_AttVal(mf,fcolatt,rank+1,0.0,NULL);


  MESH_UpdateAttributes(mesh, comm);

  idx = 0;
  while ((mv = MESH_Next_GhostVertex(mesh,&idx))) {
    MEnt_Get_AttVal(mv,vcolatt,&vcol,&rval,&pval);

    mpid = MV_MasterParID(mv);

    CHECK_EQUAL(mpid+1,vcol);
  }

  idx = 0;
  while ((mf = MESH_Next_GhostFace(mesh,&idx))) {
    MEnt_Get_AttVal(mf,fcolatt,&fcol,&rval,&pval);

    mpid = MF_MasterParID(mf);

    CHECK_EQUAL(mpid+1,fcol);
  }

  return;
}



TEST(UpdAtt2D_Weave) {

  int i, idx, ival;
  int vcol, fcol, mpid;
  int nproc, rank, status, dim;
  double rval;
  void *pval;
  Mesh_ptr mesh;
  MVertex_ptr mv;
  MFace_ptr mf;
  char filename[256];


  MSTK_Init();
  MSTK_Comm comm = MPI_COMM_WORLD;


  MPI_Comm_size(comm,&nproc);
  MPI_Comm_rank(comm,&rank);


  mesh = MESH_New(UNKNOWN_REP);

   sprintf(filename,"parallel/4proc/quad3x2.mstk.%-1d",rank);
   status = MESH_InitFromFile(mesh,filename,comm);

   CHECK(status);

   CHECK(!MESH_Num_Regions(mesh));
   CHECK(MESH_Num_Faces(mesh));


   int input_type = 0;  /* no parallel info present in meshes */
   int num_ghost_layers = 1; /* always */
   int topodim = 2;
   status = MSTK_Weave_DistributedMeshes(mesh, topodim, 
                                         num_ghost_layers, input_type, comm);

   CHECK(status);

   status = MESH_CheckTopo(mesh);
   CHECK(status);

   int DebugWait=0;
   while (DebugWait);


  MAttrib_ptr vcolatt, fcolatt;
  vcolatt = MAttrib_New(mesh,"vcolatt",INT,MVERTEX);

  idx = 0;
  while (mv = MESH_Next_OverlapVertex(mesh,&idx))
    MEnt_Set_AttVal(mv,vcolatt,rank+1,0.0,NULL);

  fcolatt = MAttrib_New(mesh,"fcolatt",INT,MFACE);

  idx = 0;
  while (mf = MESH_Next_OverlapFace(mesh,&idx))
    MEnt_Set_AttVal(mf,fcolatt,rank+1,0.0,NULL);


  MESH_UpdateAttributes(mesh,comm);

  idx = 0;
  while ((mv = MESH_Next_GhostVertex(mesh,&idx))) {
    MEnt_Get_AttVal(mv,vcolatt,&vcol,&rval,&pval);

    mpid = MV_MasterParID(mv);

    CHECK_EQUAL(mpid+1,vcol);
  }

  idx = 0;
  while ((mf = MESH_Next_GhostFace(mesh,&idx))) {
    MEnt_Get_AttVal(mf,fcolatt,&fcol,&rval,&pval);

    mpid = MF_MasterParID(mf);

    CHECK_EQUAL(mpid+1,fcol);
  }

  return;

}

}


