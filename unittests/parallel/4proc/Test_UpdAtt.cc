#include <UnitTest++.h>

#include "../../../include/MSTK.h"

SUITE(Parallel) {
TEST(UpdAtt2D_Dist) {

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
  MSTK_Set_Comm(MPI_COMM_WORLD);

  nproc = MSTK_Comm_size();
  rank = MSTK_Comm_rank();

  Mesh_ptr mesh0;
  if (rank == 0) {
    mesh0 = MESH_New(UNKNOWN_REP);

    sprintf(filename,"parallel/4proc/quad10x10.mstk");
    MESH_InitFromFile(mesh0, filename);
    
    if (MESH_Num_Regions(mesh0) > 0) {
      fprintf(stderr,"Code is for surface meshes only. Exiting...\n");
    }
    else if (MESH_Num_Faces(mesh0) > 0)
      dim = 2;
    else {
      fprintf(stderr,"Mesh is neither solid nor surface mesh. Exiting...\n");
      exit(-1);
    }
  }

//  DebugWait=1;
//  while (DebugWait);

  int ring = 1; /* One ring ghosts */
  int with_attr = 1; /* Do allow exchange of attributes */
  int method;

#if defined (_MSTK_HAVE_METIS)
  method = 0;
#elif defined (_MSTK_HAVE_ZOLTAN)
  method = 1;
#else
  fprintf(stderr,"Cannot find partitioner\n");
  int status = 0;
  CHECK(status);
#endif

  MSTK_Mesh_Distribute(mesh0, &mesh, &dim, ring, with_attr, method);

  if (rank == 0) MESH_Delete(mesh0);

  MAttrib_ptr vcolatt, fcolatt;
  vcolatt = MAttrib_New(mesh,"vcolatt",INT,MVERTEX);

  idx = 0;
  while (mv = MESH_Next_OverlapVertex(mesh,&idx))
    MEnt_Set_AttVal(mv,vcolatt,rank+1,0.0,NULL);

  fcolatt = MAttrib_New(mesh,"fcolatt",INT,MFACE);

  idx = 0;
  while (mf = MESH_Next_OverlapFace(mesh,&idx))
    MEnt_Set_AttVal(mf,fcolatt,rank+1,0.0,NULL);


  MSTK_UpdateAttr(mesh);

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
  MSTK_Set_Comm(MPI_COMM_WORLD);


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
   int topodim = 2;
   status = MSTK_Weave_DistributedMeshes(mesh, topodim, 
                                         num_ghost_layers, input_type);

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


  MSTK_UpdateAttr(mesh);

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


