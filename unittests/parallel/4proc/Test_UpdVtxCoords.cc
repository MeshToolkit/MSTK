#include <UnitTest++.h>

#include "../../../include/MSTK.h"

SUITE(Parallel) {
TEST(VertexUpdate2D) {

  int i, idx, ival;
  int nv;
  int nproc, rank, status, dim;
  double xyz[3], expxyz[3], rval, delta=0.02;
  void *pval;
  Mesh_ptr mesh;
  MVertex_ptr mv;
  char filename[256];


  MSTK_Init(MPI_COMM_WORLD);


  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);


  //  mesh = MESH_New(UNKNOWN_REP);

  //  sprintf(filename,"parallel/4proc/quad3x2.mstk.%-1d",rank);
  //  status = MESH_InitFromFile(mesh,filename);

  //  CHECK(status);

  //  CHECK(!MESH_Num_Regions(mesh));
  //  CHECK(MESH_Num_Faces(mesh));


  //  int input_type = 0;  /* no parallel info present in meshes */
  //  int num_ghost_layers = 1; /* always */
  //  status = MSTK_Weave_DistributedMeshes(mesh, num_ghost_layers, input_type, rank, nproc, MPI_COMM_WORLD);

  //  CHECK(status);


  if (rank == 0) {
    Mesh_ptr mesh0 = MESH_New(UNKNOWN_REP);

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
    mesh = mesh0;
  }

//  DebugWait=1;
//  while (DebugWait);

  int ring = 1; /* One ring ghosts */
  int with_attr = 1; /* Do allow exchange of attributes */
  int method = 0; /* Use Metis as the partitioner */
  MSTK_Mesh_Distribute(&mesh, &dim, ring, with_attr, method, rank, nproc, 
		       MPI_COMM_WORLD);



  MAttrib_ptr xyzatt;
  xyzatt = MAttrib_New(mesh,"xyzatt",VECTOR,MVERTEX,3);
  double *oxyz = (double *)malloc(3*MESH_Num_Vertices(mesh)*sizeof(double));

  idx = 0; i = 0;
  while (mv = MESH_Next_OverlapVertex(mesh,&idx)) {
    MV_Coords(mv,&(oxyz[3*i]));

    MEnt_Set_AttVal(mv,xyzatt,ival,rval,&(oxyz[3*i]));

    xyz[0] = oxyz[3*i] + delta;
    xyz[1] = oxyz[3*i+1] + delta;
    xyz[2] = oxyz[3*i+2];
    MV_Set_Coords(mv,xyz);

    i++;
  }

  fprintf(stderr,"Deformed mesh proc %-d\n",rank);

  MESH_UpdateVertexCoords(mesh,rank,nproc,MPI_COMM_WORLD);
  MSTK_UpdateAttr(mesh,rank,nproc,MPI_COMM_WORLD);

  fprintf(stderr,"Updated vertex coordinates proc %-d\n",rank);

  nv = MESH_Num_Vertices(mesh); /* includes ghost vertices */

  double *oxyz1;
  idx = 0; i = 0;
  while ((mv = MESH_Next_GhostVertex(mesh,&idx))) {
    MEnt_Get_AttVal(mv,xyzatt,&ival,&rval,&pval);
    oxyz1 = (double *) pval;
    expxyz[0] = oxyz1[0] + delta;
    expxyz[1] = oxyz1[1] + delta;
    expxyz[2] = oxyz1[2];
    
    MV_Coords(mv,xyz);
    CHECK_ARRAY_EQUAL(expxyz,xyz,3);
    if (expxyz[0] != xyz[0] || expxyz[1] != xyz[1]) {
      fprintf(stderr,"On processor %-d, Vertex %-2d (GID=%-2d) has a coordinate mismatch\n",rank,MV_ID(mv),MV_GlobalID(mv));
      fprintf(stderr,"Expected %lf %lf but got %lf %lf\n",expxyz[0],expxyz[1],xyz[0],xyz[1]);
      return;
    }
  }

  free(oxyz);
  return;
}
}
