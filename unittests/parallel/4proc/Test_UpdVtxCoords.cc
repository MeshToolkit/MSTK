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
      status = 0;
      CHECK(status);
    }
    else if (MESH_Num_Faces(mesh0) > 0)
      dim = 2;
    else {
      fprintf(stderr,"Mesh is neither solid nor surface mesh. Exiting...\n");
      status = 0;
      CHECK(status);
    }

    CHECK(MESH_CheckTopo(mesh0));
  }

  int ring = 1; /* One ring ghosts */
  int with_attr = 1; /* Do allow exchange of attributes */
  int del_inmesh = 1; /* Delete input mesh after partitioning */
  PartitionMethod method;
#if defined (_MSTK_HAVE_METIS)
  method = METIS;
#elif defined (_MSTK_HAVE_ZOLTAN)
  method = ZOLTAN_RCB;
#else
  fprintf(stderr,"Cannot find partitioner\n");
  status = 0;
  CHECK(status);
#endif

  mesh = NULL;
  status = MSTK_Mesh_Distribute(mesh0, &mesh, &dim, ring, with_attr, method, 
				del_inmesh, comm);

  CHECK(MESH_CheckTopo(mesh));

  /* if (rank == 0) MESH_Delete(mesh0); */

  MAttrib_ptr xyzatt;
  xyzatt = MAttrib_New(mesh,"xyzatt",VECTOR,MVERTEX,3);
  double *oxyz = (double *) new double[3*MESH_Num_Vertices(mesh)];

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

  //  fprintf(stderr,"Deformed mesh proc %-d\n",rank);

  MESH_UpdateVertexCoords(mesh, comm);
  MESH_UpdateAttributes(mesh, comm);

  //  fprintf(stderr,"Updated vertex coordinates proc %-d\n",rank);

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

  delete [] oxyz;
  return;
}
}
