#include <UnitTest++.h>
#include <cstdio>

#include "MSTK.h"
#include "MSTK_private.h"

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
      status = 0;
      CHECK(status);
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



TEST(GathScatAttrib_2D) {

  int i, idx, ival, status;
  int vcol, fcol, mpid;
  int nproc, rank, dim;
  double rval;
  void *pval;
  Mesh_ptr mesh;
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
    
    if (MESH_Num_Regions(mesh0) > 0)
      fprintf(stderr,"Code is for surface meshes only. Exiting...\n");
    else if (MESH_Num_Faces(mesh0) > 0)
      dim = 2;
    else {
      fprintf(stderr,"Mesh is neither solid nor surface mesh. Exiting...\n");
      status = 0;
      CHECK(status);
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



  MAttrib_ptr vertatt, edgeatt;
  vertatt = MAttrib_New(mesh,"vertatt",DOUBLE,MVERTEX);
  edgeatt = MAttrib_New(mesh,"edgeatt",DOUBLE,MEDGE);


  for (int i = 1; i < 5; i++) {
    MAttOpType op = (MAttOpType) i;

    MVertex_ptr mv;
    idx = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx)))
      MEnt_Set_AttVal(mv,vertatt,0,rank+1,NULL);

    MESH_GathScat1Attribute(mesh,vertatt,op,comm);


    /* Now use plain old MPI to determine the right answer */

    int nv = MESH_Num_Vertices(mesh);
    int *allnv = (int *) new int[nproc];
    MPI_Allgather(&nv,1,MPI_INT,allnv,1,MPI_INT,MPI_COMM_WORLD);

    int maxnv = 0;
    for (int j = 0; j < nproc; j++)
      if (allnv[j] > maxnv)
        maxnv = allnv[j];
    int *locgidlist = (int *) new int[maxnv];
    double *locvallist = (double *) new double[maxnv];
    int *globgidlist = (int *) new int[maxnv*nproc];
    double *globvallist = (double *) new double[maxnv*nproc];
    for (int i = 0; i < maxnv; i++) {
      locgidlist[i] = 0;
      locvallist[i] = 0.0;
    }
    for (int i = 0; i < maxnv*nproc; i++) {
      globgidlist[i] = 0;
      globvallist[i] = 0.0;
    }

    int j = 0;
    idx = 0; 
    while ((mv = MESH_Next_Vertex(mesh,&idx))) {
      locgidlist[j] = MV_GlobalID(mv);
      locvallist[j] = rank+1;
      j++;
    }

    MPI_Allgather(locgidlist,maxnv,MPI_INT,globgidlist,maxnv,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(locvallist,maxnv,MPI_DOUBLE,globvallist,maxnv,MPI_DOUBLE,MPI_COMM_WORLD);

    idx = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx))) {
      int idx;
      
      int ival;
      double rval;
      void *pval;

      MEnt_Get_AttVal(mv,vertatt,&ival,&rval,&pval);

      int gid = MV_GlobalID(mv);

      switch (op) {
        case ATTOP_MAX: {
          double maxval = -1;
          for (int k = 0; k < maxnv*nproc; k++) 
            if (globgidlist[k] == gid)
              if (globvallist[k] > maxval)
                maxval = globvallist[k];
          CHECK_EQUAL(maxval,rval);
          break;
        }
        case ATTOP_MIN: {
          double minval = 1e+16;
          idx = 0;
          for (int k = 0; k < maxnv*nproc; k++) 
            if (globgidlist[k] == gid)
              if (globvallist[k] < minval)
                minval = globvallist[k];          
          CHECK_EQUAL(minval,rval);
          break;
        }
        case ATTOP_SUM: case ATTOP_AVG: {
          double sum = 0.0;
          int nval = 0;
          idx = 0;
          for (int k = 0; k < maxnv*nproc; k++) 
            if (globgidlist[k] == gid) {
              sum += globvallist[k];
              nval++;
            }
          if (op == ATTOP_SUM)
            CHECK_EQUAL(sum,rval);
          else
            CHECK_CLOSE(sum/nval,rval,1.0e-12);
        }
        default: 
          break;
      }
    } /* for each vertex */
    delete [] locgidlist;
    delete [] locvallist;
    delete [] globgidlist;
    delete [] globvallist;





    MEdge_ptr me;
    idx = 0;
    while (me = MESH_Next_Edge(mesh,&idx))
      MEnt_Set_AttVal(me,edgeatt,0,rank+1,NULL);

    MESH_GathScat1Attribute(mesh,edgeatt,op,comm);


    /* Now use plain old MPI to determine the right answer */

    int ne = MESH_Num_Edges(mesh);
    int *allne = (int *) new int [nproc];
    MPI_Allgather(&ne,1,MPI_INT,allne,1,MPI_INT,MPI_COMM_WORLD);

    int maxne = 0;
    for (int j = 0; j < nproc; j++)
      if (allne[j] > maxne)
        maxne = allne[j];
    locgidlist = new int[maxne];
    locvallist = new double [maxne];
    globgidlist = new int [maxne*nproc];
    globvallist = new double [maxne*nproc];
    for (int i = 0; i < maxne; i++) {
      locgidlist[i] = 0;
      locvallist[i] = 0.0;
    }
    for (int i = 0; i < maxne*nproc; i++) {
      globgidlist[i] = 0;
      globvallist[i] = 0.0;
    }

    j = 0;
    idx = 0; 
    while ((me = MESH_Next_Edge(mesh,&idx))) {
      locgidlist[j] = ME_GlobalID(me);
      locvallist[j] = rank+1;
      j++;
    }

    MPI_Allgather(locgidlist,maxne,MPI_INT,globgidlist,maxne,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(locvallist,maxne,MPI_DOUBLE,globvallist,maxne,MPI_DOUBLE,MPI_COMM_WORLD);

    idx = 0;
    while ((me = MESH_Next_Edge(mesh,&idx))) {
      int idx;
      
      int ival;
      double rval;
      void *pval;

      MEnt_Get_AttVal(me,edgeatt,&ival,&rval,&pval);

      int gid = ME_GlobalID(me);

      switch (op) {
        case ATTOP_MAX: {
          double maxval = -1;
          for (int k = 0; k < maxne*nproc; k++) 
            if (globgidlist[k] == gid)
              if (globvallist[k] > maxval)
                maxval = globvallist[k];
          CHECK_EQUAL(maxval,rval);
          break;
        }
        case ATTOP_MIN: {
          double minval = 1e+16;
          idx = 0;
          for (int k = 0; k < maxne*nproc; k++) 
            if (globgidlist[k] == gid)
              if (globvallist[k] < minval)
                minval = globvallist[k];          
          CHECK_EQUAL(minval,rval);
          break;
        }
        case ATTOP_SUM: case ATTOP_AVG: {
          double sum = 0.0;
          int nval = 0;
          idx = 0;
          for (int k = 0; k < maxne*nproc; k++) 
            if (globgidlist[k] == gid) {
              sum += globvallist[k];
              nval++;
            }
          if (op == ATTOP_SUM)
            CHECK_EQUAL(sum,rval);
          else
            CHECK_CLOSE(sum/nval,rval,1.0e-12);
        }
        default: 
          break;
      }
    } /* for each edge */
    delete [] locgidlist;
    delete [] locvallist;
    delete [] globgidlist;
    delete [] globvallist;


  } /* For each GathScat operation type */

  return;
} /* GathScatAttrib_2D */
 



TEST(XchngEdgeAttrib_2D) {

  int i, idx, ival, status;
  int vcol, fcol, mpid;
  int nproc, rank, dim;
  double rval;
  void *pval;
  Mesh_ptr mesh;
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
    
    if (MESH_Num_Regions(mesh0) > 0)
      fprintf(stderr,"Code is for surface meshes only. Exiting...\n");
    else if (MESH_Num_Faces(mesh0) > 0)
      dim = 2;
    else {
      fprintf(stderr,"Mesh is neither solid nor surface mesh. Exiting...\n");
      status = 0;
      CHECK(status);
    }

    status = MESH_CheckTopo(mesh0);
    CHECK(status);
  }

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

  MAttrib_ptr edgeatt;
  edgeatt = MAttrib_New(mesh,"edgeatt",DOUBLE,MEDGE);

  MEdge_ptr me;
  idx = 0;
  while (me = MESH_Next_Edge(mesh,&idx))
    MEnt_Set_AttVal(me,edgeatt,0,rank+1,NULL);
  
  MESH_XchngEdgeFaceAttrib(mesh,edgeatt,comm);

  /* Mark/count all edges on the partition boundary */
#ifdef MSTK_USE_MARKERS
  int prtnmark = MSTK_GetMarker();
#else
  MAttrib_ptr prtnatt = MAttrib_New(mesh, "prtnatt", INT, MALLTYPE);
#endif

  int ne = 0;
  idx = 0; 
  while ((me = MESH_Next_Edge(mesh,&idx))) {
    int foundowned=0, foundghost=0;
    List_ptr efaces = ME_Faces(me);
    MFace_ptr vf;
    int idx2 = 0;
    while ((vf = List_Next_Entry(efaces,&idx2))) {
      if (MF_PType(vf) == PGHOST)
        foundghost = 1;
      else
        foundowned = 1;
    }
    if (!foundghost || !foundowned) continue;

#ifdef MSTK_USE_MARKERS
    MEnt_Mark(me,prtnmark);
#else
    MEnt_Set_AttVal(me, prtnatt, 1, 0.0, NULL);
#endif
    ne++;
  }

  /* Now use plain old MPI to determine the right answer */
  
  int *allne = (int *) new int[nproc];
  MPI_Allgather(&ne,1,MPI_INT,allne,1,MPI_INT,MPI_COMM_WORLD);
  
  int maxne = 0;
  for (int j = 0; j < nproc; j++)
    if (allne[j] > maxne)
      maxne = allne[j];
  int *locgidlist = new int[maxne];
  double *locvallist = new double[maxne];
  int *globgidlist = new int[maxne*nproc];
  double *globvallist = new double[maxne*nproc];
  for (int i = 0; i < maxne; i++) {
    locgidlist[i] = 0;
    locvallist[i] = 0.0;
  }
  for (int i = 0; i < maxne*nproc; i++) {
    globgidlist[i] = 0;
    globvallist[i] = 0.0;
  }
  
  int j = 0;
  idx = 0; 
  while ((me = MESH_Next_Edge(mesh,&idx))) {
    int emarked;
#ifdef MSTK_USE_MARKERS
    emarked = MEnt_IsMarked(me, prtnmark);
#else
    MEnt_Get_AttVal(me, prtnatt, &emarked, &rval, &pval);
#endif
    if (!emarked) continue;
    locgidlist[j] = ME_GlobalID(me);
    locvallist[j] = rank+1;
    j++;
  }
  
  MPI_Allgather(locgidlist,maxne,MPI_INT,globgidlist,maxne,MPI_INT,MPI_COMM_WORLD);
  MPI_Allgather(locvallist,maxne,MPI_DOUBLE,globvallist,maxne,MPI_DOUBLE,MPI_COMM_WORLD);
  
  idx = 0;
  while ((me = MESH_Next_Edge(mesh,&idx))) { 
    int emarked;
#ifdef MSTK_USE_MARKERS
    emarked = MEnt_IsMarked(me, prtnmark);
#else
    MEnt_Get_AttVal(me, prtnatt, &emarked, &rval, &pval);
#endif
    if (!emarked) continue;

    int gid = ME_GlobalID(me);
    
    MEnt_Get_AttVal(me,edgeatt,&ival,&rval,&pval);
    
    double rval_mpi=-1e+10;
    int found = 0;
    for (int p = 0; p < nproc; p++) {
      if (p == rank) continue;
      for (int k = 0; k < maxne; k++) 
        if (globgidlist[p*maxne+k] == gid) {
          rval_mpi = globvallist[p*maxne+k];
          found = 1;
          break;
        }
      if (found) break;
    }
    CHECK_EQUAL(rval_mpi,rval);
  } /* for each edge */

  delete [] locgidlist;
  delete [] locvallist;
  delete [] globgidlist;
  delete [] globvallist;

#ifdef MSTK_USE_MARKERS
  idx = 0;
  while ((me = MESH_Next_Edge(mesh,&idx)))
    MEnt_Unmark(me,prtnmark);

  MSTK_FreeMarker(prtnmark);
#else
  MAttrib_Delete(prtnatt);
#endif
    
  return;
} /* XchngEdgeAttrib_2D */
 

}


