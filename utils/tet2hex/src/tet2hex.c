/* Translator between mesh formats */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _MSTK_HAVE_MPI
#include <mpi.h>
#endif

#include "MSTK.h" 

int main(int argc, char *argv[]) {
  char infname[256], outfname[256];
  Mesh_ptr mesh;
  int len, ok;
  int build_classfn=1, partition=0, weave=0, use_geometry=0, parallel_check=0;
  int check_topo=0;
  int num_ghost_layers=0, partmethod=0;
  MshFmt inmeshformat, outmeshformat;


  if (argc < 3) {
    fprintf(stderr,"\n");
    fprintf(stderr,"usage: tet2hex infilename outfilename\n\n");
    exit(-1);
  }

#ifdef MSTK_HAVE_MPI
  MPI_Init(&argc,&argv);
#endif


  MSTK_Init();


  int rank=0, numprocs=1;
#ifdef MSTK_HAVE_MPI

  MSTK_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&numprocs);

#else
  MSTK_Comm comm = NULL;
#endif
  

  strcpy(infname,argv[argc-2]);
  strcpy(outfname,argv[argc-1]);


  len = strlen(infname);
  if (len > 5 && strncmp(&(infname[len-5]),".mstk",5) == 0)
    inmeshformat = MSTK;
  else if (len > 4 && strncmp(&(infname[len-4]),".gmv",4) == 0)
    inmeshformat = GMV;
  else if ((len > 4 && strncmp(&(infname[len-4]),".exo",4) == 0) ||
	   (len > 2 && strncmp(&(infname[len-2]),".g",2) == 0) ||
	   (len > 2 && strncmp(&(infname[len-2]),".e",2) == 0))
      inmeshformat = EXODUSII;
  else if (len > 4 && strncmp(&(infname[len-4]),".par",4) == 0) 
    inmeshformat = NEMESISI;
  else if ((len > 4 && strncmp(&(infname[len-4]),".avs",4) == 0) ||
	   (len > 4 && strncmp(&(infname[len-4]),".inp",4) == 0))
    inmeshformat = AVSUCD;
  else if (len > 4 && strncmp(&(infname[len-4]),".x3d",4) == 0)
    inmeshformat = X3D;
  else {
    fprintf(stderr,"Unrecognized input mesh format\n");
    fprintf(stderr,"Recognized input mesh formats: MSTK, GMV, ExodusII\n");
    exit(-1);
  }

  len = strlen(outfname);
  if (len > 5 && strncmp(&(outfname[len-5]),".mstk",5) == 0)
    outmeshformat = MSTK;
  else if (len > 4 && strncmp(&(outfname[len-4]),".gmv",4) == 0)
    outmeshformat = GMV;
  else if ((len > 4 && strncmp(&(outfname[len-4]),".exo",4) == 0) ||
	   (len > 2 && strncmp(&(outfname[len-2]),".g",2) == 0) ||
	   (len > 2 && strncmp(&(outfname[len-2]),".e",2) == 0))
    outmeshformat = EXODUSII;
  else if ((len > 4 && strncmp(&(outfname[len-4]),".par",4) == 0))
    outmeshformat = NEMESISI;
  else if (len > 4 && strncmp(&(outfname[len-4]),".vtk",4) == 0)
    outmeshformat = VTK;
  else if (len > 5 && strncmp(&(outfname[len-5]),".cgns",5) == 0)
    outmeshformat = CGNS;
  else if ((len > 4 && strncmp(&(outfname[len-4]),".avs",4) == 0) ||
	   (len > 4 && strncmp(&(outfname[len-4]),".inp",4) == 0))
    outmeshformat = AVSUCD;
  else if (len > 4 && strncmp(&(outfname[len-4]),".x3d",4) == 0)
    outmeshformat = X3D;
  else if (len > 3 && strncmp(&(outfname[len-3]),".dx",3) == 0)
    outmeshformat = DX;
  else {
    fprintf(stderr,"Unrecognized output mesh format\n");
    fprintf(stderr,"Recognized output mesh formats: MSTK, GMV, ExodusII, NemesisI, CGNS, VTK, AVS/UCD\n");
    exit(-1);
  }


  if (inmeshformat == MSTK) {
    mesh = MESH_New(UNKNOWN_REP);

    if (rank == 0) {

      fprintf(stderr,"Reading file in MSTK format...");

      ok = MESH_InitFromFile(mesh,infname,comm);

      fprintf(stderr,"Done\n");
    }
  }
  else {
    mesh = MESH_New(F1);

    int opts[5]={0,0,0,0,0};

    ok = 0;
    switch(inmeshformat) {
    case GMV:
      if (rank == 0)
        fprintf(stderr,"Importing mesh from GMV file...");
      ok = MESH_ImportFromFile(mesh,infname,"gmv",opts,comm);
      break;
    case EXODUSII:
      if (rank == 0)
        fprintf(stderr,"Importing mesh from ExodusII file...");
      opts[0] = 0; /* don't partition while importing - do it later */
      opts[1] = 0;
      opts[2] = 0; /* no ghost layers */  
      opts[3] = 0;
      ok = MESH_ImportFromFile(mesh,infname,"exodusii",opts,comm);
      break;
    case X3D:
      if (rank == 0)
        fprintf(stderr,"Importing mesh from X3D format...");
      ok = MESH_ImportFromFile(mesh,infname,"x3d",opts,comm);
      break;
    default:
      fprintf(stderr,"Cannot import from unrecognized format. ");
    }

    if (ok) {
      if (rank == 0)
        fprintf(stderr,"Done\n");
    }
    else {
      if (rank == 0)
        fprintf(stderr,"Failed\n");
      exit(-1);
    }
  }


  // Check that all elements are tets

  int alltets = 1;
  int idx = 0;
  MRegion_ptr mr;
  while ((mr = MESH_Next_Region(mesh, &idx))) {
    List_ptr rfaces = MR_Faces(mr);
    if (List_Num_Entries(rfaces) != 4) {  // In theory its possible to have a
                                          // 4-sided region that is not a 
                                          // tetrahedron, but we ignore this
      alltets = 0;
      break;
    }
  }


  // Call the Tet to Hex conversion routine. Each edge is split in two
  // edges at the midpoint. Each triangular face is split into three
  // quads, each of which is formed by the centroid of the face, two
  // midpoints of edges and a vertex of the triangle; Each tet is
  // split into 4 tets, each of which is formed by the centroids of
  // two adjacent faces, the mid-point of the common edge and a vertex
  // at one of the ends of the edge.

  fprintf(stderr, "Converting to Hex mesh by splitting tets...\n");
  Mesh_ptr hexmesh;
  ok = MESH_Tet2Hex(mesh, &hexmesh);
  if (!ok) {
    fprintf(stderr, "Tet->Hex conversion failed. Exiting...\n");
    exit(-1);
  } else
    fprintf(stderr, "Done\n");

  if (outmeshformat == MSTK) {
    if (rank == 0)
      fprintf(stderr,"Writing mesh to MSTK file...");
    MESH_WriteToFile(hexmesh, outfname, MESH_RepType(mesh), comm);
    fprintf(stderr,"Done\n");
  }
  else {
    switch(outmeshformat) {
    case GMV:
      if (rank == 0)
        fprintf(stderr,"Exporting mesh to GMV format...");
      ok = MESH_ExportToFile(hexmesh, outfname, "gmv", -1, NULL, NULL, comm);
      break;
    case EXODUSII:case NEMESISI:
      if (rank == 0)
        fprintf(stderr,"Exporting mesh to ExodusII/NemesisI format...");
      ok = MESH_ExportToFile(hexmesh, outfname, "exodusii", -1, NULL, NULL,
			     comm);
      break;
    case X3D:
      if (rank == 0)
        fprintf(stderr,"Exporting mesh to FLAG/X3D format...");
      ok = MESH_ExportToFile(hexmesh, outfname, "x3d", -1, NULL, NULL, comm);
      break;
    default:
      if (rank == 0)
        fprintf(stderr,"Cannot export mesh to unrecognized format. \n");      
    }

    if (rank == 0) {
      if (ok)
        fprintf(stderr,"Done\n");
      else {
        fprintf(stderr,"Failed\n");
        exit(-1);
      }
    }
  }


  MESH_Delete(mesh);
  MESH_Delete(hexmesh);
 
  // while (1);

#ifdef MSTK_HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
