/* Translator between mesh formats */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _MSTK_HAVE_MPI
#include <mpi.h>
#endif

#include "MSTK.h" 

typedef enum {MSTK,GMV,EXODUSII,CGNS,VTK,STL,AVSUCD,DX,X3D} MshFmt;

int main(int argc, char *argv[]) {
  char infname[256], outfname[256];
  Mesh_ptr mesh;
  int len, ok, build_classfn=1, partition=0;
  MshFmt inmeshformat, outmeshformat;
  

  if (argc < 3) {
    fprintf(stderr,"usage: %s <--classify=yes|no> infilename outfilename\n",argv[0]);
    exit(-1);
  }

  if (argc > 3) {
    int i;
    for (i = 1; i < argc-2; i++)
      if (strncmp(argv[i],"--classify",10) == 0) {
        if (strncmp(argv[i],"--classify=y",12) == 0) 
          build_classfn = 1;
        else if (strncmp(argv[i],"--classify=n",12) == 0)
          build_classfn = 0;
      }
      else if (strncmp(argv[i],"--partition",11) == 0) {
        if (strncmp(argv[i],"--partition=y",13) == 0)
          partition = 1;
        else if (strncmp(argv[i],"--partition=n",13) == 0) 
          partition = 0;
      }
      else
        fprintf(stderr,"Unrecognized option...Ignoring\n");
  }

  strcpy(infname,argv[argc-2]);
  strcpy(outfname,argv[argc-1]);


  len = strlen(infname);
  if (len > 5 && strncmp(&(infname[len-5]),".mstk",5) == 0) {
    inmeshformat = MSTK;
    if (argc > 3 && strncmp(argv[1],"--classify=y",12) == 0)
      build_classfn = 1;
    else
      build_classfn = 0;
  }
  else if (len > 4 && strncmp(&(infname[len-4]),".gmv",4) == 0)
    inmeshformat = GMV;
  else if ((len > 4 && strncmp(&(infname[len-4]),".exo",4) == 0) ||
	   (len > 2 && strncmp(&(infname[len-2]),".g",2) == 0) ||
	   (len > 2 && strncmp(&(infname[len-2]),".e",2) == 0) ||
           (len > 4 && strncmp(&(infname[len-4]),".par",4) == 0)) 
    inmeshformat = EXODUSII;
  else if (len > 4 && strncmp(&(infname[len-4]),".vtk",4) == 0)
    inmeshformat = VTK;
  else if (len > 5 && strncmp(&(infname[len-5]),".cgns",5) == 0)
    inmeshformat = CGNS;
  else if ((len > 4 && strncmp(&(infname[len-4]),".avs",4) == 0) ||
	   (len > 4 && strncmp(&(infname[len-4]),".inp",4) == 0))
    inmeshformat = AVSUCD;
  else if (len > 4 && strncmp(&(infname[len-4]),".x3d",4) == 0)
    inmeshformat = X3D;
  else {
    fprintf(stderr,"Unrecognized input mesh format\n");
    fprintf(stderr,"Recognized input mesh formats: MSTK, GMV, ExodusII, CGNS, VTK, AVS/UCD\n");
    exit(-1);
  }

  len = strlen(outfname);
  if (len > 5 && strncmp(&(outfname[len-5]),".mstk",5) == 0)
    outmeshformat = MSTK;
  else if (len > 4 && strncmp(&(outfname[len-4]),".gmv",4) == 0)
    outmeshformat = GMV;
  else if ((len > 4 && strncmp(&(outfname[len-4]),".exo",4) == 0) ||
	   (len > 2 && strncmp(&(outfname[len-2]),".g",2) == 0) ||
	   (len > 2 && strncmp(&(outfname[len-2]),".e",2) == 0) ||
           (len > 4 && strncmp(&(outfname[len-4]),".par",4) == 0)) 
    outmeshformat = EXODUSII;
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
    fprintf(stderr,"Recognized output mesh formats: MSTK, GMV, ExodusII, CGNS, VTK, AVS/UCD\n");
    exit(-1);
  }


#ifdef MSTK_HAVE_MPI
  MPI_Init(&argc,&argv);
#endif


  MSTK_Init();


#ifdef MSTK_HAVE_MPI

  MSTK_Set_Comm(MPI_COMM_WORLD);

  int rank = MSTK_Comm_rank();
  int numprocs = MSTK_Comm_size();

#endif
  
  if (inmeshformat == MSTK) {
    mesh = MESH_New(UNKNOWN_REP);

    fprintf(stderr,"Reading file in MSTK format...");

    ok = MESH_InitFromFile(mesh,infname);

    fprintf(stderr,"Done\n");
  }
  else {
    mesh = MESH_New(F1);

    switch(inmeshformat) {
    case GMV:
      fprintf(stderr,"Importing mesh from GMV file...");
      ok = MESH_ImportFromFile(mesh,infname,"gmv");
      break;
    case EXODUSII:
      fprintf(stderr,"Importing mesh from ExodusII file...");
      ok = MESH_ImportFromFile(mesh,infname,"exodusii");
      break;
    case CGNS:
      fprintf(stderr,"Cannot import mesh from CGNS format. ");
      break;
    case VTK:
      fprintf(stderr,"Cannot import mesh from VTK format. ");
      break;
    case AVSUCD:
      fprintf(stderr,"Cannot import mesh from AVS format. ");
      break;
    case X3D:
      fprintf(stderr,"Importing mesh from X3D format...");
      ok = MESH_ImportFromFile(mesh,infname,"x3d");
      break;
    default:
      fprintf(stderr,"Cannot import from unrecognized format. ");
    }

    if (ok)
      fprintf(stderr,"Done\n");
    else {
      fprintf(stderr,"Failed\n");
      exit(-1);
    }
  }




  if (build_classfn) {
    fprintf(stderr,"Building classification information....");

    ok = MESH_BuildClassfn(mesh);  

    if (ok)
      fprintf(stderr,"Done\n");
    else {
      fprintf(stderr,"Failed\n");
      exit(-1);
    }
  }


#ifdef MSTK_HAVE_MPI

  if (partition) {

    /* Need to make sure that we have a mesh only on processor 0 */
    /* For now we cannot repartition                             */

    int dim = 3;
    if (rank > 0) {
      if (mesh && MESH_Num_Vertices(mesh) != 0) {
        fprintf(stderr,"Mesh is already distributed - cannot repartition\n");
        MPI_Finalize();
        exit(-1);
      }
    }

    if (rank == 0) {
      fprintf(stderr,"Partitioning mesh into %d parts...",numprocs);

      if (MESH_Num_Regions(mesh))
        dim = 3;
      else if (MESH_Num_Faces(mesh))
        dim = 2;
      else {
        fprintf(stderr,"Partitioning possible only for 2D or 3D meshes\n");
        exit(-1);
      }
    }
     
    int ring = 0; /* No ghost ring of elements */
    int with_attr = 1; /* Do allow exchange of attributes */
    int method = 0; /* Use Metis as the partitioner */
        
    int ok = MSTK_Mesh_Distribute(&mesh, &dim, ring, with_attr,
                                  rank, numprocs, MPI_COMM_WORLD);

    if (rank == 0) {
      if (ok)
        fprintf(stderr,"done\n");
      else {
        fprintf(stderr,"failed\n");
        exit(-1);
      }
    }
  }

#endif /* MSTK_HAVE_MPI */



  if (outmeshformat == MSTK) {
    fprintf(stderr,"Writing mesh to MSTK file...");
    MESH_WriteToFile(mesh,outfname,MESH_RepType(mesh));
    fprintf(stderr,"Done\n");
  }
  else {
    switch(outmeshformat) {
    case GMV:
      fprintf(stderr,"Exporting mesh to GMV format...");
      ok = MESH_ExportToFile(mesh,outfname,"gmv",-1,NULL,NULL);
      break;
    case EXODUSII:
      fprintf(stderr,"Exporting mesh to ExodusII format...");
      ok = MESH_ExportToFile(mesh,outfname,"exodusii",-1,NULL,NULL);
      break;
    case CGNS:
      fprintf(stderr,"Cannot export to CGNS format. ");
      break;
    case VTK:
      fprintf(stderr,"Cannot export to VTK format. ");
      break;
    case AVSUCD:
      fprintf(stderr,"Cannot export to AVS format. ");
      break;
    case X3D:
      fprintf(stderr,"Exporting mesh to FLAG/X3D format...");
      ok = MESH_ExportToFile(mesh,outfname,"flag",-1,NULL,NULL);
      break;
    case STL:
      fprintf(stderr,"Exporting mesh to STL format...");
      ok = MESH_ExportToFile(mesh, outfname,"stl",-1,NULL,NULL);
      break;
    case DX:
      fprintf(stderr,"Exporting mesh to DX format...");
      ok = MESH_ExportToDX(mesh, outfname, 1);
      break;
    default:
      fprintf(stderr,"Cannot export mesh to unrecognized format. \n");      
    }

    if (ok)
      fprintf(stderr,"Done\n");
    else {
      fprintf(stderr,"Failed\n");
      exit(-1);
    }
  }


  MESH_Delete(mesh);

#ifdef MSTK_HAVE_MPI
  MPI_Finalize();
#endif

  return 1;
}
