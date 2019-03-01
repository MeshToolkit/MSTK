/* Translator between mesh formats */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _MSTK_HAVE_MPI
#include <mpi.h>
#endif

#include "MSTK.h"

MshFmt getFormat(char *filename) {
  int len = strlen(filename);
  if (len > 5 && strncmp(&(filename[len-5]),".mstk",5) == 0)
    return MSTK;
  else if (len > 4 && strncmp(&(filename[len-4]),".gmv",4) == 0)
    return GMV;
  else if ((len > 4 && strncmp(&(filename[len-4]),".exo",4) == 0) ||
	   (len > 2 && strncmp(&(filename[len-2]),".g",2) == 0) ||
	   (len > 2 && strncmp(&(filename[len-2]),".e",2) == 0))
    return EXODUSII;
  else if (len > 4 && strncmp(&(filename[len-4]),".par",4) == 0) 
    return NEMESISI;
  else if ((len > 4 && strncmp(&(filename[len-4]),".avs",4) == 0) ||
	   (len > 4 && strncmp(&(filename[len-4]),".inp",4) == 0))
    return AVSUCD;
  else if (len > 4 && strncmp(&(filename[len-4]),".x3d",4) == 0)
    return X3D;
  else if (len > 4 && strncmp(&(filename[len-4]),".vtk",4) == 0)
    return VTK;
  else if (len > 5 && strncmp(&(filename[len-5]),".cgns",5) == 0)
    return CGNS;
  else {
    fprintf(stderr,"Unrecognized mesh format\n");
    fprintf(stderr,"Recognized mesh formats: MSTK, GMV, ExodusII, NemesisI, AVS/UCD, VTK, CGNS\n");
    exit(-1);
  }
}

int main(int argc, char *argv[]) {
  char infname[256], outfname[256];
  int len, ok;
  int build_classfn=1, partition=-1, weave=-1, use_geometry=0, parallel_check=0;
  int check_topo=0;
  int num_ghost_layers=0, partmethod=0;
  MshFmt inmeshfmt, outmeshfmt;
  FILE *fp;

  if (argc < 3) {
    fprintf(stderr,"\n");
    fprintf(stderr,"usage: meshconvert <--classify=0|n|1|y|2> <--partition=y|1|n|0> <--partition-method=0|1|2> <--parallel-check=y|1|n|0> <--weave=y|1|n|0> <--num-ghost-layers=?> <--check-topo=y|1|n|0> infilename outfilename\n\n");
    fprintf(stderr,"partition-method = 0, METIS\n");
    fprintf(stderr,"                 = 1, ZOLTAN with GRAPH partioning\n");
    fprintf(stderr,"                 = 2, ZOLTAN with RCB partitioning\n");
    fprintf(stderr,"Choose 2 if you want to avoid partitioning models\n");
    fprintf(stderr,"with high aspect ratio along the short directions\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"weave = 0, Do not weave distributed meshes for inter-processor connectivity\n");
    fprintf(stderr,"      = 1, Weave distributed meshes for inter-processor connectivity\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"classify = 0/n, No mesh classification information derived\n");
    fprintf(stderr,"         = 1/y, Mesh entity classification derived from material IDs\n");
    fprintf(stderr,"         = 2, As in 1 with additional inference from boundary geometry\n");
    fprintf(stderr,"CLASSIFICATION: Relationship of mesh entities to geometric model/domain\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"check-topo = 0/n, No checking of topological coonsistency of classification info\n");
    fprintf(stderr,"           = 1/y, Check topological consistency of classification info\n");
    fprintf(stderr,"\n");
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

  if (rank == 0) {
    fprintf(stderr,"\nApp to convert unstructured meshes between formats\n");
    fprintf(stderr,"Contact: Rao Garimella (rao@lanl.gov)\n\n");
  }

  if (argc > 3) {
    int i;
    for (i = 1; i < argc-2; i++) {
      if (strncmp(argv[i],"--classify",10) == 0) {
        if (strncmp(argv[i]+11,"y",1) == 0 ||
            strncmp(argv[i]+11,"1",1) == 0)
          build_classfn = 1;
        else if (strncmp(argv[i]+11,"n",1) == 0 ||
                 strncmp(argv[i]+11,"0",1) == 0)
          build_classfn = 0;
        else if (strncmp(argv[i]+11,"2",1) == 0) {
          build_classfn = 1;
          use_geometry = 1;
        }
        else
          MSTK_Report("meshconvert",
                      "--classify option should be 0, 1, 2, y or n",
                      MSTK_FATAL);          
           
      }
      else if (strncmp(argv[i],"--partition=",12) == 0) {
        if (strncmp(argv[i]+12,"y",1) == 0 ||
            strncmp(argv[i]+12,"1",1) == 0)
          partition = 1;  // new method for partitioning serial file
        else if (strncmp(argv[i]+12,"o",1) == 0 ||
                 strncmp(argv[i]+12,"2",1) == 0) 
          partition = 2;  // old method for partitioning serial file
        else if (strncmp(argv[i]+12,"n",1) == 0 ||
                 strncmp(argv[i]+12,"0",1) == 0) 
          partition = 0;  // no partitioning
        else
          MSTK_Report("meshconvert",
                      "--partition option should be y, n, 1 or 0",
                      MSTK_FATAL);          
      }
      else if (strncmp(argv[i],"--partition-method",18) == 0 ||
               strncmp(argv[i],"--partition_method",18) == 0) {
        sscanf(argv[i]+19,"%d",&partmethod);
        partition = 1;
      }
      else if (strncmp(argv[i],"--parallel-check",16) == 0 ||
               strncmp(argv[i],"--parallel_check",16) == 0) {
        if (strncmp(argv[i]+17,"y",1) == 0 ||
            strncmp(argv[i]+17,"1",1) == 1)
          parallel_check = 1;
        else if (strncmp(argv[i]+17,"n",13) == 0 ||
                 strncmp(argv[i]+17,"0",13) == 0) 
          parallel_check = 0;
      }
      else if (strncmp(argv[i],"--weave",7) == 0) {
        if (strncmp(argv[i],"--weave=y",11) == 0 ||
            strncmp(argv[i],"--weave=1",11) == 0)
          weave = 1;
        else if (strncmp(argv[i],"--weave=n",11) == 0 ||
                 strncmp(argv[i],"--weave=0",11) == 0) 
          weave = 0;
        else
          MSTK_Report("meshconvert",
                      "--weave option should be y, n, 1 or 0",
                      MSTK_FATAL);          
           
      }
      else if (strncmp(argv[i],"--num_ghost_layers",18) == 0) {
        sscanf(argv[i]+19,"%d",&num_ghost_layers);
      }
      else if (strncmp(argv[i],"--check-topo",12) == 0) {
        if (strncmp(argv[i],"--check-topo=y",14) == 0 ||
            strncmp(argv[i],"--check-topo=1",14) == 0)
          check_topo=1;
      }
      else
        fprintf(stderr,"Unrecognized option...Ignoring\n");
    }
  }

  /* what format is the input mesh file in? */
  strcpy(infname,argv[argc-2]);
  inmeshfmt = getFormat(infname);
  if (inmeshfmt == EXODUSII && weave == 1) {
    inmeshfmt = NEMESISI;
  }

  /* what format should the output mesh file be in? */
  strcpy(outfname,argv[argc-1]);
  outmeshfmt = getFormat(outfname);



  /* Do we have single or multiple input files? */
  int serial_file = 0, parallel_file = 0;
  if ((fp = fopen(infname,"r"))) {
    serial_file = 1;
    fclose(fp);
  }

#ifdef MSTK_HAVE_MPI
  if (numprocs > 1) {
    /* Assume that all parallel files are of the form
       "base.ext.numprocs.rank" or "base.ext.rank"*/

    int filecount = 0;
    for (int r = 0; r < numprocs; r++) {
      char tmpfname[256];
      sprintf(tmpfname,"%s.%-d.%-d",infname,numprocs,r);
      if ((fp = fopen(tmpfname,"r"))) {
	filecount++;
	fclose(fp);
      }
    }
    parallel_file = (filecount == numprocs) ? 1 : 0;


    if (serial_file && parallel_file) {
      /* If both serial and parallel files are present, warn the user
	 if the defaults don't make sense */
    
      if (partition > 0 && weave > 0) { /* conflicting options were specified */
	MSTK_Report("meshconvert",
		    "Have serial and parallel input files and BOTH partition and weave options were specified. Pick one: --partition=1|2 will partiton the serial file and --weave=1|y will weave the parallel files",
		    MSTK_FATAL);
      } else if ((partition == -1) && (weave == -1)) {
	if (rank == 0)
	  MSTK_Report("meshconvert",
		      "Have serial and parallel input files but NEITHER partition nor weave options were specifeid. Partitioning serial file", MSTK_WARN);
	partition = 1;
      }
      
    } else if (serial_file || parallel_file) {
      if (serial_file && partition == -1)  /* serial file and partition option not set */
	partition = 1;
      else if (parallel_file && weave == -1)  /* parallel files and weave option not set */
	weave = 1;
    } else {
      MSTK_Report("meshconvert",
		  "Did not find input file",MSTK_FATAL);
    }
  }
#endif


  /* now read the mesh */

  Mesh_ptr mesh = NULL;
  int opts[5]={0,0,0,0,0};
  
  if (serial_file) {
    if (inmeshfmt == EXODUSII) {
      /* Exodus II read is special - we have a way of doing a faster reading
	 the serial file and partitioning */
      
      if (rank == 0)
	fprintf(stderr,"Importing mesh from ExodusII file...");
      opts[0] = (partition > 0) ? 1 : 0;
      opts[1] = (partition > 0) ? partition : 0;
      opts[2] = 1;
      opts[3] = partmethod;

      mesh = MESH_New(F1);
      ok = MESH_ImportFromFile(mesh,infname,"exodusii",opts,comm);

    } else {

      Mesh_ptr serial_mesh = NULL;


      if (rank == 0) {
	ok = 0;
	switch(inmeshfmt) {
	case MSTK: {
	  serial_mesh = MESH_New(UNKNOWN_REP);
	  fprintf(stderr,"Reading file in MSTK format...");
	  ok = MESH_InitFromFile(serial_mesh,infname,comm);
	  fprintf(stderr,"Done\n");
	  break;
	}
	case GMV: {
	  fprintf(stderr,"Importing mesh from GMV file...");
	  serial_mesh = MESH_New(F1);
	  ok = MESH_ImportFromFile(serial_mesh,infname,"gmv",opts,comm);
	  break;
	}
	case X3D: {
	  fprintf(stderr,"Importing mesh from X3D format...");
	  serial_mesh = MESH_New(F1);
	  ok = MESH_ImportFromFile(serial_mesh,infname,"x3d",opts,comm);
	  break;
	}
	case CGNS: case VTK: case AVSUCD: 
	  if (rank == 0)
	    MSTK_Report("meshconvert","Cannot import mesh from this format. ",
			MSTK_FATAL);
	  break;
	default:
	  if (rank == 0)
	    MSTK_Report("meshconvert","Cannot import from unrecognized format. ",
			MSTK_FATAL);
	}

	if (ok)
	  fprintf(stderr,"Done\n");
	else {
	  fprintf(stderr,"Failed\n");
	  exit(-1);
	}
      }

      if (partition > 0) {
	int ring = 1; /* 1 layer of ghost elements */
	int with_attr = 1; /* Do allow exchange of attributes */
	int del_inmesh = 1; /* Delete input mesh after partitioning */

	int dim;
	if (rank == 0)
	  dim = MESH_Num_Regions(serial_mesh) ? 3 : 2;

#ifdef MSTK_HAVE_MPI
	MPI_Bcast(&dim, 1, MPI_INT, 0, comm);
	int ok = MSTK_Mesh_Distribute(serial_mesh, &mesh, &dim, ring, with_attr,
				      partmethod, del_inmesh, comm);
#else
        MSTK_Report("meshconvert",
                    "Request for partitioning in serial run - use mpirun",
                    MSTK_FATAL);
#endif
      } else
	mesh = serial_mesh;
    }

    /* Figure out geometric classification */
    
    if (build_classfn) {  /* works correctly only for un-partitioned meshes */
      if (rank == 0) {
	fprintf(stderr,"Building classification information....");
	
	ok = MESH_BuildClassfn(mesh,use_geometry);  

	if (ok)
	  fprintf(stderr,"Done\n");
	else {
	  fprintf(stderr,"Failed\n");
	  exit(-1);
	}
      }
    }
    
    /* Check that the imported serial mesh topology is ok */

    if (check_topo) {  /* works correctly only for un-partitioned meshes */
      if (rank == 0) {
	fprintf(stderr,"Checking mesh topology....");

	ok = MESH_CheckTopo(mesh);
	
	if (ok)
	  fprintf(stderr,"Done\n");
	else {
	  fprintf(stderr,"Failed\n");
	  exit(-1);
	}
      }
    }
  } else {  /* parallel files */

    char parfilename[256];
    ok = 0;
    switch(inmeshfmt) {
    case MSTK: {
      if (rank == 0)
	fprintf(stderr,"Reading file in MSTK format...");

      sprintf(parfilename,"%s.%-d.%-d",infname,numprocs,rank);
      mesh = MESH_New(UNKNOWN_REP);
      
      ok = MESH_InitFromFile(mesh,parfilename,comm);
      
      break;
    }
    case GMV: {
      if (rank == 0)
	fprintf(stderr,"Importing mesh from GMV file...");

      sprintf(parfilename,"%s.%05d",infname,rank+1);
      mesh = MESH_New(F1);
      
      ok = MESH_ImportFromFile(mesh,parfilename,"gmv",opts,comm);

      break;
    }
    case NEMESISI: {
      if (rank == 0)
	fprintf(stderr,"Reading file in MSTK format...");

      sprintf(parfilename,"%s.%-d.%-d",infname,numprocs,rank);
      mesh = MESH_New(UNKNOWN_REP);
      
      ok = MESH_ImportFromFile(mesh,parfilename,"nenesisi",opts,comm);
      
      break;
    }
    case CGNS: {
      if (rank == 0)
	fprintf(stderr,"Cannot import mesh from CGNS format. ");
      break;
    }
    case VTK: {
      if (rank == 0)
	fprintf(stderr,"Cannot import mesh from VTK format. ");
      break;
    }
    case AVSUCD: {
      if (rank == 0)
	  fprintf(stderr,"Cannot import mesh from AVS format. ");
      break;
    }
    case X3D: {
      if (rank == 0)
	fprintf(stderr,"Importing mesh from X3D format...");

      sprintf(parfilename,"%s.%05d",infname,rank+1);
      mesh = MESH_New(F1);

      ok = MESH_ImportFromFile(mesh,parfilename,"x3d",opts,comm);

      break;
    }
    default:
      fprintf(stderr,"Cannot import from unrecognized format. ");
    }
    
    if (ok)
      fprintf(stderr,"Done\n");
    else {
      fprintf(stderr,"Failed\n");
      exit(-1);
    }

#ifdef MSTK_HAVE_MPI

    if (weave == 1) {

      int input_type = 0;
      if (inmeshfmt == NEMESISI)
	input_type = 1;
      else if (inmeshfmt == X3D)
	input_type = 2;

      int num_ghost_layers = 1;
      int dim = MESH_Num_Regions(mesh) ? 3 : 2;
      MSTK_Weave_DistributedMeshes(mesh, dim, num_ghost_layers, input_type,
				   comm);

    }
    
    if (numprocs > 1 && parallel_check == 1) {

      /* Do a parallel consistency check too */
      /* Not checking the mesh geometry here */
      
      ok = ok && MESH_Parallel_Check(mesh,comm);
      
      int allok = 0;
      MPI_Reduce(&ok,&allok,1,MPI_INT,MPI_MIN,0,comm);
      
      if (rank == 0 && allok)
	fprintf(stderr,"Parallel checks passed...\n");
    }
#endif

  }  /* if (serial_file) {...} else {...} */

    

    

  if (outmeshfmt == MSTK) {
    if (rank == 0)
      fprintf(stderr,"Writing mesh to MSTK file...");
    MESH_WriteToFile(mesh,outfname,MESH_RepType(mesh),comm);
    fprintf(stderr,"Done\n");
  }
  else {
    switch(outmeshfmt) {
    case GMV:
      if (rank == 0)
        fprintf(stderr,"Exporting mesh to GMV format...");
      ok = MESH_ExportToFile(mesh,outfname,"gmv",-1,NULL,NULL,comm);
      break;
    case EXODUSII:case NEMESISI:
      if (rank == 0)
        fprintf(stderr,"Exporting mesh to ExodusII/NemesisI format...");
      ok = MESH_ExportToFile(mesh,outfname,"exodusii",-1,NULL,NULL,comm);
      break;
    case CGNS:
      if (rank == 0)
        fprintf(stderr,"Cannot export to CGNS format. ");
      break;
    case VTK:
      if (rank == 0)
        fprintf(stderr,"Cannot export to VTK format. ");
      break;
    case AVSUCD:
      if (rank == 0)
        fprintf(stderr,"Cannot export to AVS format. ");
      break;
    case X3D:
      if (rank == 0)
        fprintf(stderr,"Exporting mesh to FLAG/X3D format...");
      ok = MESH_ExportToFile(mesh,outfname,"x3d",-1,NULL,NULL,comm);
      break;
    case STL:
      if (rank == 0)
        fprintf(stderr,"Exporting mesh to STL format...");
      ok = MESH_ExportToFile(mesh, outfname,"stl",-1,NULL,NULL,comm);
      break;
    case DX:
      if (rank == 0)
        fprintf(stderr,"Exporting mesh to DX format...");
      ok = MESH_ExportToDX(mesh, outfname, 1);
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
 
#ifdef MSTK_HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
