/* Translator between mesh formats */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _MSTK_HAVE_MPI
#include <mpi.h>
#endif

#include "MSTK.h" 

typedef enum {MSTK,GMV,EXODUSII,NEMESISI,CGNS,VTK,STL,AVSUCD,DX,X3D} MshFmt;

int main(int argc, char *argv[]) {
  char infname[256], outfname[256];
  Mesh_ptr mesh;
  int len, ok, build_classfn=1, partition=0, weave=0, use_geometry=0;
  int num_ghost_layers=0, partmethod=0;
  MshFmt inmeshformat, outmeshformat;


  if (argc < 3) {
    fprintf(stderr,"\n");
    fprintf(stderr,"usage: meshconvert <--classify=0|n|1|y|2> <--partition=y|1|n|0> <--partition-method=0|1|2> <--weave=y|1|n|0> <--num-ghost-layers=?> infilename outfilename\n\n");
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
  

  if (argc > 3) {
    int i;
    for (i = 1; i < argc-2; i++) {
      if (strncmp(argv[i],"--classify",10) == 0) {
        if (strncmp(argv[i],"--classify=y",12) == 0 ||
            strncmp(argv[i],"--classify=1",12) == 0)
          build_classfn = 1;
        else if (strncmp(argv[i],"--classify=n",12) == 0 ||
                 strncmp(argv[i],"--classify=0",12) == 0)
          build_classfn = 0;
        else if (strncmp(argv[i],"--classify=2",12) == 0) {
          build_classfn = 1;
          use_geometry = 1;
        }
        else
          MSTK_Report("meshconvert",
                      "--classify option should be 0, 1, 2, y or n",
                      MSTK_FATAL);          
           
      }
      else if (strncmp(argv[i],"--partition",11) == 0 &&
               strncmp(argv[i],"--partition-method",18) != 0) {
        if (strncmp(argv[i],"--partition=y",13) == 0 ||
            strncmp(argv[i],"--partition=1",13) == 1)
          partition = 1;
        else if (strncmp(argv[i],"--partition=n",13) == 0 ||
                 strncmp(argv[i],"--partition=0",13) == 1) 
          partition = 0;
        else 
          MSTK_Report("meshconvert",
                      "--partition option should be y, n, 1 or 0",
                      MSTK_FATAL);          
      }
      else if (strncmp(argv[i],"--partition-method",18) == 0) {
        sscanf(argv[i]+19,"%d",&partmethod);
        partition = 1;
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
      else
        fprintf(stderr,"Unrecognized option...Ignoring\n");
    }

    /* If running on multiple processors, assume that you either partition
       a serial mesh or weave together distributed meshes */

    if (numprocs > 1 && weave == 0 && partition == 0) 
      partition = 1; 
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
    fprintf(stderr,"Recognized input mesh formats: MSTK, GMV, ExodusII, NemesisI, AVS/UCD\n");
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
      fprintf(stderr,"Importing mesh from GMV file...");
      ok = MESH_ImportFromFile(mesh,infname,"gmv",opts,comm);
      break;
    case EXODUSII:
      fprintf(stderr,"Importing mesh from ExodusII file...");
      opts[0] = 0; /* don't partition while importing - do it later */
      opts[1] = 0;
      opts[2] = 0; /* no ghost layers */  
      opts[3] = 0;
      ok = MESH_ImportFromFile(mesh,infname,"exodusii",opts,comm);
      break;
    case NEMESISI:
      fprintf(stderr,"Importing mesh from NemesisI file...");
      opts[0] = weave ? 1 : 0; 
      opts[1] = num_ghost_layers; /* no ghost layers */      
      ok = MESH_ImportFromFile(mesh,infname,"nemesisi",opts,comm);      
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
      ok = MESH_ImportFromFile(mesh,infname,"x3d",opts,comm);
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

    ok = MESH_BuildClassfn(mesh,use_geometry);  

    if (ok)
      fprintf(stderr,"Done\n");
    else {
      fprintf(stderr,"Failed\n");
      exit(-1);
    }
  }

  /* Check that the imported mesh is ok */

  ok = MESH_CheckTopo(mesh);


#ifdef MSTK_HAVE_MPI

  if (partition) {

    /* Need to make sure that we have a mesh only on processor 0 */
    /* For now we cannot repartition                             */

    int prepartitioned = 0, mesh_present = 0;

    int dim = 3;
    if (rank > 0) {
      if (mesh && MESH_Num_Vertices(mesh) != 0)
        mesh_present = 1;
    }
    
    MPI_Reduce(&mesh_present,&prepartitioned,1,MPI_INT,MPI_MAX,0,comm);
    MPI_Bcast(&prepartitioned,1,MPI_INT,0,comm);

    if (!prepartitioned) {
    
      Mesh_ptr mesh0=NULL;
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
        
        mesh0 = mesh;
        mesh = NULL;
      }
      
      int ring = 0; /* No ghost ring of elements */
      int with_attr = 1; /* Do allow exchange of attributes */
      
      int ok = MSTK_Mesh_Distribute(mesh0, &mesh, &dim, ring, with_attr,
                                    partmethod, comm);
      
      if (rank == 0) {
        if (ok)
          fprintf(stderr,"done\n");
        else {
          fprintf(stderr,"failed\n");
          exit(-1);
        }
      }

      if (rank == 0)
        MESH_Delete(mesh0);

    }    
  }

#endif /* MSTK_HAVE_MPI */


#ifdef MSTK_HAVE_MPI
  if (numprocs > 1) {

    /* Do a parallel consistency check too */

    ok = ok && MESH_Parallel_Check(mesh,comm);
    
  }
#endif

  /* Not checking the mesh geometry here */


  if (outmeshformat == MSTK) {
    fprintf(stderr,"Writing mesh to MSTK file...");
    MESH_WriteToFile(mesh,outfname,MESH_RepType(mesh),comm);
    fprintf(stderr,"Done\n");
  }
  else {
    switch(outmeshformat) {
    case GMV:
      fprintf(stderr,"Exporting mesh to GMV format...");
      ok = MESH_ExportToFile(mesh,outfname,"gmv",-1,NULL,NULL,comm);
      break;
    case EXODUSII:case NEMESISI:
      fprintf(stderr,"Exporting mesh to ExodusII format...");
      ok = MESH_ExportToFile(mesh,outfname,"exodusii",-1,NULL,NULL,comm);
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
      ok = MESH_ExportToFile(mesh,outfname,"x3d",-1,NULL,NULL,comm);
      break;
    case STL:
      fprintf(stderr,"Exporting mesh to STL format...");
      ok = MESH_ExportToFile(mesh, outfname,"stl",-1,NULL,NULL,comm);
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
 
  // while (1);

#ifdef MSTK_HAVE_MPI
  MPI_Finalize();
#endif

  return 1;
}
