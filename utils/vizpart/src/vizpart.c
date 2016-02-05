#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef MSTK_HAVE_MPI
#include "mpi.h"
#endif

#include "MSTK.h"


int main(int argc, char **argv) {
  int len, ok;
  int i, j, id, idx, ncells, file_found;
  int num_parts, imethod;
  int nproc, rank;  
  Mesh_ptr mesh;
  char basename[256], meshname[256], gmvfilename[256], method[256];


  if (argc == 1) {
    fprintf(stderr,"Usage: %s --num=num_parts --method=Metis/Zoltan (default=Metis) --file=name.mstk\n",argv[0]);
    exit(-1);    
  }

  file_found = 0;
  num_parts = 4;
  strcpy(method,"Metis");
  for (i = 1; i < argc; i++) {
    if (strncmp(argv[i],"--file=",7) == 0) {
      sscanf(argv[i]+7,"%s",meshname);
      file_found = 1;
    }
    else if (strncmp(argv[i],"--num=",6) == 0)
      sscanf(argv[i]+6,"%d",&num_parts);
    else if (strncmp(argv[i],"--method=",9) == 0) 
      sscanf(argv[i]+9,"%s",method);
    else if (strncmp(argv[i],"--help",6) == 0) {
      fprintf(stderr,"Usage: %s --num=num_parts --method=Metis/Zoltan (default=Metis) --file=name.mstk\n",argv[0]);
      exit(-1);    
    }      
  }

  if (strncasecmp(method,"metis",5) == 0)
    imethod = 0;
  else if (strncasecmp(method,"zoltan",6) == 0)
    imethod = 1;
  else {
    fprintf(stderr,"vizpart: Partitioning method not recognized\n");
    exit(-1);
  }

  if (!file_found) {
    fprintf(stderr,"Must specify input filename using --file argument\n");
    exit(-1);
  }


  MPI_Init(&argc,&argv);
      
  MSTK_Init();

  MSTK_Comm comm = MPI_COMM_WORLD;

  MPI_Comm_size(comm,&nproc);
  MPI_Comm_rank(comm,&rank);

  if (imethod == 1 && nproc != num_parts) {
    fprintf(stderr,"Zoltan based partitioner: Number of processors must equal requested number of partition\n");
    exit(-1);
  }

  
  strcpy(basename,meshname);
  len = strlen(meshname);
  if (len > 5 && strncmp(&(meshname[len-5]),".mstk",5) == 0)
    basename[len-5] = '\0';
  else
    strcat(meshname,".mstk");

  mesh = MESH_New(UNKNOWN_REP);

  if (rank == 0) {

    ok = MESH_InitFromFile(mesh,meshname,comm);
    if (!ok) {
      fprintf(stderr,"Cannot file input file %s\n\n\n",meshname);
      exit(-1);
    }

    if ( (ncells = MESH_Num_Regions(mesh)) > 0) {
      printf("3D mesh with %d regions\n", ncells);
    }
    else if ( (ncells = MESH_Num_Faces(mesh)) > 0)
      printf("2D mesh with %d faces\n", ncells);
    else {
      fprintf(stderr,"Mesh is neither solid nor surface mesh. Exiting...\n");
      exit(-1);
    }
  }

  Mesh_ptr *submeshes = (Mesh_ptr*) malloc(num_parts*sizeof(Mesh_ptr));
  int *part;
  MESH_Get_Partitioning(mesh, imethod, &part, comm);

  if(rank == 0) {

    int del_inmesh = 0;
    MESH_Partition(mesh, num_parts, part, submeshes);

    idx = 0;
    if(MESH_Num_Regions(mesh)) {
      MAttrib_ptr region_part = MAttrib_New(mesh, "part", INT, MREGION);
      MRegion_ptr mr;
      while ((mr = MESH_Next_Region(mesh, &idx))) {
        id = MR_ID(mr);
	MEnt_Set_AttVal(mr,region_part,part[id-1],0,NULL);
      }
    }
    else {
      MAttrib_ptr face_part = MAttrib_New(mesh, "part", INT, MFACE);
      MFace_ptr mf;
      while ((mf = MESH_Next_Face(mesh, &idx))) {
        id = MF_ID(mf);
	MEnt_Set_AttVal(mf,face_part,part[idx-1],0,NULL);
      }
    }

    sprintf(gmvfilename,"%s_part.gmv",basename);
    MESH_ExportToGMV(mesh,gmvfilename,0,NULL,NULL,comm);
    for( i = 0; i < num_parts; i++) {
      sprintf(gmvfilename,"%s_part.%04d.%04d.gmv",basename,i,0);
      MESH_ExportToGMV(submeshes[i],gmvfilename,0,NULL,NULL,comm);
    }

    for (i = 0; i < num_parts; i++)
      MESH_Delete(submeshes[i]);
  }


  free(submeshes);
  MESH_Delete(mesh);
  free(part);

  MPI_Finalize();
  return 1;
}


  
