#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"


#include "MSTK.h"



int main(int argc, char *argv[]) {
  Mesh_ptr refmesh, mesh;
  int num, rank, ncells;
  char filename[256], basename[256], gmvfilename[256], sname[256];
  int i, j, k, idx, idx1, idx2, nfv, nsteps, dim;
  int converged, mkid, mkid2, ival;
  void *pval;
  double vxyz[3], fxyz[10][3], rval;
  double maxrho;
  MFace_ptr mf, vf, ff;
  MVertex_ptr mv, fv;
  List_ptr fverts, vfaces, ffaces;
  MAttrib_ptr rhoatt;
  int DebugWait=1;
  FILE *fp;
  
  
  if(argc == 1) {
    printf("Pseudo-advection of shock front in parallel\n");
    printf("usage: mpirun -np 4 %s meshname.mstk\n",argv[0]);
    exit(-1);
  }

  
  /*  MPI_Init(&argc,&argv); */

  MSTK_Init();

  MPI_Comm_size(MPI_COMM_WORLD,&num);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  

  fprintf(stderr,"On processor %d out of %d processors\n",rank,num);

  strcpy(filename,argv[1]);
  strcpy(basename,argv[1]);
  basename[strlen(basename)-5] = '\0';


  mesh = MESH_New(UNKNOWN_REP);
  if (rank == 0) {
    MESH_InitFromFile(mesh, filename);
    
    if ( (ncells = MESH_Num_Regions(mesh)) > 0) {
      printf("3D meshes\n");
    }
    else if ( (ncells = MESH_Num_Faces(mesh)) > 0)
      printf("2D meshes with %d faces\n", ncells);
    else {
      fprintf(stderr,"Mesh is neither solid nor surface mesh. Exiting...\n");
      exit(-1);
    }
  }

  DebugWait=0;
  while (DebugWait);
  Mesh_ptr *submeshes = (Mesh_ptr*) MSTK_malloc(num*sizeof(Mesh_ptr));
  int *part;
  MESH_Get_Partition(mesh, num, 1, rank, MPI_COMM_WORLD, &part);

  if(rank == 0) {
    MSTK_Mesh_Partition(mesh, num, part, 0, 0, submeshes);
    MAttrib_ptr g2latt = MAttrib_New(mesh,"Global2Local",POINTER,MALLTYPE);
    idx = 0;
    if(MESH_Num_Regions(mesh)) {
      MAttrib_ptr region_part = MAttrib_New(mesh, "part", INT, MREGION);
      MRegion_ptr mr;
      while (mr = MESH_Next_Region(mesh, &idx)) 
	MEnt_Set_AttVal(mr,region_part,part[idx-1],0,NULL);
    }
    else {
      MAttrib_ptr face_part = MAttrib_New(mesh, "part", INT, MFACE);
      MFace_ptr mf;
      while (mf = MESH_Next_Face(mesh, &idx)) 
	MEnt_Set_AttVal(mf,face_part,part[idx-1],0,NULL);
    }
    if(rank == 0) { 
      sprintf(gmvfilename,"%s_zoltan.gmv",basename,i,0);
      MESH_ExportToGMV(mesh,gmvfilename,0,NULL,NULL);
      for( i = 0; i < num; i++) {
	sprintf(gmvfilename,"%s_zoltan.%04d.%04d.gmv",basename,i,0);
	MESH_ExportToGMV(submeshes[i],gmvfilename,0,NULL,NULL);
	printf("here on %d  processor %d\n",i, rank);
      }
      
    }
  }
    printf("here on processor %d\n",rank);
    /*
      MSTK_free(mesh);
      for (i = 0; i < num; i++)
      MSTK_free(submeshes[i]);
    */
    /*free(part);*/
    MPI_Finalize();
    
    return 0;
}
  
