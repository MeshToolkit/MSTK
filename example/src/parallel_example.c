#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

#include "MSTK.h"



int main(int argc, char *argv[]) {
  Mesh_ptr refmesh, mesh, mymesh;
  int num, rank;
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
    printf("usage: mpirun -np %s 4 meshname.mstk\n",argv[0]);
    exit(-1);
  }

  
  MPI_Init(&argc,&argv);

  MSTK_Init(MPI_COMM_WORLD);

  MPI_Comm_size(MPI_COMM_WORLD,&num);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  

  fprintf(stderr,"On processor %d out of %d processors\n",rank,num);

  strcpy(filename,argv[1]);
  strcpy(basename,argv[1]);
  basename[strlen(basename)-5] = '\0';


  mesh = MESH_New(UNKNOWN_REP);

  if (rank == 0) {
    MESH_InitFromFile(mesh, filename);
    
    if (MESH_Num_Regions(mesh) > 0) {
      fprintf(stderr,"Code is for surface meshes only. Exiting...\n");
    }
    else if (MESH_Num_Faces(mesh) > 0)
      dim = 2;
    else {
      fprintf(stderr,"Mesh is neither solid nor surface mesh. Exiting...\n");
      exit(-1);
    }
    mymesh = mesh;
  }

  DebugWait=0;
  while (DebugWait);

  int ring = 1; /* One ring ghosts */
  int with_attr = 1; /* Do allow exchange of attributes */
  int method = 0; /* Use Metis as the partitioner */
  MSTK_Mesh_Distribute(&mymesh, &dim, ring, with_attr, method, 
                       rank, num, MPI_COMM_WORLD);

  /* Can just do this too */

  /*
  MSTK_Mesh_Read_Distribute(&mymesh,filename,2,1,1,rank,num,MPI_COMM_WORLD);
  */


  /* Simulate a diagonal shock front moving across the mesh */

  rhoatt = MAttrib_New(mymesh,"density",DOUBLE,MFACE);
  

  mkid = MSTK_GetMarker();
  mkid2 = MSTK_GetMarker();


  /* Find the cell(s) connected to the upper left corner vertex */

  idx = 0; 
  while ((mv = MESH_Next_Vertex(mymesh,&idx))) {
    MV_Coords(mv,vxyz);

    if (vxyz[0] == 0.0 && vxyz[1] == 1.0) {

      vfaces = MV_Faces(mv);
      idx1 = 0;
      while ((vf = List_Next_Entry(vfaces,&idx1)))
	MEnt_Set_AttVal(vf,rhoatt,0,10.0,NULL);
      List_Delete(vfaces);

      break;
    }
  }



  sprintf(gmvfilename,"%s.gmv.%04d.%04d",basename,rank,0);

  MESH_ExportToGMV(mymesh,gmvfilename,0,NULL,NULL);




  /* Update attributes across processors */
  
  fprintf(stderr,"Updating attributes on proc %-d...",rank);
  MSTK_UpdateAttr(mymesh, rank, num, MPI_COMM_WORLD);
  fprintf(stderr,"done\n");


  sprintf(sname,"out%-3d",rank);
  fp = fopen(sname,"w");
  
  /* Propagate the density value through the mesh */

  nsteps = 10;
  for (i = 0; i < nsteps; i++) {

    if (rank == 0) fprintf(stderr,"Step %-d\n",i+1);
    fprintf(fp,"i = %-d\n",i+1);

    idx = 0; 
    while ((mf = MESH_Next_Face(mymesh,&idx))) {
      if (MEnt_PType(mf) == PGHOST) continue;

      ffaces = List_New(10);
      
      fverts = MF_Vertices(mf,1,0);
      idx1 = 0;
      while ((fv = List_Next_Entry(fverts,&idx1))) {
	
	vfaces = MV_Faces(fv);
	idx2 = 0;
	while ((vf = List_Next_Entry(vfaces,&idx2))) {
	  if (!MEnt_IsMarked(vf,mkid)) {
	    List_Add(ffaces,vf);
	    MEnt_Mark(vf,mkid);
	  }
	}
	List_Delete(vfaces);
	
      }
      List_Delete(fverts);
      
      List_Unmark(ffaces,mkid);

      maxrho = 0.0;
      idx1 = 0;
      while ((ff = List_Next_Entry(ffaces,&idx1))) {
	MEnt_Get_AttVal(ff,rhoatt,&ival,&rval,&pval);
	if (rval > maxrho)
	  maxrho = rval;
      }


      fprintf(fp,"MF %-d ffaces: ",MF_ID(mf));
      idx1 = 0;
      while ((ff = List_Next_Entry(ffaces,&idx1))) {
	fprintf(fp," %-d",MF_ID(ff));
	if (MEnt_PType(ff) == PGHOST)
	  fprintf(fp,"G");
	MEnt_Get_AttVal(ff,rhoatt,&ival,&rval,&pval);
	fprintf(fp," (%-2.0lf)",rval);
      }
      fprintf(fp,"\n");


      if (maxrho > 0.0)
	MEnt_Mark(mf,mkid2);
      List_Delete(ffaces);

    } /* while (mf = MESH_Next_Face.... */

    idx = 0;
    while ((mf = MESH_Next_Face(mymesh,&idx))) {
      if (MEnt_IsMarked(mf,mkid2)) {
	MEnt_Set_AttVal(mf,rhoatt,0,10,NULL);
	MEnt_Unmark(mf,mkid2);
      }
    }


    sprintf(gmvfilename,"%s.gmv.%04d.%04d",basename,rank,i+1);

    fprintf(stderr,"Exported %s on proc %-d...",gmvfilename,rank);
    MESH_ExportToGMV(mymesh,gmvfilename,0,NULL,NULL);
    fprintf(stderr,"done\n");

    fprintf(stderr,"Updating attributes on proc %-d...",rank);
    MSTK_UpdateAttr(mymesh, rank, num, MPI_COMM_WORLD);
    fprintf(stderr,"done\n");

    fprintf(fp,"\n\n");

  } /* for (i = 0; i < nsteps; i++) */


  MPI_Finalize();
  return 0;
}
