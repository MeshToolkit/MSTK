#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
    printf("serial partition with metis\n");
    printf("usage: %s -n 4 meshname.mstk\n",argv[0]);
    exit(-1);
  }

  
  strcpy(filename,argv[3]);
  strcpy(basename,argv[3]);
  basename[strlen(basename)-5] = '\0';

  mesh = MESH_New(UNKNOWN_REP);
  
  MESH_InitFromFile(mesh, filename);
  int ncells = MESH_Num_Regions(mesh);
  if (ncells > 0) {
    printf("3D mesh\n");
  }
  else if ((ncells = MESH_Num_Faces(mesh)) > 0)
    printf("2D mesh\n");
  else {
    fprintf(stderr,"Mesh is neither solid nor surface mesh. Exiting...\n");
    exit(-1);
  }
  num = atoi(argv[2]);


  Mesh_ptr *submeshes = (Mesh_ptr*) MSTK_malloc(num*sizeof(Mesh_ptr));

  int *part;
  MESH_Get_Partition(mesh, num, 0, 0, NULL, &part);
  MSTK_Mesh_Partition(mesh, num, part, 0, 0, submeshes);
  
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
  sprintf(gmvfilename,"%s_metis.gmv",basename,i,0);
  MESH_ExportToGMV(mesh,gmvfilename,0,NULL,NULL);
  MESH_Delete(mesh);
  for( i = 0; i < num; i++) {
    sprintf(gmvfilename,"%s_metis.%04d.%04d.gmv",basename,i,0);
    MESH_ExportToGMV(submeshes[i],gmvfilename,0,NULL,NULL);
    MESH_Delete(submeshes[i]);
  }	
  /*  free(part);*/
  MSTK_free(submeshes);
  return 0;
}
