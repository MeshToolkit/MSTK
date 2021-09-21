/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif


int MESH_ExportToOFF(Mesh_ptr mesh, const char *filename) {
  FILE		        *fp;
  char                  mesg[256];


  if (MESH_Num_Regions(mesh) || !MESH_Num_Faces(mesh)) {
    MSTK_Report("MESH_ExportToOFF", "Only surface meshes can be exported to OFF format", MSTK_ERROR);
    return 0;
  }

  if (!(fp = fopen(filename,"w"))) {
    sprintf(mesg,"Couldn't open output file %s\n",filename);
    MSTK_Report("MESH_ExportToSTL",mesg,MSTK_ERROR);
    return 0;
  }

  /* Write header */

  fprintf(fp,"OFF\n");  // Apparently optional; turn off if it causes problems
  fprintf(fp,"%d %d %d\n",MESH_Num_Vertices(mesh), MESH_Num_Faces(mesh), MESH_Num_Edges(mesh));


  /* Write vertex info */
  
  int idx = 0;
  MVertex_ptr mv;
  while ((mv = MESH_Next_Vertex(mesh, &idx))) {
    double xyz[3];
    MV_Coords(mv, xyz);
    fprintf(fp,"%20.12lf %20.12lf %20.12lf\n",xyz[0],xyz[1],xyz[2]);
  }


  /* Write face vertex info */
  MFace_ptr mf;
  idx = 0;
  while ((mf = MESH_Next_Face(mesh, &idx))) {
    List_ptr fverts = MF_Vertices(mf,1,0);
    fprintf(fp,"%d ",List_Num_Entries(fverts));
    for (int i = 0; i < List_Num_Entries(fverts); i++)
      fprintf(fp," %d",MV_ID(List_Entry(fverts,i)));
    fprintf(fp,"\n");
    List_Delete(fverts);
  }

  return 1;
}

#ifdef __cplusplus
}
#endif

  
  
