#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif


int MESH_ExportToSTL(Mesh_ptr mesh, const char *filename) {
  List_ptr		fverts;
  MFace_ptr	        face;
  char                  mesg[256];
  int                   i, j, nfv, idx;
  double		fxyz[MAXPV2][3], normal[3], len, vec1[3], vec2[3];
  FILE		        *fp;

  
  if (!(fp = fopen(filename,"w"))) {
    sprintf(mesg,"Couldn't open output file %s\n",filename);
    MSTK_Report("MESH_ExportToSTL",mesg,ERROR);
    return 0;
  }

  if (MESH_Num_Regions(mesh)) {
    MSTK_Report("MESH_ExportToSTL","Cannot export solid meshes\n",ERROR);
    return 0;
  }

  idx = 0;
  while ((face = MESH_Next_Face(mesh,&idx))) {
    if (MF_Num_Edges(face) != 3) {
      MSTK_Report("MESH_ExportToSTL","Can only export triangular meshes to STL\n",ERROR);
      MSTK_Report("MESH_ExportToSTL","If you need polygonal meshes to be exported\n",MESG);
      MSTK_Report("MESH_ExportToSTL","as a triangulated surface mesh to STL, \n",MESG);
      MSTK_Report("MESH_ExportToSTL","please send a mail to rao@lanl.gov\n",MESG);
      return 0;
    }
  }
  

  /* Opening line for STL file */
  fprintf(fp,"solid\n");

  idx = 0;
  while ((face = MESH_Next_Face(mesh,&idx))) {

    fverts = MF_Vertices(face,1,0);
    nfv = List_Num_Entries(fverts);

    for (i = 0; i < nfv; i++)
      MV_Coords(List_Entry(fverts,i),fxyz[i]);
      
    {
      vec1[0] = fxyz[1][0]-fxyz[0][0];
      vec1[1] = fxyz[1][1]-fxyz[0][1];
      vec1[2] = fxyz[1][2]-fxyz[0][2];
      
      vec2[0] = fxyz[2][0]-fxyz[0][0];
      vec2[1] = fxyz[2][1]-fxyz[0][1];
      vec2[2] = fxyz[2][2]-fxyz[0][2];
      
      normal[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
      normal[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
      normal[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
      
      len = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
      
      normal[0] /= len; normal[1] /= len; normal[2] /= len;
    }

    fprintf(fp,"facet normal %12.8lf %12.8lf %12.8lf\n",normal[0],normal[1],normal[2]);
    fprintf(fp,"  outer loop\n");
    for (j = 0; j < 3; j++) 
      fprintf(fp,"    vertex %20.10lf %20.10lf %20.10lf\n",fxyz[j][0],fxyz[j][1],fxyz[j][2]);
    fprintf(fp,"endloop");
    fprintf(fp,"endfacet");
  }

  fprintf(fp,"endsolid");

  fclose(fp);

  return 1;
}

#ifdef __cplusplus
  }
#endif

