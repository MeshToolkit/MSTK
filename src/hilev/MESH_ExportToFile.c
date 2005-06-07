#include <stdio.h>
#include <stdlib.h>
#include <MSTK.h>

#ifdef __cplusplus
extern "C" {
#endif

int MESH_ExportToFile(Mesh_ptr mesh, const char *filename, const char *format,
		      const int natt, const char **attnames, int *opts) {

  if (strncmp(format,"gmv",3) == 0) {
    return MESH_ExportToGMV(mesh,filename,natt,attnames,opts);
  }
  /*  else (strncmp(format,"avs",3) == 0) {
    return MESH_ImportFromUCD(mesh,filename,natt,attnames);
    } */
  /*  else (strncmp(format,"dx",2) == 0) {
    return MESH_ExportToDX(mesh,filename,natt,attnames);
    } */
  else {
    MSTK_Report("MESH_ImportFromFile","Unsupported import format",ERROR);
    return 0;
  }

}

#ifdef __cplusplus
}
#endif
