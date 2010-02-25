#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <MSTK.h>

#ifdef __cplusplus
extern "C" {
#endif

int MESH_ImportFromFile(Mesh_ptr mesh, const char *filename, const char *format) {

  if (strncmp(format,"gmv",3) == 0) {
    return MESH_ImportFromGMV(mesh,filename);
  }
  /*  else (strncmp(format,"avs",3) == 0) {
    return MESH_ImportFromUCD(mesh,filename);
    } */
  else {
    MSTK_Report("MESH_ImportFromFile","Unsupported import format",ERROR);
    return 0;
  }

}

#ifdef __cplusplus
}
#endif
