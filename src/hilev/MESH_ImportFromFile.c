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
  else if (strncmp(format,"exo",3) == 0) {
#ifdef ENABLE_ExodusII
    return MESH_ImportFromExodusII(mesh,filename);
#else
    MSTK_Report("MESH_ImportFromFile","Exodus II file support not built in",MSTK_ERROR);
#endif
  } 
  else {
    MSTK_Report("MESH_ImportFromFile","Unsupported import format",MSTK_ERROR);
    return 0;
  }

}

#ifdef __cplusplus
}
#endif
