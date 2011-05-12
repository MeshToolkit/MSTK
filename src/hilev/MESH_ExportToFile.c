#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <MSTK.h>

#ifdef __cplusplus
extern "C" {
#endif

int MESH_ExportToFile(Mesh_ptr mesh, const char *filename, const char *format,
		      const int natt, const char **attnames, int *opts) {

  if (strncmp(format,"gmv",3) == 0) {
    return MESH_ExportToGMV(mesh,filename,natt,attnames,opts);
  }
  else if (strncmp(format,"flag",4) == 0) {
    return MESH_ExportToFLAGX3D(mesh,filename,natt,attnames,opts);
  }
  else if (strncmp(format,"exodusii",8) == 0) {
#ifdef ENABLE_ExodusII
    return MESH_ExportToExodusII(mesh,filename,natt,attnames,opts);
#else
    MSTK_Report("MESH_ExportFromFile","Exodus II file support not built in",ERROR);
#endif
  } 
  else {
    MSTK_Report("MESH_ExportFromFile","Unsupported export format",ERROR);
    return 0;
  }

}

#ifdef __cplusplus
}
#endif
