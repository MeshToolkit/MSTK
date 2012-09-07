#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <MSTK.h>

#ifdef __cplusplus
extern "C" {
#endif

int MESH_ExportToFile(Mesh_ptr mesh, const char *filename, const char *format,
		      const int natt, const char **attnames, const int *opts, MSTK_Comm comm) {

  if (strncmp(format,"mstk",4) == 0) {
    return MESH_WriteToFile(mesh,filename,MESH_RepType(mesh),comm);
  }
  else if (strncmp(format,"gmv",3) == 0) {
    return MESH_ExportToGMV(mesh,filename,natt,attnames,opts,comm);
  }
  else if (strncmp(format,"x3d",4) == 0) {
    return MESH_ExportToFLAGX3D(mesh,filename,natt,attnames,opts,comm);
  }
  else if (strncmp(format,"exodusii",8) == 0) {
#ifdef ENABLE_ExodusII
    return MESH_ExportToExodusII(mesh,filename,natt,attnames,opts,comm);
#else
    MSTK_Report("MESH_ExportFromFile","Exodus II file support not built in",MSTK_ERROR);
#endif
  } 
  else if (strncmp(format, "stl", 3) == 0) {
    return MESH_ExportToSTL(mesh,filename);
  }
  else if (strncmp(format, "dx", 3) == 0) {
    return MESH_ExportToDX(mesh,filename,1); /* binary export */
  }


  MSTK_Report("MESH_ExportFromFile","Unsupported export format",MSTK_ERROR);
  return 0;

}

#ifdef __cplusplus
}
#endif
