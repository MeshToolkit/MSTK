#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <MSTK.h>

#ifdef __cplusplus
extern "C" {
#endif
  /* Function to export MSTK mesh to MSTK and other formats */

  /* if natt = 0, all attributes are written out 
     if natt = -1, no attributes are written out 
     if natt > 0, only attributes specified in attnames are written out 

     opts is an array of flags that controls how mesh is exported and is
     specific to each format. See the specific MESH_ExportTo"Format" file
 
     If you want default options, call the functions as

     MESH_ExportToFile(mymesh,myfile,myformat,0,NULL,NULL,NULL);
  */

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
  else if (strncmp(format,"exo",3) == 0) {
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
