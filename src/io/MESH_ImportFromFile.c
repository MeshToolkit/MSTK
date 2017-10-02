#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <MSTK.h>

#ifdef __cplusplus
extern "C" {
#endif

  int MESH_ImportFromFile(Mesh_ptr mesh, const char *infilename, const char *informat, int *opts, MSTK_Comm comm) {

    char format[16];
    
    if (informat)
      strcpy(format,informat);
    else {
      char *ext;
      if (infilename[0] == '.') {
        ext = strchr(infilename,'.');
        ++ext;
        ext = strchr(ext,'.');
      }
      else
        ext = strchr(infilename,'.');
      if (ext != NULL) {
        ++ext;                 /* go past the . */
        if (strcmp(ext,"mstk") == 0) 
          strcpy(format,"mstk");
        else if (strcmp(ext,"gmv") == 0)
          strcpy(format,"gmv");
        else if (strcmp(ext,"exo") == 0)
          strcpy(format,"exodusii");
        else if (strcmp(ext,"par") == 0)
          strcpy(format,"nemesisi");
        else {
          fprintf(stderr,"Unknown file format\n");
          return -1;
        }
      }
      else {
        fprintf(stderr,"Unknown file format\n");
        return -1;
      }
    }


    if (strcmp(format,"mstk") == 0)
      return MESH_InitFromFile(mesh,infilename,comm);
    else {
      if (MESH_RepType(mesh) == UNKNOWN_REP)
        MESH_SetRepType(mesh,F1);
      
      if (strcmp(format,"gmv") == 0)
        return MESH_ImportFromGMV(mesh,infilename,comm);
      else if (strcmp(format,"exodusii") == 0) {
#ifdef ENABLE_ExodusII
        return MESH_ImportFromExodusII(mesh,infilename,opts,comm);
#else
        MSTK_Report("MESH_ImportFromFile","Exodus II file support not built in",MSTK_ERROR);
#endif
      } 
      else if (strcmp(format,"nemesisi") == 0) {
#ifdef ENABLE_ExodusII
        return MESH_ImportFromNemesisI(mesh,infilename,opts,comm);
#else
        MSTK_Report("MESH_ImportFromFile","Exodus II/Nemesis I file support not built in",MSTK_ERROR);
#endif
      } 
      else if (strcmp(format,"x3d") == 0) {
        int rank = 0;
        int numprocs = 1;
        return MESH_ImportFromFLAGX3D(mesh,infilename,comm);
      }
    }
    
    MSTK_Report("MESH_ImportFromFile","Unsupported import format",MSTK_ERROR);
    return 0;
    
  }

#ifdef __cplusplus
}
#endif
