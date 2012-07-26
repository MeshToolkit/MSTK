#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "MSTK.h"

int MESH_ExportToDXBin(Mesh_ptr, const char *fname);

int MESH_ExportToDX(Mesh_ptr mesh, const char *fname, int binary) {

  return MESH_ExportToDXBin(mesh, fname);

}

