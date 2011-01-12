#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "MSTK.h"

int MESH_ExportToDX(Mesh_ptr mesh, const char *fname, int binary) {

  MESH_ExportToDXBin(mesh, fname);

}

