/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "MSTK.h"

int MESH_ExportToDXBin(Mesh_ptr, const char *fname);

int MESH_ExportToDX(Mesh_ptr mesh, const char *fname, int binary) {

  return MESH_ExportToDXBin(mesh, fname);

}

