/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>

#include "MSTK.h"

/* Get the partitioning of mesh by reading a coloring file -
   Doesn't actually do anything to the mesh */

#ifdef __cplusplus
extern "C" {
#endif


int MESH_PartitionWithColoringFile(Mesh_ptr mesh, int nparts, int **part) {
  int ncells;
  FILE* fid;
  size_t count;

  ncells = MESH_Num_Regions(mesh);
  *part = (int *) malloc(ncells*sizeof(int));

  fid = fopen("coloring.bin", "rb");
  if (fid == NULL) {
    fprintf(stderr,"Nonexistent coloring file \"coloring.bin\"\n");
    exit(-1);
  }

  count = fread(*part, sizeof(int), (size_t) ncells, fid);
  fprintf(stderr,"read %zu of %i region colors from coloring.bin\n", count, ncells);
  if (count != (size_t) ncells) {
    fprintf(stderr,"Error reading coloring file \"coloring.bin\"\n");
    exit(-1);
  }
  fclose(fid);
  return 1;
}
