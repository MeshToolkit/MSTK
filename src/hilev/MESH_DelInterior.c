/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Routine to remove interior of a mesh and only leave a boundary mesh.
   For 3D meshes, only a surface mesh is left
   For 2D meshes, only the boundary edge mesh is left
   For 1D meshes, only the bounding vertices are left
   Does not do anything if the mesh has only vertices
*/

int MESH_DelInterior(Mesh_ptr mesh) {
  int i, dim=-1, nr, nf, ne, nv;
  MVertex_ptr vertex;
  MEdge_ptr edge;
  MFace_ptr face;
  MRegion_ptr region;

  nr = MESH_Num_Regions(mesh);
  nf = MESH_Num_Faces(mesh);
  ne = MESH_Num_Edges(mesh);
  nv = MESH_Num_Vertices(mesh);

  if (nr)
    dim = 3;
  else if (nf)
    dim = 2;
  else if (ne)
    dim = 1;
  else {
    /* Only vertices */
    return 0;
  }

  /* With the current storage scheme for mesh entities or more generally lists,
     it is more efficient to delete elements starting from the back */
  
  for (i = nr-1; i >= 0; i--) { 
    region = MESH_Region(mesh,i);    
    MR_Delete(region,0);
  }

  for (i = nf-1; i >= 0; i--) { 
    face = MESH_Face(mesh,i);
    if (MF_GEntDim(face) == dim)
      MF_Delete(face,0);
  }

  for (i = ne-1; i >= 0; i--) { 
    edge = MESH_Edge(mesh,i);
    if (ME_GEntDim(edge) == dim)
      ME_Delete(edge,0);
  }

  for (i = nv-1; i >= 0; i--) { 
    vertex = MESH_Vertex(mesh,i);
    if (MV_GEntDim(vertex) == dim)
      MV_Delete(vertex,0);
  }

  return 1;
}


#ifdef __cplusplus
}
#endif
