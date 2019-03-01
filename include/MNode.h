/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#ifndef _H_MNode
#define _H_MNode

#ifdef __cplusplus
extern "C" {
#endif


  typedef struct MNode {
    int id;
    int floating;
    double coord[3];
    Mesh_ptr mesh;
    MTopoEntity_ptr mentity;
  } MNode, *MNode_ptr;
    
  int MN_ID(MNode_ptr mnode);
  int MN_Is_Floating(MNode_ptr mnode);
  void MN_SetCoords(MNode_ptr mnode);
  void MN_Coords(MNode_ptr mnode, double *xyz);
  int MN_MTopoDim(MNode_ptr mnode);
  MTopoEntity_ptr MN_MTopoEntity(MNode_ptr mnode);

#ifdef __cplusplus
}
#endif
