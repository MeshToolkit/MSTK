/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <UnitTest++.h>
#include <cstdio>
#include <vector>

#include "MSTK.h"


TEST(Mesh_Attributes) {

  MSTK_Init();

  /* Read a 2x2x2 mesh of hexes with full classification in MSTK format */

  Mesh_ptr mesh = MESH_New(UNKNOWN_REP);
  int ok = MESH_InitFromFile(mesh,"serial/reghex3D.mstk",NULL);
  CHECK_EQUAL(ok, 1);

  // Put different types of attributes on different entities
  char attname[256];
  MAttrib_ptr intatt = MAttrib_New(mesh, "intatt", INT, MALLTYPE);
  MAttrib_Get_Name(intatt, attname);
  CHECK_EQUAL("intatt", attname);
  CHECK_EQUAL(INT, MAttrib_Get_Type(intatt));
  CHECK_EQUAL(MALLTYPE, MAttrib_Get_EntDim(intatt));

  MAttrib_ptr dblatt = MAttrib_New(mesh, "dblatt", DOUBLE, MVERTEX);
  MAttrib_Get_Name(dblatt, attname);
  CHECK_EQUAL("dblatt", attname);
  CHECK_EQUAL(DOUBLE, MAttrib_Get_Type(dblatt));
  CHECK_EQUAL(MVERTEX, MAttrib_Get_EntDim(dblatt));
  
  MAttrib_ptr ptratt = MAttrib_New(mesh, "ptratt", POINTER, MREGION);
  MAttrib_Get_Name(ptratt, attname);
  CHECK_EQUAL("ptratt", attname);
  CHECK_EQUAL(POINTER, MAttrib_Get_Type(ptratt));
  CHECK_EQUAL(MREGION, MAttrib_Get_EntDim(ptratt));

  MAttrib_ptr vecatt = MAttrib_New(mesh, "vecatt", VECTOR, MEDGE, 3);
  MAttrib_Get_Name(vecatt, attname);
  CHECK_EQUAL("vecatt", attname);
  CHECK_EQUAL(VECTOR, MAttrib_Get_Type(vecatt));
  CHECK_EQUAL(MEDGE, MAttrib_Get_EntDim(vecatt));

  MAttrib_ptr tnsratt = MAttrib_New(mesh, "tnsratt", TENSOR, MFACE, 6);
  MAttrib_Get_Name(tnsratt, attname);
  CHECK_EQUAL("tnsratt", attname);
  CHECK_EQUAL(TENSOR, MAttrib_Get_Type(tnsratt));
  CHECK_EQUAL(MFACE, MAttrib_Get_EntDim(tnsratt));


  // Set 

  int idx = 0;
  MVertex_ptr mv;
  while ((mv = MESH_Next_Vertex(mesh, &idx)))
    if (MV_ID(mv)%2) MEnt_Set_IntAttVal(mv, intatt, 23);

  idx = 0;
  MEdge_ptr me;
  while ((me = MESH_Next_Edge(mesh, &idx)))
    MEnt_Set_IntAttVal(me, intatt, 929);

  double xyz[3];
  idx = 0;
  while ((mv = MESH_Next_Vertex(mesh, &idx))) {
    MV_Coords(mv, xyz);
    MEnt_Set_DblAttVal(mv, dblatt, (xyz[0]+xyz[1]+xyz[2]));
  }

  MRegion_ptr mr;
  idx = 0;
  while ((mr = MESH_Next_Region(mesh, &idx))) {
    List_ptr rfaces = MR_Faces(mr);
    MEnt_Set_PtrAttVal(mr, ptratt, List_Entry(rfaces, 1));
  }

  idx = 0;
  while ((me = MESH_Next_Edge(mesh, &idx))) {
    double vxyz0[3];
    double vxyz1[3];
    double vec[3];
    MV_Coords(ME_Vertex(me, 0), vxyz0);
    MV_Coords(ME_Vertex(me, 1), vxyz1);
    int i;
    for (i = 0; i < 3; i++)
      vec[i] = vxyz1[i]-vxyz0[i];
    MEnt_Set_VecAttVal(me, vecatt, vec);
  }

  idx = 0;
  MFace_ptr mf;
  while ((mf = MESH_Next_Face(mesh, &idx))) {
    double fxyz[MAXPV2][3];
    int nfv;

    MF_Coords(mf, &nfv, fxyz);

    double center[3] = {0.0, 0.0, 0.0};
    int i, j;
    for (i = 0; i < nfv; i++) {
      for (j = 0; j < 3; j++) {
        center[j] += fxyz[i][j];
      }
    }
    for (j = 0; j < 3; j++)
      center[j] /= nfv;

    double tencomp[6];
    for (i = 0; i < 3; i++)
      tencomp[i] = center[3-i-1];
    for (i = 0; i < 3; i++)
      tencomp[3+i] = center[i];

    MEnt_Set_TnsrAttVal(mf, tnsratt, tencomp);
  }

  // Retrieve and check

  idx = 0;
  while ((mv = MESH_Next_Vertex(mesh, &idx)))
    if (MV_ID(mv)%2) CHECK_EQUAL(23, MEnt_Get_IntAttVal(mv, intatt));

  idx = 0;
  while ((me = MESH_Next_Edge(mesh, &idx)))
    CHECK_EQUAL(929, MEnt_Get_IntAttVal(me, intatt));

  
  idx = 0;
  while ((mv = MESH_Next_Vertex(mesh, &idx))) {
    MV_Coords(mv, xyz);
    double val = MEnt_Get_DblAttVal(mv, dblatt);
    CHECK_EQUAL((xyz[0]+xyz[1]+xyz[2]), val);
  }

  while ((mr = MESH_Next_Region(mesh, &idx))) {
    List_ptr rfaces = MR_Faces(mr);
    MFace_ptr rf = MEnt_Get_PtrAttVal(mr, ptratt);
    CHECK_EQUAL(List_Entry(rfaces, 1), rf);
  }

  idx = 0;
  while ((me = MESH_Next_Edge(mesh, &idx))) {
    double *vec, evec[3];
    ME_Vec(me, evec);
    vec = MEnt_Get_VecAttVal(me, vecatt);
    CHECK_ARRAY_EQUAL(evec, vec, 3);
  }

  idx = 0;
  while ((mf = MESH_Next_Face(mesh, &idx))) {
    double fxyz[MAXPV2][3];
    int nfv;

    MF_Coords(mf, &nfv, fxyz);

    double center[3] = {0.0, 0.0, 0.0};
    int i, j;
    for (i = 0; i < nfv; i++) {
      for (j = 0; j < 3; j++) {
        center[j] += fxyz[i][j];
      }
    }
    for (j = 0; j < 3; j++)
      center[j] /= nfv;

    double *tencomp = MEnt_Get_TnsrAttVal(mf, tnsratt);    
    for (i = 0; i < 3; i++)
      CHECK_EQUAL(center[3-i-1], tencomp[i]);
    for (i = 0; i < 3; i++)
      CHECK_EQUAL(center[i], tencomp[3+i]);
  }

  MESH_Delete(mesh);
}


