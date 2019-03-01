/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <UnitTest++.h>

#include "MSTK.h"

// Test if we can split an edge along with connected regions in a tet mesh

TEST(ME_Split_SimplexMesh) 
{
  int idx, idx2, i, j, nv, ner, ok, nsplit;
  Mesh_ptr mesh;
  MVertex_ptr rv, vsplit, ev[2], cmnrv, opprv[2];
  MEdge_ptr me;
  MRegion_ptr vr, er0, er1;
  List_ptr eregions, vregions, rverts0, rverts1;
  double exyz[2][3], cxyz[3];

  MSTK_Init();

  mesh = MESH_New(UNKNOWN_REP);
  ok = MESH_InitFromFile(mesh,"serial/twotets.mstk",NULL);
  CHECK_EQUAL(ok,1);

  CHECK(MESH_Num_Vertices(mesh) > 0);
  
  idx = 0; nsplit = 0;
  while ((me = MESH_Next_Edge(mesh,&idx))) {

    /* Store regions of face before it gets deleted */
    eregions = ME_Regions(me);
    ner = List_Num_Entries(eregions);
    if (ner < 2) continue;
      
    for (i = 0; i < 2; i++) {
      ev[i] = ME_Vertex(me,i);
      MV_Coords(ev[i],exyz[i]);
    }
        
    cxyz[0] = cxyz[1] = cxyz[2] = 0.0;
    for (i = 0; i < 2; i++) {
      for (j = 0; j < 3; j++)
        cxyz[j] += exyz[i][j];
    }
    for (j = 0; j < 3; j++)
      cxyz[j] /= 2;
    
    cmnrv = 0;
    er0 = List_Entry(eregions,0);
    er1 = List_Entry(eregions,1);
    rverts0 = MR_Vertices(er0);
    rverts1 = MR_Vertices(er1);
    idx2 = 0;
    while ((rv = List_Next_Entry(rverts0,&idx2))) {
      if (rv != ev[0] && rv != ev[1]) {
        if (List_Contains(rverts1,rv))
          cmnrv = rv;
        else
          opprv[0] = rv;
      }
    }
    List_Delete(rverts0);

    idx2 = 0;
    while ((rv = List_Next_Entry(rverts1,&idx2))) {
      if (rv != ev[0] && rv != ev[1] && rv != cmnrv)
        opprv[1] = rv;
    }
    List_Delete(rverts1);

    List_Delete(eregions);

    CHECK(cmnrv);
    CHECK(opprv[0]);
    CHECK(opprv[1]);

    vsplit = ME_Split_SimplexMesh(me,cxyz);
    CHECK(vsplit);
    nsplit++;

    vregions = MV_Regions(vsplit);
    CHECK_EQUAL(2*ner,List_Num_Entries(vregions));

    idx2 = 0;
    while ((vr = List_Next_Entry(vregions,&idx2))) {
      rverts0 = MR_Vertices(vr);
      CHECK_EQUAL(4,List_Num_Entries(rverts0));

      for (i = 0, nv = 0; i < 4; i++) {
        rv = List_Entry(rverts0,i);
        if (rv == vsplit)
          nv++;
        else if (rv == ev[0] || rv == ev[1])
          nv++;
        else if (rv == cmnrv)
          nv++;
        else if (rv == opprv[0] || rv == opprv[1])
          nv++;
      }
      List_Delete(rverts0);

      /* Did we find all the four vertices we expected? */
      CHECK_EQUAL(4,nv);

    }
    List_Delete(vregions);

    if (nsplit == 1)
      break;
  }

}
