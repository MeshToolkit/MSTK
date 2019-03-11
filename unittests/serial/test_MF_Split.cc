/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <UnitTest++.h>

#include "MSTK.h"

// Test if we can split a face correctly in a general 3D mesh

TEST(MF_Split) 
{
  int idx, i, j, n, nfv, ok;
  Mesh_ptr mesh;
  MVertex_ptr vsplit;
  MFace_ptr mf;
  MRegion_ptr fr;
  List_ptr fregions, vfaces;
  double fxyz[MAXPV2][3], cxyz[3];

  MSTK_Init();

  mesh = MESH_New(UNKNOWN_REP);
  ok = MESH_InitFromFile(mesh,"serial/reghex3D.mstk",NULL);
  CHECK_EQUAL(ok,1);

  CHECK(MESH_Num_Vertices(mesh) > 0);
  
  idx = 0; n = 0;
  while ((mf = MESH_Next_Face(mesh,&idx))) {
    if (MF_GEntDim(mf) == 3) {

      /* Store regions of face before it gets deleted */
      fregions = MF_Regions(mf);
      
      MF_Coords(mf,&nfv,fxyz);
      cxyz[0] = cxyz[1] = cxyz[2] = 0.0;
      for (i = 0; i < nfv; i++) {
        for (j = 0; j < 3; j++)
          cxyz[j] += fxyz[i][j];
      }
      for (j = 0; j < 3; j++)
        cxyz[j] /= nfv;

      vsplit = MF_Split(mf,cxyz);
      CHECK(vsplit);

      vfaces = MV_Faces(vsplit);
      CHECK_EQUAL(nfv,List_Num_Entries(vfaces));

      for (i = 0; i < nfv; i++) {
        MFace_ptr vf = List_Entry(vfaces,i);
        idx = 0;
        while ((fr = List_Next_Entry(fregions,&idx)))
          CHECK(MR_UsesEntity(fr,vf,MFACE));
      }
      List_Delete(fregions);

      n++;
      if (n > 2)
        break;
    }
  }

  MESH_Delete(mesh);
}
