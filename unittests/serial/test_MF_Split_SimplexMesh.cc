/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <UnitTest++.h>

#include "MSTK.h"

// Test if we can split a face along with connected regions in a tet mesh

TEST(MF_Split_SimplexMesh) 
{
  int idx, idx2, i, j, nv, nfv, ok, nsplit;
  Mesh_ptr mesh;
  MVertex_ptr rv, vsplit, fv[3], opprv[2];
  MFace_ptr mf;
  MRegion_ptr vr, fr;
  List_ptr fverts, fregions, vregions, rverts;
  double fxyz[3][3], cxyz[3];

  MSTK_Init();

  mesh = MESH_New(UNKNOWN_REP);
  ok = MESH_InitFromFile(mesh,"serial/twotets.mstk",NULL);
  CHECK_EQUAL(ok,1);

  CHECK(MESH_Num_Vertices(mesh) > 0);
  
  idx = 0; nsplit = 0;
  while ((mf = MESH_Next_Face(mesh,&idx))) {
    if (MF_GEntDim(mf) == 3) {

      /* Store regions of face before it gets deleted */
      fregions = MF_Regions(mf);
      
      fverts = MF_Vertices(mf,1,0);
      nfv = List_Num_Entries(fverts);
      CHECK_EQUAL(3,nfv);

      for (i = 0; i < 3; i++) {
        fv[i] = List_Entry(fverts,i);
        MV_Coords(fv[i],fxyz[i]);
      }
        
      cxyz[0] = cxyz[1] = cxyz[2] = 0.0;
      for (i = 0; i < nfv; i++) {
        for (j = 0; j < 3; j++)
          cxyz[j] += fxyz[i][j];
      }
      for (j = 0; j < 3; j++)
        cxyz[j] /= nfv;

      for (i = 0; i < 2; i++) {
        opprv[i] = 0;
        fr = List_Entry(fregions,i);
        rverts = MR_Vertices(fr);
        idx2 = 0;
        while ((rv = List_Next_Entry(rverts,&idx2))) {
          if (!List_Contains(fverts,rv)) {
            opprv[i] = rv;
            break;
          }
        }
        List_Delete(rverts);
      }

      vsplit = MF_Split_SimplexMesh(mf,cxyz);
      CHECK(vsplit);
      nsplit++;

      vregions = MV_Regions(vsplit);
      CHECK_EQUAL(6,List_Num_Entries(vregions));

      idx2 = 0;
      while ((vr = List_Next_Entry(vregions,&idx2))) {
        rverts = MR_Vertices(vr);
        CHECK_EQUAL(4,List_Num_Entries(rverts));

        for (i = 0, nv = 0; i < 4; i++) {
          rv = List_Entry(rverts,i);
          if (rv == vsplit)
            nv++;
          else if (List_Contains(fverts,rv))
            nv++;
          else {
            /* fourth vertex of child region - should be opprv[0] or
               opprv[1] */
            CHECK((rv == opprv[0]) || (rv == opprv[1]));
            if ((rv == opprv[0]) || (rv == opprv[1]))
              nv++;
          }
        }
        List_Delete(rverts);

        /* Did we find all the four vertices we expected? */
        CHECK_EQUAL(4,nv);
      }
      List_Delete(vregions);

      List_Delete(fverts);
      List_Delete(fregions);
    }
    if (nsplit == 5) 
      break;
  }

  MESH_Delete(mesh);
}
