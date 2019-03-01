#include <UnitTest++.h>
#include <iostream>
#include <vector>
#include "omp.h"

#include "MSTK.h"

// Test to make sure entity marking works correctly in multi-threaded runs
// Sadly, I am not able to make this fail with gcc 5.4.0 on Ubuntu 16.04 even
// without any mutex locks around the code that fetches markers. But the test
// does pass with multiple threads.


#ifdef MSTK_USE_MARKERS

TEST(Mesh_Markers) {
 
  MSTK_Init();

  /* Read a 2x2x2 mesh of hexes with full classification in MSTK format */

  Mesh_ptr mesh = MESH_New(F1);
  int ok = MESH_ImportFromFile(mesh,"serial/degenpoly3.exo",NULL,NULL,NULL);
  CHECK_EQUAL(ok, 1);

  /* Calculate get mesh vertices via regions and tag them once and only
   * once using entity marking */
  
  int nv = MESH_Num_Vertices(mesh);
  int *vcounts1 = new int[nv];
  int *vcountsN = new int[nv];
  for (int i = 0; i < nv; ++i) {
    vcounts1[i] = 0;
    vcountsN[i] = 0;
  }

  int maxthreads = 8;
  omp_set_num_threads(maxthreads);
  std::cerr << "Max threads is " << maxthreads << "\n";
  
  int mkid = MSTK_GetMarker();
  CHECK(mkid != -1);

  int nr = MESH_Num_Regions(mesh);
  
#pragma omp parallel default(none)              \
  shared(mesh, mkid, nr, vcounts1, vcountsN, maxthreads)
  {
    int nthreads = omp_get_num_threads();
    CHECK_EQUAL(maxthreads, nthreads);

#pragma omp for
    for (int ir = 0; ir < nr; ir++) {
      MRegion_ptr mr = MESH_Region(mesh, ir);

      // Even though we can use MR_Vertices, we will gather vertices
      // of the region explicitly here so as to get multiple markers
      // from different threads

      int mkid2 = MSTK_GetMarker();  // use to build vertices of regions
      int mkid3 = MSTK_GetMarker();  // use to build edges of regions
      List_ptr rvertices = List_New(10);
      List_ptr redges = List_New(10);
      
      List_ptr rfaces = MR_Faces(mr);
      int idx2 = 0;
      MFace_ptr rf;
      while ((rf = List_Next_Entry(rfaces, &idx2))) {
        List_ptr fedges = MF_Edges(rf, 1, 0);
        int idx3 = 0;
        MEdge_ptr fe;
        while ((fe = List_Next_Entry(fedges, &idx3))) {
          if (!MEnt_IsMarked(fe, mkid2)) {
            List_Add(redges, fe);
            MEnt_Mark(fe, mkid2);
          }

          for (int j = 0; j < 2; j++) {
            MVertex_ptr ev = ME_Vertex(fe, j);
            if (!MEnt_IsMarked(ev, mkid2)) {
              List_Add(rvertices, ev);
              MEnt_Mark(ev, mkid3);
            }          
          }
        }
        List_Delete(fedges);
      }
      List_Delete(rfaces);
        
      idx2 = 0;
      MVertex_ptr rv;
      while ((rv = List_Next_Entry(rvertices, &idx2))) {
        int vid = MV_ID(rv);
        vcountsN[vid-1]++;
        if (!MEnt_IsMarked(rv, mkid)) {
          vcounts1[vid-1]++;
          MEnt_Mark(rv, mkid);
        }
      }

      List_Unmark(redges, mkid2);
      MSTK_FreeMarker(mkid2);
      List_Delete(redges);

      List_Unmark(rvertices, mkid3);
      MSTK_FreeMarker(mkid3);
      List_Delete(rvertices);
    }
  }
    
  for (int i = 0; i < nv; i++) {
    CHECK_EQUAL(1, vcounts1[i]);
    CHECK(vcountsN[i] > 0);
  }
  
  int idx = 0;
  MVertex_ptr mv;
  while ((mv = MESH_Next_Vertex(mesh, &idx))) {
    MEnt_Unmark(mv, mkid);
  }
  MSTK_FreeMarker(mkid);

}

#endif
