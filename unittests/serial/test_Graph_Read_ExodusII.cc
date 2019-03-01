#include <UnitTest++.h>

#include <cstdlib>

#include "MSTK.h"
#include "MSTK_private.h"

// Test if we can read the Exodus II file to construct only the adjacency graph
// Do this by comparing against the adjacency graph we get from a full mesh

// Assume that if we write out a mesh in Exodus II format and read it
// back in the element numbering will not change

TEST(READ_GRAPH_2D_SIMPLE) {
  
  MSTK_Init();

  Mesh_ptr mesh = MESH_New(F1);

  mesh = MESH_Gen_Structured(0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 2, 2, 0);

  MESH_ExportToFile(mesh, "struct2-temp.exo", "exo", 0, NULL, NULL, NULL);

  // Element graph of mesh in memory
  int nelems1;
  int *adjoff1, *adjelems1;
  MESH_GetElementGraph(mesh, &nelems1, &adjoff1, &adjelems1);

  // Sort the adjacent element list for each element in increasing order
  for (int i = 0; i < nelems1; i++) {
    int off = adjoff1[i];
    int nadj = adjoff1[i+1] - adjoff1[i];
    qsort(adjelems1 + off, nadj, sizeof(int), compareINT);
  }

  int ndim;
  int nelems2;
  int *adjoff2, *adjelems2;
  ExodusII_GetElementGraph("serial/mesh5x5-skewed.exo", &ndim, &nelems2, &adjoff2,
			   &adjelems2, false, NULL);
  ExodusII_GetElementGraph("struct2-temp.exo", &ndim, &nelems2, &adjoff2,
			   &adjelems2, false, NULL);
  for (int i = 0; i < nelems2; i++) {
    int off = adjoff2[i];
    int nadj = adjoff2[i+1] - adjoff2[i];
    qsort(adjelems2 + off, nadj, sizeof(int), compareINT);
  }

  CHECK_EQUAL(nelems1, nelems2);
  
  for (int i = 0; i < nelems1; i++)
    CHECK_EQUAL(adjoff1[i], adjoff2[i]);

  int nentries = adjoff1[nelems1];
  for (int i = 0; i < nentries; i++)
    CHECK_EQUAL(adjelems1[i], adjelems2[i]);
}

TEST(READ_GRAPH_3D_SIMPLE) {

  MSTK_Init();

  Mesh_ptr mesh = MESH_New(F1);

  mesh = MESH_Gen_Structured(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3);

  MESH_ExportToFile(mesh, "struct3-temp.exo", "exo", 0, NULL, NULL, NULL);

  // Element graph of mesh in memory
  int nelems1;
  int *adjoff1, *adjelems1;
  MESH_GetElementGraph(mesh, &nelems1, &adjoff1, &adjelems1);

  // Sort the adjacent element list for each element in increasing order
  for (int i = 0; i < nelems1; i++) {
    int off = adjoff1[i];
    int nadj = adjoff1[i+1] - adjoff1[i];
    qsort(adjelems1 + off, nadj, sizeof(int), compareINT);
  }

  int ndim;
  int nelems2;
  int *adjoff2, *adjelems2;
  double (*elemcens)[3];
  ExodusII_GetElementGraph("struct3-temp.exo", &ndim, &nelems2, &adjoff2,
			   &adjelems2, true, &elemcens);
  for (int i = 0; i < nelems2; i++) {
    int off = adjoff2[i];
    int nadj = adjoff2[i+1] - adjoff2[i];
    qsort(adjelems2 + off, nadj, sizeof(int), compareINT);
  }

  CHECK_EQUAL(nelems1, nelems2);
  
  for (int i = 0; i < nelems1; i++)
    CHECK_EQUAL(adjoff1[i], adjoff2[i]);

  int nentries = adjoff1[nelems1];
  for (int i = 0; i < nentries; i++)
    CHECK_EQUAL(adjelems1[i], adjelems2[i]);
}

    
// Assume that if we write out a mesh in Exodus II format and read it
// back in the element numbering will not change

TEST(READ_GRAPH_2D_LARGE) {

  MSTK_Init();

  Mesh_ptr mesh = MESH_New(F1);

  MESH_ImportFromFile(mesh, "struct2-temp.exo", "exodusii", NULL, NULL);

  // Element graph of mesh in memory
  int nelems1;
  int *adjoff1, *adjelems1;
  MESH_GetElementGraph(mesh, &nelems1, &adjoff1, &adjelems1);

  // Sort the adjacent element list for each element in increasing order
  for (int i = 0; i < nelems1; i++) {
    int off = adjoff1[i];
    int nadj = adjoff1[i+1] - adjoff1[i];
    qsort(adjelems1 + off, nadj, sizeof(int), compareINT);
  }

  int ndim;
  int nelems2;
  int *adjoff2, *adjelems2;
  ExodusII_GetElementGraph("struct2-temp.exo", &ndim, &nelems2, &adjoff2,
			   &adjelems2, false, NULL);
  for (int i = 0; i < nelems2; i++) {
    int off = adjoff2[i];
    int nadj = adjoff2[i+1] - adjoff2[i];
    qsort(adjelems2 + off, nadj, sizeof(int), compareINT);
  }

  CHECK_EQUAL(nelems1, nelems2);
  
  for (int i = 0; i < nelems1; i++)
    CHECK_EQUAL(adjoff1[i], adjoff2[i]);

  int nentries = adjoff1[nelems1];
  for (int i = 0; i < nentries; i++)
    CHECK_EQUAL(adjelems1[i], adjelems2[i]);
}

TEST(READ_GRAPH_3D_LARGE) {

  MSTK_Init();

  Mesh_ptr mesh = MESH_New(F1);

  MESH_ImportFromFile(mesh, "struct3-temp.exo", "exodusii", NULL, NULL);

  // Element graph of mesh in memory
  int nelems1;
  int *adjoff1, *adjelems1;
  MESH_GetElementGraph(mesh, &nelems1, &adjoff1, &adjelems1);

  // Sort the adjacent element list for each element in increasing order
  for (int i = 0; i < nelems1; i++) {
    int off = adjoff1[i];
    int nadj = adjoff1[i+1] - adjoff1[i];
    qsort(adjelems1 + off, nadj, sizeof(int), compareINT);
  }

  int ndim;
  int nelems2;
  int *adjoff2, *adjelems2;
  ExodusII_GetElementGraph("struct3-temp.exo", &ndim, &nelems2, &adjoff2,
			   &adjelems2, false, NULL);
  for (int i = 0; i < nelems2; i++) {
    int off = adjoff2[i];
    int nadj = adjoff2[i+1] - adjoff2[i];
    qsort(adjelems2 + off, nadj, sizeof(int), compareINT);
  }

  CHECK_EQUAL(nelems1, nelems2);
  
  for (int i = 0; i < nelems1; i++)
    CHECK_EQUAL(adjoff1[i], adjoff2[i]);

  int nentries = adjoff1[nelems1];
  for (int i = 0; i < nentries; i++)
    CHECK_EQUAL(adjelems1[i], adjelems2[i]);
}

