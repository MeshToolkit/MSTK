#include <UnitTest++.h>

#include "../../include/MSTK.h"

// Test if we can remove a face correctly from a region (check list of
// remaining faces and their directions)

TEST(MF_Rem_Edge_F1F3)
{
  int ok, i, j, nfe, nv=5;
  Mesh_ptr mesh;
  MFace_ptr mf;
  MVertex_ptr verts[5];
  List_ptr felist;
  double xyz[5][3] = {{0.0,0.0,0.0},
                      {1.0,0.0,0.0},
                      {1.5,0.5,0.0},
                      {1.0,1.0,0.0},
                      {0.0,1.0,0.0}};
  int evidx[5][2] = {{0,1},{2,1},{2,3},{4,3},{4,0}};
  int fedirs[5] = {1,0,1,0,1};
  MEdge_ptr fedges[5];


  MSTK_Init();

  mesh = MESH_New(F1);

  
  for (i = 0; i < nv; i++) {
    verts[i] = MV_New(mesh);
    MV_Set_Coords(verts[i],xyz[i]);
  }

  for (i = 0; i < nv; i++) {
    fedges[i] = ME_New(mesh);
    ME_Set_Vertex(fedges[i],0,verts[evidx[i][0]]);
    ME_Set_Vertex(fedges[i],1,verts[evidx[i][1]]);
  }

  mf = MF_New(mesh);
  MF_Set_Edges(mf,nv,fedges,fedirs);

  // Remove edges from face one by one and see if 
  // the leftover edges have the right IDs and dirs

  for (i = 0; i < nv; i++) {
    MF_Rem_Edge(mf,fedges[i]);
    
    felist = MF_Edges(mf,1,0);
    nfe = List_Num_Entries(felist);

    CHECK_EQUAL(nv-(i+1),nfe);  // Right number of edges?

    for (j = 0; j < nfe; j++) {
      CHECK_EQUAL(fedges[i+1+j],List_Entry(felist,j));
      CHECK_EQUAL(fedirs[i+1+j],MF_EdgeDir_i(mf,j));
    }                  

    List_Delete(felist);
  }

}
