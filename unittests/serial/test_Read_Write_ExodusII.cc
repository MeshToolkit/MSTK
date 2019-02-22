#include <UnitTest++.h>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "../../include/MSTK.h"
#include "../../include/MSTK_private.h"

// Test if we can correctly read a polygonal exodus mesh 
// Then test if we can correctly write it and read it back

TEST(Read_Write_ExodusII_Poly2) 
{
  int idx, ok, nfe;
  Mesh_ptr mesh, mesh2;
  MFace_ptr mf;
  MVertex_ptr mv;
  int num_tri, num_quad, num_penta;

  MSTK_Init();

  mesh = MESH_New(UNKNOWN_REP);
  ok = MESH_ImportFromFile(mesh,"serial/poly2.exo",NULL,NULL,NULL);
  CHECK_EQUAL(ok,1);

  CHECK(MESH_Num_Vertices(mesh) > 0);
  
  idx = 0; num_tri = 0; num_quad = 0; num_penta = 0;
  while ((mf = MESH_Next_Face(mesh,&idx))) {
    int ne = MF_Num_Edges(mf);
    if (ne == 3)
      num_tri++;
    else if (ne == 4)
      num_quad++;
    else if (ne == 5)
      num_penta++;
  }

  CHECK_EQUAL(num_tri,1);
  CHECK_EQUAL(num_quad,1);
  CHECK_EQUAL(num_penta,1);

  ok = MESH_CheckTopo(mesh);
  CHECK_EQUAL(ok, 1);

  ok = MESH_ExportToFile(mesh,"./poly2-tmp.exo",NULL,0,NULL,NULL,NULL);

  CHECK_EQUAL(ok,1);

  mesh2 = MESH_New(UNKNOWN_REP);
  ok = MESH_ImportFromFile(mesh2,"./poly2-tmp.exo",NULL,NULL,NULL);

  idx = 0; num_tri = 0; num_quad = 0; num_penta = 0;
  while ((mf = MESH_Next_Face(mesh2,&idx))) {
    int ne = MF_Num_Edges(mf);
    if (ne == 3)
      num_tri++;
    else if (ne == 4)
      num_quad++;
    else if (ne == 5)
      num_penta++;
  }

  CHECK_EQUAL(num_tri,1);
  CHECK_EQUAL(num_quad,1);
  CHECK_EQUAL(num_penta,1);

  ok = MESH_CheckTopo(mesh2);
  CHECK_EQUAL(ok, 1);

  MESH_Delete(mesh);
  MESH_Delete(mesh2);
}




TEST(Write_Read_ExodusII_HexMesh) {

  MFace_ptr mf;
  MVertex_ptr mv;
  MSet_ptr mset;

  MSTK_Init();

  /* Read a 2x2x2 mesh of hexes with full classification in MSTK format */ 

  Mesh_ptr mesh = MESH_New(UNKNOWN_REP);
  int ok = MESH_InitFromFile(mesh,"serial/reghex3D.mstk",NULL);
  CHECK_EQUAL(ok, 1);

  ok = MESH_CheckTopo(mesh);
  CHECK_EQUAL(ok, 1);

  /* Reclassify one layer of hexes as belonging to a different ID region */

  int elblock_numregs[2] = {6,2};
  int elblock_regids[2][7] = {{3,4,5,6,7,8},
                               {1,2,0,0,0,0,0}};
  
  for (int i = 0; i < 2; i++) {
    MRegion_ptr mr = MESH_Region(mesh,i);
    MR_Set_GEntID(mr,2);
  }

  /* Reclassify the faces between region 1 and region 2 as being on an
     internal surface */

  int idx = 0;
  while ((mf = MESH_Next_Face(mesh,&idx))) {

    if (MF_GEntDim(mf) == 3) {
      List_ptr fregs = MF_Regions(mf);
      if (MR_GEntID(List_Entry(fregs,0)) != MR_GEntID(List_Entry(fregs,1))) {
        MF_Set_GEntDim(mf,2);
        MF_Set_GEntID(mf,7);
      }
    }

  }


  /* Create two element sets */
  int elementset_regids[2][4];

  for (int i = 0; i < 2; i++) {
    char elementsetname[256];
    sprintf(elementsetname,"elemset_%-d",i+1);

    mset = MSet_New(mesh,elementsetname,MREGION);

    for (int j = 0; j < 4; j++) {
      MRegion_ptr mr = MESH_Region(mesh,i*4+j);
      MSet_Add(mset,mr);
      elementset_regids[i][j] = MR_ID(mr);
    }
  }

  /* Create 2 sidesets - 1 external (top surface) and one internal */

  idx = 0;
  double maxz = -1e8;
  while ((mv = MESH_Next_Vertex(mesh,&idx))) {
    double vxyz[3];
    MV_Coords(mv,vxyz);
    if (vxyz[2] > maxz)
      maxz = vxyz[2];
  }

  mset = MSet_New(mesh,"sideset_1",MFACE);
  idx = 0;
  while ((mf = MESH_Next_Face(mesh,&idx))) {
    double fxyz[MAXPV2][3];
    int nfv;
    int zmatch = 1;
    MF_Coords(mf,&nfv,fxyz);
    for (int i = 0; i < nfv; i++) {
      if (fxyz[i][2] != maxz) {
        zmatch=0;
        break;
      }
    }
    if (zmatch) MSet_Add(mset,mf);
  }

  mset = MSet_New(mesh,"sideset_2",MFACE);
    
  idx = 0;
  while ((mf = MESH_Next_Face(mesh,&idx))) {
    if (MF_GEntDim(mf) == 2 && MF_GEntID(mf) == 7) {
      MSet_Add(mset,mf);
    }
  }

  /* Create 2 nodesets - one composed of nodes on faces and the other 
   composed of nodes on edges */

  int nodeset_numverts[2] = {12,6};
  int nodeset_vertids[2][12] = {{2,4,6,8,10,12,16,18,20,22,24,26},
                                {5,11,13,15,17,23,0,0,0,0,0,0}};
  
  for (int i = 0; i < 2; i++) {
    char nodesetname[256];
    sprintf(nodesetname,"nodeset_%-d",11*(i+1));
    
    mset = MSet_New(mesh,nodesetname,MVERTEX);

    idx = 0;
    int j = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx))) {
      if (MV_GEntDim(mv) == i+1) {
        MSet_Add(mset,mv);
        nodeset_vertids[i][j] = MV_ID(mv);
        j++;
      }
    }
  }

  /* Now export to an Exodus II file */

  ok = MESH_ExportToFile(mesh,"temp.exo","exodusii",0,NULL,NULL,NULL);


  /* Now create another mesh and import this file back */

  Mesh_ptr mesh2 = MESH_New(F1);

  ok = MESH_ImportFromFile(mesh2,"temp.exo","exodusii",NULL,NULL);
  CHECK_EQUAL(ok, 1);

  ok = MESH_CheckTopo(mesh2);
  CHECK_EQUAL(ok, 1);

  /* Now verify that we retrieved all the model regions (element
     blocks) as expected */

  int nr = MESH_Num_Regions(mesh2);
  CHECK_EQUAL(nr,MESH_Num_Regions(mesh));

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < elblock_numregs[i]; j++) {
      int rid = elblock_regids[i][j];
      MRegion_ptr mr = MESH_RegionFromID(mesh2,rid);
      CHECK(mr);
      CHECK_EQUAL(MR_GEntID(mr),i+1);
    }
  }

  /* Verify that we retrieved element sets as expected */

  for (int i = 0; i < 2; i++) {
    char elementsetname[256];
    sprintf(elementsetname,"elemset_%-d",i+1);

    mset = MESH_MSetByName(mesh2,elementsetname);
    CHECK(mset);

    for (int j = 0; j < 4; j++) {
      MRegion_ptr mr = (MRegion_ptr) MSet_Entry(mset,j);
      int rid = MR_ID(mr);
      int found = 0;
      for (int k = 0; k < 4; k++) {
        if (elementset_regids[i][k] == rid) {
          found = 1;
          break;
        }
      }
      CHECK_EQUAL(found,1);
    }
  }


  /* Verify that we retrieved side sets as expected */

  // Verify that sideset 1 has faces on the top surface

  mset = MESH_MSetByName(mesh2,"sideset_1");
  CHECK(mset);
  CHECK_EQUAL(4,MSet_Num_Entries(mset));
  for (int i = 0; i < 4; i++) {
    double fxyz[MAXPV2][3];
    int nfv, zmatch=1;
 
    mf = MSet_Entry(mset,i);
    MF_Coords(mf,&nfv,fxyz);
    
    for (int j = 0; j < nfv; j++) {
      if (fxyz[i][2] != maxz) {
        zmatch = 0;
        break;
      }
    }
    CHECK(zmatch);
  }

  // Verify that sideset 2 has only interior faces

  mset = MESH_MSetByName(mesh2,"sideset_2");
  CHECK(mset);
  CHECK_EQUAL(4,MSet_Num_Entries(mset));
  for (int i = 0; i < 4; i++) {
    List_ptr fregs = MF_Regions(MSet_Entry(mset,i));
    CHECK_EQUAL(2,List_Num_Entries(fregs));
    List_Delete(fregs);
  }  


  /* Verify that we retrieved node sets as expected */
  /* We cannot ensure that the node IDs in the original 
     MSTK mesh are the same as those in the mesh imported
     Exodus II files - so we will only match coordinates */

  for (int i = 0; i < 2; i++) {
    char nodesetname[256];
    sprintf(nodesetname,"nodeset_%-d",11*(i+1));

    mset = MESH_MSetByName(mesh2,nodesetname);
    CHECK(mset);

    for (int j = 0; j < nodeset_numverts[i]; j++) {
      mv = (MVertex_ptr) MSet_Entry(mset,j);
      double xyz[3];
      MV_Coords(mv,xyz);
      int found = 0;
      for (int k = 0; k < nodeset_numverts[i]; k++) {
        MVertex_ptr mv2 = MESH_VertexFromID(mesh,nodeset_vertids[i][k]);
        double xyz2[3];
        MV_Coords(mv2,xyz2);
        double dist2 = ((xyz[0]-xyz2[0])*(xyz[0]-xyz2[0]) +
                        (xyz[1]-xyz2[1])*(xyz[1]-xyz2[1]) +
                        (xyz[2]-xyz2[2])*(xyz[2]-xyz2[2]));
        if (dist2 < 1.0e-16) {
          found = 1;
          break;
        }
      }
      CHECK_EQUAL(found,1);
    }
  }

  MESH_Delete(mesh);
  MESH_Delete(mesh2);
}


// Read/write a layered mesh with (1) a mixture of polyhedral and
// standard elements abutting each other in different blocks -
// standard elements in a higher block number than the polyhedral
// elements (2) disjoint groups of polyhedral elements in a single
// element block (3) degenerate and zero volume polyhedral
// elements. All of these conditions stress the part of the code that
// determines region-face directions for polyhedral elements

TEST(Write_Read_ExodusII_DegeneratePoly3) {
  int ok;
  Mesh_ptr mesh1, mesh2;

  MSTK_Init();

  mesh1 = MESH_New(UNKNOWN_REP);
  ok = MESH_ImportFromFile(mesh1,"serial/degenpoly3.exo",NULL,NULL,NULL);
  CHECK_EQUAL(ok,1);

  int use_geometry = 1;
  ok &= MESH_BuildClassfn(mesh1, use_geometry);
  ok &= MESH_CheckTopo(mesh1);
  CHECK_EQUAL(ok, 1);


  ok = MESH_ExportToFile(mesh1,"./degenpoly3-tmp.exo",NULL,0,NULL,NULL,NULL);
  CHECK_EQUAL(ok,1);


  mesh2 = MESH_New(UNKNOWN_REP);
  ok = MESH_ImportFromFile(mesh2,"./degenpoly3-tmp.exo",NULL,NULL,NULL);
  CHECK_EQUAL(ok, 1);

  ok = MESH_CheckTopo(mesh2);
  CHECK_EQUAL(ok, 1);

  MESH_Delete(mesh1);
  MESH_Delete(mesh2);
}


// Read/write a multi-block polyhedral mesh where the blocks are disjoint

TEST(Write_Read_ExodusII_MultiBlockPoly3) {
  int ok;
  Mesh_ptr mesh1, mesh2;

  MSTK_Init();

  mesh1 = MESH_New(UNKNOWN_REP);
  ok = MESH_ImportFromFile(mesh1,"serial/multiblock-poly3.exo",NULL,NULL,NULL);
  CHECK_EQUAL(ok,1);

  int use_geometry = 1;
  ok &= MESH_BuildClassfn(mesh1, use_geometry);
  ok &= MESH_CheckTopo(mesh1);
  CHECK_EQUAL(ok, 1);


  ok = MESH_ExportToFile(mesh1,"./multiblock-poly3-tmp.exo",NULL,0,NULL,NULL,NULL);
  CHECK_EQUAL(ok,1);


  mesh2 = MESH_New(UNKNOWN_REP);
  ok = MESH_ImportFromFile(mesh2,"./multiblock-poly3-tmp.exo",NULL,NULL,NULL);
  CHECK_EQUAL(ok, 1);

  ok = MESH_CheckTopo(mesh2);
  CHECK_EQUAL(ok, 1);

  MESH_Delete(mesh1);
  MESH_Delete(mesh2);
}


TEST(Write_Read_ExodusII_Variables) {
  int i, idx;
  
  MSTK_Init();

  /* Read a 2x2x2 mesh of hexes with full classification in MSTK format */ 

  Mesh_ptr mesh = MESH_New(UNKNOWN_REP);
  int ok = MESH_InitFromFile(mesh,"serial/reghex3D.mstk",NULL);
  CHECK_EQUAL(ok, 1);

  ok = MESH_CheckTopo(mesh);
  CHECK_EQUAL(ok, 1);

  /* Reclassify one layer of hexes as belonging to a different ID region */

  int elblock_numregs[2] = {6,2};
  int elblock_regids[2][7] = {{3,4,5,6,7,8},
                               {1,2,0,0,0,0,0}};
  
  MRegion_ptr mr;
  for (int i = 0; i < 2; i++) {
    mr = MESH_Region(mesh,i);
    MR_Set_GEntID(mr,2);
  }

  int nr = MESH_Num_Regions(mesh);
  int nv = MESH_Num_Vertices(mesh);

  /* Create a scalar and a vector attribute each on elements and nodes */

  MAttrib_ptr elvalatt_out = MAttrib_New(mesh,"elval",DOUBLE,MREGION);
  MAttrib_ptr elvecatt_out = MAttrib_New(mesh,"elvec",VECTOR,MREGION,3);

  double *elval_out = (double *) new double[nr];
  double **elvec_out = (double **) new double *[nr];
  for (int k = 0; k < nr; k++)
    elvec_out[k] = new double[3];

  idx = 0; i = 0;
  while ((mr = MESH_Next_Region(mesh,&idx))) {
    elval_out[i] = 2.5*MR_ID(mr);
    elvec_out[i][0] = elval_out[i];
    elvec_out[i][1] = elval_out[i]+1;
    elvec_out[i][2] = elval_out[i]+2;

    MEnt_Set_AttVal(mr,elvalatt_out,0,elval_out[i],NULL);
    MEnt_Set_AttVal(mr,elvecatt_out,0,0.0,elvec_out[i]);
    i++;
  }


  double *ndval_out = (double *) new double[nv];
  double **ndvec_out = (double **) new double *[nv];
  for (int k = 0; k < nv; k++) 
    ndvec_out[k] = (double *) new double[3];

  MAttrib_ptr ndvalatt_out = MAttrib_New(mesh,"ndval",DOUBLE,MVERTEX);
  MAttrib_ptr ndvecatt_out = MAttrib_New(mesh,"ndvec",VECTOR,MVERTEX,3);

  MVertex_ptr mv;
  idx = 0; i = 0;
  while ((mv = MESH_Next_Vertex(mesh,&idx))) {
    ndval_out[i] = 2.5*MV_ID(mv);
    ndvec_out[i][0] = ndval_out[i];
    ndvec_out[i][1] = ndval_out[i]+1;
    ndvec_out[i][2] = ndval_out[i]+2;

    MEnt_Set_AttVal(mv,ndvalatt_out,0,ndval_out[i],NULL);
    MEnt_Set_AttVal(mv,ndvecatt_out,0,0.0,ndvec_out[i]);
    i++;
  }

  /* Now export to an Exodus II file */

  ok = MESH_ExportToFile(mesh,"temp.exo","exodusii",0,NULL,NULL,NULL);
  CHECK_EQUAL(ok, 1);


  /* Now create another mesh and import this file back */

  Mesh_ptr mesh2 = MESH_New(F1);

  ok = MESH_ImportFromFile(mesh2,"temp.exo","exodusii",NULL,NULL);
  CHECK_EQUAL(ok, 1);

  ok = MESH_CheckTopo(mesh2);
  CHECK_EQUAL(ok, 1);

  /* Now verify that we retrieved all the attributes on elements -
     because of the way elements are written and read back into Exodus
     II files (by element block) the ordering stays the same */

  MAttrib_ptr elvalatt_in = MESH_AttribByName(mesh2,"elval");
  CHECK(elvalatt_in);
  CHECK_EQUAL(DOUBLE,MAttrib_Get_Type(elvalatt_in));
  CHECK_EQUAL(MREGION,MAttrib_Get_EntDim(elvalatt_in));

  MAttrib_ptr elvecatt_in = MESH_AttribByName(mesh2,"elvec");
  CHECK(elvecatt_in);
  CHECK_EQUAL(VECTOR,MAttrib_Get_Type(elvecatt_in));
  CHECK_EQUAL(MREGION,MAttrib_Get_EntDim(elvecatt_in));
  CHECK_EQUAL(3,MAttrib_Get_NumComps(elvecatt_in));

 
  idx = 0; i = 0;
  while ((mr = MESH_Next_Region(mesh2,&idx))) {
    int ival;
    double rval;
    void *pval;

    MEnt_Get_AttVal(mr,elvalatt_in,&ival,&rval,&pval);
    CHECK_EQUAL(elval_out[i],rval);

    double vec[3];
    MEnt_Get_AttVal(mr,elvecatt_in,&ival,&rval,&pval);
    CHECK_ARRAY_EQUAL(elvec_out[i],(double *)pval,3);
    i++;
  }


  /* Now verify that we retrieved all the attributes on nodes -
   * because of the way nodes are written to Exodus II and read back,
   * the ordering of nodes may not be the same and we have to figure
   * out a map between the nodes of the original mesh and the
   * re-imported mesh */

  MAttrib_ptr ndvalatt_in = MESH_AttribByName(mesh2,"ndval");
  CHECK(ndvalatt_in);
  CHECK_EQUAL(DOUBLE,MAttrib_Get_Type(ndvalatt_in));
  CHECK_EQUAL(MVERTEX,MAttrib_Get_EntDim(ndvalatt_in));

  MAttrib_ptr ndvecatt_in = MESH_AttribByName(mesh2,"ndvec");
  CHECK(ndvecatt_in);
  CHECK_EQUAL(VECTOR,MAttrib_Get_Type(ndvecatt_in));
  CHECK_EQUAL(MVERTEX,MAttrib_Get_EntDim(ndvecatt_in));
  CHECK_EQUAL(3,MAttrib_Get_NumComps(ndvecatt_in));

  double *node_map = new double[nv];
  for (i = 0; i < nv; i++) node_map[i] = -1;
  idx = 0; i = 0;
  while ((mv = MESH_Next_Vertex(mesh, &idx))) {
    double coord[3];
    MV_Coords(mv, coord);
    
    int idx2 = 0, k = 0, found = 0;
    MVertex_ptr mv2;
    while ((mv2 = MESH_Next_Vertex(mesh2, &idx2))) {
      double coord2[3];
      MV_Coords(mv2, coord2);

      double dist2 = ((coord[0]-coord2[0])*(coord[0]-coord2[0]) +
                      (coord[1]-coord2[1])*(coord[1]-coord2[1]) +
                      (coord[2]-coord2[2])*(coord[2]-coord2[2]));
      if (dist2 < 1.0e-12) {
        found = 1;
        break;
      }
      k++;
    }
    CHECK(found);
    node_map[k] = i;
    i++;
  }
 
  idx = 0; i = 0;
  while ((mv = MESH_Next_Vertex(mesh2,&idx))) {
    int ival;
    double rval;
    void *pval;

    int k = node_map[i];
    MEnt_Get_AttVal(mv,ndvalatt_in,&ival,&rval,&pval);
    CHECK_EQUAL(ndval_out[k],rval);

    MEnt_Get_AttVal(mv,ndvecatt_in,&ival,&rval,&pval);
    CHECK_ARRAY_EQUAL(ndvec_out[k],(double *)pval,3);
    i++;
  }

  delete [] elval_out;
  for (int k = 0; k < nr; k++)
    delete [] elvec_out[k];
  delete [] elvec_out;
  delete [] ndval_out;
  for (int k = 0; k < nv; k++)
    delete [] ndvec_out[k];
  delete [] ndvec_out;

  MESH_Delete(mesh);
  MESH_Delete(mesh2);
}


TEST(Write_Read_ExodusII_Partial_2D) {

  MSTK_Init();

  /* Read a 5x5 mesh of quads */

  Mesh_ptr mesh = MESH_New(F1);
  int ok = MESH_ReadExodusII_Serial(mesh,"parallel/4proc/mesh5x5-skewed.exo",0);
  CHECK_EQUAL(ok, 1);

  int nelems_total = MESH_Num_Faces(mesh);

  /* Create a subset of elements to do a partial read */

  int *temp_ptr = new int;
  srand((long)temp_ptr);
  delete temp_ptr;

  for (int i = 0; i < 3; i++) {
    int nelems = rand() % nelems_total + 1;
    
    int *elems = new int[nelems];
    for (int j = 0; j < nelems; j++) {  // generate nelems unique IDs
      int done = 0;
      while (!done) {
	int id = rand() % nelems_total + 1;
	int found = 0;
	for (int k = 0; k < j; k++)
	  if (elems[k] == id) {
	    found = 1;
	    break;
	  }
	if (!found) {
	  done = 1;  // found a unique ID
	  elems[j] = id;
	}
      }
    }
	
    Mesh_ptr mesh2 = MESH_New(F1);
    ok = MESH_ReadExodusII_Partial(mesh2,"parallel/4proc/mesh5x5-skewed.exo",0,
				    nelems, elems);
    CHECK_EQUAL(ok, 1);
    MESH_Delete(mesh2);
    delete [] elems;
  }

  MESH_Delete(mesh);
}

TEST(Write_Read_ExodusII_Partial_3D) {

  MSTK_Init();

  /* Read a 5x5 mesh of quads */

  Mesh_ptr mesh = MESH_New(F1);
  int ok = MESH_ReadExodusII_Serial(mesh,"parallel/8proc/hex_3x3x3_ss.exo",0);
  CHECK_EQUAL(ok, 1);

  int nelems_total = MESH_Num_Regions(mesh);

  /* Create a subset of elements to do a partial read */

  int *temp_ptr = new int;
  srand((long)temp_ptr);
  delete temp_ptr;

  for (int i = 0; i < 3; i++) {
    int nelems = rand() % nelems_total + 1;
    
    int *elems = new int[nelems];
    for (int j = 0; j < nelems; j++) {  // generate nelems unique IDs
      int done = 0;
      while (!done) {
	int id = rand() % nelems_total + 1;
	int found = 0;
	for (int k = 0; k < j; k++)
	  if (elems[k] == id) {
	    found = 1;
	    break;
	  }
	if (!found) {
	  done = 1;  // found a unique ID
	  elems[j] = id;
	}
      }
    }
	
    Mesh_ptr mesh2 = MESH_New(F1);
    ok = MESH_ReadExodusII_Partial(mesh2,"parallel/8proc/hex_3x3x3_ss.exo",0,
				    nelems, elems);
    CHECK_EQUAL(ok, 1);
    MESH_Delete(mesh2);
    delete [] elems;
  }

  MESH_Delete(mesh);
}
