#include <iostream>
#include <fstream>

#include <UnitTest++.h>
#include <exodusII.h>

#include "MSTK.h"

// Test if we can correctly read a polygonal exodus mesh 
// Then test if we can correctly write it and read it back

int LocalFaceNumber(MFace_ptr mf, MRegion_ptr mr);
int mstk2exo_face_map[5][6]={{0,0,0,0,0,0},
			    {4,1,2,3,0,0}, /* TET */
			    {0,0,0,0,0,0}, /* PYRAMID, no support in exo*/
		            {4,5,1,2,3,0}, /* PRISM, must verify nums */
		            {5,6,1,2,3,4}};/* HEX */

TEST(Verify_Exoatt) 
{
  int idx, ok;
  Mesh_ptr mesh, mesh2;
  MFace_ptr mf;
  MRegion_ptr mr;
  MVertex_ptr mv;
  char attribname[64] = "myvars";
  int cellsetid = 10001, sidesetid = 909;

  MSTK_Init();

  // Import an exodus II mesh

  mesh = MESH_New(UNKNOWN_REP);
  ok = MESH_ImportFromFile(mesh,"reghex3D.exo","exodusii",NULL,NULL);
  CHECK_EQUAL(ok,1);

  int nr = MESH_Num_Regions(mesh);
  CHECK(nr);                    // assume this is a hex mesh of a box
  
  // Write out an attribute file for the regions 
  // The attribute value is 1.0 for half of the regions and 7.7 for the others

  int nr2 = (int) (nr/2.0);

  std::filebuf fb1;
  fb1.open("attribs.txt",std::ios::out);
  std::ostream os1(&fb1);
  os1 << "1\n";  // Number of attributes
  os1 << attribname << " CELL -1\n";  // attribute with attribname on all 
                                     // cells of the mesh (hence -1)
  
  for (int i = 0; i < nr2; i++) 
    os1 << "1.0\n";
  for (int i = nr2; i < nr; i++)
    os1 << "7.7\n";

  fb1.close();

  // Write out entity sets file

  std::filebuf fb2;
  fb2.open("sets.txt",std::ios::out);
  std::ostream os2(&fb2);
  os2 << "2\n";  // Two sets to be read
  
  // write out cell set containing the first nr2 cells

  os2 << cellsetid << " CELL " << nr2 << "\n";
  for (int i = 0; i < nr2; i++)
    os2 << i+1 << "\n";

  // write out faces of the top surface

  double maxz = -1.0e8;
  idx = 0;
  while ((mv = MESH_Next_Vertex(mesh,&idx))) {
    double xyz[3];
    MV_Coords(mv,xyz);
    if (xyz[2] > maxz)
      maxz = xyz[2];
  }

  List_ptr topfaces = List_New(10);
  idx = 0;
  while ((mf = MESH_Next_Face(mesh,&idx))) {
    int nfv;
    double fxyz[4][3];
    int zmatch=1;

    MF_Coords(mf,&nfv,fxyz);
    for (int i = 0; i < nfv; i++)
      if (fxyz[i][2] != maxz) {
        zmatch = 0;
        break;
      }
    if (zmatch) List_Add(topfaces,mf);
  }

  os2 << sidesetid << " FACE " << List_Num_Entries(topfaces) << "\n";
  idx = 0;
  while ((mf = List_Next_Entry(topfaces,&idx))) {
    List_ptr fverts = MF_Vertices(mf,1,0);
    int nfv = List_Num_Entries(fverts);
    os2 << nfv;
    for (int i = 0; i < nfv; i++)
      os2 << " " << MV_ID(List_Entry(fverts,i));
    os2 << "\n";
    List_Delete(fverts);
  }
  List_Delete(topfaces);

  fb2.close();
  
  

  // Run the exoatt command and augment input.exo with attributes and entity sets
  system("../exoatt --attfile=attribs.txt --setfile=sets.txt reghex3D.exo output.exo"); 



  // Next read the output file from exoatt and verify that the mesh has the 
  // attributes and entity sets

  MSTK_Init();

  mesh2 = MESH_New(UNKNOWN_REP);
  ok = MESH_ImportFromFile(mesh2,"output.exo","exodusii",NULL,NULL);
  CHECK_EQUAL(ok,1);


  // Verify the cell set

  std::stringstream cellsetname;
  cellsetname << "elemset_" << cellsetid;
  MSet_ptr mset = MESH_MSetByName(mesh2,cellsetname.str().c_str());
  CHECK(mset);
  CHECK_EQUAL(nr2,MSet_Num_Entries(mset));

  // Verify the side set

  std::stringstream sidesetname;
  sidesetname << "sideset_" << sidesetid;
  mset = MESH_MSetByName(mesh2,sidesetname.str().c_str());
  CHECK(mset);
  idx = 0;
  while ((mf = (MFace_ptr) MSet_Next_Entry(mset,&idx))) {
    int nfv;
    double fxyz[4][3];
    int zmatch=1;

    MF_Coords(mf,&nfv,fxyz);

    for (int i = 0; i < 4; i++)
      if (fxyz[i][2] != maxz) {
        zmatch = 0;
        break;
      }

    CHECK(zmatch);
  }

  // Since we don't import entity variables into the mesh, we have to
  // read them directly from the Exodus II file to verify

  int exoid=0, status;
  int comp_ws = sizeof(double), io_ws = 0;
  float version;

  exoid = ex_open("output.exo", EX_READ, &comp_ws, &io_ws, &version);
  CHECK(exoid);

  int nvars;
#ifdef EXODUS_6_DEPRECATED
  status = ex_get_var_param(exoid, "e", &nvars);
#else
  status = ex_get_variable_param(exoid, EX_ELEM_BLOCK, &nvars);
#endif
  CHECK_EQUAL(1,nvars);

  char var_name[256];
#ifdef EXODUS_6_DEPRECATED
  status = ex_get_var_name(exoid, "e", 1, var_name);
#else
  status = ex_get_variable_name(exoid, EX_ELEM_BLOCK, 1, var_name);
#endif
  CHECK_EQUAL(attribname,var_name);

  double var;
  for (int i = 0; i < nr2; i++) {
#ifdef EXODUS_6_DEPRECATED
    status = ex_get_elem_var_time(exoid, 1, i+1, 1, 1, &var);
#else
    status = ex_get_var_time(exoid, EX_ELEM_BLOCK, 1, i+1, 1, 1, &var);
#endif
    CHECK_EQUAL(1.0,var);
  }

  for (int i = nr2; i < nr; i++) {
#ifdef EXODUS_6_DEPRECATED
    status = ex_get_elem_var_time(exoid, 1, i+1, 1, 1, &var);
#else
    status = ex_get_var_time(exoid, EX_ELEM_BLOCK, 1, i+1, 1, 1, &var);
#endif
    CHECK_EQUAL(7.7,var);
  }

  ex_close(exoid);
  
}
      

  /* Local ID of face in region according to MSTK conventions - Zero based */

  int LocalFaceNumber(MFace_ptr mf, MRegion_ptr mr) {
    MRType mrtype = MR_ElementType(mr);

    switch (mrtype) {
    case TET: case PYRAMID: case PRISM: case HEX: {
      int i, j;
      int allfound = 0;

      List_ptr rverts = MR_Vertices(mr);
      List_ptr fverts = MF_Vertices(mf,1,0);

      int nrf = MSTK_nrf_template[mrtype];
      for (i = 0; i < nrf; i++) {
        int nrfv = MSTK_rfv_template[mrtype][i][0];

        allfound = 1;
        for (j = 0; j < nrfv; j++) {
          MVertex_ptr rv = List_Entry(rverts,MSTK_rfv_template[mrtype][i][j+1]);
          if (!List_Contains(fverts,rv)) {
            allfound = 0;
            break;
          }
        }
        if (allfound)
          break;
      }
      
      List_Delete(rverts);
      List_Delete(fverts);

      if (allfound)
        return i;
      else {
        MSTK_Report("MF_LocalID_in_Region","Face not found in region",MSTK_ERROR);
        return -1;
      }
      break;
    }
    case POLYHED: case RUNKNOWN: {
      List_ptr rfaces = MR_Faces(mr);
      int i, found = 0;
      for (i = 0; i < List_Num_Entries(rfaces); i++) {
        if (mf == List_Entry(rfaces,i)) { /* will this work for reduced representations? */
          found = 1;
          break;
        }
      }
      List_Delete(rfaces);
      if (found)
        return i;
      else {
        MSTK_Report("MF_LocalID_in_Region","Face not found in region",MSTK_ERROR);
        return -1;
      }
      break;
    }
    case RDELETED: default:
      return -1;
    }
  }
	    

