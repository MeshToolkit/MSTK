#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Create a new hexmesh from a tet mesh by splitting edges, faces
     and regions. This will allocate a new mesh pointer and assign it
     to hexmesh

     Konstantin Lipnikov, lipnikov@lanl.gov
     10/14/2005     
*/

  int MESH_Tet2Hex(Mesh_ptr tetmesh, Mesh_ptr *hexmesh) { 

    int          i, j, k, n, nV, nR, idv, idf, ide, idr;
    double       xyz[3], xyz_lst[100][3];
    
    MVertex_ptr  ptr_v, ptr_v_new, ptr_v_face[3];
    MVertex_ptr  old_ptr_v_regn[4], ptr_v_regn[15], ptr_v_tet[8];
    MEdge_ptr    ptr_e;
    MFace_ptr    ptr_f;
    MRegion_ptr  ptr_r, ptr_r_new;
    
    List_ptr     lst_v;
    
    MAttrib_ptr  atr_v, atr_e, atr_f, atr_r;
    
    int          ival;
    double       dval;
    void        *vval;
    
    static int tet_hex[4][8] = {{0, 4, 13, 5, 6, 12, 14, 11},
				{1, 7, 13, 4, 8, 10, 14, 12},
				{2, 5, 13, 7, 9, 11, 14, 10},
				{3, 6, 11, 9, 8, 12, 14, 10}};
    
    nV = MESH_Num_Vertices(tetmesh);
    nR = MESH_Num_Regions(tetmesh);
    
    *hexmesh = MESH_New(MESH_RepType(tetmesh));
    
    
    /* copy vertices */

    atr_v = MAttrib_New(tetmesh, "vertex", POINTER, MVERTEX);

    idv = 0;
    while( (ptr_v = MESH_Next_Vertex(tetmesh, &idv)) ) {
      MV_Coords(ptr_v, xyz);
      
      ptr_v_new = MV_New(*hexmesh);
      MV_Set_Coords(ptr_v_new, xyz);
      MV_Set_GEntDim(ptr_v_new,MV_GEntDim(ptr_v));
      MV_Set_GEntID(ptr_v_new,MV_GEntID(ptr_v));

      MEnt_Set_AttVal(ptr_v, atr_v, 0, 0, (void *)ptr_v_new);
    }


    /* add new edge vertices */
    atr_e = MAttrib_New(tetmesh, "edge", POINTER, MEDGE);

    ide = 0;
    while( (ptr_e = MESH_Next_Edge(tetmesh, &ide)) ) {
      for( i=0; i<2; i++ ) {
        ptr_v = ME_Vertex(ptr_e, i);
        MV_Coords(ptr_v, xyz_lst[i]);
      }

      for( i=0; i<3; i++ ) xyz[i] = (xyz_lst[0][i] + xyz_lst[1][i]) / 2;

      ptr_v_new = MV_New(*hexmesh);
      MV_Set_Coords(ptr_v_new, xyz);
      MV_Set_GEntDim(ptr_v_new,ME_GEntDim(ptr_e));
      MV_Set_GEntID(ptr_v_new,ME_GEntID(ptr_e));

      MEnt_Set_AttVal(ptr_e, atr_e, 0, 0, (void*)ptr_v_new);
    }

  
    /* add new face vertices */
    atr_f = MAttrib_New(tetmesh, "face", POINTER, MFACE);

    idf = 0;
    while( (ptr_f = MESH_Next_Face(tetmesh, &idf)) ) {
      MF_Coords(ptr_f, &n, xyz_lst);
      if( n != 3 ) {
	MSTK_Report("MESH_Tet2Hex","Input mesh is not a tetrahedral mesh",MSTK_WARN);
	return 0;
      }

      for( i=0; i<3; i++ ) xyz[i] = (xyz_lst[0][i] + xyz_lst[1][i] + xyz_lst[2][i]) / 3;

      ptr_v_new = MV_New(*hexmesh);
      MV_Set_Coords(ptr_v_new, xyz);
      MV_Set_GEntDim(ptr_v_new,MF_GEntDim(ptr_f));
      MV_Set_GEntID(ptr_v_new,MF_GEntID(ptr_f));


      MEnt_Set_AttVal(ptr_f, atr_f, 0, 0, (void*)ptr_v_new);
    }


    /* add new region vertex */
    atr_r = MAttrib_New(tetmesh, "region", POINTER, MREGION);

    idr = 0;
    while( (ptr_r = MESH_Next_Region(tetmesh, &idr)) ) {
      MR_Coords(ptr_r, &n, xyz_lst);

      for( i=0; i<3; i++ ) xyz[i] = 0;

      for( k=0; k<n; k++ ) 
	for( i=0; i<3; i++ ) xyz[i] += xyz_lst[k][i];

      for( i=0; i<3; i++ ) xyz[i] /= n;

      ptr_v_new = MV_New(*hexmesh);
      MV_Set_Coords(ptr_v_new, xyz);
      MV_Set_GEntDim(ptr_v_new,MR_GEntDim(ptr_r));
      MV_Set_GEntID(ptr_v_new,MR_GEntID(ptr_r));

      MEnt_Set_AttVal(ptr_r, atr_r, 0, 0, (void*)ptr_v_new);
    }


    /* convert tets to hexes */
    idr = 0;
    while( (ptr_r = MESH_Next_Region(tetmesh, &idr)) ) {
      lst_v = MR_Vertices(ptr_r);
      n = List_Num_Entries(lst_v);

      for( i=0; i<4; i++ ) {
	old_ptr_v_regn[i] = List_Entry(lst_v, i);
	MEnt_Get_AttVal(old_ptr_v_regn[i], atr_v, 0, 0, &(ptr_v_regn[i]));
      }

      for( i=0;   i<4; i++ )
	for( j=i+1; j<4; j++ ) {
	  ptr_e = MVs_CommonEdge(old_ptr_v_regn[i], old_ptr_v_regn[j]);

	  MEnt_Get_AttVal(ptr_e, atr_e, &ival, &dval, &vval);
	  ptr_v_regn[n] = (MVertex_ptr)(vval);
	  MV_Coords(ptr_v_regn[n], xyz);
	  n++;
	}

      for( i=0; i<4; i++ ) {
        for( j=0; j<3; j++ ) 
	  ptr_v_face[j] = old_ptr_v_regn[(i + j + 1) % 4];
        ptr_f = MVs_CommonFace(3, ptr_v_face);

        MEnt_Get_AttVal(ptr_f, atr_f, &ival, &dval, &vval);
        ptr_v_regn[n] = (MVertex_ptr)(vval);
        n++;
      }

      MEnt_Get_AttVal(ptr_r, atr_r, &ival, &dval, &vval);
      ptr_v_regn[n] = (MVertex_ptr)(vval);
      n++;

      for( k=0; k<4; k++ ) {
        for( i=0; i<8; i++ ) ptr_v_tet[i] = ptr_v_regn[tet_hex[k][i]];

	ptr_r_new = MR_New(*hexmesh);

        MR_Set_Vertices(ptr_r_new, 8, ptr_v_tet, 0, NULL);

	MR_Set_GEntDim(ptr_r_new,MR_GEntDim(ptr_r));
	MR_Set_GEntID(ptr_r_new,MR_GEntID(ptr_r));
      }

      List_Delete(lst_v);
    }



    /* remove attributes on edges of tet mesh */
    idv = 0;
    while( (ptr_v = MESH_Next_Vertex(tetmesh, &idv)) ) {
      MEnt_Rem_AttVal(ptr_v, atr_v);
    }
    MAttrib_Delete(atr_v);

    ide = 0;
    while( (ptr_e = MESH_Next_Edge(tetmesh, &ide)) ) {
      MEnt_Rem_AttVal(ptr_e, atr_e);
    }
    MAttrib_Delete(atr_e);

  
    /* remove attributes on faces of tet mesh */
    idf = 0;
    while( (ptr_f = MESH_Next_Face(tetmesh, &idf)) ) {
      MEnt_Rem_AttVal(ptr_f, atr_f);
    }
    MAttrib_Delete(atr_f);

    /* remove attributes on regions of tet mesh */
    idr = 0;
    while( (ptr_r = MESH_Next_Region(tetmesh, &idr)) ) {
      MEnt_Rem_AttVal(ptr_r, atr_r);
    }
    MAttrib_Delete(atr_r);

  
    return 1;
  }

#ifdef _cplusplus
}
#endif
