#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Create a new quadmesh from a tri mesh by splitting edges and faces.
     This will allocate a new mesh pointer and assign it to quadmesh

     Rao Garimella
     Adapted from MESH_Tet2Hex written by Konstantin Lipnikov
*/

  int MESH_Tri2Quad(Mesh_ptr trimesh, Mesh_ptr *quadmesh) { 

    int          i, j, k, n, nV, nF, idv, ide, idf;
    double       xyz[3], xyz_lst[100][3];
    
    MVertex_ptr  ptr_v, ptr_v_new, ptr_v_face[7];
    MVertex_ptr  old_ptr_v_face[4], ptr_v_quad[8];
    MEdge_ptr    ptr_e;
    MFace_ptr    ptr_f, ptr_f_new;
    
    List_ptr     lst_v, lst_e;
    
    MAttrib_ptr  atr_v, atr_e, atr_f;
    
    int          ival;
    double       dval;
    void        *vval;
    
    static int tri_quad[3][4] = {{0, 3, 6, 5},
				 {3, 1, 4, 6},
				 {5, 6, 4, 2}};

    
    nV = MESH_Num_Vertices(trimesh);
    nF = MESH_Num_Faces(trimesh);
    
    *quadmesh = MESH_New(MESH_RepType(trimesh));
    
    
    /* copy vertices */

    atr_v = MAttrib_New(trimesh, "vertex", POINTER, MVERTEX);

    idv = 0;
    while( (ptr_v = MESH_Next_Vertex(trimesh, &idv)) ) {
      MV_Coords(ptr_v, xyz);
      
      ptr_v_new = MV_New(*quadmesh);
      MV_Set_Coords(ptr_v_new, xyz);
      MV_Set_GEntDim(ptr_v_new,MV_GEntDim(ptr_v));
      MV_Set_GEntID(ptr_v_new,MV_GEntID(ptr_v));

      MEnt_Set_AttVal(ptr_v, atr_v, 0, 0, (void *)ptr_v_new);
    }


    /* add new edge vertices */
    atr_e = MAttrib_New(trimesh, "edge", POINTER, MEDGE);

    ide = 0;
    while( (ptr_e = MESH_Next_Edge(trimesh, &ide)) ) {
      for( i=0; i<2; i++ ) {
        ptr_v = ME_Vertex(ptr_e, i);
        MV_Coords(ptr_v, xyz_lst[i]);
      }

      for( i=0; i<3; i++ ) xyz[i] = (xyz_lst[0][i] + xyz_lst[1][i]) / 2;

      ptr_v_new = MV_New(*quadmesh);
      MV_Set_Coords(ptr_v_new, xyz);
      MV_Set_GEntDim(ptr_v_new,ME_GEntDim(ptr_e));
      MV_Set_GEntID(ptr_v_new,ME_GEntID(ptr_e));

      MEnt_Set_AttVal(ptr_e, atr_e, 0, 0, (void*)ptr_v_new);
    }

  
    /* add new face vertices */
    atr_f = MAttrib_New(trimesh, "face", POINTER, MFACE);

    idf = 0;
    while( (ptr_f = MESH_Next_Face(trimesh, &idf)) ) {
      MF_Coords(ptr_f, &n, xyz_lst);
      if( n != 3 ) {
	MSTK_Report("MESH_Tri2Quad","Input mesh is not a triangle mesh",MSTK_WARN);
	return 0;
      }

      for( i=0; i<3; i++ ) xyz[i] = (xyz_lst[0][i] + xyz_lst[1][i] + xyz_lst[2][i]) / 3;

      ptr_v_new = MV_New(*quadmesh);
      MV_Set_Coords(ptr_v_new, xyz);
      MV_Set_GEntDim(ptr_v_new,MF_GEntDim(ptr_f));
      MV_Set_GEntID(ptr_v_new,MF_GEntID(ptr_f));


      MEnt_Set_AttVal(ptr_f, atr_f, 0, 0, (void*)ptr_v_new);
    }



    /* convert tris to quads */
    idf = 0;
    while( (ptr_f = MESH_Next_Face(trimesh, &idf)) ) {
      lst_v = MF_Vertices(ptr_f, 1, 0);
      n = List_Num_Entries(lst_v);
      lst_e = MF_Edges(ptr_f, 1, List_Entry(lst_v, 0));

      for( i=0; i<3; i++ ) {
	old_ptr_v_face[i] = List_Entry(lst_v, i);
	MEnt_Get_AttVal(old_ptr_v_face[i], atr_v, 0, 0, &(ptr_v_face[i]));
      }
      
      for( i=0;   i<3; i++ ) {
	ptr_e = List_Entry(lst_e, i);
	
	MEnt_Get_AttVal(ptr_e, atr_e, &ival, &dval, &vval);
	ptr_v_face[n] = (MVertex_ptr)(vval);
	MV_Coords(ptr_v_face[n], xyz);
	n++;
      }
      
      MEnt_Get_AttVal(ptr_f, atr_f, &ival, &dval, &vval);
      ptr_v_face[n] = (MVertex_ptr)(vval);
      n++;
      
      for( k=0; k<3; k++ ) {
        for( i=0; i<4; i++ ) ptr_v_quad[i] = ptr_v_face[tri_quad[k][i]];

	ptr_f_new = MF_New(*quadmesh);

        MF_Set_Vertices(ptr_f_new, 4, ptr_v_quad);

	MF_Set_GEntDim(ptr_f_new, MF_GEntDim(ptr_f));
	MF_Set_GEntID(ptr_f_new, MF_GEntID(ptr_f));
      }

      List_Delete(lst_v);
      List_Delete(lst_e);
    }



    /* remove attributes on edges of tri mesh */
    idv = 0;
    while( (ptr_v = MESH_Next_Vertex(trimesh, &idv)) ) {
      MEnt_Rem_AttVal(ptr_v, atr_v);
    }
    MAttrib_Delete(atr_v);

    ide = 0;
    while( (ptr_e = MESH_Next_Edge(trimesh, &ide)) ) {
      MEnt_Rem_AttVal(ptr_e, atr_e);
    }
    MAttrib_Delete(atr_e);

  
    /* remove attributes on faces of tri mesh */
    idf = 0;
    while( (ptr_f = MESH_Next_Face(trimesh, &idf)) ) {
      MEnt_Rem_AttVal(ptr_f, atr_f);
    }
    MAttrib_Delete(atr_f);

    return 1;
  }

#ifdef _cplusplus
}
#endif
