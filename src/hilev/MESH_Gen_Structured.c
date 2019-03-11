/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Mesh.h"
#include "MSTK.h"
#include "MSTK_private.h"


/* Generate 2D and 3D structured meshes */
/* lx, ly, lz - lower limits in X, Y and Z directions */
/* ux, uy, uz - upper limits in X, Y and Z directions */
/* nx, ny, nz - number of cells in X, Y and Z directions */
/* SPECIFY nz=0 FOR 2D MESHES */

/* Mesh entities will have correct classification against a square or
   cubical geometric domain */

/*

  Index directions for classification templates

  k   j
  |  /
  | /
  |/___ i


  Model vertex, edge and face enumeration for classification templates 


         MODEL                   MODEL                  MODEL
         VERTICES                EDGES                  FACES

     7 ____________ 8          ______7_____           ____________  
      /|          /|          /|          /|         /|      2   /| 
     / |         / |       12/ |8      11/ | 	    / |  4      / | 
   5/___________/6 |        /_____3_____/  |6	   /___________/  | 
    |  |        |  |        |  |        |  | 	   |  |        | 5| 
    |  |________|__|        |  |_____5__|__| 	   |6 |_1______|__| 
    |  /3       |  /4      4|  /        |  / 	   |  /        |  / 
    | /         | /         | /9       2| /10	   | /      3  | /  
    |/__________|/          |/__________|/   	   |/__________|/   
   1             2                1
                                                   
                                                    Front  - F1
						    Back   - F2
						    Bottom - F3
						    Top    - F4
						    Left   - F6
						    Right  - F5

*/

Mesh_ptr MESH_Gen_Structured(double lx, double ly, double lz,
			     double ux, double uy, double uz,
			     int nx, int ny, int nz) {
  int mesh_type;
  int i, j, k, ii, jj, kk, gid, gdim;
  double xyz[3], dx, dy, dz;
  MVertex_ptr mv, rverts[8], fverts[4], everts[2];
  MEdge_ptr me;
  MFace_ptr mf;
  MRegion_ptr mr;
  Mesh_ptr mesh;
  char outfile[256];
  int vgid_tmpl[3][3][3] = {{{1,4,5},{9,6,12},{3,8,7}},{{1,1,3},{3,1,4},{5,2,7}},{{2,2,6},{10,5,11},{4,6,8}}};
  int vgdim_tmpl[3][3][3]= {{{0,1,0},{1,2,1}, {0,1,0}},{{1,2,1},{2,3,2},{1,2,1}},{{0,1,0},{1,2,1},{0,1,0}}};
  int egdim_tmpl[3][3] = {{1,2,1},{2,3,2},{1,2,1}};
  int egid_tmpl2[3][3] = {{4,6,8},{1,1,2},{2,5,6}};  /* Y direction edges (iterating over i,k) */
  int egid_tmpl1[3][3] = {{9,6,12},{3,1,4},{10,5,11}}; /* Z direction edges (iterating over i,j)*/
  int egid_tmpl0[3][3] = {{1,1,3},{3,1,4},{5,2,7}}; /* X direction edges (iterating over j,k) */
  int fgdim_tmpl[3] = {2,3,2};
  int fgid_tmpl0[3] = {6,1,5};
  int fgid_tmpl1[3] = {1,1,2};
  int fgid_tmpl2[3] = {3,1,4};


  MSTK_Init();
  mesh = MESH_New(F1);

  mesh_type = (nz > 0) ? 3 : 2; /* solid */
  if (ny == 0)
    MSTK_Report("MESH_Gen_Structured",
		"Cannot generate wire meshes", MSTK_FATAL);
  
  if (mesh_type == 2) {

    /* 2D Mesh */

    MVertex_ptr **verts;

    dx = (ux-lx)/nx;
    dy = (uy-ly)/ny;

    verts = (MVertex_ptr **) malloc((nx+1)*sizeof(MVertex_ptr *));
    for (i = 0; i < nx+1; i++)
      verts[i] = (MVertex_ptr *) malloc((ny+1)*sizeof(MVertex_ptr));
 
    xyz[2] = 0.0;
    for (j = 0; j < ny+1; j++) {
      xyz[1] = (j == ny) ? uy : ly + j*dy;

      for (i = 0; i < nx+1; i++) {
        xyz[0] = (i == nx) ? ux : lx + i*dx;

        mv = MV_New(mesh);
        MV_Set_Coords(mv,xyz);

        if (i == 0) {
          if (j == 0) {
            MV_Set_GEntDim(mv,0);
            MV_Set_GEntID(mv,1);
          }
          else if (j == ny) {
            MV_Set_GEntDim(mv,0);
            MV_Set_GEntID(mv,4);	  
          }
          else {
            MV_Set_GEntDim(mv,1);
            MV_Set_GEntID(mv,4);
          }
        }
        else if (i == nx) {
          if (j == 0) {
            MV_Set_GEntDim(mv,0);
            MV_Set_GEntID(mv,2);
          }
          else if (j == ny) {
            MV_Set_GEntDim(mv,0);
            MV_Set_GEntID(mv,3);	  
          }
          else {
            MV_Set_GEntDim(mv,1);
            MV_Set_GEntID(mv,2);
          }
        }
        else {
          if (j == 0) {
            MV_Set_GEntDim(mv,1);
            MV_Set_GEntID(mv,1);
          }
          else if (j == ny) {
            MV_Set_GEntDim(mv,1);
            MV_Set_GEntID(mv,3);
          }
          else {
            MV_Set_GEntDim(mv,2);
            MV_Set_GEntID(mv,1);
          }
        }

        verts[i][j] = mv;
      }
    }


    for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {

        MVertex_ptr v0, v1;
        List_ptr fedges[4];
        int dir[4];

        mf = MF_New(mesh);
      
        /* edge 0 */
        v0 = verts[i][j];
        v1 = verts[i+1][j];
        fedges[0] = MVs_CommonEdge(v0,v1);
        if (fedges[0])
          dir[0] = (ME_Vertex(fedges[0],0) == v0) ? 1 : 0;
        else {
          me = ME_New(mesh);
	
          ME_Set_Vertex(me,0,v0);
          ME_Set_Vertex(me,1,v1);
	
          if (j == 0) {
            ME_Set_GEntDim(me,1);
            ME_Set_GEntID(me,1);
          }
          else {
            ME_Set_GEntDim(me,2);
            ME_Set_GEntID(me,1);
          }
	
          fedges[0] = me;
          dir[0] = 1;
        }
      
      
        /* edge 1 */
        v0 = verts[i+1][j];
        v1 = verts[i+1][j+1];
        fedges[1] = MVs_CommonEdge(v0,v1);
        if (fedges[1])
          dir[1] = (ME_Vertex(fedges[1],0) == v0) ? 1 : 0;
        else {
          me = ME_New(mesh);
	
          ME_Set_Vertex(me,0,v0);
          ME_Set_Vertex(me,1,v1);
	
          if (i+1 == nx) {
            ME_Set_GEntDim(me,1);
            ME_Set_GEntID(me,2);
          }
          else {
            ME_Set_GEntDim(me,2);
            ME_Set_GEntID(me,1);
          }
	
          fedges[1] = me;
          dir[1] = 1;
        }
      
      
        /* edge 2 */
        v0 = verts[i+1][j+1];
        v1 = verts[i][j+1];
        fedges[2] = MVs_CommonEdge(v0,v1);
        if (fedges[2])
          dir[2] = (ME_Vertex(fedges[2],0) == v0) ? 1 : 0;
        else {
          me = ME_New(mesh);
	
          ME_Set_Vertex(me,0,v0);
          ME_Set_Vertex(me,1,v1);
	
          if (j+1 == ny) {
            ME_Set_GEntDim(me,1);
            ME_Set_GEntID(me,3);
          }
          else {
            ME_Set_GEntDim(me,2);
            ME_Set_GEntID(me,1);
          }
	
          fedges[2] = me;
          dir[2] = 1;
        }
      
      
        /* edge 3 */
        v0 = verts[i][j+1];
        v1 = verts[i][j];
        fedges[3] = MVs_CommonEdge(v0,v1);
        if (fedges[3])
          dir[3] = (ME_Vertex(fedges[3],0) == v0) ? 1 : 0;
        else {
          me = ME_New(mesh);
	
          ME_Set_Vertex(me,0,v0);
          ME_Set_Vertex(me,1,v1);
	
          if (i == 0) {
            ME_Set_GEntDim(me,1);
            ME_Set_GEntID(me,4);
          }
          else {
            ME_Set_GEntDim(me,2);
            ME_Set_GEntID(me,1);
          }
	
          fedges[3] = me;
          dir[3] = 1;
        }


        MF_Set_Edges(mf,4,fedges,dir);

        MF_Set_GEntDim(mf,2);
        MF_Set_GEntID(mf,1);
      }
    }
   

    for (i = 0; i < nx; i++)
      free(verts[i]);
    free(verts);

  }
  else {

    /* 3D Mesh */

    MVertex_ptr ***verts;

    dx = (ux-lx)/nx;
    dy = (uy-ly)/ny;
    dz = (uz-lz)/nz;

    verts = (MVertex_ptr ***) malloc((nx+1)*sizeof(MVertex_ptr **));
    for (j = 0; j < nx+1; j++) {
      verts[j] = (MVertex_ptr **) malloc((ny+1)*sizeof(MVertex_ptr *)); 
      for (k = 0; k < ny+1; k++) 
        verts[j][k] = (MVertex_ptr *) malloc((nz+1)*sizeof(MVertex_ptr));
    }

    for (k = 0; k < nz+1; k++) {
      xyz[2] = (k == nz) ? uz : lz + k*dz;
      kk =  (k%nz) ? 1 : (k ? 2 : 0);

      for (j = 0; j < ny+1; j++) {
        xyz[1] = (j == ny) ? uy : ly + j*dy;      
        jj = (j%ny) ? 1 : (j ? 2 : 0);

        for (i = 0; i < nx+1; i++) {
          xyz[0] = (i == nx) ? ux : lx + i*dx;
          ii = (i%nx) ? 1 : (i ? 2 : 0);
	
          mv = MV_New(mesh);
          MV_Set_Coords(mv,xyz);	
          verts[i][j][k] = mv;

          gdim  = vgdim_tmpl[ii][jj][kk];
          MV_Set_GEntDim(mv,gdim);

          gid = vgid_tmpl[ii][jj][kk];
          MV_Set_GEntID(mv,gid);
        }
      }
    }


    /* Create the edges explicitly to get the classification right */
    for (i = 0; i < nx+1; i++) {
      for (j = 0; j < ny+1; j++) {
        for (k = 0; k < nz; k++) {
          me = ME_New(mesh);

          everts[0] = verts[i][j][k];
          everts[1] = verts[i][j][k+1];
          ME_Set_Vertex(me,0,everts[0]);
          ME_Set_Vertex(me,1,everts[1]);

          ii = (i%nx) ? 1 : (i ? 2 : 0);
          jj = (j%ny) ? 1 : (j ? 2 : 0);
          gdim = egdim_tmpl[ii][jj];
          gid = egid_tmpl2[ii][jj];

          ME_Set_GEntDim(me,gdim);
          ME_Set_GEntID(me,gid);
        }
      }
    }
	
    for (i = 0; i < nx+1; i++) {
      for (k = 0; k < nz+1; k++) {
        for (j = 0; j < ny; j++) {
          me = ME_New(mesh);

          everts[0] = verts[i][j][k];
          everts[1] = verts[i][j+1][k];
          ME_Set_Vertex(me,0,everts[0]);
          ME_Set_Vertex(me,1,everts[1]);

          ii = (i%nx) ? 1 : (i ? 2 : 0);
          kk = (k%nz) ? 1 : (k ? 2 : 0);
          gdim = egdim_tmpl[ii][kk];
          gid = egid_tmpl1[ii][kk];

          ME_Set_GEntDim(me,gdim);
          ME_Set_GEntID(me,gid);
        }
      }
    }
	
    for (j = 0; j < ny+1; j++) {
      for (k = 0; k < nz+1; k++) {
        for (i = 0; i < nx; i++) {
          me = ME_New(mesh);

          everts[0] = verts[i][j][k];
          everts[1] = verts[i+1][j][k];
          ME_Set_Vertex(me,0,everts[0]);
          ME_Set_Vertex(me,1,everts[1]);

          jj = (j%ny) ? 1 : (j ? 2 : 0);
          kk = (k%nz) ? 1 : (k ? 2 : 0);
          gdim = egdim_tmpl[jj][kk];
          gid = egid_tmpl0[jj][kk];

          ME_Set_GEntDim(me,gdim);
          ME_Set_GEntID(me,gid);
        }
      }
    }
	


    /* Create the faces explicitly to get the classification right */
    for (i = 0; i < nx+1; i++) {
      for (j = 0; j < ny; j++) {
        for (k = 0; k < nz; k++) {
          mf = MF_New(mesh);

          fverts[0] = verts[i][j][k];
          fverts[1] = verts[i][j+1][k];
          fverts[2] = verts[i][j+1][k+1];
          fverts[3] = verts[i][j][k+1];
          MF_Set_Vertices(mf,4,fverts);

          ii = (i%nx) ? 1 : (i ? 2 : 0);
          gdim = fgdim_tmpl[ii];
          gid = fgid_tmpl0[ii];

          MF_Set_GEntDim(mf,gdim);
          MF_Set_GEntID(mf,gid);
        }
      }
    }
	
    for (j = 0; j < ny+1; j++) {
      for (i = 0; i < nx; i++) {
        for (k = 0; k < nz; k++) {
          mf = MF_New(mesh);

          fverts[0] = verts[i][j][k];
          fverts[1] = verts[i+1][j][k];
          fverts[2] = verts[i+1][j][k+1];
          fverts[3] = verts[i][j][k+1];
          MF_Set_Vertices(mf,4,fverts);

          jj = (j%ny) ? 1 : (j ? 2 : 0);
          gdim = fgdim_tmpl[jj];
          gid = fgid_tmpl1[jj];

          MF_Set_GEntDim(mf,gdim);
          MF_Set_GEntID(mf,gid);
        }
      }
    }
	
    for (k = 0; k < nz+1; k++) {
      for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
          mf = MF_New(mesh);

          fverts[0] = verts[i][j][k];
          fverts[1] = verts[i+1][j][k];
          fverts[2] = verts[i+1][j+1][k];
          fverts[3] = verts[i][j+1][k];
          MF_Set_Vertices(mf,4,fverts);

          kk = (k%nz) ? 1 : (k ? 2 : 0);
          gdim = fgdim_tmpl[kk];
          gid = fgid_tmpl2[kk];

          MF_Set_GEntDim(mf,gdim);
          MF_Set_GEntID(mf,gid);
        }
      }
    }
	

    /* Not the most efficient way but the easiest to code */

    for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
        for (k = 0; k < nz; k++) {
          mr = MR_New(mesh);
          MR_Set_GEntID(mr,1);
	
          rverts[0] = verts[i][j][k];       rverts[1] = verts[i+1][j][k]; 
          rverts[2] = verts[i+1][j+1][k];   rverts[3] = verts[i][j+1][k];
          rverts[4] = verts[i][j][k+1];     rverts[5] = verts[i+1][j][k+1]; 
          rverts[6] = verts[i+1][j+1][k+1]; rverts[7] = verts[i][j+1][k+1];

          MR_Set_Vertices(mr, 8, rverts, 6, NULL);
        }
      }
    }

    for (i = 0; i < nx+1; i++) {
      for (j = 0; j < ny+1; j++)
        free(verts[i][j]);
      free(verts[i]);
    }
    free(verts);

  }

  return mesh;
  
}  /* MESH_Gen_Structured */
  
