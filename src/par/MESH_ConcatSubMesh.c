/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#define _H_Mesh_Private

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "Mesh.h"
#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif



  /* 
     This function concatenates submesh into mesh, based on global ID

     NOTE: DESIGNED TO ADD ONLY ONE LAYER OF ELEMENTS FROM THE
     SUBMESHES - IT WOULD BE GOOD TO GENERALIZE IT (03/07/2019)

     Author(s): Duo Wang, Rao Garimella
  */

  int MESH_ConcatSubMesh_Face(Mesh_ptr mesh, int num, Mesh_ptr *submeshes);
  int MESH_ConcatSubMesh_Region(Mesh_ptr mesh, int num, Mesh_ptr *submeshes);


  int MESH_ConcatSubMesh(Mesh_ptr mesh, int topodim, int num, Mesh_ptr *submeshes) {
    
    if (topodim == 3)
      MESH_ConcatSubMesh_Region(mesh, num, submeshes);
    else if (topodim == 2) 
      MESH_ConcatSubMesh_Face(mesh, num, submeshes);
    else {
      MSTK_Report("MESH_ConcatSubMesh()","only send volume or surface mesh",MSTK_ERROR);
      exit(-1);
    }
    return 1;
  }

  int MESH_ConcatSubMesh_Face(Mesh_ptr mesh, int num, Mesh_ptr *submeshes) {
    int nfv, nfe, i, j, k, ival;
    MVertex_ptr mv, new_mv, sub_mv;
    MEdge_ptr me, new_me, sub_me;
    MFace_ptr new_mf, sub_mf;
    List_ptr mfverts, mfedges;
    int add_face, idx, global_id, iloc, *loc;
    double coor[3], rval;
    void *pval;
    Mesh_ptr submesh;

    List_ptr parbndry_verts = List_New(10);
    List_ptr parbndry_edges = List_New(10);

    MEdge_ptr *fedges = (MEdge_ptr *) malloc(MAXPV2*sizeof(MEdge_ptr));
    int *fedirs = (int *) malloc(MAXPV2*sizeof(int));

    MAttrib_ptr parbndryatt = MAttrib_New(mesh, "on_parbndry", INT, MVERTEX);
    
    /* collect edges and vertices on the partition boundary */
    int num_parbndry_edges = 0;
    idx = 0;
    while ((me = MESH_Next_Edge(mesh,&idx))) 
      if (ME_PType(me) != PINTERIOR) {
        List_Add(parbndry_edges,me);
        num_parbndry_edges++;
      }
    int num_parbndry_verts = 0;
    idx = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx)))
      if (MV_PType(mv) != PINTERIOR) {
        List_Add(parbndry_verts,mv);
        MEnt_Set_AttVal(mv, parbndryatt, 1, 0.0, NULL);
        num_parbndry_verts++;
      }
    /* sort based on global ID */
    List_Sort(parbndry_edges,num_parbndry_edges,sizeof(MEdge_ptr),compareGlobalID);
    List_Sort(parbndry_verts,num_parbndry_verts,sizeof(MVertex_ptr),compareGlobalID);

    int *parbndry_vert_gids = (int *) malloc(num_parbndry_verts*sizeof(int));
    int *parbndry_edge_gids = (int *) malloc(num_parbndry_edges*sizeof(int));

    /* store them in array for binary search */
    for (i = 0; i < num_parbndry_edges; i++) {
      me = List_Entry(parbndry_edges,i);
      parbndry_edge_gids[i] = ME_GlobalID(me);
    }
    for (i = 0; i < num_parbndry_verts; i++) {
      mv = List_Entry(parbndry_verts,i);
      parbndry_vert_gids[i] = MV_GlobalID(mv);
    }

    
    /* Make list of new edges and vertices which will be updated
       with each mesh that is concatenated */
    int max_vnew = 0, max_enew = 0;
    for (i = 0; i < num; i++) {
      max_vnew += MESH_Num_Vertices(submeshes[i]);
      max_enew += MESH_Num_Edges(submeshes[i]);
    }

    int num_new_verts = 0, num_new_edges = 0; 
    int *new_vert_gids = (int *) malloc(max_vnew*sizeof(int));
    int *new_edge_gids = (int *) malloc(max_enew*sizeof(int));

    List_ptr new_verts = List_New(max_vnew);
    List_ptr new_edges = List_New(max_enew);


    /* Now process each mesh and add a layer of ghost elements from
       each of them to the main partition */
    
    for (i = 0; i < num; i++) {
      submesh = submeshes[i];
    
      MAttrib_ptr vidatt = MAttrib_New(submesh, "tempvid", POINTER, MVERTEX);
      MAttrib_ptr eidatt = MAttrib_New(submesh, "tempeid", POINTER, MEDGE);

      idx = 0;
      while ((sub_mf = MESH_Next_Face(submesh, &idx))) {
        add_face = 0;
      
        /* Find matching vertices between the submesh and main mesh */
      
        mfverts = MF_Vertices(sub_mf,1,0);
        nfv = List_Num_Entries(mfverts);
        for (j = 0; j < nfv; j++) {
          sub_mv = List_Entry(mfverts,j);
        
          /* Does the vertex have a known counterpart on the partition
           * boundary of the main mesh? */
          MEnt_Get_AttVal(sub_mv, vidatt, &ival, &rval, &mv);

          if (mv) {
            int on_parbndry=0;
            MEnt_Get_AttVal(mv, parbndryatt, &on_parbndry, &rval, &pval);
            if (on_parbndry)
              add_face = 1; 
          } else {
        
            /* Does the global ID of this vertex of the sub mesh face
             * match the global ID of a partition boundary vertex in
             * the main mesh? */
            
            global_id = MV_GlobalID(sub_mv);
            loc = (int *) bsearch(&global_id, parbndry_vert_gids, num_parbndry_verts, sizeof(int),
                                  compareINT);
            if (loc) {  /* found a match */
              add_face = 1; 
              iloc = loc - parbndry_vert_gids;
              mv = List_Entry(parbndry_verts,iloc); 
              /* here set the ghost vertex property, only necessary when the input submeshes are not consistent */
              if (MV_PType(mv) == PGHOST && MV_PType(sub_mv) != PGHOST) {
                MV_Set_GEntDim(mv,MV_GEntDim(sub_mv));
                MV_Set_GEntID(mv,MV_GEntID(sub_mv));
              }
              
              MEnt_Set_AttVal(sub_mv, vidatt, 0, 0.0, mv);
            }
          }
        }
        List_Delete(mfverts);
      
        /* Find matching edges between the submesh and main mesh */
      
        mfedges = MF_Edges(sub_mf,1,0);
        nfe = List_Num_Entries(mfedges);
        for (j = 0; j < nfe; j++) {
          sub_me = List_Entry(mfedges,j);

          /* Does the edge have a known counterpart on the partition
           * boundary of the main mesh */
          MEnt_Get_AttVal(sub_me, eidatt, &ival, &rval, &me);

          if (!me) {
            /* Does the global ID of this edge of the sub mesh face
             * match the global ID of a partition boundary edge in the
             * main mesh? */
            
            global_id = ME_GlobalID(sub_me);
            loc = (int *) bsearch(&global_id, parbndry_edge_gids, num_parbndry_edges, sizeof(int),
                                  compareINT);
            if (loc) {
              iloc = loc - parbndry_edge_gids;
              me = List_Entry(parbndry_edges,iloc); 
              /* here set the ghost edge property, only necessary when the input submeshes are not consistent */
              if (ME_PType(me) == PGHOST && ME_PType(sub_me) != PGHOST) {
                ME_Set_GEntDim(me,ME_GEntDim(sub_me));
                ME_Set_GEntID(me,ME_GEntID(sub_me));
              }
              
	      MEnt_Set_AttVal(sub_me, eidatt, 0, 0.0, me);
            }
          }
        }
        
        if (!add_face) {
          List_Delete(mfedges);
          continue;
        }

        new_mf = MF_New(mesh); /* add face */
        MF_Set_GEntDim(new_mf,MF_GEntDim(sub_mf));
        MF_Set_GEntID(new_mf,MF_GEntID(sub_mf));
        MF_Set_PType(new_mf,PGHOST);
        MF_Set_MasterParID(new_mf,MF_MasterParID(sub_mf));
        MF_Set_GlobalID(new_mf,MF_GlobalID(sub_mf));
      
        nfe = List_Num_Entries(mfedges);
        for (j = 0; j < nfe; j++) {
          sub_me = List_Entry(mfedges,j);
          global_id = ME_GlobalID(sub_me);
          fedirs[j] = MF_EdgeDir_i(sub_mf,j) == 1 ? 1 : 0;

          new_me = NULL;	  
	  MEnt_Get_AttVal(sub_me, eidatt, &ival, &rval, &new_me);

          if (!new_me) {
            /* search in the ghost layer if another edge with
             * this global ID has been added */
            loc = (int *) bsearch(&global_id, new_edge_gids, num_new_edges,
                                  sizeof(int), compareINT);
            if (loc) {
              iloc = loc - new_edge_gids;
              new_me = List_Entry(new_edges, iloc);
              MEnt_Set_AttVal(sub_me, eidatt, 0, 0.0, new_me);
            }
          }

          if (new_me) {
            if (MV_GlobalID(ME_Vertex(new_me,0)) != MV_GlobalID(ME_Vertex(sub_me,0)))
              fedirs[j] = 1 - fedirs[j];  /* if the edge dir is not the same, reverse the edge dir */
          } else  {  /* add a new edge to main mesh */
            
            new_me = ME_New(mesh);
            ME_Set_GEntDim(new_me,ME_GEntDim(sub_me));
            ME_Set_GEntID(new_me,ME_GEntID(sub_me));
            ME_Set_PType(new_me,PGHOST);
            ME_Set_MasterParID(new_me,ME_MasterParID(sub_me));
            ME_Set_GlobalID(new_me,ME_GlobalID(sub_me));
	  
            MEnt_Set_AttVal(sub_me, eidatt, 0, 0.0, new_me);
	    List_Add(new_edges, new_me);
          
            for (k = 0; k < 2; k++) {
              sub_mv = ME_Vertex(sub_me,k);
              global_id = MV_GlobalID(sub_mv);

              new_mv = NULL;
              MEnt_Get_AttVal(sub_mv, vidatt, &ival, &rval, &new_mv);
	      if (!new_mv) {
		/* search in the ghost layer if another vertex with
                 * this global ID has been added */
                loc = (int *) bsearch(&global_id, new_vert_gids, num_new_verts,
                                      sizeof(int), compareINT);
                if (loc) {
                  iloc = loc - new_vert_gids;
                  new_mv = List_Entry(new_verts, iloc);
                  MEnt_Set_AttVal(sub_mv, vidatt, 0, 0.0, new_mv);
                }
              }

              if (!new_mv) {  /* add a new vertex to main mesh */
                new_mv = MV_New(mesh);
                MV_Set_GEntDim(new_mv,MV_GEntDim(sub_mv));
                MV_Set_GEntID(new_mv,MV_GEntID(sub_mv));
                MV_Set_PType(new_mv,PGHOST);
                MV_Set_MasterParID(new_mv,MV_MasterParID(sub_mv));
                MV_Set_GlobalID(new_mv,MV_GlobalID(sub_mv));
                MV_Coords(sub_mv,coor);
                MV_Set_Coords(new_mv,coor);
	      
                MEnt_Set_AttVal(sub_mv, vidatt, 0, 0.0, new_mv);
		List_Add(new_verts, new_mv);
              }
              ME_Set_Vertex(new_me,k,new_mv);  /* set edge-vertex */
            }
          }								
          fedges[j] = new_me;
        }
        MF_Set_Edges(new_mf,nfe,fedges,fedirs); /* set face-edge */

        List_Delete(mfedges);
      }

      idx = 0;
      while ((sub_mv = MESH_Next_Vertex(submesh, &idx)))
	MEnt_Rem_AttVal(sub_mv, vidatt);
      MAttrib_Delete(vidatt);
      idx = 0;
      while ((sub_me = MESH_Next_Edge(submesh, &idx)))
	MEnt_Rem_AttVal(sub_me, eidatt);
      MAttrib_Delete(eidatt);

      /* Sort the added entity lists by GlobalID */
      num_new_edges = List_Num_Entries(new_edges);
      List_Sort(new_edges, num_new_edges, sizeof(MEdge_ptr), compareGlobalID);
      for (j = 0; j < num_new_edges; j++)
        new_edge_gids[j] = ME_GlobalID(List_Entry(new_edges, j));

      num_new_verts = List_Num_Entries(new_verts);
      List_Sort(new_verts, num_new_verts, sizeof(MVertex_ptr), compareGlobalID);
      for (j = 0; j < num_new_verts; j++)
        new_vert_gids[j] = MV_GlobalID(List_Entry(new_verts, j));
    }

    idx = 0;
    while ((mv = List_Next_Entry(parbndry_verts, &idx)))
      MEnt_Rem_AttVal(mv, parbndryatt);
    MAttrib_Delete(parbndryatt);
    
    List_Delete(parbndry_edges);
    List_Delete(parbndry_verts);
    List_Delete(new_edges);
    List_Delete(new_verts);

    free(parbndry_vert_gids);
    free(parbndry_edge_gids);
    free(new_vert_gids);
    free(new_edge_gids);

    free(fedges);
    free(fedirs);

    return 1;
  }

  /* right now assume there are no overlapped regions */

  int MESH_ConcatSubMesh_Region(Mesh_ptr mesh, int num, Mesh_ptr *submeshes) {
    int nrf, nre, nrv, nfe, i, j, k, num_parbndry_verts, num_parbndry_edges, num_parbndry_faces, ival;
    MVertex_ptr mv, new_mv, sub_mv;
    MEdge_ptr me, new_me, sub_me;
    MFace_ptr mf, new_mf, sub_mf;
    MRegion_ptr new_mr, sub_mr;
    List_ptr mrfaces, mredges, mrverts, mfedges;
    int add_region, idx, global_id, iloc, *loc;
    double coor[3], rval;
    void *pval;
    Mesh_ptr submesh;

    List_ptr parbndry_verts = List_New(10);        
    List_ptr parbndry_edges = List_New(10);
    List_ptr parbndry_faces = List_New(10);

    MFace_ptr *rfaces = (MFace_ptr *) malloc(MAXPF3*sizeof(MFace_ptr));
    int *rfdirs = (int *) malloc(MAXPF3*sizeof(int));
    MEdge_ptr *fedges = (MEdge_ptr *) malloc(MAXPV2*sizeof(MEdge_ptr));
    int *fedirs = (int *) malloc(MAXPV2*sizeof(int));

    MAttrib_ptr parbndryatt = MAttrib_New(mesh, "on_parbndry", INT, MVERTEX);
    
    /* collect faces, edges and vertices on the partition boundary */
    idx = 0; num_parbndry_faces = 0;
    while ((mf = MESH_Next_Face(mesh,&idx))) 
      if (MF_PType(mf) != PINTERIOR) {
        List_Add(parbndry_faces,mf);
        num_parbndry_faces++;
      }
    idx = 0; num_parbndry_edges = 0;
    while ((me = MESH_Next_Edge(mesh,&idx))) 
      if (ME_PType(me) != PINTERIOR) {
        List_Add(parbndry_edges,me);
        num_parbndry_edges++;
      }
    idx = 0; num_parbndry_verts = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx)))
      if (MV_PType(mv) != PINTERIOR) {
        List_Add(parbndry_verts,mv);
        MEnt_Set_AttVal(mv, parbndryatt, 1, 0.0, NULL);
        num_parbndry_verts++;
      }
    
    /* sort based on global ID */
    List_Sort(parbndry_faces,num_parbndry_faces,sizeof(MFace_ptr),compareGlobalID);
    List_Sort(parbndry_edges,num_parbndry_edges,sizeof(MEdge_ptr),compareGlobalID);
    List_Sort(parbndry_verts,num_parbndry_verts,sizeof(MVertex_ptr),compareGlobalID);

    int *parbndry_vert_gids = (int *)malloc(num_parbndry_verts*sizeof(int));
    int *parbndry_edge_gids = (int *)malloc(num_parbndry_edges*sizeof(int));
    int *parbndry_face_gids = (int *)malloc(num_parbndry_faces*sizeof(int));

    /* store them in array for binary search */
    for (i = 0; i < num_parbndry_faces; i++) {
      mf = List_Entry(parbndry_faces,i);
      parbndry_face_gids[i] = MF_GlobalID(mf);
    }
    for (i = 0; i < num_parbndry_edges; i++) {
      me = List_Entry(parbndry_edges,i);
      parbndry_edge_gids[i] = ME_GlobalID(me);
    }
    for (i = 0; i < num_parbndry_verts; i++) {
      mv = List_Entry(parbndry_verts,i);
      parbndry_vert_gids[i] = MV_GlobalID(mv);
    }

   /* Make list of new edges and vertices which will be updated
       with each mesh that is concatenated */
    int max_vnew = 0, max_enew = 0, max_fnew = 0;
    for (i = 0; i < num; i++) {
      max_vnew += MESH_Num_Vertices(submeshes[i]);
      max_enew += MESH_Num_Edges(submeshes[i]);
      max_fnew += MESH_Num_Faces(submeshes[i]);
    }

    int num_new_verts = 0, num_new_edges = 0, num_new_faces = 0; 
    int *new_vert_gids = (int *) malloc(max_vnew*sizeof(int));
    int *new_edge_gids = (int *) malloc(max_enew*sizeof(int));
    int *new_face_gids = (int *) malloc(max_fnew*sizeof(int));

    List_ptr new_verts = List_New(max_vnew);
    List_ptr new_edges = List_New(max_enew);
    List_ptr new_faces = List_New(max_fnew);
  

    /* Now process each mesh and add a layer of ghost elements from
       each of them to the main partition */
       
    for (i = 0; i < num; i++) {
      submesh = submeshes[i];

      MAttrib_ptr vidatt = MAttrib_New(submesh, "tempvid", POINTER, MVERTEX);
      MAttrib_ptr eidatt = MAttrib_New(submesh, "tempeid", POINTER, MEDGE);
      MAttrib_ptr fidatt = MAttrib_New(submesh, "tempfid", POINTER, MFACE);

      idx = 0;
      while ((sub_mr = MESH_Next_Region(submesh, &idx))) {
        add_region = 0;

        /* Find matching vertices between submesh and main mesh */

        mrverts = MR_Vertices(sub_mr);
        nrv = List_Num_Entries(mrverts);
        for (j = 0; j < nrv; j++) {
          sub_mv = List_Entry(mrverts,j);

          MEnt_Get_AttVal(sub_mv, &vidatt, &ival, &rval, &mv);

          if (mv) {
            int on_parbndry=0;
            MEnt_Get_AttVal(mv, &parbndryatt, &on_parbndry, &rval, &pval);
            if (on_parbndry)
              add_region = 1; 
          } else {

            /* Does the global ID of this vertex of the sub mesh region
             * match the global ID of a boundary vertex in the main
             * mesh? */
            
            global_id = MV_GlobalID(sub_mv);
            loc = (int *) bsearch(&global_id, parbndry_vert_gids, num_parbndry_verts, sizeof(int),
                                  compareINT);
            if (loc) {
              add_region = 1; 
              iloc = loc - parbndry_vert_gids;
              mv = List_Entry(parbndry_verts,iloc); 
              /* here set the ghost vertex property, only necessary when the input submeshes are not consistent */
              if(MV_PType(mv) == PGHOST && MV_PType(sub_mv) != PGHOST) {
                MV_Set_GEntDim(mv,MV_GEntDim(sub_mv));
                MV_Set_GEntID(mv,MV_GEntID(sub_mv));
              }
              
              MEnt_Set_AttVal(sub_mv, vidatt, 0, 0.0, mv);
            }
          }
        }
        List_Delete(mrverts);

        /* Find matching edges between submesh and main mesh */

        mredges = MR_Edges(sub_mr);
        nre = List_Num_Entries(mredges);
        for (j = 0; j < nre; j++) {
          sub_me = List_Entry(mredges,j);
          
          /* Does the edge already have a counterpart in the main mesh? */
          MEnt_Get_AttVal(sub_me, eidatt, &ival, &rval, &me);

          if (!me) {
            /* Does the global ID of this edge of the sub mesh region
             * match the global ID of a boundary edge in the main
             * mesh? */
            
            global_id = ME_GlobalID(sub_me);
            loc = (int *) bsearch(&global_id, parbndry_edge_gids, num_parbndry_edges, sizeof(int),
                                  compareINT);
            if (loc) {
              add_region = 1; 
              iloc = loc - parbndry_edge_gids;
              me = List_Entry(parbndry_edges,iloc); 
              /* here set the ghost edge property, only necessary when the input submeshes are not consistent */
              if(ME_PType(me) == PGHOST && ME_PType(sub_me) != PGHOST) {
                ME_Set_GEntDim(me,ME_GEntDim(sub_me));
                ME_Set_GEntID(me,ME_GEntID(sub_me));
              }

              MEnt_Set_AttVal(sub_me, eidatt, 0, 0.0, me);
            }
          }
        }
        List_Delete(mredges);
          
        /* Find matching faces between submesh and main mesh */

        mrfaces = MR_Faces(sub_mr);
        nrf = List_Num_Entries(mrfaces);
        for (j = 0; j < nrf; j++) {
          sub_mf = List_Entry(mrfaces,j);

          MEnt_Get_AttVal(sub_mf, fidatt, &ival, &rval, &mf);

          if (!mf) {
            /* Does the global ID of this face of the sub mesh region
             * match the global ID of a boundary face in the main
             * mesh? */
            
            global_id = MF_GlobalID(sub_mf);
            loc = (int *) bsearch(&global_id, parbndry_face_gids, num_parbndry_faces, sizeof(int),
                                  compareINT);
            if (loc) {
              iloc = loc - parbndry_face_gids;
              mf = List_Entry(parbndry_faces,iloc); 
              /* here set the ghost edge property, only necessary when the input submeshes are not consistent */
              if (MF_PType(mf) == PGHOST && MF_PType(sub_mf) != PGHOST) {
                MF_Set_GEntDim(mf,MF_GEntDim(sub_mf));
                MF_Set_GEntID(mf,MF_GEntID(sub_mf));
              }

              MEnt_Set_AttVal(sub_mf, fidatt, 0, 0.0, mf);
            }
          }
        }

        if (!add_region) {
          List_Delete(mrfaces);
          continue;
        }
        
        new_mr = MR_New(mesh);                  /* add region */
        MR_Set_GEntDim(new_mr,MR_GEntDim(sub_mr));
        MR_Set_GEntID(new_mr,MR_GEntID(sub_mr));
        MR_Set_PType(new_mr,PGHOST);
        MR_Set_MasterParID(new_mr,MR_MasterParID(sub_mr));
        MR_Set_GlobalID(new_mr,MR_GlobalID(sub_mr));
	
        nrf = List_Num_Entries(mrfaces);
        int i2;
        for(i2 = 0; i2 < nrf; i2++) {
          sub_mf = List_Entry(mrfaces,i2);
          global_id = MF_GlobalID(sub_mf);
          rfdirs[i2] = MR_FaceDir_i(sub_mr,i2) == 1 ? 1 : 0;

          new_mf = NULL;
          MEnt_Get_AttVal(sub_mf, fidatt, &ival, &rval, &new_mf);

          if (!new_mf) {
            /* search in the ghost layer if another face with
             * this global ID has been added */
            loc = (int *) bsearch(&global_id, new_face_gids, num_new_faces,
                                  sizeof(int), compareINT);
            if (loc) {
              iloc = loc - new_face_gids;
              new_mf = List_Entry(new_faces, iloc);
              MEnt_Set_AttVal(sub_mf, fidatt, 0, 0.0, new_mf);
            }
          }

          if (new_mf) {
            List_ptr mfverts = MF_Vertices(sub_mf,1,0);
            int fvgid0[2];
            fvgid0[0] = MF_GlobalID(List_Entry(mfverts,0));
            fvgid0[1] = MF_GlobalID(List_Entry(mfverts,1));
            List_Delete(mfverts);

            mfverts = MF_Vertices(new_mf,1,0);
            int nfv = List_Num_Entries(mfverts);
            int fvgid1[MAXPV2];
            for (j = 0; j < nfv; j++)
              fvgid1[j] = MF_GlobalID(List_Entry(mfverts,j));
            List_Delete(mfverts);

            for (j = 0; j < nfv; j++) {
              if (fvgid1[j] == fvgid0[0]) {
                if (fvgid1[(j+nfv-1)%nfv] == fvgid0[1]) /* reverse dir */
                  rfdirs[i2] = !rfdirs[i2];
                break;
              }
            }                  
          }
          else {  /* add a new face to main mesh */
            new_mf = MF_New(mesh); /* add face */
            MF_Set_GEntDim(new_mf,MF_GEntDim(sub_mf));
            MF_Set_GEntID(new_mf,MF_GEntID(sub_mf));
            MF_Set_PType(new_mf,PGHOST);
            MF_Set_MasterParID(new_mf,MF_MasterParID(sub_mf));
            MF_Set_GlobalID(new_mf,MF_GlobalID(sub_mf));
	    
            MEnt_Set_AttVal(sub_mf, fidatt, 0, 0.0, new_mf);
            List_Add(new_faces, new_mf);
	    
            mfedges = MF_Edges(sub_mf,1,0);
            nfe = List_Num_Entries(mfedges);
            for(j = 0; j < nfe; j++) {
              sub_me = List_Entry(mfedges,j);
              global_id = ME_GlobalID(sub_me);
              
              fedirs[j] = MF_EdgeDir_i(sub_mf,j) == 1 ? 1 : 0;

              new_me = NULL;
              MEnt_Get_AttVal(sub_me, eidatt, &ival, &rval, &new_me);

              if (!new_me) {
                /* search in the ghost layer if another edge with
                 * this global ID has been added */
                loc = (int *) bsearch(&global_id, new_edge_gids, num_new_edges,
                                      sizeof(int), compareINT);
                if (loc) {
                  iloc = loc - new_edge_gids;
                  new_me = List_Entry(new_edges, iloc);
                  MEnt_Set_AttVal(sub_me, eidatt, 0, 0.0, new_me);
                }
              }
              
              if (new_me) {
                if(MV_GlobalID(ME_Vertex(new_me,0)) != MV_GlobalID(ME_Vertex(sub_me,0)))
                  fedirs[j] = 1 - fedirs[j];  /* if the edge dir is not the same, reverse the edge dir */
	      
              } else {  /* add a new edge to main mesh */
                new_me = ME_New(mesh);      /* add new edge and copy information */
                ME_Set_GEntDim(new_me,ME_GEntDim(sub_me));
                ME_Set_GEntID(new_me,ME_GEntID(sub_me));
                ME_Set_PType(new_me,PGHOST);
                ME_Set_MasterParID(new_me,ME_MasterParID(sub_me));
                ME_Set_GlobalID(new_me,ME_GlobalID(sub_me));
		
                MEnt_Set_AttVal(sub_me, eidatt, 0, 0.0, new_me);
                List_Add(new_edges, new_me);

                for(k = 0; k < 2; k++) {
                  sub_mv = ME_Vertex(sub_me,k);
                  global_id = MV_GlobalID(sub_mv);

                  new_mv = NULL;
                  MEnt_Get_AttVal(sub_mv, vidatt, &ival, &rval, &new_mv);

                  if (!new_mv) {
                    /* search in the ghost layer if another vertex with
                     * this global ID has been added */
                    loc = (int *) bsearch(&global_id, new_vert_gids, num_new_verts,
                                          sizeof(int), compareINT);
                    if (loc) {
                      iloc = loc - new_vert_gids;
                      new_mv = List_Entry(new_verts, iloc);
                      MEnt_Set_AttVal(sub_mv, vidatt, 0, 0.0, new_mv);
                    }
                  }
              
                  if (!new_mv) {  /* add new vertex to main mesh */
                    new_mv = MV_New(mesh);  /* add new vertex and copy information */
                    MV_Set_GEntDim(new_mv,MV_GEntDim(sub_mv));
                    MV_Set_GEntID(new_mv,MV_GEntID(sub_mv));
                    MV_Set_PType(new_mv,PGHOST);
                    MV_Set_MasterParID(new_mv,MV_MasterParID(sub_mv));
                    MV_Set_GlobalID(new_mv,MV_GlobalID(sub_mv));
                    MV_Coords(sub_mv,coor);
                    MV_Set_Coords(new_mv,coor);
		    
                    MEnt_Set_AttVal(sub_mv, vidatt, 0, 0.0, new_mv);
                    List_Add(new_verts, new_mv);
                  }
                  ME_Set_Vertex(new_me,k,new_mv);  /* set edge-vertex */
                }
              }							
              fedges[j] = new_me;
            }
            MF_Set_Edges(new_mf,nfe,fedges,fedirs); /* set face-edge */
            List_Delete(mfedges);
          }
          rfaces[i2] = new_mf;
        }
        MR_Set_Faces(new_mr,nrf,rfaces,rfdirs); /* set region-face */

        List_Delete(mrfaces);
      }

      idx = 0;
      while ((sub_mv = MESH_Next_Vertex(submesh, &idx)))
	MEnt_Rem_AttVal(sub_mv, vidatt);
      MAttrib_Delete(vidatt);
      idx = 0;
      while ((sub_me = MESH_Next_Edge(submesh, &idx)))
	MEnt_Rem_AttVal(sub_me, eidatt);
      MAttrib_Delete(eidatt);
      idx = 0;
      while ((sub_mf = MESH_Next_Face(submesh, &idx)))
	MEnt_Rem_AttVal(sub_mf, fidatt);
      MAttrib_Delete(fidatt);

      /* Sort the added entity lists by GlobalID */
      num_new_faces = List_Num_Entries(new_faces);
      List_Sort(new_faces, num_new_faces, sizeof(MFace_ptr), compareGlobalID);
      for (j = 0; j < num_new_faces; j++)
        new_face_gids[j] = MF_GlobalID(List_Entry(new_faces, j));

      num_new_edges = List_Num_Entries(new_edges);
      List_Sort(new_edges, num_new_edges, sizeof(MEdge_ptr), compareGlobalID);
      for (j = 0; j < num_new_edges; j++)
        new_edge_gids[j] = ME_GlobalID(List_Entry(new_edges, j));

      num_new_verts = List_Num_Entries(new_verts);
      List_Sort(new_verts, num_new_verts, sizeof(MVertex_ptr), compareGlobalID);
      for (j = 0; j < num_new_verts; j++)
        new_vert_gids[j] = MV_GlobalID(List_Entry(new_verts, j));
      
    }      

    idx = 0;
    while ((mv = List_Next_Entry(parbndry_verts, &idx)))
      MEnt_Rem_AttVal(mv, parbndryatt);
    MAttrib_Delete(parbndryatt);
    
    List_Delete(parbndry_faces);
    List_Delete(parbndry_edges);
    List_Delete(parbndry_verts);
    List_Delete(new_faces);
    List_Delete(new_edges);
    List_Delete(new_verts);

    free(parbndry_vert_gids);
    free(parbndry_edge_gids);
    free(parbndry_face_gids);
    free(new_face_gids);
    free(new_edge_gids);
    free(new_vert_gids);
    free(fedges);
    free(fedirs);
    free(rfaces);
    free(rfdirs);

    return 1;
  }

#ifdef __cplusplus
}
#endif

