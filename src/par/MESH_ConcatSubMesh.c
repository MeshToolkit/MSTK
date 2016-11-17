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

     NOTE: SEEMS TO BE DESIGNED TO ADD ONLY ONE LAYER OF ELEMENTS IN
     THE SUBMESHES - IT WOULD BE GOOD TO GENERALIZE IT

     Author(s): Duo Wang, Rao Garimella
  */

  int MESH_ConcatSubMesh_Face(Mesh_ptr mesh, int num, Mesh_ptr *submeshes);
  int MESH_ConcatSubMesh_Region(Mesh_ptr mesh, int num, Mesh_ptr *submeshes);


  int MESH_ConcatSubMesh(Mesh_ptr mesh, int topodim, int num, Mesh_ptr *submeshes) {
    
    fprintf(stderr,"In MESH_ConcatSubMesh....................\n");
    
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

  /* check if l has the same entity as m, based on global ID */
  MEntity_ptr entity_on_list(MEntity_ptr m, List_ptr *l) {
    int i, num = List_Num_Entries(*l);
    for(i = 0; i < num; i++)
      if(MEnt_GlobalID(m) == MEnt_GlobalID(List_Entry(*l,i)))
        return List_Entry(*l,i);
    return NULL;
  }
  
  int MESH_ConcatSubMesh_Face(Mesh_ptr mesh, int num, Mesh_ptr *submeshes) {
    int nfv, nfe, i, j, k, nbv, nbe;
    MVertex_ptr mv, new_mv, sub_mv;
    MEdge_ptr me, new_me, sub_me, *fedges;
    MFace_ptr new_mf, sub_mf;
    List_ptr mfverts, mfedges;
    List_ptr added_edges, added_verts;
    List_ptr boundary_edges, boundary_verts;
    int add_face, idx, *fedirs, global_id, iloc, *loc;
    double coor[3], rval;
    void *pval;
    Mesh_ptr submesh;
    int *MV_to_list_id, *ME_to_list_id, *MV_global_id, *ME_global_id;  

    added_edges = List_New(10);          
    added_verts = List_New(10);
  
    boundary_verts = List_New(10);
    boundary_edges = List_New(10);

    fedges = (MEdge_ptr *) malloc(MAXPV2*sizeof(MEdge_ptr));
    fedirs = (int *) malloc(MAXPV2*sizeof(int));

    /* collect boundary edges and vertices */
    idx = 0; nbe = 0;
    while ((me = MESH_Next_Edge(mesh,&idx))) 
      if (ME_PType(me) != PINTERIOR) {
        List_Add(boundary_edges,me);
        nbe++;
      }
    idx = 0; nbv = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx)))
      if (MV_PType(mv) != PINTERIOR) {
        List_Add(boundary_verts,mv);
        nbv++;
      }
    /* sort based on global ID */
    List_Sort(boundary_edges,nbe,sizeof(MEdge_ptr),compareGlobalID);
    List_Sort(boundary_verts,nbv,sizeof(MVertex_ptr),compareGlobalID);

    MV_global_id = (int *) malloc(nbv*sizeof(int));
    ME_global_id = (int *) malloc(nbe*sizeof(int));

    /* store them in array for binary search */
    for (i = 0; i < nbe; i++) {
      me = List_Entry(boundary_edges,i);
      ME_global_id[i] = ME_GlobalID(me);
    }
    for (i = 0; i < nbv; i++) {
      mv = List_Entry(boundary_verts,i);
      MV_global_id[i] = MV_GlobalID(mv);
    }

  
    for (i = 0; i < num; i++) {
      submesh = submeshes[i];
    
#ifdef MSTK_USE_MARKERS  
      int mkvid = MSTK_GetMarker(); /* mark the vertices in submeshes that is in mesh or already added into mesh */
      int mkeid = MSTK_GetMarker(); /* mark the edges in submeshes that is in mesh or already added into mesh */
    
      int mkvid2 = MSTK_GetMarker(); /* mark the vertices on mesh boundary, used to decide whether to add a face  */
      int mkeid2 = MSTK_GetMarker(); /* mark the edges on mesh boundary, used to decide whether to add a face  */
      List_ptr marked_edges = List_New(10);          
      List_ptr marked_verts = List_New(10);
#else
      MAttrib_ptr vidatt = MAttrib_New(submesh, "tempvid", INT, MVERTEX);
      MAttrib_ptr vidatt2 = MAttrib_New(submesh, "tempvid2", INT, MVERTEX);
      MAttrib_ptr eidatt = MAttrib_New(submesh, "tempeid", INT, MEDGE);
      MAttrib_ptr eidatt2 = MAttrib_New(submesh, "tempeid2", INT, MEDGE);
#endif

      MV_to_list_id = (int *)malloc(MESH_Num_Vertices(submesh)*sizeof(int));
      ME_to_list_id = (int *)malloc(MESH_Num_Edges(submesh)*sizeof(int));

      idx = 0;
      while ((sub_mf = MESH_Next_Face(submesh, &idx))) {
        add_face = 0;
      
        /* Find matching vertices between the submesh and main mesh */
      
        mfverts = MF_Vertices(sub_mf,1,0);
        nfv = List_Num_Entries(mfverts);
        for (j = 0; j < nfv; j++) {
          sub_mv = List_Entry(mfverts,j);
        
          /* Does the vertex have a known counterpart in the main mesh? */
          int vmarked;
#ifdef MSTK_USE_MARKERS
          vmarked = MEnt_IsMarked(sub_mv,mkvid2);
#else
          MEnt_Get_AttVal(sub_mv, vidatt2, &vmarked, &rval, &pval);
#endif
          if (vmarked) {
            add_face = 1; 
          } else {
        
            /* Does the global ID of this vertex of the sub mesh face
             * match the global ID of a boundary vertex in the main
             * mesh? */
            
            global_id = MV_GlobalID(sub_mv);
            loc = (int *) bsearch(&global_id, MV_global_id, nbv, sizeof(int),
                                  compareINT);
            if (loc) {  /* found a match */
              add_face = 1; 
              iloc = loc - MV_global_id;
              mv = List_Entry(boundary_verts,iloc); 
              /* here set the ghost vertex property, only necessary when the input submeshes are not consistent */
              if (MV_PType(mv) == PGHOST && MV_PType(sub_mv) != PGHOST) {
                MV_Set_GEntDim(mv,MV_GEntDim(sub_mv));
                MV_Set_GEntID(mv,MV_GEntID(sub_mv));
              }
              
              MV_to_list_id[MV_ID(sub_mv)-1] = MV_ID(mv)-1;
#ifdef MSTK_USE_MARKERS 
              MEnt_Mark(sub_mv,mkvid);
              MEnt_Mark(sub_mv,mkvid2);
              List_Add(marked_verts,sub_mv);
#else
              MEnt_Set_AttVal(sub_mv, vidatt, 1, 0.0, NULL);
              MEnt_Set_AttVal(sub_mv, vidatt2, 1, 0.0, NULL);
#endif
            }
          }
        }
      
        /* Find matching edges between the submesh and main mesh */
      
        mfedges = MF_Edges(sub_mf,1,0);
        nfe = List_Num_Entries(mfedges);
        for (j = 0; j < nfe; j++) {
          sub_me = List_Entry(mfedges,j);
        
          /* Does the edge have a known counterpart in the main mesh? */
          int emarked;
#ifdef MSTK_USE_MARKERS
          emarked = MEnt_IsMarked(sub_me,mkeid2);
#else
          MEnt_Get_AttVal(sub_me, eidatt2, &emarked, &rval, &pval);
#endif
          if (emarked) {
            add_face = 1; 
          } else {
        
            /* Does the global ID of this vertex of the sub mesh face
             * match the global ID of a boundary vertex in the main
             * mesh? */
            
            global_id = ME_GlobalID(sub_me);
            loc = (int *) bsearch(&global_id, ME_global_id, nbe, sizeof(int),
                                  compareINT);
            if (loc) {
              add_face = 1; 
              iloc = loc - ME_global_id;
              me = List_Entry(boundary_edges,iloc); 
              /* here set the ghost edge property, only necessary when the input submeshes are not consistent */
              if (ME_PType(me) == PGHOST && ME_PType(sub_me) != PGHOST) {
                ME_Set_GEntDim(me,ME_GEntDim(sub_me));
                ME_Set_GEntID(me,ME_GEntID(sub_me));
              }
              ME_to_list_id[ME_ID(sub_me)-1] = ME_ID(me)-1;
#ifdef MSTK_USE_MARKERS
              MEnt_Mark(sub_me,mkeid);
              MEnt_Mark(sub_me,mkeid2);
              List_Add(marked_edges,sub_me); 
#else
              MEnt_Set_AttVal(sub_me, eidatt, 1, 0.0, NULL);
              MEnt_Set_AttVal(sub_me, eidatt2, 1, 0.0, NULL);
#endif
            }
          }
        }
        List_Delete(mfverts);
        
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
          fedirs[j] = MF_EdgeDir_i(sub_mf,j) == 1 ? 1 : 0;
          new_me = NULL;
          int emarked;
#ifdef MSTK_USE_MARKERS
          emarked = MEnt_IsMarked(sub_me,mkeid);
#else
          MEnt_Get_AttVal(sub_me, eidatt, &emarked, &rval, &pval);
#endif
          if (emarked) /* has matching edge in the main mesh */
            new_me = MESH_Edge(mesh,ME_to_list_id[ME_ID(sub_me)-1]); 
          else  /* check if edge has already been added to main mesh */
            new_me = (MEdge_ptr)entity_on_list(sub_me,&added_edges); /* check if it is already added */
          if (new_me)
            if (MV_GlobalID(ME_Vertex(new_me,0)) != MV_GlobalID(ME_Vertex(sub_me,0)))
              fedirs[j] = 1 - fedirs[j];  /* if the edge dir is not the same, reverse the edge dir */
	
          if (!new_me)  {  /* add a new edge to main mesh */
            new_me = ME_New(mesh);
            ME_Set_GEntDim(new_me,ME_GEntDim(sub_me));
            ME_Set_GEntID(new_me,ME_GEntID(sub_me));
            ME_Set_PType(new_me,PGHOST);
            ME_Set_MasterParID(new_me,ME_MasterParID(sub_me));
            ME_Set_GlobalID(new_me,ME_GlobalID(sub_me));
	  
            ME_to_list_id[ME_ID(sub_me)-1] = ME_ID(new_me)-1;
            List_Add(added_edges,new_me);
          
            /* Mark edge as having been added to the main mesh */
#ifdef MSTK_USE_MARKERS
            MEnt_Mark(sub_me,mkeid);
            List_Add(marked_edges,sub_me);
#else
            MEnt_Set_AttVal(sub_me, eidatt, 1, 0.0, NULL);
#endif
          
            for (k = 0; k < 2; k++) {
              sub_mv = ME_Vertex(sub_me,k);
              new_mv = NULL;
            
              int vmarked;
#ifdef MSTK_USE_MARKERS
              vmarked = MEnt_IsMarked(sub_mv,mkvid);
#else
              MEnt_Get_AttVal(sub_mv, vidatt, &vmarked, &rval, &pval);
#endif
              if (vmarked) /* has matching vertex in main mesh */
                new_mv = MESH_Vertex(mesh,MV_to_list_id[MV_ID(sub_mv)-1]);
              else /* check if edge has already been added to main mesh */
                new_mv = (MVertex_ptr)entity_on_list(sub_mv,&added_verts);
	    
              if (!new_mv) {  /* add a new vertex to main mesh */
                new_mv = MV_New(mesh);
                MV_Set_GEntDim(new_mv,MV_GEntDim(sub_mv));
                MV_Set_GEntID(new_mv,MV_GEntID(sub_mv));
                MV_Set_PType(new_mv,PGHOST);
                MV_Set_MasterParID(new_mv,MV_MasterParID(sub_mv));
                MV_Set_GlobalID(new_mv,MV_GlobalID(sub_mv));
                MV_Coords(sub_mv,coor);
                MV_Set_Coords(new_mv,coor);
	      
                MV_to_list_id[MV_ID(sub_mv)-1] = MV_ID(new_mv)-1;
                List_Add(added_verts,new_mv);
              
                /* Mark the vertex as having been added to the main mesh */
#ifdef MSTK_USE_MARKERS
                MEnt_Mark(sub_mv,mkvid);
                List_Add(marked_verts,sub_mv);
#else
                MEnt_Set_AttVal(sub_mv, vidatt, 1, 0.0, NULL);
#endif
              }
              ME_Set_Vertex(new_me,k,new_mv);  /* set edge-vertex */
            }
          }								
          fedges[j] = new_me;
        }
        MF_Set_Edges(new_mf,nfe,fedges,fedirs); /* set face-edge */

        List_Delete(mfedges);
      }

#ifdef MSTK_USE_MARKERS
      List_Unmark(marked_edges,mkeid);
      List_Unmark(marked_verts,mkvid);
      List_Unmark(marked_edges,mkeid2);
      List_Unmark(marked_verts,mkvid2);
      List_Delete(marked_edges);
      List_Delete(marked_verts);
      MSTK_FreeMarker(mkeid);
      MSTK_FreeMarker(mkeid2);
      MSTK_FreeMarker(mkvid);
      MSTK_FreeMarker(mkvid2);
#else
      MAttrib_Delete(vidatt);
      MAttrib_Delete(vidatt2);
      MAttrib_Delete(eidatt);
      MAttrib_Delete(eidatt2);
#endif

      free(MV_to_list_id);
      free(ME_to_list_id);
    }


    List_Delete(boundary_edges);
    List_Delete(boundary_verts);
    List_Delete(added_edges);
    List_Delete(added_verts);

    free(MV_global_id);
    free(ME_global_id);
    free(fedges);
    free(fedirs);
    return 1;
  }

  /* right now assume there are no overlapped regions */

  int MESH_ConcatSubMesh_Region(Mesh_ptr mesh, int num, Mesh_ptr *submeshes) {
    int nrf, nre, nrv, nfe, i, j, k, nbv, nbe, nbf;
    MVertex_ptr mv, new_mv, sub_mv;
    MEdge_ptr me, new_me, sub_me, *fedges;
    MFace_ptr mf, new_mf, sub_mf, *rfaces;
    MRegion_ptr new_mr, sub_mr;
    List_ptr mrfaces, mredges, mrverts, mfedges;
    List_ptr added_faces, added_edges, added_verts;
    List_ptr boundary_faces, boundary_edges, boundary_verts;
    int add_region, idx, *rfdirs, *fedirs, global_id, iloc, *loc;
    double coor[3], rval;
    void *pval;
    Mesh_ptr submesh;
    int *MF_to_list_id, *MV_to_list_id, *ME_to_list_id, *MV_global_id, *ME_global_id, *MF_global_id;  

    added_faces = List_New(10);  boundary_verts = List_New(10);        
    added_edges = List_New(10);  boundary_edges = List_New(10);
    added_verts = List_New(10);  boundary_faces = List_New(10);

    rfaces = (MFace_ptr *) malloc(MAXPF3*sizeof(MFace_ptr));
    rfdirs = (int *) malloc(MAXPF3*sizeof(int));
    fedges = (MEdge_ptr *) malloc(MAXPV2*sizeof(MEdge_ptr));
    fedirs = (int *) malloc(MAXPV2*sizeof(int));

    /* collect boundary faces, edges and vertices */
    idx = 0; nbf = 0;
    while ((mf = MESH_Next_Face(mesh,&idx))) 
      if (MF_PType(mf) != PINTERIOR) {
        List_Add(boundary_faces,mf);
        nbf++;
      }
    idx = 0; nbe = 0;
    while ((me = MESH_Next_Edge(mesh,&idx))) 
      if (ME_PType(me) != PINTERIOR) {
        List_Add(boundary_edges,me);
        nbe++;
      }
    idx = 0; nbv = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx)))
      if (MV_PType(mv) != PINTERIOR) {
        List_Add(boundary_verts,mv);
        nbv++;
      }
    /* sort based on global ID */
    List_Sort(boundary_faces,nbf,sizeof(MFace_ptr),compareGlobalID);
    List_Sort(boundary_edges,nbe,sizeof(MEdge_ptr),compareGlobalID);
    List_Sort(boundary_verts,nbv,sizeof(MVertex_ptr),compareGlobalID);

    MV_global_id = (int *)malloc(nbv*sizeof(int));
    ME_global_id = (int *)malloc(nbe*sizeof(int));
    MF_global_id = (int *)malloc(nbf*sizeof(int));

    /* store them in array for binary search */
    for (i = 0; i < nbf; i++) {
      mf = List_Entry(boundary_faces,i);
      MF_global_id[i] = MF_GlobalID(mf);
    }
    for (i = 0; i < nbe; i++) {
      me = List_Entry(boundary_edges,i);
      ME_global_id[i] = ME_GlobalID(me);
    }
    for (i = 0; i < nbv; i++) {
      mv = List_Entry(boundary_verts,i);
      MV_global_id[i] = MV_GlobalID(mv);
    }

    for (i = 0; i < num; i++) {
      submesh = submeshes[i];

#ifdef MSTK_USE_MARKERS  
      int mkvid = MSTK_GetMarker(); /* mark the vertices in submeshes that is in mesh or already added into mesh */
      int mkeid = MSTK_GetMarker(); /* mark the edges in submeshes that is in mesh or already added into mesh */
      int mkfid = MSTK_GetMarker(); /* mark the faces in submeshes that is in mesh or already added into mesh */
    
      int mkvid2 = MSTK_GetMarker(); /* mark the vertices on mesh boundary, used to decide whether to add a face  */
      int mkeid2 = MSTK_GetMarker(); /* mark the edges on mesh boundary, used to decide whether to add a face  */
      int mkfid2 = MSTK_GetMarker(); /* mark the faces on mesh boundary, used to decide whether to add a face  */
      List_ptr marked_faces = List_New(10);
      List_ptr marked_edges = List_New(10);
      List_ptr marked_verts = List_New(10);
#else
      MAttrib_ptr vidatt = MAttrib_New(submesh, "tempvid", INT, MVERTEX);
      MAttrib_ptr vidatt2 = MAttrib_New(submesh, "tempvid2", INT, MVERTEX);
      MAttrib_ptr eidatt = MAttrib_New(submesh, "tempeid", INT, MEDGE);
      MAttrib_ptr eidatt2 = MAttrib_New(submesh, "tempeid2", INT, MEDGE);
      MAttrib_ptr fidatt = MAttrib_New(submesh, "tempfid", INT, MFACE);
      MAttrib_ptr fidatt2 = MAttrib_New(submesh, "tempfid2", INT, MFACE);
#endif

      MV_to_list_id = (int *)malloc(MESH_Num_Vertices(submesh)*sizeof(int));
      ME_to_list_id = (int *)malloc(MESH_Num_Edges(submesh)*sizeof(int));
      MF_to_list_id = (int *)malloc(MESH_Num_Faces(submesh)*sizeof(int));

      idx = 0;
      while ((sub_mr = MESH_Next_Region(submesh, &idx))) {
        add_region = 0;

        /* Find matching faces between submesh and main mesh */

        mrfaces = MR_Faces(sub_mr);
        nrf = List_Num_Entries(mrfaces);
        for (j = 0; j < nrf; j++) {
          sub_mf = List_Entry(mrfaces,j);

          /* Does the face already have a counterpart in the main mesh? */
          int fmarked;
#ifdef MSTK_USE_MARKERS
          fmarked = MEnt_IsMarked(sub_mf,mkfid2);
#else
          MEnt_Get_AttVal(sub_mf, fidatt2, &fmarked, &rval, &pval);
#endif
          if (fmarked) {
            add_region = 1; 
          } else {

            /* Does the global ID of this face of the sub mesh region
             * match the global ID of a boundary face in the main
             * mesh? */
            
            global_id = MF_GlobalID(sub_mf);
            loc = (int *) bsearch(&global_id, MF_global_id, nbf, sizeof(int),
                                  compareINT);
            if (loc) {
              add_region = 1; 
              iloc = loc - MF_global_id;
              mf = List_Entry(boundary_faces,iloc); 
              /* here set the ghost edge property, only necessary when the input submeshes are not consistent */
              if (MF_PType(mf) == PGHOST && MF_PType(sub_mf) != PGHOST) {
                MF_Set_GEntDim(mf,MF_GEntDim(sub_mf));
                MF_Set_GEntID(mf,MF_GEntID(sub_mf));
              }
              MF_to_list_id[MF_ID(sub_mf)-1] = MF_ID(mf)-1;
#ifdef MSTK_USE_MARKERS
              MEnt_Mark(sub_mf,mkfid);
              MEnt_Mark(sub_mf,mkfid2);
              List_Add(marked_faces,sub_mf); 
#else
              MEnt_Set_AttVal(sub_mf, fidatt, 1, 0.0, NULL);
              MEnt_Set_AttVal(sub_mf, fidatt2, 1, 0.0, NULL);
#endif
            }
          }
        }

        /* Find matching edges between submesh and main mesh */

        mredges = MR_Edges(sub_mr);
        nre = List_Num_Entries(mredges);
        for (j = 0; j < nre; j++) {
          sub_me = List_Entry(mredges,j);
          
          /* Does the edge already have a counterpart in the main mesh? */
          int emarked;
#ifdef MSTK_USE_MARKERS
          emarked = MEnt_IsMarked(sub_me,mkeid2);
#else
          MEnt_Get_AttVal(sub_me, eidatt2, &emarked, &rval, &pval);
#endif
          if (emarked) {
            add_region = 1; 
          } else {

            /* Does the global ID of this edge of the sub mesh region
             * match the global ID of a boundary edge in the main
             * mesh? */
            
            global_id = ME_GlobalID(sub_me);
            loc = (int *)bsearch(&global_id,
                                 ME_global_id,
                                 nbe,
                                 sizeof(int),
                                 compareINT);
            if (loc) {
              add_region = 1; 
              iloc = loc - ME_global_id;
              me = List_Entry(boundary_edges,iloc); 
              /* here set the ghost edge property, only necessary when the input submeshes are not consistent */
              if(ME_PType(me) == PGHOST && ME_PType(sub_me) != PGHOST) {
                ME_Set_GEntDim(me,ME_GEntDim(sub_me));
                ME_Set_GEntID(me,ME_GEntID(sub_me));
              }
              ME_to_list_id[ME_ID(sub_me)-1] = ME_ID(me)-1;
#ifdef MSTK_USE_MARKERS
              MEnt_Mark(sub_me,mkeid);
              MEnt_Mark(sub_me,mkeid2);
              List_Add(marked_edges,sub_me); 
#else
              MEnt_Set_AttVal(sub_me, eidatt, 1, 0.0, NULL);
              MEnt_Set_AttVal(sub_me, eidatt2, 1, 0.0, NULL);
#endif
            }
          }
        }
          
        /* Find matching vertices between submesh and main mesh */

        mrverts = MR_Vertices(sub_mr);
        nrv = List_Num_Entries(mrverts);
        for (j = 0; j < nrv; j++) {
          sub_mv = List_Entry(mrverts,j);

          /* Does this vertex have a counterpart in the main mesh? */
          int vmarked;
#ifdef MSTK_USE_MARKERS         
          vmarked = MEnt_IsMarked(sub_mv,mkvid2);
#else
          MEnt_Get_AttVal(sub_mv, &vidatt2, &vmarked, &rval, &pval);
#endif
          if (vmarked) {
            add_region = 1; 
          } else {

            /* Does the global ID of this vertex of the sub mesh region
             * match the global ID of a boundary vertex in the main
             * mesh? */
            
            global_id = MV_GlobalID(sub_mv);
            loc = (int *)bsearch(&global_id,
                                 MV_global_id,
                                 nbv,
                                 sizeof(int),
                                 compareINT);
            if (loc) {
              add_region = 1; 
              iloc = loc - MV_global_id;
              mv = List_Entry(boundary_verts,iloc); 
              /* here set the ghost vertex property, only necessary when the input submeshes are not consistent */
              if(MV_PType(mv) == PGHOST && MV_PType(sub_mv) != PGHOST) {
                MV_Set_GEntDim(mv,MV_GEntDim(sub_mv));
                MV_Set_GEntID(mv,MV_GEntID(sub_mv));
              }
              
              MV_to_list_id[MV_ID(sub_mv)-1] = MV_ID(mv)-1;
#ifdef MSTK_USE_MARKERS 
              MEnt_Mark(sub_mv,mkvid);
              MEnt_Mark(sub_mv,mkvid2);
              List_Add(marked_verts,sub_mv);
#else
              MEnt_Set_AttVal(sub_mv, vidatt, 1, 0.0, NULL);
              MEnt_Set_AttVal(sub_mv, vidatt2, 1, 0.0, NULL);
#endif
            }
          }
        }
        List_Delete(mrverts);
        List_Delete(mredges);
	
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
          rfdirs[i2] = MR_FaceDir_i(sub_mr,i2) == 1 ? 1 : 0;
          new_mf = NULL;
          int fmarked;
#ifdef MSTK_USE_MARKERS
          fmarked = MEnt_IsMarked(sub_mf,mkfid);
#else
          MEnt_Get_AttVal(sub_mf, fidatt, &fmarked, &rval, &pval);
#endif
          if (fmarked)  /* has matching face in main mesh */
            new_mf = MESH_Face(mesh,MF_to_list_id[MF_ID(sub_mf)-1]); 
          else  /* check if face has already been added to main mesh */
            new_mf = (MFace_ptr)entity_on_list(sub_mf,&added_faces); /* check if it is already added */
	  
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
	    
            MF_to_list_id[MF_ID(sub_mf)-1] = MF_ID(new_mf)-1;
            List_Add(added_faces,new_mf);
#ifdef MSTK_USE_MARKERS
            MEnt_Mark(sub_mf,mkfid);
            List_Add(marked_faces,sub_mf);
#else
            MEnt_Set_AttVal(sub_mf, fidatt, 1, 0.0, NULL);
#endif
	    
            mfedges = MF_Edges(sub_mf,1,0);
            nfe = List_Num_Entries(mfedges);
            for(j = 0; j < nfe; j++) {
              sub_me = List_Entry(mfedges,j);
              fedirs[j] = MF_EdgeDir_i(sub_mf,j) == 1 ? 1 : 0;
              new_me = NULL;
              int emarked;
#ifdef MSTK_USE_MARKERS
              emarked = MEnt_IsMarked(sub_me,mkeid);
#else
              MEnt_Get_AttVal(sub_me, eidatt, &emarked, &rval, &pval);
#endif
              if (emarked)  /* has a matching edge in main mesh */
                new_me = MESH_Edge(mesh,ME_to_list_id[ME_ID(sub_me)-1]); 
              else   /* check if edge has already been added to main mesh */
                new_me = (MEdge_ptr)entity_on_list(sub_me,&added_edges); /* check if it is already added */
              if (new_me)
                if(MV_GlobalID(ME_Vertex(new_me,0)) != MV_GlobalID(ME_Vertex(sub_me,0)))
                  fedirs[j] = 1 - fedirs[j];  /* if the edge dir is not the same, reverse the edge dir */
	      
              if (!new_me)  {  /* add a new edge to main mesh */
                new_me = ME_New(mesh);      /* add new edge and copy information */
                ME_Set_GEntDim(new_me,ME_GEntDim(sub_me));
                ME_Set_GEntID(new_me,ME_GEntID(sub_me));
                ME_Set_PType(new_me,PGHOST);
                ME_Set_MasterParID(new_me,ME_MasterParID(sub_me));
                ME_Set_GlobalID(new_me,ME_GlobalID(sub_me));
		
                ME_to_list_id[ME_ID(sub_me)-1] = ME_ID(new_me)-1;
                List_Add(added_edges,new_me);
#ifdef MSTK_USE_MARKERS
                MEnt_Mark(sub_me,mkeid);
                List_Add(marked_edges,sub_me);
#else
                MEnt_Set_AttVal(sub_me, eidatt, 1, 0.0, NULL);
#endif
                for(k = 0; k < 2; k++) {
                  sub_mv = ME_Vertex(sub_me,k);
                  new_mv = NULL;

                  int vmarked;
#ifdef MSTK_USE_MARKERS
                  vmarked = MEnt_IsMarked(sub_mv,mkvid);
#else
                  MEnt_Get_AttVal(sub_mv, vidatt, &vmarked, &rval, &pval);
#endif
                  if (vmarked)  /* has matching vertex in main mesh */
                    new_mv = MESH_Vertex(mesh,MV_to_list_id[MV_ID(sub_mv)-1]);
                  else  /* check if vertex has already been added to main mesh */
                    new_mv = (MVertex_ptr)entity_on_list(sub_mv,&added_verts);
		  
                  if (!new_mv) {  /* add new vertex to main mesh */
                    new_mv = MV_New(mesh);  /* add new vertex and copy information */
                    MV_Set_GEntDim(new_mv,MV_GEntDim(sub_mv));
                    MV_Set_GEntID(new_mv,MV_GEntID(sub_mv));
                    MV_Set_PType(new_mv,PGHOST);
                    MV_Set_MasterParID(new_mv,MV_MasterParID(sub_mv));
                    MV_Set_GlobalID(new_mv,MV_GlobalID(sub_mv));
                    MV_Coords(sub_mv,coor);
                    MV_Set_Coords(new_mv,coor);
		    
                    MV_to_list_id[MV_ID(sub_mv)-1] = MV_ID(new_mv)-1;
                    List_Add(added_verts,new_mv);
#ifdef MSTK_USE_MARKERS
                    MEnt_Mark(sub_mv,mkvid);
                    List_Add(marked_verts,sub_mv);
#else
                    MEnt_Set_AttVal(sub_mv, vidatt, 1, 0.0, NULL);
#endif
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

#ifdef MSTK_USE_MARKERS
      List_Unmark(marked_faces,mkfid);
      List_Unmark(marked_faces,mkfid2);
      List_Delete(marked_faces);
      List_Unmark(marked_edges,mkeid);
      List_Unmark(marked_edges,mkeid2);
      List_Delete(marked_edges);
      List_Unmark(marked_verts,mkvid);
      List_Unmark(marked_verts,mkvid2);
      List_Delete(marked_verts);
      MSTK_FreeMarker(mkfid);
      MSTK_FreeMarker(mkeid);
      MSTK_FreeMarker(mkvid);
      MSTK_FreeMarker(mkfid2);
      MSTK_FreeMarker(mkeid2);
      MSTK_FreeMarker(mkvid2);
#else
      MAttrib_Delete(vidatt);
      MAttrib_Delete(vidatt2);
      MAttrib_Delete(eidatt);
      MAttrib_Delete(eidatt2);
      MAttrib_Delete(fidatt);
      MAttrib_Delete(fidatt2);
#endif

      free(MV_to_list_id);
      free(ME_to_list_id);
      free(MF_to_list_id);
    }      
      

    List_Delete(boundary_faces);
    List_Delete(boundary_edges);
    List_Delete(boundary_verts);
    List_Delete(added_faces);
    List_Delete(added_edges);
    List_Delete(added_verts);

    free(MV_global_id);
    free(ME_global_id);
    free(MF_global_id);
    free(fedges);
    free(fedirs);
    free(rfaces);
    free(rfdirs);

    return 1;
  }

#ifdef __cplusplus
}
#endif

