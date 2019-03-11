/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* 
     This function adds ghost entities and labels POVERLAP entities

     ring: ghost cell size, 
     0: only ghost vertex on processor boundary
     1: 1-ring of ghost face(surface mesh) or region(volume mesh)

     must call MESH_BuildPBoundary() first

  */

  int MESH_Vol_AddGhost_FN(Mesh_ptr mesh, Mesh_ptr submesh, int part_no, int ring);
  int MESH_Surf_AddGhost_FN(Mesh_ptr mesh, Mesh_ptr submesh, int part_no, int ring);
  int MESH_Vol_AddGhost_R4(Mesh_ptr mesh, Mesh_ptr submesh, int part_no, int ring);
  int MESH_Vol_AddGhost_R1R2(Mesh_ptr mesh, Mesh_ptr submesh, int part_no, int ring);
  int MESH_Surf_AddGhost_R1R2R4(Mesh_ptr mesh, Mesh_ptr submesh, int part_no, int ring);

  static int (*MESH_Vol_AddGhost_jmp[MSTK_MAXREP])(Mesh_ptr mesh, 
						   Mesh_ptr submesh, 
						   int part_no, int ring) =
  {MESH_Vol_AddGhost_FN, MESH_Vol_AddGhost_FN, MESH_Vol_AddGhost_R1R2,
   MESH_Vol_AddGhost_R1R2, MESH_Vol_AddGhost_R4};

  static int (*MESH_Surf_AddGhost_jmp[MSTK_MAXREP])(Mesh_ptr mesh, 
						    Mesh_ptr submesh, 
						    int part_no, int ring) =
  {MESH_Surf_AddGhost_FN, MESH_Surf_AddGhost_FN, MESH_Surf_AddGhost_R1R2R4,
   MESH_Surf_AddGhost_R1R2R4, MESH_Surf_AddGhost_R1R2R4};


int MESH_AddGhost(Mesh_ptr mesh, Mesh_ptr submesh, int part_no, int ring) {
  int nf, nr;
  RepType rtype;


  /* basic mesh information */
  nf = MESH_Num_Faces(submesh);
  nr = MESH_Num_Regions(submesh);
  rtype = MESH_RepType(submesh);

  if (ring > 0) {

    if (nr) {
      (*MESH_Vol_AddGhost_jmp[rtype])(mesh,submesh,part_no,ring);
      MESH_Build_GhostLists(submesh,3);
    }
    else if(nf) {
      (*MESH_Surf_AddGhost_jmp[rtype])(mesh,submesh,part_no,ring);
      MESH_Build_GhostLists(submesh,2);
    }
    else {
      MSTK_Report("MESH_AddGhost()","This is not a valid mstk file",MSTK_ERROR);
      exit(-1);
    }

  }


  /* Mark all PBoundary entities as POVERLAP or PGHOST */

  int idx;
  
  MVertex_ptr lmv;
  idx = 0;
  while((lmv = MESH_Next_Vertex(submesh,&idx))) {
    if (MV_OnParBoundary(lmv)) {
      if(MV_MasterParID(lmv) == part_no)
        MV_Set_PType(lmv,POVERLAP);
      else
        MV_Set_PType(lmv,PGHOST);
    }
  }
  
  MEdge_ptr lme;
  idx = 0;
  while((lme = MESH_Next_Edge(submesh,&idx))) {
    if (ME_OnParBoundary(lme)) {
      if(MV_MasterParID(lme) == part_no)
        ME_Set_PType(lme,POVERLAP);
      else
        ME_Set_PType(lme,PGHOST);
    }
  }
  
  if (nr) {
    MFace_ptr lmf;
    idx = 0;
    while((lmf = MESH_Next_Face(submesh,&idx))) {      
      if (MF_OnParBoundary(lmf)) {
        if(MF_MasterParID(lmf) == part_no)
          MF_Set_PType(lmf,POVERLAP);
        else
          MF_Set_PType(lmf,PGHOST);
      }
    }
  }
  

  return 1;
}



int MESH_Surf_AddGhost_FN(Mesh_ptr mesh, Mesh_ptr submesh, int part_no, int ring) {
  int i, j, k, l, ir, idx;
  int nvf, nfv, nfe, idx2, found;
  double xyz[3];
  MVertex_ptr gmv, gmv2, lmv, lmv2;
  MEdge_ptr lme, gme;
  MFace_ptr lmf, gmf;
  List_ptr vfaces, fverts, fedges;
  List_ptr bverts, bverts2, ovfaces, gbverts, gbverts2;
  List_ptr lmvlist, lmelist, lmflist;
  MAttrib_ptr g2latt, l2gatt;
  double rval;
  void *pval;

  g2latt = MESH_AttribByName(mesh,"Global2Local");
  l2gatt = MESH_AttribByName(submesh,"Local2Global");


  /* Mark the list of global elements in this submesh */

#ifdef MSTK_USE_MARKERS
  int mkfid = MSTK_GetMarker();  
#else
  MAttrib_ptr gfidatt = MAttrib_New(mesh, "temp_gf_id", INT, MFACE);
  MAttrib_ptr lfidatt = MAttrib_New(submesh, "temp_lf_id", INT, MFACE);
#endif
  idx = 0;
  while((lmf = MESH_Next_Face(submesh,&idx))) {
    MEnt_Get_AttVal(lmf,l2gatt,0,0,&gmf);
#ifdef MSTK_USE_MARKERS
    MEnt_Mark(gmf,mkfid);
#else
    MEnt_Set_AttVal(gmf, gfidatt, 1, 0.0, NULL);
#endif
  }

  /* Collect all the PBOUNDARY vertices of submesh */
#ifdef MSTK_USE_MARKERS
  int mkvid = MSTK_GetMarker();
#else
  MAttrib_ptr lvidatt = MAttrib_New(mesh, "temp_lv_id", INT, MVERTEX);
  MAttrib_ptr gvidatt = MAttrib_New(submesh, "temp_gv_id", INT, MVERTEX);
#endif

  bverts = List_New(10);
  gbverts = List_New(10);
  idx = 0;
  while ((lmv = MESH_Next_Vertex(submesh,&idx))) {
    if (MV_OnParBoundary(lmv)) {
      List_Add(bverts,lmv);

      MEnt_Get_AttVal(lmv,l2gatt,0,0,&gmv);
      List_Add(gbverts,gmv);

#ifdef MSTK_USE_MARKERS
      MEnt_Mark(lmv,mkvid);
      MEnt_Mark(gmv,mkvid);
#else
      MEnt_Set_AttVal(lmv, lvidatt, 1, 0.0, NULL);
      MEnt_Set_AttVal(gmv, gvidatt, 1, 0.0, NULL);
#endif
    }
  }

  /* Mark 'N' rings of elements going in from the boundary as overlap elements */


  ovfaces = List_New(List_Num_Entries(bverts));

  for (ir = 0; ir < ring; ir++) {

    bverts2 = List_New(List_Num_Entries(bverts));

    idx = 0;
    while ((lmv = List_Next_Entry(bverts,&idx))) {

      vfaces = MV_Faces(lmv);
      nvf = List_Num_Entries(vfaces);

      for (j = 0; j < nvf; j++) {
        lmf = List_Entry(vfaces,j);
        
        int fmarked;
#ifdef MSTK_USE_MARKERS 
        fmarked = MEnt_IsMarked(lmf,mkfid);
#else
        MEnt_Get_AttVal(lmf, lfidatt, &fmarked, &rval, &pval);
#endif
        if (!fmarked) {
          List_Add(ovfaces,lmf);
          MF_Set_PType(lmf,POVERLAP);          

#ifdef MSTK_USE_MARKERS
          MEnt_Mark(lmf,mkfid);
#else
          MEnt_Set_AttVal(lmf, lfidatt, 1, 0.0, NULL);
#endif

          fedges = MF_Edges(lmf,1,0);
          nfe = List_Num_Entries(fedges);

          for (k = 0; k < nfe; k++) {
            lme = List_Entry(fedges,k);
            if (!ME_OnParBoundary(lme))
              ME_Set_PType(lme,POVERLAP);

            for (l = 0; l < 2; l++) {
              lmv2 = ME_Vertex(lme,l);
              int vmarked;
#ifdef MSTK_USE_MARKERS
              vmarked = MEnt_IsMarked(lmv2,mkvid);
#else
              MEnt_Get_AttVal(lmv2, lvidatt, &vmarked, &rval, &pval);
#endif
              if (!vmarked) {
                if (!ME_OnParBoundary(lmv2))
                  MV_Set_PType(lmv2,POVERLAP);
                List_Add(bverts2,lmv2);
#ifdef MSTK_USE_MARKERS
                MEnt_Mark(lmv2,mkvid);
#else
                MEnt_Set_AttVal(lmv2, lvidatt, 1, 0.0, NULL);
#endif
              }
            } /* l */

          } /* k */
          List_Delete(fedges);
        }
      } /* j */
      List_Delete(vfaces);

    }
#ifdef MSTK_USE_MARKERS
    List_Unmark(bverts,mkvid);
#else
    idx = 0;
    while ((lmv = List_Next_Entry(bverts, &idx)))
      MEnt_Set_AttVal(lmv, lvidatt, 0, 0.0, NULL);
#endif
    List_Delete(bverts);
    
    bverts = bverts2;

  } /* for (i = 0; i < ring; i++) */
  
#ifdef MSTK_USE_MARKERS
  List_Unmark(bverts,mkvid);
#else
  idx = 0;
  while ((lmv = List_Next_Entry(bverts, &idx)))
    MEnt_Set_AttVal(lmv, lvidatt, 0, 0.0, NULL);
#endif
  List_Delete(bverts);

#ifdef MSTK_USE_MARKERS
  List_Unmark(ovfaces,mkfid);
#else
  idx = 0;
  while ((lmf = List_Next_Entry(ovfaces, &idx)))
    MEnt_Set_AttVal(lmf, lfidatt, 0, 0.0, NULL);
#endif
  List_Delete(ovfaces);


  /* Duplicate 'N' rings of elements around submesh from global mesh */

  MEdge_ptr *lfedges = (MEdge_ptr *) malloc(MAXPV2*sizeof(MEdge_ptr));
  int *lfedirs = (int *) malloc(MAXPV2*sizeof(int));

  for (ir = 0; ir < ring; ir++) {
    
    gbverts2 = List_New(List_Num_Entries(gbverts));

    idx = 0; 
    while ((gmv = List_Next_Entry(gbverts,&idx))) {

      vfaces = MV_Faces(gmv);
      nvf = List_Num_Entries(vfaces);

      for (j = 0; j < nvf; j++) {
	gmf = List_Entry(vfaces,j);
        int fmarked;
#ifdef MSTK_USE_MARKERS
        fmarked = MEnt_IsMarked(gmf,mkfid);
#else
        MEnt_Get_AttVal(gmf, gfidatt, &fmarked, &rval, &pval);
#endif
	if (!fmarked) { /* not overlapping with submesh */

          /* Found a ghost face that has not yet been processed */
          /* Add to list of ghost faces and duplicate in submesh */

#ifdef MSTK_USE_MARKERS
	  MEnt_Mark(gmf,mkfid);
#else
          MEnt_Set_AttVal(gmf, gfidatt, 1, 0.0, NULL);
#endif

          fedges = MF_Edges(gmf,1,0);
          nfe = List_Num_Entries(fedges);
          
          for (k = 0; k < nfe; k++) {
            gme = List_Entry(fedges,k);
            MEnt_Get_AttVal(gme,g2latt,0,0,&lmelist);
            
            if (!lmelist)
              MSTK_Report("MESH_AddGhost","No local edge list found with global edge list",MSTK_FATAL);
            
            idx2 = 0;
            found = 0;
            while ((lme = List_Next_Entry(lmelist,&idx2))) {
              if (ME_Mesh(lme) == submesh) {
                found = 1;
                break;
              }
            }
            
            if (!found) { 
              lme = ME_New(submesh);
              ME_Set_GEntID(lme,ME_GEntID(gme));
              ME_Set_GEntDim(lme,ME_GEntDim(gme));
              ME_Set_PType(lme,PGHOST);

              for (l = 0; l < 2; l++) {
                gmv2 = ME_Vertex(gme,l);
                MEnt_Get_AttVal(gmv2,g2latt,0,0,&lmvlist);

                idx2 = 0;
                found = 0;
                while ((lmv = List_Next_Entry(lmvlist,&idx2))) {
                  if (MV_Mesh(lmv) == submesh) {
                    found = 1;
                    break;
                  }
                }
                
                if (!found) {
                  lmv = MV_New(submesh);
                  MV_Coords(gmv2,xyz);
                  MV_Set_Coords(lmv,xyz);
                  MV_Set_GEntID(lmv,MV_GEntID(gmv2));
                  MV_Set_GEntDim(lmv,MV_GEntDim(gmv2));
                  MV_Set_PType(lmv,PGHOST);
                  MV_Set_GlobalID(lmv,MV_GlobalID(gmv2));
                  MV_Set_MasterParID(lmv,MV_MasterParID(gmv2));
                  MEnt_Set_AttVal(lmv,l2gatt,0,0,gmv2);
                  List_Add(lmvlist,lmv);
                }

                ME_Set_Vertex(lme,l,lmv);

                /* Also, if this is an unmarked global vertex then it
                   represent an outer vertex of the layer of ghost
                   elements. Add it to the gbverts2 list to fetch the
                   next layer of ghost elements */
                    
                int vmarked;
#ifdef MSTK_USE_MARKERS
                vmarked = MEnt_IsMarked(gmv2,mkvid);
#else
                MEnt_Get_AttVal(gmv2, gvidatt, &vmarked, &rval, &pval);
#endif
                if (!vmarked) {
                  List_Add(gbverts2,gmv2);

#ifdef MSTK_USE_MARKERS
                  MEnt_Mark(gmv2,mkvid);
#else
                  MEnt_Set_AttVal(gmv2, gvidatt, 1, 0.0, NULL);
#endif
                }
              } /* l */

              ME_Set_GlobalID(lme,ME_GlobalID(gme));
              ME_Set_MasterParID(lme,ME_MasterParID(gme));
              MEnt_Set_AttVal(lme,l2gatt,0,0,gme);
              List_Add(lmelist,lme);
            }

            lfedges[k] = lme;
            lfedirs[k] = MF_EdgeDir_i(gmf,k);	
          } /* k */
          
          List_Delete(fedges);          
          
          lmf = MF_New(submesh);
          MF_Set_GEntID(lmf,MF_GEntID(gmf));
          MF_Set_GEntDim(lmf,MF_GEntDim(gmf));
          MF_Set_PType(lmf,PGHOST);
          MF_Set_Edges(lmf,nfe,lfedges,lfedirs);
          MF_Set_GlobalID(lmf,MF_GlobalID(gmf));
          MF_Set_MasterParID(lmf,MF_MasterParID(gmf));
          
          MEnt_Set_AttVal(lmf,l2gatt,0,0,gmf);
          MEnt_Get_AttVal(gmf,g2latt,0,0,&lmflist); /* do we need to add link for ghost faces? */
          if (!lmflist)
            MSTK_Report("MESH_AddGhost","Local face list not found with global face",MSTK_FATAL);
          List_Add(lmflist,lmf);
        }

      } /* j */
      List_Delete(vfaces);

    }

#ifdef MSTK_USE_MARKERS
    List_Unmark(gbverts,mkvid);    
#else
    idx = 0;
    while ((gmv = List_Next_Entry(gbverts, &idx)))
      MEnt_Rem_AttVal(gmv, gvidatt);
#endif
    List_Delete(gbverts);
    
#ifdef MSTK_USE_MARKERS
    List_Mark(gbverts2,mkvid);
#else
    idx = 0;
    while ((gmv = List_Next_Entry(gbverts2, &idx)))
      MEnt_Set_AttVal(gmv, gvidatt, 1, 0.0, NULL);
#endif
    gbverts = gbverts2;
  }

#ifdef MSTK_USE_MARKERS
  List_Unmark(gbverts,mkvid);
#else
    idx = 0;
    while ((gmv = List_Next_Entry(gbverts, &idx)))
      MEnt_Rem_AttVal(gmv, gvidatt);
#endif
  List_Delete(gbverts);

  free(lfedges);
  free(lfedirs);


#ifdef MSTK_USE_MARKERS
  /* Unmark all the global faces corresponding to this submesh and its
     ghost layer */

  idx = 0;
  while ((lmf = MESH_Next_Face(submesh,&idx))) {
    MEnt_Get_AttVal(lmf,l2gatt,0,0,&gmf);
    if (gmf)
      MEnt_Unmark(gmf,mkfid);
  }

  MSTK_FreeMarker(mkvid);
  MSTK_FreeMarker(mkfid);

#else

  MAttrib_Delete(gfidatt);
  MAttrib_Delete(lfidatt);
  MAttrib_Delete(gvidatt);
  MAttrib_Delete(lvidatt);

#endif

  return 1;
}


int MESH_Vol_AddGhost_FN(Mesh_ptr mesh, Mesh_ptr submesh, int part_no, int ring) {
  int i, j, k, l, ir, idx;
  int nvr, nvf, nrf, nre, nrv, nfv, nfe, idx2, found;
  double xyz[3];
  MVertex_ptr gmv, gmv2, lmv, lmv2;
  MEdge_ptr lme, gme;
  MFace_ptr lmf, gmf;
  MRegion_ptr lmr, gmr;
  List_ptr vregions, rfaces, fedges, fverts, redges, rverts;
  List_ptr bverts, bverts2, gghregs, gbverts, gbverts2, ovregions;
  List_ptr lmvlist, lmelist, lmflist, lmrlist;
  MAttrib_ptr g2latt, l2gatt;
  double rval;
  void *pval;

  g2latt = MESH_AttribByName(mesh,"Global2Local");
  l2gatt = MESH_AttribByName(submesh,"Local2Global");


  /* Mark the list of global entities in this submesh */

#ifdef MSTK_USE_MARKERS
  int mkrid = MSTK_GetMarker();  
#else
  MAttrib_ptr gridatt = MAttrib_New(mesh, "temp_gr_id", INT, MREGION); 
  MAttrib_ptr lridatt = MAttrib_New(submesh, "temp_lr_id", INT, MREGION);
#endif

  idx = 0;
  while ((lmr = MESH_Next_Region(submesh,&idx))) {
    MEnt_Get_AttVal(lmr,l2gatt,0,0,&gmr);
#ifdef MSTK_USE_MARKERS
    MEnt_Mark(gmr,mkrid);
#else
    MEnt_Set_AttVal(gmr,gridatt,1,0.0,NULL);
#endif
  }

  /* Collect all the PBOUNDARY vertices of submesh */

#ifdef MSTK_USE_MARKERS
  int mkvid = MSTK_GetMarker();
#else
  MAttrib_ptr gvidatt = MAttrib_New(mesh, "temp_gv_id", INT, MVERTEX);
  MAttrib_ptr lvidatt = MAttrib_New(submesh, "temp_lv_id", INT, MVERTEX);
#endif
  bverts = List_New(10);
  gbverts = List_New(10);
  idx = 0;
  while ((lmv = MESH_Next_Vertex(submesh,&idx))) {
    if (MV_OnParBoundary(lmv)) {
      List_Add(bverts,lmv);

      MEnt_Get_AttVal(lmv,l2gatt,0,0,&gmv);
      List_Add(gbverts,gmv);

#ifdef MSTK_USE_MARKERS
      MEnt_Mark(lmv,mkvid);
      MEnt_Mark(gmv,mkvid);
#else
      MEnt_Set_AttVal(lmv, lvidatt, 1, 0.0, NULL);
      MEnt_Set_AttVal(gmv, gvidatt, 1, 0.0, NULL);
#endif
    }
  }

 
  /* Mark 'N' rings of elements going in from the boundary as overlap elements */

  ovregions = List_New(List_Num_Entries(bverts));

  for (ir = 0; ir < ring; ir++) {

    bverts2 = List_New(List_Num_Entries(bverts));

    idx = 0;
    while ((lmv = List_Next_Entry(bverts,&idx))) {

      vregions = MV_Regions(lmv);
      nvr = List_Num_Entries(vregions);
      for (i = 0; i < nvr; i++) {
        lmr = List_Entry(vregions,i);

        int rmarked;
#ifdef MSTK_USE_MARKERS
        rmarked = MEnt_IsMarked(lmr,mkrid);
#else
        MEnt_Get_AttVal(lmr, lridatt, &rmarked, &rval, &pval);
#endif
        if (!rmarked) {
          List_Add(ovregions,lmr);
          MR_Set_PType(lmr,POVERLAP);
#ifdef MSTK_USE_MARKERS
          MEnt_Mark(lmr,mkrid);
#else
          MEnt_Set_AttVal(lmr, lridatt, 1, 0.0, NULL);
#endif
          rfaces = MR_Faces(lmr);
          nrf = List_Num_Entries(rfaces);
          for (j = 0; j < nrf; j++) {
            lmf = List_Entry(rfaces,j);
            if (!MF_OnParBoundary(lmf))
              MF_Set_PType(lmf,POVERLAP);

            fedges = MF_Edges(lmf,1,0);
            nfe = List_Num_Entries(fedges);
            for (k = 0; k < nfe; k++) {
              lme = List_Entry(fedges,k);
              if (!ME_OnParBoundary(lme))
                ME_Set_PType(lme,POVERLAP);

              int l;
              for (l = 0; l < 2; l++) {
                lmv2 = ME_Vertex(lme,l);
                int vmarked;
#ifdef MSTK_USE_MARKERS
                vmarked = MEnt_IsMarked(lmv2,mkvid);
#else
                MEnt_Get_AttVal(lmv2, lvidatt, &vmarked, &rval, &pval);
#endif
                if (!vmarked) {
                  if (!MV_OnParBoundary(lmv2))
                    MV_Set_PType(lmv2,POVERLAP);
                  List_Add(bverts2,lmv2);
#ifdef MSTK_USE_MARKERS
                  MEnt_Mark(lmv2,mkvid);
#else
                  MEnt_Set_AttVal(lmv2, lvidatt, 1, 0.0, NULL);
#endif
                }
              }

            }
            List_Delete(fedges);
          }

          List_Delete(rfaces);
        }

      }
      List_Delete(vregions);

    }
#ifdef MSTK_USE_MARKERS
    List_Unmark(bverts,mkvid);
#else
    idx = 0;
    while ((lmv = List_Next_Entry(bverts,&idx)))
      MEnt_Set_AttVal(lmv, lvidatt, 0, 0.0, NULL);
#endif
    List_Delete(bverts);
    
    bverts = bverts2;
    
  } /* for (i = 0; i < ring; i++) */

#ifdef MSTK_USE_MARKERS
  List_Unmark(bverts,mkvid);
#else
  idx = 0;
  while ((lmv = List_Next_Entry(bverts,&idx)))
    MEnt_Set_AttVal(lmv, lvidatt, 0, 0.0, NULL);
#endif
  List_Delete(bverts);

#ifdef MSTK_USE_MARKERS
  List_Unmark(ovregions,mkrid);
#else
    idx = 0;
    while ((lmr = List_Next_Entry(ovregions,&idx)))
      MEnt_Set_AttVal(lmr, lridatt, 0, 0.0, NULL);
#endif
  List_Delete(ovregions);


  /* Duplicate 'N' rings of elements around submesh from global mesh */

  MFace_ptr *lrfaces = (MFace_ptr *) malloc(MAXPF3*sizeof(MFace_ptr));
  int *lrfdirs = (int *) malloc(MAXPF3*sizeof(int));
  MEdge_ptr *lfedges = (MEdge_ptr *) malloc(MAXPV2*sizeof(MEdge_ptr));
  int *lfedirs = (int *) malloc(MAXPV2*sizeof(int));

  for (ir = 0; ir < ring; ir++) {
    
    gbverts2 = List_New(List_Num_Entries(gbverts));

    idx = 0; 
    while ((gmv = List_Next_Entry(gbverts,&idx))) {

      vregions = MV_Regions(gmv);
      nvr = List_Num_Entries(vregions);

      for (i = 0; i < nvr; i++) {
	gmr = List_Entry(vregions,i);
        int rmarked;
#ifdef MSTK_USE_MARKERS
        rmarked = MEnt_IsMarked(gmr,mkrid);
#else
        MEnt_Get_AttVal(gmr, gridatt, &rmarked, &rval, &pval);
#endif
	if (!rmarked) { /* not overlapping with submesh */

          /* Found a ghost region that has not yet been processed */
          /* Add to list of ghost regions and duplicate in submesh */

#ifdef MSTK_USE_MARKERS
	  MEnt_Mark(gmr,mkrid);
#else
          MEnt_Set_AttVal(gmr, gridatt, 1, 0.0, NULL);
#endif

          rfaces = MR_Faces(gmr);
          nrf = List_Num_Entries(rfaces);

          for (j = 0; j < nrf; j++) {
            gmf = List_Entry(rfaces,j);
            MEnt_Get_AttVal(gmf,g2latt,0,0,&lmflist);
            
            if (!lmflist)
              MSTK_Report("MESH_AddGhost",
                          "No local face list found with global face",
                          MSTK_FATAL);
            
            idx2 = 0;
            found = 0;
            while ((lmf = List_Next_Entry(lmflist,&idx2))) {
              if (MF_Mesh(lmf) == submesh) {
                found = 1;
                break;
              }
            }
            
            if (!found) {
              lmf = MF_New(submesh);
              MF_Set_GEntID(lmf,MF_GEntID(gmf));
              MF_Set_GEntDim(lmf,MF_GEntDim(gmf));
              MF_Set_PType(lmf,PGHOST);
              
              fedges = MF_Edges(gmf,1,0);
              nfe = List_Num_Entries(fedges);
              for (k = 0; k < nfe; k++) {
                gme = List_Entry(fedges,k);
                MEnt_Get_AttVal(gme,g2latt,0,0,&lmelist);
                
                idx2 = 0;
                found = 0;
                while ((lme = List_Next_Entry(lmelist,&idx2))) {
                  if (ME_Mesh(lme) == submesh) {
                    found = 1;
                    break;
                  }
                }
                
                if (!found) {
                  lme = ME_New(submesh);
                  ME_Set_GEntID(lme,ME_GEntID(gme));
                  ME_Set_GEntDim(lme,ME_GEntDim(gme));
                  ME_Set_PType(lme,PGHOST);
                  
                  for (l = 0; l < 2; l++) {
                    gmv2 = ME_Vertex(gme,l);
                    MEnt_Get_AttVal(gmv2,g2latt,0,0,&lmvlist);
                    
                    found = 0;
                    idx2 = 0;
                    while ((lmv = List_Next_Entry(lmvlist,&idx2))) {
                      if (MV_Mesh(lmv) == submesh) {
                        found = 1;
                        break;
                      }
                    }
                    
                    if (!found) {
                      lmv = MV_New(submesh);
                      MV_Coords(gmv2,xyz);
                      MV_Set_Coords(lmv,xyz);
                      MV_Set_GEntID(lmv,MV_GEntID(gmv2));
                      MV_Set_GEntDim(lmv,MV_GEntDim(gmv2));
                      MV_Set_PType(lmv,PGHOST);
                      MV_Set_GlobalID(lmv,MV_GlobalID(gmv2));
                      MV_Set_MasterParID(lmv,MV_MasterParID(gmv2));
                      MEnt_Set_AttVal(lmv,l2gatt,0,0,gmv2);
                      List_Add(lmvlist,lmv);
                    }
                    
                    ME_Set_Vertex(lme,l,lmv);

                    /* Also, if this is an unmarked global vertex then
                       it represent an outer vertex of the layer of
                       ghost elements. Add it to the gbverts2 list to
                       fetch the next layer of ghost elements */
                    
                    int vmarked;
#ifdef MSTK_USE_MARKERS                    
                    vmarked = MEnt_IsMarked(gmv2,mkvid);
#else
                    MEnt_Get_AttVal(gmv2, gvidatt, &vmarked, &rval, &pval);
#endif
                    if (!vmarked) {
                      List_Add(gbverts2,gmv2);
#ifdef MSTK_USE_MARKERS
                      MEnt_Mark(gmv2,mkvid);
#else
                      MEnt_Set_AttVal(gmv2, gvidatt, 1, 0.0, NULL);
#endif
                    }
                    
                  }
                  ME_Set_GlobalID(lme,ME_GlobalID(gme));
                  ME_Set_MasterParID(lme,ME_MasterParID(gme));
                  MEnt_Set_AttVal(lme,l2gatt,0,0,gme);
                  List_Add(lmelist,lme);
                }
                
                lfedges[k] = lme;
                lfedirs[k] = MF_EdgeDir_i(gmf,k);
              } /* k */
              List_Delete(fedges);

              MF_Set_Edges(lmf,nfe,lfedges,lfedirs);
              MF_Set_GlobalID(lmf,MF_GlobalID(gmf));
              MF_Set_MasterParID(lmf,MF_MasterParID(gmf));
              MEnt_Set_AttVal(lmf,l2gatt,0,0,gmf);
              List_Add(lmflist,lmf);
            }

            lrfaces[j] = lmf;
            lrfdirs[j] = MR_FaceDir_i(gmr,j);

          } /* j */
          
          List_Delete(rfaces);
          
          lmr = MR_New(submesh);
          MR_Set_GEntID(lmr,MR_GEntID(gmr));
          MR_Set_GEntDim(lmr,MR_GEntDim(gmr));
          MR_Set_PType(lmr,PGHOST);
          MR_Set_Faces(lmr,nrf,lrfaces,lrfdirs);
          MR_Set_GlobalID(lmr,MR_GlobalID(gmr));
          MR_Set_MasterParID(lmr,MR_MasterParID(gmr));
          
          MEnt_Set_AttVal(lmr,l2gatt,0,0,gmr);
          MEnt_Get_AttVal(gmr,g2latt,0,0,&lmrlist);
          if (!lmrlist)
            MSTK_Report("MESH_AddGhost","Local region list not found with global region",MSTK_FATAL);
          List_Add(lmrlist,lmr);          
	}

      } /* i */
      List_Delete(vregions);

    }

#ifdef MSTK_USE_MARKERS
    List_Unmark(gbverts,mkvid);
#else
    idx = 0;
    while ((gmv = List_Next_Entry(gbverts,&idx)))
      MEnt_Set_AttVal(gmv, gvidatt, 0, 0.0, NULL);
#endif
    List_Delete(gbverts);
    
#ifdef MSTK_USE_MARKERS
    List_Mark(gbverts2,mkvid);
#else
    idx = 0;
    while ((gmv = List_Next_Entry(gbverts2,&idx)))
      MEnt_Set_AttVal(gmv, gvidatt, 1, 0.0, NULL);
#endif
    gbverts = gbverts2;
  } /* ir */

#ifdef MSTK_USE_MARKERS
  List_Unmark(gbverts,mkvid);
#else
  idx = 0;
  while ((gmv = List_Next_Entry(gbverts,&idx)))
    MEnt_Set_AttVal(gmv, gvidatt, 0, 0.0, NULL);
#endif
  List_Delete(gbverts);

  free(lfedges);
  free(lfedirs);
  free(lrfaces);
  free(lrfdirs);

  /* Unmark all the global faces corresponding to this submesh and its
     ghost layer */

#ifdef MSTK_USE_MARKERS

  idx = 0;
  while ((lmr = MESH_Next_Region(submesh,&idx))) {
    MEnt_Get_AttVal(lmr,l2gatt,0,0,&gmr);
    if (gmr)
      MEnt_Unmark(gmr,mkrid);
  }

  MSTK_FreeMarker(mkvid);
  MSTK_FreeMarker(mkrid);

#else

  MAttrib_Delete(lvidatt);
  MAttrib_Delete(gvidatt);
  MAttrib_Delete(lridatt);
  MAttrib_Delete(gridatt);

#endif

  return 1;
}




int MESH_Surf_AddGhost_R1R2R4(Mesh_ptr mesh, Mesh_ptr submesh, int part_no, int ring) {
  MSTK_Report("MESH_Surf_AddGhost_R1R2R4","Not Implemented",MSTK_FATAL);
  return 0;
}

int MESH_Vol_AddGhost_R1R2(Mesh_ptr mesh, Mesh_ptr submesh, int part_no, int ring) {
  MSTK_Report("MESH_Vol_AddGhost_R1R2","Not Implemented",MSTK_FATAL);
  return 0;
}
  
int MESH_Vol_AddGhost_R4(Mesh_ptr mesh, Mesh_ptr submesh, int part_no, int ring) {
  MSTK_Report("MESH_Vol_AddGhost_R4","Not Implemented",MSTK_FATAL);
  return 0;
}
  
#ifdef __cplusplus
}
#endif

