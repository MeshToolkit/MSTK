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
  int i, j, k, l, ir, idx, mkvid, mkfid;
  int nvf, nfv, nfe, idx2, found;
  double xyz[3];
  MVertex_ptr gmv, gmv2, lmv, lmv2;
  MEdge_ptr lme, gme;
  MFace_ptr lmf, gmf;
  List_ptr vfaces, fverts, fedges;
  List_ptr bverts, bverts2, ovfaces, gbverts, gbverts2;
  List_ptr lmvlist, lmelist, lmflist;
  MAttrib_ptr g2latt, l2gatt;

  g2latt = MESH_AttribByName(mesh,"Global2Local");
  l2gatt = MESH_AttribByName(submesh,"Local2Global");


  /* Mark the list of global elements in this submesh */

  mkfid = MSTK_GetMarker();  
  idx = 0;
  while((lmf = MESH_Next_Face(submesh,&idx))) {
    MEnt_Get_AttVal(lmf,l2gatt,0,0,&gmf);
    MEnt_Mark(gmf,mkfid);
  }

  /* Collect all the PBOUNDARY vertices of submesh */
  mkvid = MSTK_GetMarker();
  bverts = List_New(10);
  gbverts = List_New(10);
  idx = 0;
  while ((lmv = MESH_Next_Vertex(submesh,&idx))) {
    if (MV_OnParBoundary(lmv)) {
      List_Add(bverts,lmv);
      MEnt_Mark(lmv,mkvid);

      MEnt_Get_AttVal(lmv,l2gatt,0,0,&gmv);
      List_Add(gbverts,gmv);
      MEnt_Mark(gmv,mkvid);
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
        
        if (!MEnt_IsMarked(lmf,mkfid)) {        
          List_Add(ovfaces,lmf);
          MEnt_Mark(lmf,mkfid);
          MF_Set_PType(lmf,POVERLAP);          

          fedges = MF_Edges(lmf,1,0);
          nfe = List_Num_Entries(fedges);

          for (k = 0; k < nfe; k++) {
            lme = List_Entry(fedges,k);
            if (!ME_OnParBoundary(lme))
              ME_Set_PType(lme,POVERLAP);

            for (l = 0; l < 2; l++) {
              lmv2 = ME_Vertex(lme,l);
              if (!MEnt_IsMarked(lmv2,mkvid)) {
                if (!ME_OnParBoundary(lmv2))
                  MV_Set_PType(lmv2,POVERLAP);
                List_Add(bverts2,lmv2);
                MEnt_Mark(lmv2,mkvid);
              }
            } /* l */

          } /* k */
          List_Delete(fedges);
        }
      } /* j */
      List_Delete(vfaces);

    }
    List_Unmark(bverts,mkvid);
    List_Delete(bverts);
    
    bverts = bverts2;

  } /* for (i = 0; i < ring; i++) */
  
  List_Unmark(bverts,mkvid);
  List_Delete(bverts);

  List_Unmark(ovfaces,mkfid);
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

	if (!MEnt_IsMarked(gmf,mkfid)) { /* not overlapping with submesh */

          /* Found a ghost face that has not yet been processed */
          /* Add to list of ghost faces and duplicate in submesh */

	  MEnt_Mark(gmf,mkfid);

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
                    
                if (!MEnt_IsMarked(gmv2,mkvid)) {
                  List_Add(gbverts2,gmv2);
                  MEnt_Mark(gmv2,mkvid);
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

    List_Unmark(gbverts,mkvid);    
    List_Delete(gbverts);
    
    List_Mark(gbverts2,mkvid);
    gbverts = gbverts2;
  }
  List_Unmark(gbverts,mkvid);
  List_Delete(gbverts);

  free(lfedges);
  free(lfedirs);


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

  return 1;
}


int MESH_Vol_AddGhost_FN(Mesh_ptr mesh, Mesh_ptr submesh, int part_no, int ring) {
  int i, j, k, l, ir, idx, mkvid, mkrid;
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

  g2latt = MESH_AttribByName(mesh,"Global2Local");
  l2gatt = MESH_AttribByName(submesh,"Local2Global");


  /* Mark the list of global entities in this submesh */

  mkrid = MSTK_GetMarker();  
  idx = 0;
  while((lmr = MESH_Next_Region(submesh,&idx))) {
    MEnt_Get_AttVal(lmr,l2gatt,0,0,&gmr);
    MEnt_Mark(gmr,mkrid);
  }

  /* Collect all the PBOUNDARY vertices of submesh */

  mkvid = MSTK_GetMarker();
  bverts = List_New(10);
  gbverts = List_New(10);
  idx = 0;
  while ((lmv = MESH_Next_Vertex(submesh,&idx))) {
    if (MV_OnParBoundary(lmv)) {
      List_Add(bverts,lmv);
      MEnt_Mark(lmv,mkvid);

      MEnt_Get_AttVal(lmv,l2gatt,0,0,&gmv);
      List_Add(gbverts,gmv);
      MEnt_Mark(gmv,mkvid);
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

        if (!MEnt_IsMarked(lmr,mkrid)) {
          List_Add(ovregions,lmr);
          MEnt_Mark(lmr,mkrid);
          MR_Set_PType(lmr,POVERLAP);

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
                if (!MEnt_IsMarked(lmv2,mkvid)) {
                  if (!MV_OnParBoundary(lmv2))
                    MV_Set_PType(lmv2,POVERLAP);
                  List_Add(bverts2,lmv2);
                  MEnt_Mark(lmv2,mkvid);
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
    List_Unmark(bverts,mkvid);
    List_Delete(bverts);

    bverts = bverts2;

  } /* for (i = 0; i < ring; i++) */

  List_Unmark(bverts,mkvid);
  List_Delete(bverts);

  List_Unmark(ovregions,mkrid);
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

	if (!MEnt_IsMarked(gmr,mkrid)) { /* not overlapping with submesh */

          /* Found a ghost region that has not yet been processed */
          /* Add to list of ghost regions and duplicate in submesh */

	  MEnt_Mark(gmr,mkrid);

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
                    
                    if (!MEnt_IsMarked(gmv2,mkvid)) {
                      List_Add(gbverts2,gmv2);
                      MEnt_Mark(gmv2,mkvid);
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

    List_Unmark(gbverts,mkvid);
    List_Delete(gbverts);
    
    List_Mark(gbverts2,mkvid);
    gbverts = gbverts2;
  } /* ir */
  List_Unmark(gbverts,mkvid);
  List_Delete(gbverts);
  MSTK_FreeMarker(mkvid);

  free(lfedges);
  free(lfedirs);
  free(lrfaces);
  free(lrfdirs);

  /* Unmark all the global faces corresponding to this submesh and its
     ghost layer */

  idx = 0;
  while ((lmr = MESH_Next_Region(submesh,&idx))) {
    MEnt_Get_AttVal(lmr,l2gatt,0,0,&gmr);
    if (gmr)
      MEnt_Unmark(gmr,mkrid);
  }


  MSTK_FreeMarker(mkvid);
  MSTK_FreeMarker(mkrid);

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

