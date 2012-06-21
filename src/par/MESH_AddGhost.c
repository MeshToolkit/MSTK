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

     must call PMESH_BuildPBoundary() first

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
  int nf, nr, idx;
  MVertex_ptr lmv;
  RepType rtype;

  /* if 0 ring, only mark processor boundary vertex */
  if (ring == 0) {
    idx = 0;
    while((lmv = MESH_Next_Vertex(submesh,&idx))) {
      if(MV_PType(lmv) == PBOUNDARY) {
	if(MV_MasterParID(lmv) == part_no)
	  MV_Set_PType(lmv,POVERLAP);
	else
	  MV_Set_PType(lmv,PGHOST);
      }
    }
    return 1;
  }

  /* basic mesh information */
  nf = MESH_Num_Faces(submesh);
  nr = MESH_Num_Regions(submesh);
  rtype = MESH_RepType(submesh);

  if (nr)
    (*MESH_Vol_AddGhost_jmp[rtype])(mesh,submesh,part_no,ring);
  else if(nf) 
    (*MESH_Surf_AddGhost_jmp[rtype])(mesh,submesh,part_no,ring);
  else {
    MSTK_Report("MESH_AddGhost()","This is not a valid mstk file",MSTK_ERROR);
    exit(-1);
  }

  MESH_Build_GhostLists(submesh);

  return 1;
}



int MESH_Surf_AddGhost_FN(Mesh_ptr mesh, Mesh_ptr submesh, int part_no, int ring) {
  int i, j, k, idx, vertex_id, edge_id, face_id, mkvid, mkvid2, mkfid;
  int nvf, nfv, nfe, idx2, found;
  double xyz[3];
  MVertex_ptr gmv, gmv2, lmv, lmv2;
  MEdge_ptr lme, gme;
  MFace_ptr lmf, gmf;
  List_ptr vfaces, fverts, fedges, gghfaces, gbverts, gbverts2;
  List_ptr lmvlist, lmelist, lmflist;
  MAttrib_ptr g2latt, l2gatt;

  g2latt = MESH_AttribByName(mesh,"Global2Local");
  l2gatt = MESH_AttribByName(submesh,"Local2Global");


  /* Mark the list of global entities in this submesh */

  mkvid = MSTK_GetMarker();
  idx = 0;
  while((lmv = MESH_Next_Vertex(submesh,&idx))) {
    MEnt_Get_AttVal(lmv,l2gatt,0,0,&gmv);
    MEnt_Mark(gmv,mkvid);
  }


  mkfid = MSTK_GetMarker();  
  idx = 0;
  while((lmf = MESH_Next_Face(submesh,&idx))) {
    MEnt_Get_AttVal(lmf,l2gatt,0,0,&gmf);
    MEnt_Mark(gmf,mkfid);
  }


  /* Loop through PBOUNDARY vertices of submesh. Mark any connected
     faces as POVERLAP. Also mark any subentities of these faces not
     marked as PBOUNDARY as POVERLAP */

  /* THIS CODE WON'T WORK FOR MORE THAN 1 LAYER OF GHOST ELEMENTS */

  idx = 0;
  while ((lmv = MESH_Next_Vertex(submesh,&idx))) {
    if (MV_PType(lmv) != PBOUNDARY) continue;

    vfaces = MV_Faces(lmv);
    nvf = List_Num_Entries(vfaces);
    for (i = 0; i < nvf; i++) {
      lmf = List_Entry(vfaces,i);
      MF_Set_PType(lmf,POVERLAP);

      fedges = MF_Edges(lmf,1,0);
      nfe = List_Num_Entries(fedges);
      for (j = 0; j < nfe; j++) {
	lme = List_Entry(fedges,j);
	if (ME_PType(lme) != PBOUNDARY) 
	  ME_Set_PType(lme,POVERLAP);
      }
      List_Delete(fedges);

      fverts = MF_Vertices(lmf,1,0);
      nfv = nfe;
      for (j = 0; j < nfv; j++) {
	lmv2 = List_Entry(fverts,j);
	if (MV_PType(lmv2) != PBOUNDARY)
	  MV_Set_PType(lmv2,POVERLAP);
      }
      List_Delete(fverts);
    }
    List_Delete(vfaces);
  }


  /* Get a 'N' rings of elements around submesh from global mesh */
  if (ring > 1) {
    MSTK_Report("MESH_AddGhost","Code above is not setup for more than one layer of ghost elements",MSTK_FATAL);
  }

  /* Simultaneously make a list of global vertices on partition boundary */
  
  mkvid2 = MSTK_GetMarker();

  gbverts = List_New(10);
  idx = 0;
  while ((lmv = MESH_Next_Vertex(submesh,&idx))) {
    if (MV_PType(lmv) == PBOUNDARY) {
      MEnt_Get_AttVal(lmv,l2gatt,0,0,&gmv);
      List_Add(gbverts,gmv);
      MEnt_Mark(gmv,mkvid2);
    }
  }


  gghfaces = List_New(List_Num_Entries(gbverts)); /* global ghost faces */

  for (i = 0; i < ring; i++) {
    
    gbverts2 = List_New(List_Num_Entries(gbverts));

    idx = 0; 
    while ((gmv = List_Next_Entry(gbverts,&idx))) {

      vfaces = MV_Faces(gmv);
      nvf = List_Num_Entries(vfaces);

      for (j = 0; j < nvf; j++) {
	gmf = List_Entry(vfaces,j);

	if (!MEnt_IsMarked(gmf,mkfid)) {
	  List_Add(gghfaces,gmf);
	  MEnt_Mark(gmf,mkfid);

	  fverts = MF_Vertices(gmf,1,0);
	  nfv = List_Num_Entries(fverts);
	  for (k = 0; k < nfv; k++) {
	    gmv2 = List_Entry(fverts,k);
	    if (!MEnt_IsMarked(gmv2,mkvid2)) {
	      List_Add(gbverts2,gmv2);
	      MEnt_Mark(gmv2,mkvid2);
	    }
	  }
	  List_Delete(fverts);
	}
      }
      List_Delete(vfaces);

    }

    List_Unmark(gbverts,mkvid2);
    List_Delete(gbverts);
    
    List_Mark(gbverts2,mkvid2);
    gbverts = gbverts2;
  }
  List_Unmark(gbverts,mkvid2);
  List_Delete(gbverts);
  MSTK_FreeMarker(mkvid2);

  List_Unmark(gghfaces,mkfid);
  
 
  
  /* Duplicate these elements and their subentities as ghost entities
     for the submesh */

  MEdge_ptr *lfedges = (MEdge_ptr *) malloc(MAXPV2*sizeof(MEdge_ptr));
  int *lfedirs = (int *) malloc(MAXPV2*sizeof(int));

  idx = 0;
  while ((gmf = List_Next_Entry(gghfaces,&idx))) {

    fverts = MF_Vertices(gmf,1,0);
    nfv = List_Num_Entries(fverts);

    for (i = 0; i < nfv; i++) {
      gmv = List_Entry(fverts,i);
      MEnt_Get_AttVal(gmv,g2latt,0,0,&lmvlist);

      if (!lmvlist)
	MSTK_Report("MESH_AddGhost","No local vertex list found with global vertex",MSTK_FATAL);

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
	MV_Coords(gmv,xyz);
	MV_Set_Coords(lmv,xyz);
	MV_Set_GEntID(lmv,MV_GEntID(gmv));
	MV_Set_GEntDim(lmv,MV_GEntDim(gmv));
	MV_Set_PType(lmv,PGHOST);
	MV_Set_GlobalID(lmv,MV_GlobalID(gmv));
	MV_Set_MasterParID(lmv,MV_MasterParID(gmv));
	MEnt_Set_AttVal(lmv,l2gatt,0,0,gmv);
	List_Add(lmvlist,lmv);
      }
    }

    List_Delete(fverts);

    fedges = MF_Edges(gmf,1,0);
    nfe = List_Num_Entries(fedges);

    for (i = 0; i < nfe; i++) {
      gme = List_Entry(fedges,i);
      MEnt_Get_AttVal(gme,g2latt,0,0,&lmelist);

      if (!lmelist)
	MSTK_Report("MESH_AddGhost","No local edge list found with gloval edge list",MSTK_FATAL);

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
	for (j = 0; j < 2; j++) {
	  gmv = ME_Vertex(gme,j);
	  MEnt_Get_AttVal(gmv,g2latt,0,0,&lmvlist);

	  idx2 = 0;
	  found = 0;
	  while ((lmv = List_Next_Entry(lmvlist,&idx2))) {
	    if (MV_Mesh(lmv) == submesh) {
	      found = 1;
	      break;
	    }
	  }

	  if (!lmv)
	    MSTK_Report("MESH_Surf_AddGhost_FN","Missing ghost vertex",MSTK_ERROR);

	  ME_Set_Vertex(lme,j,lmv);
	}
	ME_Set_GlobalID(lme,ME_GlobalID(gme));
	ME_Set_MasterParID(lme,ME_MasterParID(gme));
	MEnt_Set_AttVal(lme,l2gatt,0,0,gme);
	List_Add(lmelist,lme);
      }
      lfedges[i] = lme;
      lfedirs[i] = MF_EdgeDir_i(gmf,i);	
    }

    List_Delete(fedges);


    lmf = MF_New(submesh);
    MF_Set_GEntID(lmf,MF_GEntID(gmf));
    MF_Set_GEntDim(lmf,MF_GEntDim(gmf));
    MF_Set_PType(lmf,PGHOST);
    MF_Set_Edges(lmf,nfe,lfedges,lfedirs);
    MF_Set_GlobalID(lmf,MF_GlobalID(gmf));
    MF_Set_MasterParID(lmf,MF_MasterParID(gmf));

    MEnt_Set_AttVal(lmf,l2gatt,0,0,gmf);
    MEnt_Get_AttVal(gmf,g2latt,0,0,&lmflist);
    if (!lmflist)
      MSTK_Report("MESH_AddGhost","Local face list not found with global face",MSTK_FATAL);
    List_Add(lmflist,lmf);
  }
  List_Delete(gghfaces);

  free(lfedges);
  free(lfedirs);


  /* Mark all entities marked as PBOUNDARY as POVERLAP or
     PGHOST */

  idx = 0;
  while ((lmv = MESH_Next_Vertex(submesh,&idx))) {
    if (MV_PType(lmv) == PBOUNDARY) {
      if (MV_MasterParID(lmv) == part_no)
	MV_Set_PType(lmv,POVERLAP);
      else
	MV_Set_PType(lmv,PGHOST);
    }
  }

  idx = 0;
  while ((lme = MESH_Next_Edge(submesh,&idx))) {
    if (ME_PType(lme) == PBOUNDARY) {
      if (ME_MasterParID(lme) == part_no)
	ME_Set_PType(lme,POVERLAP);
      else
	ME_Set_PType(lme,PGHOST);
    }
  }



  idx = 0;
  while ((lmv = MESH_Next_Vertex(submesh,&idx))) {
    MEnt_Get_AttVal(lmv,l2gatt,0,0,&gmv);
    if (gmv)
      MEnt_Unmark(gmv,mkvid);
  }

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
  int i, j, k, idx, vertex_id, edge_id, face_id, region_id, mkvid, mkvid2, mkrid;
  int nvr, nvf, nrf, nre, nrv, nfv, nfe, idx2, found;
  double xyz[3];
  MVertex_ptr gmv, gmv2, lmv, lmv2;
  MEdge_ptr lme, gme;
  MFace_ptr lmf, gmf;
  MRegion_ptr lmr, gmr;
  List_ptr vregions, rfaces, fedges, fverts, redges, rverts;
  List_ptr gghregs, gbverts, gbverts2, lmvlist, lmelist, lmflist, lmrlist;
  MAttrib_ptr g2latt, l2gatt;

  g2latt = MESH_AttribByName(mesh,"Global2Local");
  l2gatt = MESH_AttribByName(submesh,"Local2Global");


  /* Mark the list of global entities in this submesh */

  mkvid = MSTK_GetMarker();
  idx = 0;
  while((lmv = MESH_Next_Vertex(submesh,&idx))) {
    MEnt_Get_AttVal(lmv,l2gatt,0,0,&gmv);
    MEnt_Mark(gmv,mkvid);
  }

  mkrid = MSTK_GetMarker();  
  idx = 0;
  while((lmr = MESH_Next_Region(submesh,&idx))) {
    MEnt_Get_AttVal(lmr,l2gatt,0,0,&gmr);
    MEnt_Mark(gmr,mkrid);
  }


  /* Loop through PBOUNDARY vertices of submesh. Mark any connected
     regions as POVERLAP. Also mark any subentities of these regions not
     marked as PBOUNDARY as POVERLAP */

  idx = 0;
  while ((lmv = MESH_Next_Vertex(submesh,&idx))) {
    if (MV_PType(lmv) != PBOUNDARY) continue;

    vregions = MV_Regions(lmv);
    nvr = List_Num_Entries(vregions);
    for (i = 0; i < nvr; i++) {
      lmr = List_Entry(vregions,i);
      MR_Set_PType(lmr,POVERLAP);

      rfaces = MR_Faces(lmr);
      nrf = List_Num_Entries(rfaces);
      for (j = 0; j < nrf; j++) {
	lmf = List_Entry(rfaces,j);
	if (MF_PType(lmf) != PBOUNDARY) 
	  MF_Set_PType(lmf,POVERLAP);
      }
      List_Delete(rfaces);

      redges = MR_Edges(lmr);
      nre = List_Num_Entries(redges);
      for (j = 0; j < nre; j++) {
	lme = List_Entry(redges,j);
	if (ME_PType(lmr) != PBOUNDARY)
	  MR_Set_PType(lmr,POVERLAP);
      }
      List_Delete(redges);

      rverts = MR_Vertices(lmr);
      nrv = List_Num_Entries(rverts);
      for (j = 0; j < nrv; j++) {
	lmv2 = List_Entry(rverts,j);
	if (MV_PType(lmv2) != PBOUNDARY)
	  MV_Set_PType(lmv2,POVERLAP);
      }
      List_Delete(rverts);
    }
    List_Delete(vregions);
  }


  /* Get a 'N' rings of elements around submesh from global mesh */

  /* Simultaneously make a list of global vertices on partition boundary */
  
  mkvid2 = MSTK_GetMarker();

  gbverts = List_New(10);
  idx = 0;
  while ((lmv = MESH_Next_Vertex(submesh,&idx))) {
    if (MV_PType(lmv) == PBOUNDARY) {
      MEnt_Get_AttVal(lmv,l2gatt,0,0,&gmv);
      List_Add(gbverts,gmv);
      MEnt_Mark(gmv,mkvid2);
    }
  }


  gghregs = List_New(List_Num_Entries(gbverts)); /* global ghost regions */

  for (i = 0; i < ring; i++) {
    
    gbverts2 = List_New(List_Num_Entries(gbverts));

    idx = 0; 
    while ((gmv = List_Next_Entry(gbverts,&idx))) {

      vregions = MV_Regions(gmv);
      nvr = List_Num_Entries(vregions);

      for (j = 0; j < nvr; j++) {
	gmr = List_Entry(vregions,j);

	if (!MEnt_IsMarked(gmr,mkrid)) {
	  List_Add(gghregs,gmr);
	  MEnt_Mark(gmr,mkrid);

	  rverts = MR_Vertices(gmr);
	  nrv = List_Num_Entries(rverts);
	  for (k = 0; k < nrv; k++) {
	    gmv2 = List_Entry(rverts,k);
	    if (!MEnt_IsMarked(gmv2,mkvid2)) {
	      List_Add(gbverts2,gmv2);
	      MEnt_Mark(gmv2,mkvid2);
	    }
	  }
	  List_Delete(rverts);
	}
      }
      List_Delete(vregions);

    }

    List_Unmark(gbverts,mkvid2);
    List_Delete(gbverts);
    
    List_Mark(gbverts2,mkvid2);
    gbverts = gbverts2;
  }
  List_Unmark(gbverts,mkvid2);
  List_Delete(gbverts);
  MSTK_FreeMarker(mkvid2);

  List_Unmark(gghregs,mkrid);
  
 
  
  /* Duplicate these elements and their subentities as ghost entities
     for the submesh */

  MFace_ptr *lrfaces = (MFace_ptr *) malloc(MAXPF3*sizeof(MFace_ptr));
  int *lrfdirs = (int *) malloc(MAXPF3*sizeof(int));
  MEdge_ptr *lfedges = (MEdge_ptr *) malloc(MAXPV2*sizeof(MEdge_ptr));
  int *lfedirs = (int *) malloc(MAXPV2*sizeof(int));

  idx = 0;
  while ((gmr = List_Next_Entry(gghregs,&idx))) {

    rverts = MR_Vertices(gmr);
    nrv = List_Num_Entries(rverts);

    for (i = 0; i < nrv; i++) {
      gmv = List_Entry(rverts,i);
      MEnt_Get_AttVal(gmv,g2latt,0,0,&lmvlist);

      if (!lmvlist)
	MSTK_Report("MESH_AddGhost","No local vertex list with global vertexx",MSTK_FATAL);

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
	MV_Coords(gmv,xyz);
	MV_Set_Coords(lmv,xyz);
	MV_Set_GEntID(lmv,MV_GEntID(gmv));
	MV_Set_GEntDim(lmv,MV_GEntDim(gmv));
	MV_Set_PType(lmv,PGHOST);
	MV_Set_GlobalID(lmv,MV_GlobalID(gmv));
	MV_Set_MasterParID(lmv,MV_MasterParID(gmv));
	MEnt_Set_AttVal(lmv,l2gatt,0,0,gmv);
	List_Add(lmvlist,lmv);
      }
    }
    List_Delete(rverts);

    redges = MR_Edges(gmr);
    nre = List_Num_Entries(redges);

    for (i = 0; i < nre; i++) {
      gme = List_Entry(redges,i);
      MEnt_Get_AttVal(gme,g2latt,0,0,&lmelist);
      
      if (!lmelist)
	MSTK_Report("MESH_AddGhost","No local edge list found with global edge",MSTK_FATAL);

      idx2 = 0;
      found = 0;
      while ((lme = List_Next_Entry(lmelist,&idx2))) {
	if (MV_Mesh(lme) == submesh) {
	  found = 1;
	  break;
	}
      }

      if (!found) {
	lme = ME_New(submesh);
	ME_Set_GEntID(lme,ME_GEntID(gme));
	ME_Set_GEntDim(lme,ME_GEntDim(gme));
	ME_Set_PType(lme,PGHOST);
	for (j = 0; j < 2; j++) {
	  gmv = ME_Vertex(gme,j);
	  MEnt_Get_AttVal(gmv,g2latt,0,0,&lmvlist);

	  found = 0;
	  idx2 = 0;
	  while ((lmv = List_Next_Entry(lmvlist,&idx2))) {
	    if (MV_Mesh(lmv) == submesh) {
	      found = 1;
	      break;
	    }
	  }

	  if (!lmv)
	    MSTK_Report("MESH_Surf_AddGhost_FN","Missing ghost vertex",MSTK_ERROR);

	  ME_Set_Vertex(lme,j,lmv);
	}
	ME_Set_GlobalID(lme,ME_GlobalID(gme));
	ME_Set_MasterParID(lme,ME_MasterParID(gme));
	MEnt_Set_AttVal(lme,l2gatt,0,0,gme);
	List_Add(lmelist,lme);
      }
    }
    List_Delete(redges);

    rfaces = MR_Faces(gmr);
    nrf = List_Num_Entries(rfaces);

    for (i = 0; i < nrf; i++) {
      gmf = List_Entry(rfaces,i);
      MEnt_Get_AttVal(gmf,g2latt,0,0,&lmflist);
      
      if (!lmflist)
	MSTK_Report("MESH_AddGhost","No local face list found with global face",MSTK_FATAL);

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
	for (j = 0; j < nfe; j++) {
	  gme = List_Entry(fedges,j);
	  MEnt_Get_AttVal(gme,g2latt,0,0,&lmelist);

	  idx2 = 0;
	  found = 0;
	  while ((lme = List_Next_Entry(lmelist,&idx2))) {
	    if (ME_Mesh(lme) == submesh) {
	      found = 1;
	      break;
	    }
	  }

	  if (!lme) 
	    MSTK_Report("MESH_Surf_AddGhost_FN","Cannot find local edge of global edge",MSTK_FATAL);
	  
	  lfedges[j] = lme;
	  lfedirs[j] = MF_EdgeDir_i(gmf,j);
	}
	MF_Set_Edges(lmf,nfe,lfedges,lfedirs);
	MF_Set_GlobalID(lmf,MF_GlobalID(gmf));
	MF_Set_MasterParID(lmf,MF_MasterParID(gmf));
	MEnt_Set_AttVal(lmf,l2gatt,0,0,gmf);
	List_Add(lmflist,lmf);
      }
      lrfaces[i] = lmf;
      lrfdirs[i] = MR_FaceDir_i(gmr,i);
    }


    lmr = MR_New(submesh);
    MR_Set_GEntID(lmr,MR_GEntID(gmf));
    MR_Set_GEntDim(lmr,MR_GEntDim(gmf));
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
  List_Delete(gghregs);

  free(lfedges);
  free(lfedirs);
  free(lrfaces);
  free(lrfdirs);


  /* Mark all entities marked as PBOUNDARY as POVERLAP or
     PGHOST */

  idx = 0;
  while ((lmv = MESH_Next_Vertex(submesh,&idx))) {
    if (MV_PType(lmv) == PBOUNDARY) {
      if (MV_MasterParID(lmv) == part_no)
	MV_Set_PType(lmv,POVERLAP);
      else
	MV_Set_PType(lmv,PGHOST);
    }
  }

  idx = 0;
  while ((lme = MESH_Next_Edge(submesh,&idx))) {
    if (ME_PType(lme) == PBOUNDARY) {
      if (ME_MasterParID(lme) == part_no)
	ME_Set_PType(lme,POVERLAP);
      else
	ME_Set_PType(lme,PGHOST);
    }
  }

  idx = 0;
  while ((lmf = MESH_Next_Face(submesh,&idx))) {
    if (MF_PType(lmf) == PBOUNDARY) {
      if (MF_MasterParID(lmf) == part_no)
	MF_Set_PType(lmf,POVERLAP);
      else
	MF_Set_PType(lmf,PGHOST);
    }
  }



  /* Remove the connections between local and global entities */

  idx = 0;
  while ((lmv = MESH_Next_Vertex(submesh,&idx))) {
    MEnt_Get_AttVal(lmv,l2gatt,0,0,&gmv);
    if (gmv)
      MEnt_Unmark(gmv,mkvid);
  }

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
}

int MESH_Vol_AddGhost_R1R2(Mesh_ptr mesh, Mesh_ptr submesh, int part_no, int ring) {
  MSTK_Report("MESH_Vol_AddGhost_R1R2","Not Implemented",MSTK_FATAL);
}
  
int MESH_Vol_AddGhost_R4(Mesh_ptr mesh, Mesh_ptr submesh, int part_no, int ring) {
  MSTK_Report("MESH_Vol_AddGhost_R4","Not Implemented",MSTK_FATAL);
}
  
#ifdef __cplusplus
}
#endif

