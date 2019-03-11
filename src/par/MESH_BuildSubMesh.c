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
     This function builds a submesh based on overlap elements 
     Author(s): Duo Wang, Rao Garimella
  */

int MESH_BuildSubMesh_Face(Mesh_ptr mesh, Mesh_ptr submesh);
int MESH_BuildSubMesh_Region(Mesh_ptr mesh, Mesh_ptr submesh);


  int MESH_BuildSubMesh(Mesh_ptr mesh, int topodim, Mesh_ptr submesh) {
  if (topodim == 3)
    MESH_BuildSubMesh_Region(mesh, submesh);
  else if(topodim == 2) 
    MESH_BuildSubMesh_Face(mesh, submesh);
  else {
    MSTK_Report("MESH_BuildSubMesh()","only send volume or surface mesh",MSTK_ERROR);
    exit(-1);
  }
  return 1;
}


  /* 
     Send 1-ring Faces to neighbor processors, and receive them 
     First update the parallel adjancy information, 
  */
int MESH_BuildSubMesh_Face(Mesh_ptr mesh, Mesh_ptr submesh) {
  int nfe, j, k;
  MVertex_ptr mv, new_mv;
  MEdge_ptr me, new_me, *fedges;
  MFace_ptr mf, new_mf;
  List_ptr mfedges;
  List_ptr faces, edges, verts;
  int idx, *fedirs;
  double coor[3];
  int *MV_to_list_id, *ME_to_list_id;
  int nv = MESH_Num_Vertices(mesh);
  int ne = MESH_Num_Edges(mesh);
  int nf = MESH_Num_Faces(mesh);

#ifdef MSTK_USE_MARKERS
  int mkvid = MSTK_GetMarker();
  int mkeid = MSTK_GetMarker();
  int mkfid = MSTK_GetMarker();
#else
  MAttrib_ptr mkatt = MAttrib_New(mesh, "mkatt", INT, MALLTYPE);
  double rval;
  void *pval;
#endif

  faces = List_New(10);
  edges = List_New(10);
  verts = List_New(10);

  MV_to_list_id = (int *)malloc((nv)*sizeof(int));
  ME_to_list_id = (int *)malloc((ne)*sizeof(int));
  /* build the 1-ring outside layer send mesh */
  idx = 0;

  fedges = (MEdge_ptr *) malloc(MAXPV2*sizeof(MEdge_ptr));
  fedirs = (int *) malloc(MAXPV2*sizeof(int));
  while( (mf = MESH_Next_Face(mesh, &idx)) ) {
    if(MF_PType(mf) != POVERLAP) continue;

    int f_inlist;
#ifdef MSTK_USE_MARKERS
    f_inlist = MEnt_IsMarked(mf,mkfid);
#else
    MEnt_Get_AttVal(mf, mkatt, &f_inlist, &rval, &pval);
#endif

    if (f_inlist) continue;
    List_Add(faces,mf);

#ifdef MSTK_USE_MARKERS
    MEnt_Mark(mf,mkfid);
#else
    MEnt_Set_AttVal(mf, mkatt, 1, 0.0, NULL);
#endif
    new_mf = MF_New(submesh);  /* add new face and copy information */
    MF_Set_GEntDim(new_mf,MF_GEntDim(mf));
    MF_Set_GEntID(new_mf,MF_GEntID(mf));
    MF_Set_PType(new_mf,MF_PType(mf));
    MF_Set_MasterParID(new_mf,MF_MasterParID(mf));
    MF_Set_GlobalID(new_mf,MF_GlobalID(mf));

    mfedges = MF_Edges(mf,1,0);
    nfe = List_Num_Entries(mfedges);
    for(j = 0; j < nfe; j++) {
      me = List_Entry(mfedges,j);
      
      int e_inlist;
#ifdef MSTK_USE_MARKERS
      e_inlist = MEnt_IsMarked(me,mkeid);
#else
      MEnt_Get_AttVal(me, mkatt, &e_inlist, &rval, &pval);
#endif
      if (e_inlist)
        new_me = MESH_Edge(submesh,ME_to_list_id[ME_ID(me)-1]);
      else {
	new_me = ME_New(submesh);    /* add new edge and copy information */
	ME_Set_GEntDim(new_me,ME_GEntDim(me));
	ME_Set_GEntID(new_me,ME_GEntID(me));
	ME_Set_PType(new_me,ME_PType(me));
	ME_Set_MasterParID(new_me,ME_MasterParID(me));
	ME_Set_GlobalID(new_me,ME_GlobalID(me));

	ME_to_list_id[ME_ID(me)-1] = ME_ID(new_me)-1;
	List_Add(edges,me);
#ifdef MSTK_USE_MARKERS
	MEnt_Mark(me,mkeid);
#else
        MEnt_Set_AttVal(me, mkatt, 1, 0.0, NULL);
#endif
 	for(k = 0; k < 2; k++) {
	  mv = ME_Vertex(me,k);

          int v_inlist;
#ifdef MSTK_USE_MARKERS
	  v_inlist = MEnt_IsMarked(mv,mkvid);
#else
	  MEnt_Get_AttVal(mv, mkatt, &v_inlist, &rval, &pval);
#endif
          if (v_inlist)
            new_mv = MESH_Vertex(submesh,MV_to_list_id[MV_ID(mv)-1]);
	  else {
	    new_mv = MV_New(submesh);  /* add new vertex and copy information */
	    MV_Set_GEntDim(new_mv,MV_GEntDim(mv));
	    MV_Set_GEntID(new_mv,MV_GEntID(mv));
	    MV_Set_PType(new_mv,MV_PType(mv));
	    MV_Set_MasterParID(new_mv,MV_MasterParID(mv));
	    MV_Set_GlobalID(new_mv,MV_GlobalID(mv));
	    MV_Coords(mv,coor);
	    MV_Set_Coords(new_mv,coor);

	    MV_to_list_id[MV_ID(mv)-1] = MV_ID(new_mv)-1;
	    List_Add(verts,mv);
#ifdef MSTK_USE_MARKERS
	    MEnt_Mark(mv,mkvid);
#else
            MEnt_Set_AttVal(mv, mkatt, 1, 0.0, NULL);
#endif
	  }
	  ME_Set_Vertex(new_me,k,new_mv);  /* set edge-vertex */
	}
      }
      
      fedges[j] = new_me;
      fedirs[j] = MF_EdgeDir_i(mf,j) == 1 ? 1 : 0;
    }
    MF_Set_Edges(new_mf,nfe,fedges,fedirs); /* set face-edge */
    List_Delete(mfedges);
  }

#ifdef MSTK_USE_MARKERS
  List_Unmark(faces,mkfid);
  List_Unmark(edges,mkeid);
  List_Unmark(verts,mkvid);
  MSTK_FreeMarker(mkfid);
  MSTK_FreeMarker(mkeid);
  MSTK_FreeMarker(mkvid);
#else
  MAttrib_Delete(mkatt);
#endif

  List_Delete(faces);
  List_Delete(edges);
  List_Delete(verts);
  free(MV_to_list_id);
  free(ME_to_list_id);
  free(fedges);
  free(fedirs);
  return 1;
}

  /* right now assume there are no overlapped regions */

int MESH_BuildSubMesh_Region(Mesh_ptr mesh, Mesh_ptr submesh) {
  int nrf, nfe, i, j, k;
  MVertex_ptr mv, new_mv;
  MEdge_ptr me, new_me, *fedges;
  MFace_ptr mf, new_mf, *rfaces;
  MRegion_ptr mr, new_mr;
  List_ptr mrfaces, mfedges;
  List_ptr faces, edges, verts;
  int idx, *fedirs, *rfdirs;
  int *MV_to_list_id, *ME_to_list_id, *MF_to_list_id;
  double coor[3];
  int nv = MESH_Num_Vertices(mesh);
  int ne = MESH_Num_Edges(mesh);
  int nf = MESH_Num_Faces(mesh);

#ifdef MSTK_USE_MARKERS
  int mkvid = MSTK_GetMarker();
  int mkeid = MSTK_GetMarker();
  int mkfid = MSTK_GetMarker();
#else
  MAttrib_ptr mkatt = MAttrib_New(mesh, "mkatt", INT, MALLTYPE);
  double rval;
  void *pval;
#endif

  faces = List_New(10);
  edges = List_New(10);
  verts = List_New(10);

  MV_to_list_id = (int *)malloc((nv)*sizeof(int));
  ME_to_list_id = (int *)malloc((ne)*sizeof(int));
  MF_to_list_id = (int *)malloc((nf)*sizeof(int));
  /* build the 1-ring outside layer send mesh */
  rfaces = (MFace_ptr *) malloc(MAXPF3*sizeof(MFace_ptr));
  rfdirs = (int *) malloc(MAXPF3*sizeof(int));

  fedges = (MEdge_ptr *) malloc(MAXPV2*sizeof(MEdge_ptr));
  fedirs = (int *) malloc(MAXPV2*sizeof(int));
  idx = 0;  
  while( (mr = MESH_Next_Region(mesh, &idx)) ) {
    if(MR_PType(mr) != POVERLAP) continue;
    new_mr = MR_New(submesh);
    MR_Set_GEntDim(new_mr,MR_GEntDim(mr));
    MR_Set_GEntID(new_mr,MR_GEntID(mr));
    MR_Set_PType(new_mr,MR_PType(mr));
    MR_Set_MasterParID(new_mr,MR_MasterParID(mr));
    MR_Set_GlobalID(new_mr,MR_GlobalID(mr));

    mrfaces = MR_Faces(mr);
    nrf = List_Num_Entries(mrfaces);
    for(i = 0; i < nrf; i++) {
      mf = List_Entry(mrfaces,i);

      int f_inlist;
#ifdef MSTK_USE_MARKERS
      f_inlist = MEnt_IsMarked(mf,mkfid);
#else
      MEnt_Get_AttVal(mf, mkatt, &f_inlist, &rval, &pval);
#endif
      if (f_inlist)
        new_mf = MESH_Face(submesh,MF_to_list_id[MF_ID(mf)-1]);
      else {
	new_mf = MF_New(submesh);  /* add new face and copy information */
	MF_Set_GEntDim(new_mf,MF_GEntDim(mf));
	MF_Set_GEntID(new_mf,MF_GEntID(mf));
	MF_Set_PType(new_mf,MF_PType(mf));
	MF_Set_MasterParID(new_mf,MF_MasterParID(mf));
	MF_Set_GlobalID(new_mf,MF_GlobalID(mf));
	
	List_Add(faces,mf);
#ifdef MSTK_USE_MARKERS
	MEnt_Mark(mf,mkfid);
#else
        MEnt_Set_AttVal(mf, mkatt, 1, 0.0, NULL);
#endif
	MF_to_list_id[MF_ID(mf)-1] = MF_ID(new_mf)-1;
	
	mfedges = MF_Edges(mf,1,0);
	nfe = List_Num_Entries(mfedges);
	for(j = 0; j < nfe; j++) {
	  me = List_Entry(mfedges,j);

          int e_inlist;
#ifdef MSTK_USE_MARKERS
          e_inlist = MEnt_IsMarked(me,mkeid);
#else
          MEnt_Get_AttVal(me, mkatt, &e_inlist, &rval, &pval);
#endif

	  if (e_inlist) 
            new_me = MESH_Edge(submesh,ME_to_list_id[ME_ID(me)-1]);
	  else {
	    new_me = ME_New(submesh);    /* add new edge and copy information */
	    ME_Set_GEntDim(new_me,ME_GEntDim(me));
	    ME_Set_GEntID(new_me,ME_GEntID(me));
	    ME_Set_PType(new_me,ME_PType(me));
	    ME_Set_MasterParID(new_me,ME_MasterParID(me));
	    ME_Set_GlobalID(new_me,ME_GlobalID(me));

	    ME_to_list_id[ME_ID(me)-1] = ME_ID(new_me)-1;
	    List_Add(edges,me);

#ifdef MSTK_USE_MARKERS
	    MEnt_Mark(me,mkeid);
#else
            MEnt_Set_AttVal(me, mkatt, 1, 0.0, NULL);
#endif
	    for(k = 0; k < 2; k++) {
	      mv = ME_Vertex(me,k);

              int v_inlist;
#ifdef MSTK_USE_MARKERS
              v_inlist = MEnt_IsMarked(mv,mkvid);
#else
              MEnt_Get_AttVal(mv, mkatt, &v_inlist, &rval, &pval);
#endif
	      if (v_inlist)
                new_mv = MESH_Vertex(submesh,MV_to_list_id[MV_ID(mv)-1]);
	      else {
		new_mv = MV_New(submesh);  /* add new vertex and copy information */
		MV_Set_GEntDim(new_mv,MV_GEntDim(mv));
		MV_Set_GEntID(new_mv,MV_GEntID(mv));
		MV_Set_PType(new_mv,MV_PType(mv));
		MV_Set_MasterParID(new_mv,MV_MasterParID(mv));
		MV_Set_GlobalID(new_mv,MV_GlobalID(mv));
		MV_Coords(mv,coor);
		MV_Set_Coords(new_mv,coor);
		
		MV_to_list_id[MV_ID(mv)-1] = MV_ID(new_mv)-1;
		List_Add(verts,mv);

#ifdef MSTK_USE_MARKERS
		MEnt_Mark(mv,mkvid);
#else
                MEnt_Set_AttVal(mv, mkatt, 1, 0.0, NULL);
#endif
	      }
	      ME_Set_Vertex(new_me,k,new_mv);  /* set edge-vertex */
	    }
	  }
	  fedges[j] = new_me;
	  fedirs[j] = MF_EdgeDir_i(mf,j) == 1 ? 1 : 0;
	}
	MF_Set_Edges(new_mf,nfe,fedges,fedirs); /* set face-edge */
	List_Delete(mfedges);
      }
      rfaces[i] = new_mf;
      rfdirs[i] = MR_FaceDir_i(mr,i) == 1 ? 1 : 0;
    }
    MR_Set_Faces(new_mr,nrf,rfaces,rfdirs); 
    List_Delete(mrfaces);
  }

#ifdef MSTK_USE_MARKERS
  List_Unmark(faces,mkfid);
  List_Unmark(edges,mkeid);
  List_Unmark(verts,mkvid);
  MSTK_FreeMarker(mkfid);
  MSTK_FreeMarker(mkeid);
  MSTK_FreeMarker(mkvid);
#else
  MAttrib_Delete(mkatt);
#endif

  List_Delete(faces);
  List_Delete(edges);
  List_Delete(verts);
  free(MV_to_list_id);
  free(ME_to_list_id);
  free(MF_to_list_id);
  free(fedges);
  free(fedirs);
  free(rfaces);
  free(rfdirs);
  return 1;
}

#ifdef __cplusplus
}
#endif

