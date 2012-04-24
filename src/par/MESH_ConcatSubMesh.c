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
     This function concatenates submesh into mesh, based on vertex global ID
     Author(s): Duo Wang, Rao Garimella
  */

int MESH_ConcatSubMesh_Face(Mesh_ptr mesh, int num, Mesh_ptr *submeshes);
int MESH_ConcatSubMesh_Region(Mesh_ptr mesh, int num, Mesh_ptr *submeshes);


int MESH_ConcatSubMesh(Mesh_ptr mesh, int num, Mesh_ptr *submeshes) {
  int nf, nr;
  RepType rtype;
  /* basic mesh information */
  rtype = MESH_RepType(mesh);
  nf = MESH_Num_Faces(mesh);
  nr = MESH_Num_Regions(mesh);
  if (nr)
    MESH_ConcatSubMesh_Region(mesh, num, submeshes);
  else if(nf) 
    MESH_ConcatSubMesh_Face(mesh, num, submeshes);
  else {
    MSTK_Report("MESH_ConcatSubMesh()","only send volume or surface mesh",MSTK_ERROR);
    exit(-1);
  }
  return 1;
}


  /* 
     Concatenate submeshes into mesh 
  */
int MESH_ConcatSubMesh_Face(Mesh_ptr mesh, int num, Mesh_ptr *submeshes) {
  int max_nv, max_ne, max_nf, nfv, nfe, i, j, k, nbv, nbe;
  MVertex_ptr mv, new_mv, sub_mv;
  MEdge_ptr me, new_me, sub_me, *fedges;
  MFace_ptr mf, new_mf, sub_mf;
  List_ptr mfverts, mfedges;
  List_ptr faces, edges, verts;
  List_ptr boundary_edges, boundary_verts;
  int add_face, idx, *fedirs, global_id, iloc, *loc;
  int mkvid, mkeid, mkvid2, mkeid2, mkfid;
  double coor[3];
  Mesh_ptr submesh;
  int *MV_to_list_id, *ME_to_list_id, *MV_global_id, *ME_global_id;  
  /* get maximum verts, edges and faces number */
  max_nv = 0; max_ne = 0; max_nf = 0;
  for(i = 0; i < num; i++) {
    submesh =submeshes[i];
    if(max_nv < MESH_Num_Vertices(submesh))
      max_nv = MESH_Num_Vertices(submesh);
    if(max_ne < MESH_Num_Edges(submesh))
      max_ne = MESH_Num_Edges(submesh);
    if(max_nf < MESH_Num_Faces(submesh))
      max_nf = MESH_Num_Faces(submesh);
  }
  
  mkvid = MSTK_GetMarker(); /* mark the vertices in submeshes that is in mesh or already added into mesh */
  mkeid = MSTK_GetMarker(); /* mark the edges in submeshes that is in mesh or already added into mesh */

  mkvid2 = MSTK_GetMarker(); /* mark the vertices in submeshes that is in mesh, used to decide whether to add a face  */
  mkeid2 = MSTK_GetMarker(); /* mark the edges in submeshes that is in mesh, used to decide whether to add a face  */

  mkfid = MSTK_GetMarker();
  faces = List_New(10);
  edges = List_New(10);          
  verts = List_New(10);
  
  boundary_verts = List_New(10);
  boundary_edges = List_New(10);

  MV_to_list_id = (int *)MSTK_malloc(num*max_nv*sizeof(int));
  ME_to_list_id = (int *)MSTK_malloc(num*max_ne*sizeof(int));

  
  fedges = (MEdge_ptr *) malloc(MAXPV2*sizeof(MEdge_ptr));
  fedirs = (int *) malloc(MAXPV2*sizeof(int));

  /* collect boundary edges and vertices */
  idx = 0; nbe = 0;
  while(me = MESH_Next_Edge(mesh,&idx)) 
    if(ME_PType(me)) {
      List_Add(boundary_edges,me);
      nbe++;
    }
  idx = 0; nbv = 0;
  while(mv = MESH_Next_Vertex(mesh,&idx))
    if(MV_PType(mv)) {
      List_Add(boundary_verts,mv);
      nbv++;
    }
  /* sort based on global ID */
  List_Sort(boundary_edges,nbe,sizeof(MEdge_ptr),compareGlobalID);
  List_Sort(boundary_verts,nbv,sizeof(MVertex_ptr),compareGlobalID);

  MV_global_id = (int *)MSTK_malloc(nbv*sizeof(int));
  ME_global_id = (int *)MSTK_malloc(nbe*sizeof(int));

  /* store them in array for binary search */
  idx = 0;
  for(i = 0; i < List_Num_Entries(boundary_edges); i++) {
    me = List_Entry(boundary_edges,i);
    ME_global_id[idx++] = ME_GlobalID(me);
  }
  idx = 0;
  for(i = 0; i < List_Num_Entries(boundary_verts); i++) {
    mv = List_Entry(boundary_verts,i);
    MV_global_id[idx++] = MV_GlobalID(mv);
  }

  /*
  printf("max_nv %d, max_ne %d\n",max_nv,max_ne);
  printf("nbe %d, nbv %d\n",nbe,nbv);
  printf("global id: \n");
  for(i = 0; i < List_Num_Entries(boundary_verts); i++)
    printf(" %d ", MV_GlobalID(List_Entry(boundary_verts,i)));
  printf("\n");
  */
  
  for(i = 0; i < num; i++) {
    submesh = submeshes[i];
    idx = 0;
    while (sub_mf = MESH_Next_Face(submesh, &idx)) {
      add_face = 0;
      /* pre store vertices from submesh that is already in mesh */
      mfverts = MF_Vertices(sub_mf,1,0);
      nfv = List_Num_Entries(mfverts);
      for(j = 0; j < nfv; j++) {
	sub_mv = List_Entry(mfverts,j);
	if(MEnt_IsMarked(sub_mv,mkvid)) {
	  if(MEnt_IsMarked(sub_mv,mkvid2))
	    add_face = 1; 
	  continue;
	}
	global_id = MV_GlobalID(sub_mv);
	loc = (int *)bsearch(&global_id,
			     MV_global_id,
			     nbv,
			     sizeof(int),
			     compareINT);
	if(loc) {
	  add_face = 1; 
	  iloc = loc - MV_global_id;
	  MV_to_list_id[i*max_nv+MV_ID(sub_mv)-1] = MV_ID(List_Entry(boundary_verts,iloc))-1;
	  List_Add(verts,sub_mv); 
	  MEnt_Mark(sub_mv,mkvid);
	  MEnt_Mark(sub_mv,mkvid2);
	}
      }
      /* pre store edges from submesh that is already in mesh */
      mfedges = MF_Edges(sub_mf,1,0);
      nfe = List_Num_Entries(mfedges);
      for(j = 0; j < nfe; j++) {
	sub_me = List_Entry(mfedges,j);
	if(MEnt_IsMarked(sub_me,mkeid)) {
	  if(MEnt_IsMarked(sub_me,mkeid2))
	    add_face = 1; 
	  continue;
	}
	global_id = ME_GlobalID(sub_me);
	loc = (int *)bsearch(&global_id,
			     ME_global_id,
			     nbe,
			     sizeof(int),
			     compareINT);
	if(loc) {
	  add_face = 1; 
	  iloc = loc - ME_global_id;
	  ME_to_list_id[i*max_ne+ME_ID(sub_me)-1] = ME_ID(List_Entry(boundary_edges,iloc))-1;
	  List_Add(edges,sub_me); 
	  MEnt_Mark(sub_me,mkeid);
	  MEnt_Mark(sub_me,mkeid2);
	}
      }
      List_Delete(mfverts);
      List_Delete(mfedges);
      if(!add_face) { continue;}
      
      new_mf = MF_New(mesh); /* add face */
      MF_Set_GEntDim(new_mf,MF_GEntDim(sub_mf));
      MF_Set_GEntID(new_mf,MF_GEntID(sub_mf));
      MF_Set_PType(new_mf,PGHOST);
      MF_Set_MasterParID(new_mf,MF_MasterParID(sub_mf));
      MF_Set_GlobalID(new_mf,MF_GlobalID(sub_mf));
      
      mfedges = MF_Edges(sub_mf,1,0);
      nfe = List_Num_Entries(mfedges);
      for(j = 0; j < nfe; j++) {
	sub_me = List_Entry(mfedges,j);
	fedirs[j] = MF_EdgeDir_i(sub_mf,j) == 1 ? 1 : 0;
	if(MEnt_IsMarked(sub_me,mkeid)) {
	  new_me = MESH_Edge(mesh,ME_to_list_id[i*max_ne+ME_ID(sub_me)-1]);
	  if(MV_GlobalID(ME_Vertex(new_me,0)) != MV_GlobalID(ME_Vertex(sub_me,0)))
	    fedirs[j] = 1 - fedirs[j];  /* if the edge dir is not the same, reverse the edge dir */
	}
	else {
	  new_me = ME_New(mesh);    /* add new edge and copy information */
	  ME_Set_GEntDim(new_me,ME_GEntDim(sub_me));
	  ME_Set_GEntID(new_me,ME_GEntID(sub_me));
	  ME_Set_PType(new_me,PGHOST);
	  ME_Set_MasterParID(new_me,ME_MasterParID(sub_me));
	  ME_Set_GlobalID(new_me,ME_GlobalID(sub_me));
	  
	  ME_to_list_id[i*max_ne+ME_ID(sub_me)-1] = ME_ID(new_me)-1;
	  List_Add(edges,sub_me);
	  MEnt_Mark(sub_me,mkeid);
	  for(k = 0; k < 2; k++) {
	    sub_mv = ME_Vertex(sub_me,k);
	    if(MEnt_IsMarked(mv,mkvid)) new_mv = MESH_Vertex(mesh,MV_to_list_id[i*max_nv+MV_ID(sub_mv)-1]);
	    else {
	      new_mv = MV_New(mesh);  /* add new vertex and copy information */
	      MV_Set_GEntDim(new_mv,MV_GEntDim(sub_mv));
	      MV_Set_GEntID(new_mv,MV_GEntID(sub_mv));
	      MV_Set_PType(new_mv,PGHOST);
	      MV_Set_MasterParID(new_mv,MV_MasterParID(sub_mv));
	      MV_Set_GlobalID(new_mv,MV_GlobalID(sub_mv));
	      MV_Coords(sub_mv,coor);
	      MV_Set_Coords(new_mv,coor);
	      
	      MV_to_list_id[i*max_nv+MV_ID(sub_mv)-1] = MV_ID(new_mv)-1;
	      List_Add(verts,sub_mv);
	      MEnt_Mark(sub_mv,mkvid);
	    }
	    ME_Set_Vertex(new_me,k,new_mv);  /* set edge-vertex */
	  }
	}								
	fedges[j] = new_me;
      }
      MF_Set_Edges(new_mf,nfe,fedges,fedirs); /* set face-edge */
      List_Delete(mfedges);
    }
  }

  List_Unmark(faces,mkfid);
  List_Unmark(edges,mkeid);
  List_Unmark(faces,mkvid);
  List_Unmark(edges,mkeid2);
  List_Unmark(verts,mkvid2);
  List_Delete(boundary_edges);
  List_Delete(boundary_verts);
  List_Delete(faces);
  List_Delete(edges);
  List_Delete(verts);
  MSTK_free(MV_to_list_id);
  MSTK_free(ME_to_list_id);
  MSTK_free(MV_global_id);
  MSTK_free(ME_global_id);
  MSTK_free(fedges);
  MSTK_free(fedirs);
  return 1;
}

  /* right now assume there are no overlapped regions */

int MESH_ConcatSubMesh_Region(Mesh_ptr mesh, int num, Mesh_ptr *submeshes) {
  return 1;
}

#ifdef __cplusplus
}
#endif

