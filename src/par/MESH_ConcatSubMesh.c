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
     Send 1-ring Faces to neighbor processors, and receive them 
     First update the parallel adjancy information, 
  */
int MESH_ConcatSubMesh_Face(Mesh_ptr mesh, int num, Mesh_ptr *submeshes) {
  int max_nv, max_ne, max_nf, nfv, nfe, i, j, k, nbv, nbe;
  MVertex_ptr mv, new_mv, sub_mv, *loc_mv;
  MEdge_ptr me, new_me, sub_me, *fedges, *loc_me;
  MFace_ptr mf, new_mf, sub_mf;
  List_ptr mfverts, mfedges;
  List_ptr faces, edges, verts;
  List_ptr boundary_edges, boundary_verts;
  int add_face, idx, *fedirs;
  int mkvid, mkeid, mkfid;
  double coor[3];
  Mesh_ptr submesh;
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
  
  mkvid = MSTK_GetMarker();
  mkeid = MSTK_GetMarker();
  mkfid = MSTK_GetMarker();
  faces = List_New(10);
  edges = List_New(10);          
  verts = List_New(10);
  
  boundary_verts = List_New(10);
  boundary_edges = List_New(10);

  int *MV_to_list_id = (int *)MSTK_malloc(num*max_nv*sizeof(int));
  int *ME_to_list_id = (int *)MSTK_malloc(num*max_ne*sizeof(int));
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
  
  List_Sort(boundary_edges,nbe,sizeof(MEdge_ptr),compareGlobalID);
  List_Sort(boundary_verts,nbv,sizeof(MVertex_ptr),compareGlobalID);
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
	printf("sub_mv global id %d\n", MV_GlobalID(sub_mv));
	if(MEnt_IsMarked(sub_mv,mkvid)) {add_face = 1; continue;} /* if already checked */
	loc_mv = (MVertex_ptr *)bsearch(&sub_mv,
					boundary_verts,
					nbv,
					sizeof(MVertex_ptr),
					compareGlobalID);
	if(loc_mv) {
	  printf("found vertex global id %d\n",MV_GlobalID(*loc_mv));
	  add_face = 1; 
	  MV_to_list_id[i*max_nv+MV_ID(sub_mv)-1] = MV_ID(*loc_mv)-1;
	  List_Add(verts,sub_mv); 
	  MEnt_Mark(sub_mv,mkvid);
	}
      }
      /* pre store edges from submesh that is already in mesh */
      mfedges = MF_Edges(sub_mf,1,0);
      nfe = List_Num_Entries(mfedges);
      for(i = 0; i < nfe; i++) {
	sub_me = List_Entry(mfedges,i);
	if(MEnt_IsMarked(sub_me,mkeid)) {add_face = 1; continue;} /* if already checked */
	loc_me = (MVertex_ptr *)bsearch(&sub_me,
					boundary_edges,
					nbe,
					sizeof(MEdge_ptr),
					compareGlobalID);
	if(loc_me) {
	  add_face = 1; 
	  MV_to_list_id[i*max_ne+MV_ID(sub_me)-1] = ME_ID(*loc_me)-1;
	  List_Add(edges,sub_me); 
	MEnt_Mark(sub_me,mkeid);
	}
      }
      List_Delete(mfverts);
      List_Delete(mfedges);
      if(!add_face) continue;
      
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
	if(MEnt_IsMarked(sub_me,mkeid)) new_me = MESH_Edge(mesh,ME_to_list_id[i*max_ne+ME_ID(sub_me)-1]);
	else {
	  new_me = ME_New(mesh);    /* add new edge and copy information */
	  ME_Set_GEntDim(new_me,ME_GEntDim(sub_me));
	  ME_Set_GEntID(new_me,ME_GEntID(sub_me));
	  ME_Set_PType(new_me,PGHOST);
	  ME_Set_MasterParID(new_me,ME_MasterParID(sub_me));
	  ME_Set_GlobalID(new_me,ME_GlobalID(sub_me));
	  
	  ME_to_list_id[ME_ID(sub_me)-1] = ME_ID(new_me)-1;
	  List_Add(edges,sub_me);
	  MEnt_Mark(sub_me,mkeid);
	  for(k = 0; k < 2; k++) {
	    sub_mv = ME_Vertex(sub_me,k);
	    if(MEnt_IsMarked(mv,mkvid)) new_mv = MESH_Vertex(mesh,MV_to_list_id[i*max_nv+MV_ID(sub_mv)-1]);
	    else {
	      new_mv = MV_New(mesh);  /* add new vertex and copy information */
	      MV_Set_GEntDim(new_mv,MV_GEntDim(mv));
	      MV_Set_GEntID(new_mv,MV_GEntID(mv));
	      MV_Set_PType(new_mv,PGHOST);
	      MV_Set_MasterParID(new_mv,MV_MasterParID(mv));
	      MV_Set_GlobalID(new_mv,MV_GlobalID(mv));
	      MV_Coords(mv,coor);
	      MV_Set_Coords(new_mv,coor);
	      
	      MV_to_list_id[MV_ID(sub_mv)-1] = MV_ID(new_mv)-1;
	      List_Add(verts,sub_mv);
	      MEnt_Mark(sub_mv,mkvid);
	    }
	    ME_Set_Vertex(new_me,k,new_mv);  /* set edge-vertex */
	  }
	}
	fedges[j] = new_me;
	fedirs[j] = MF_EdgeDir_i(sub_mf,j) == 1 ? 1 : 0;
      }
      MF_Set_Edges(new_mf,nfe,fedges,fedirs); /* set face-edge */
      List_Delete(mfedges);
    }
  }

  List_Unmark(faces,mkfid);
  List_Unmark(edges,mkeid);
  List_Unmark(verts,mkvid);
  List_Delete(boundary_edges);
  List_Delete(boundary_verts);
  List_Delete(faces);
  List_Delete(faces);
  List_Delete(edges);
  List_Delete(verts);
  MSTK_free(MV_to_list_id);
  MSTK_free(ME_to_list_id);
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

