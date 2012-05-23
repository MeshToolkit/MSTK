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

     Author(s): Duo Wang, Rao Garimella
  */

int MESH_ConcatSubMesh_Face(Mesh_ptr mesh, int num, Mesh_ptr *submeshes);
int MESH_ConcatSubMesh_Region(Mesh_ptr mesh, int num, Mesh_ptr *submeshes);


int MESH_ConcatSubMesh(Mesh_ptr mesh, int num, Mesh_ptr *submeshes) {
  int nf, nr;
  RepType rtype;
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

  /* check if l has the same entity as m, based on global ID */
MEntity_ptr entity_on_list(MEntity_ptr m, List_ptr *l) {
  int i, num = List_Num_Entries(*l);
  for(i = 0; i < num; i++)
    if(MEnt_GlobalID(m) == MEnt_GlobalID(List_Entry(*l,i)))
      return List_Entry(*l,i);
  return NULL;
}

int MESH_ConcatSubMesh_Face(Mesh_ptr mesh, int num, Mesh_ptr *submeshes) {
  int max_nv, max_ne, max_nf, nfv, nfe, i, j, k, nbv, nbe;
  MVertex_ptr mv, new_mv, sub_mv;
  MEdge_ptr me, new_me, sub_me, *fedges;
  MFace_ptr new_mf, sub_mf;
  List_ptr mfverts, mfedges;
  List_ptr edges, verts, added_edges, added_verts;
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

  mkvid2 = MSTK_GetMarker(); /* mark the vertices on mesh boundary, used to decide whether to add a face  */
  mkeid2 = MSTK_GetMarker(); /* mark the edges in on mesh boundary, used to decide whether to add a face  */

  mkfid = MSTK_GetMarker();
  edges = List_New(10);          
  verts = List_New(10);
  added_edges = List_New(10);          
  added_verts = List_New(10);
  
  boundary_verts = List_New(10);
  boundary_edges = List_New(10);

  MV_to_list_id = (int *)MSTK_malloc(num*max_nv*sizeof(int));
  ME_to_list_id = (int *)MSTK_malloc(num*max_ne*sizeof(int));

  
  fedges = (MEdge_ptr *) malloc(MAXPV2*sizeof(MEdge_ptr));
  fedirs = (int *) malloc(MAXPV2*sizeof(int));

  /* collect boundary edges and vertices */
  idx = 0; nbe = 0;
  while( (me = MESH_Next_Edge(mesh,&idx)) ) 
    if(ME_PType(me)) {
      List_Add(boundary_edges,me);
      nbe++;
    }
  idx = 0; nbv = 0;
  while( (mv = MESH_Next_Vertex(mesh,&idx)) )
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

  
  for(i = 0; i < num; i++) {
    submesh = submeshes[i];
    idx = 0;
    while ( (sub_mf = MESH_Next_Face(submesh, &idx)) ) {
      add_face = 0;
      /* pre store vertices from submesh that is already in mesh */
      mfverts = MF_Vertices(sub_mf,1,0);
      nfv = List_Num_Entries(mfverts);
      for(j = 0; j < nfv; j++) {
	sub_mv = List_Entry(mfverts,j);
	if(MEnt_IsMarked(sub_mv,mkvid2)) {
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
	  mv = List_Entry(boundary_verts,iloc); 
	  /* here set the ghost vertex property, only necessary when the input submeshes are not consistent */
	  if(MV_PType(mv) == PGHOST && MV_PType(sub_mv) != PGHOST) {
	    MV_Set_GEntDim(mv,MV_GEntDim(sub_mv));
	    MV_Set_GEntID(mv,MV_GEntID(sub_mv));
	  }

	  MV_to_list_id[i*max_nv+MV_ID(sub_mv)-1] = MV_ID(mv)-1;
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
	if(MEnt_IsMarked(sub_me,mkeid2)) {
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
	  me = List_Entry(boundary_edges,iloc); 
	  /* here set the ghost edge property, only necessary when the input submeshes are not consistent */
	  if(ME_PType(me) == PGHOST && ME_PType(sub_me) != PGHOST) {
	    ME_Set_GEntDim(me,ME_GEntDim(sub_me));
	    ME_Set_GEntID(me,ME_GEntID(sub_me));
	  }
	  ME_to_list_id[i*max_ne+ME_ID(sub_me)-1] = ME_ID(me)-1;
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
	new_me = NULL;
	if(MEnt_IsMarked(sub_me,mkeid)) /* first check if it is already in mesh */
	  new_me = MESH_Edge(mesh,ME_to_list_id[i*max_ne+ME_ID(sub_me)-1]); 
	else 
	  new_me = (MEdge_ptr)entity_on_list(sub_me,&added_edges); /* check if it is already added */
	if(new_me)
	  if(MV_GlobalID(ME_Vertex(new_me,0)) != MV_GlobalID(ME_Vertex(sub_me,0)))
	    fedirs[j] = 1 - fedirs[j];  /* if the edge dir is not the same, reverse the edge dir */
	
	if(!new_me)  {                 /* if this is really a new edge */
	  new_me = ME_New(mesh);      /* add new edge and copy information */
	  ME_Set_GEntDim(new_me,ME_GEntDim(sub_me));
	  ME_Set_GEntID(new_me,ME_GEntID(sub_me));
	  ME_Set_PType(new_me,PGHOST);
	  ME_Set_MasterParID(new_me,ME_MasterParID(sub_me));
	  ME_Set_GlobalID(new_me,ME_GlobalID(sub_me));
	  
	  ME_to_list_id[i*max_ne+ME_ID(sub_me)-1] = ME_ID(new_me)-1;
	  List_Add(edges,sub_me);
	  List_Add(added_edges,new_me);
	  MEnt_Mark(sub_me,mkeid);
	  for(k = 0; k < 2; k++) {
	    sub_mv = ME_Vertex(sub_me,k);
	    new_mv = NULL;
	    if(MEnt_IsMarked(sub_mv,mkvid)) 
	      new_mv = MESH_Vertex(mesh,MV_to_list_id[i*max_nv+MV_ID(sub_mv)-1]);
	    else
	      new_mv = (MVertex_ptr)entity_on_list(sub_mv,&added_verts);
	    
	    if(!new_mv) {
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
	      List_Add(added_verts,new_mv);
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

  List_Unmark(edges,mkeid);
  List_Unmark(verts,mkvid);
  List_Unmark(edges,mkeid2);
  List_Unmark(verts,mkvid2);
  List_Delete(boundary_edges);
  List_Delete(boundary_verts);
  List_Delete(added_edges);
  List_Delete(added_verts);

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
  int max_nv, max_ne, max_nf, max_nr, nrf, nre, nrv, nfe, h, i, j, k, nbv, nbe, nbf;
  MVertex_ptr mv, new_mv, sub_mv;
  MEdge_ptr me, new_me, sub_me, *fedges;
  MFace_ptr mf, new_mf, sub_mf, *rfaces;
  MRegion_ptr new_mr, sub_mr;
  List_ptr mrfaces, mredges, mrverts, mfedges;
  List_ptr faces, edges, verts, added_faces, added_edges, added_verts;
  List_ptr boundary_faces, boundary_edges, boundary_verts;
  int add_region, idx, *rfdirs, *fedirs, global_id, iloc, *loc;
  int mkvid, mkeid, mkfid, mkvid2, mkeid2, mkfid2;
  double coor[3];
  Mesh_ptr submesh;
  int *MF_to_list_id, *MV_to_list_id, *ME_to_list_id, *MV_global_id, *ME_global_id, *MF_global_id;  
  /* get maximum verts, edges faces and regions number */
  max_nv = 0; max_ne = 0; max_nf = 0; max_nr = 0;
  for(i = 0; i < num; i++) {
    submesh = submeshes[i];
    if(max_nv < MESH_Num_Vertices(submesh))
      max_nv = MESH_Num_Vertices(submesh);
    if(max_ne < MESH_Num_Edges(submesh))
      max_ne = MESH_Num_Edges(submesh);
    if(max_nf < MESH_Num_Faces(submesh))
      max_nf = MESH_Num_Faces(submesh);
    if(max_nr < MESH_Num_Regions(submesh))
      max_nr = MESH_Num_Regions(submesh);
  }
  
  mkvid = MSTK_GetMarker(); /* mark the vertices in submeshes that is in mesh or already added into mesh */
  mkeid = MSTK_GetMarker(); /* mark the edges in submeshes that is in mesh or already added into mesh */
  mkfid = MSTK_GetMarker(); /* mark the faces in submeshes that is in mesh or already added into mesh */

  mkvid2 = MSTK_GetMarker(); /* mark the vertices on mesh boundary, used to decide whether to add a face  */
  mkeid2 = MSTK_GetMarker(); /* mark the edges on mesh boundary, used to decide whether to add a face  */
  mkfid2 = MSTK_GetMarker(); /* mark the faces on mesh boundary, used to decide whether to add a face  */

  faces = List_New(10);   added_faces = List_New(10);  boundary_verts = List_New(10);        
  edges = List_New(10);   added_edges = List_New(10);  boundary_edges = List_New(10);
  verts = List_New(10);   added_verts = List_New(10);  boundary_faces = List_New(10);

  MV_to_list_id = (int *)MSTK_malloc(num*max_nv*sizeof(int));
  ME_to_list_id = (int *)MSTK_malloc(num*max_ne*sizeof(int));
  MF_to_list_id = (int *)MSTK_malloc(num*max_nf*sizeof(int));

  rfaces = (MFace_ptr *) malloc(MAXPF3*sizeof(MFace_ptr));
  rfdirs = (int *) malloc(MAXPF3*sizeof(int));
  fedges = (MEdge_ptr *) malloc(MAXPV2*sizeof(MEdge_ptr));
  fedirs = (int *) malloc(MAXPV2*sizeof(int));

  /* collect boundary faces, edges and vertices */
  idx = 0; nbf = 0;
  while( (mf = MESH_Next_Face(mesh,&idx)) ) 
    if(MF_PType(mf)) {
      List_Add(boundary_faces,mf);
      nbf++;
    }
  idx = 0; nbe = 0;
  while( (me = MESH_Next_Edge(mesh,&idx)) ) 
    if(ME_PType(me)) {
      List_Add(boundary_edges,me);
      nbe++;
    }
  idx = 0; nbv = 0;
  while( (mv = MESH_Next_Vertex(mesh,&idx)) )
    if(MV_PType(mv)) {
      List_Add(boundary_verts,mv);
      nbv++;
    }
  /* sort based on global ID */
  List_Sort(boundary_faces,nbf,sizeof(MFace_ptr),compareGlobalID);
  List_Sort(boundary_edges,nbe,sizeof(MEdge_ptr),compareGlobalID);
  List_Sort(boundary_verts,nbv,sizeof(MVertex_ptr),compareGlobalID);

  MV_global_id = (int *)MSTK_malloc(nbv*sizeof(int));
  ME_global_id = (int *)MSTK_malloc(nbe*sizeof(int));
  MF_global_id = (int *)MSTK_malloc(nbf*sizeof(int));

  /* store them in array for binary search */
  idx = 0;
  for(i = 0; i < List_Num_Entries(boundary_faces); i++) {
    mf = List_Entry(boundary_faces,i);
    MF_global_id[idx++] = MF_GlobalID(mf);
  }
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

  for(h = 0; h < num; h++) {
    submesh = submeshes[h];
    idx = 0;
    while ( (sub_mr = MESH_Next_Region(submesh, &idx)) ) {
	add_region = 0;
	/* pre store faces from submesh that is already in mesh */
	mrfaces = MR_Faces(sub_mr);
	nrf = List_Num_Entries(mrfaces);
	for(j = 0; j < nrf; j++) {
	  sub_mf = List_Entry(mrfaces,j);
	  if(MEnt_IsMarked(sub_mf,mkfid2)) {
	    add_region = 1; 
	    continue;
	  }
	  global_id = MF_GlobalID(sub_mf);
	  loc = (int *)bsearch(&global_id,
			       MF_global_id,
			       nbf,
			       sizeof(int),
			       compareINT);
	  if(loc) {
	    add_region = 1; 
	    iloc = loc - MF_global_id;
	    mf = List_Entry(boundary_faces,iloc); 
	    /* here set the ghost edge property, only necessary when the input submeshes are not consistent */
	    if(MF_PType(mf) == PGHOST && MF_PType(sub_mf) != PGHOST) {
	      MF_Set_GEntDim(mf,MF_GEntDim(sub_mf));
	      MF_Set_GEntID(mf,MF_GEntID(sub_mf));
	    }
	    MF_to_list_id[h*max_nf+MF_ID(sub_mf)-1] = MF_ID(mf)-1;
	    List_Add(faces,sub_mf); 
	    MEnt_Mark(sub_mf,mkfid);
	    MEnt_Mark(sub_mf,mkfid2);
	  }
	}
	/* pre store edges from submesh that is already in mesh */
	mredges = MR_Edges(sub_mr);
	nre = List_Num_Entries(mredges);
	for(j = 0; j < nre; j++) {
	  sub_me = List_Entry(mredges,j);
	  if(MEnt_IsMarked(sub_me,mkeid2)) {
	    add_region = 1; 
	    continue;
	  }
	  global_id = ME_GlobalID(sub_me);
	  loc = (int *)bsearch(&global_id,
			       ME_global_id,
			       nbe,
			       sizeof(int),
			       compareINT);
	  if(loc) {
	    add_region = 1; 
	    iloc = loc - ME_global_id;
	    me = List_Entry(boundary_edges,iloc); 
	    /* here set the ghost edge property, only necessary when the input submeshes are not consistent */
	    if(ME_PType(me) == PGHOST && ME_PType(sub_me) != PGHOST) {
	      ME_Set_GEntDim(me,ME_GEntDim(sub_me));
	      ME_Set_GEntID(me,ME_GEntID(sub_me));
	    }
	    ME_to_list_id[h*max_ne+ME_ID(sub_me)-1] = ME_ID(me)-1;
	    List_Add(edges,sub_me); 
	    MEnt_Mark(sub_me,mkeid);
	    MEnt_Mark(sub_me,mkeid2);
	  }
	}
	/* pre store vertices from submesh that is already in mesh */
	mrverts = MR_Vertices(sub_mr);
	nrv = List_Num_Entries(mrverts);
	for(j = 0; j < nrv; j++) {
	  sub_mv = List_Entry(mrverts,j);
	  if(MEnt_IsMarked(sub_mv,mkvid2)) {
	    add_region = 1; 
	    continue;
	  }
	  global_id = MV_GlobalID(sub_mv);
	  loc = (int *)bsearch(&global_id,
			       MV_global_id,
			       nbv,
			       sizeof(int),
			       compareINT);
	  if(loc) {
	    add_region = 1; 
	    iloc = loc - MV_global_id;
	    mv = List_Entry(boundary_verts,iloc); 
	    /* here set the ghost vertex property, only necessary when the input submeshes are not consistent */
	    if(MV_PType(mv) == PGHOST && MV_PType(sub_mv) != PGHOST) {
	      MV_Set_GEntDim(mv,MV_GEntDim(sub_mv));
	      MV_Set_GEntID(mv,MV_GEntID(sub_mv));
	    }

	    MV_to_list_id[h*max_nv+MV_ID(sub_mv)-1] = MV_ID(mv)-1;
	    List_Add(verts,sub_mv); 
	    MEnt_Mark(sub_mv,mkvid);
	    MEnt_Mark(sub_mv,mkvid2);

	  }
	}
	
	List_Delete(mrverts);
	List_Delete(mredges);
	List_Delete(mrfaces);
	
	if(!add_region) {continue;}
        
	new_mr = MR_New(mesh);                  /* add region */
	MR_Set_GEntDim(new_mr,MR_GEntDim(sub_mr));
	MR_Set_GEntID(new_mr,MR_GEntID(sub_mr));
	MR_Set_PType(new_mr,PGHOST);
	MR_Set_MasterParID(new_mr,MR_MasterParID(sub_mr));
	MR_Set_GlobalID(new_mr,MR_GlobalID(sub_mr));
	
	mrfaces = MR_Faces(sub_mr);
	nrf = List_Num_Entries(mrfaces);
	for(i = 0; i < nrf; i++) {
	  sub_mf = List_Entry(mrfaces,i);
	  rfdirs[i] = MR_FaceDir_i(sub_mr,i) == 1 ? 1 : 0;
	  new_mf = NULL;
	  if(MEnt_IsMarked(sub_mf,mkfid)) /* first check if it is already in mesh */
	    new_mf = MESH_Face(mesh,MF_to_list_id[h*max_nf+MF_ID(sub_mf)-1]); 
	  else 
	    new_mf = (MFace_ptr)entity_on_list(sub_mf,&added_faces); /* check if it is already added */
	  
	  if(!new_mf) {
	    new_mf = MF_New(mesh); /* add face */
	    MF_Set_GEntDim(new_mf,MF_GEntDim(sub_mf));
	    MF_Set_GEntID(new_mf,MF_GEntID(sub_mf));
	    MF_Set_PType(new_mf,PGHOST);
	    MF_Set_MasterParID(new_mf,MF_MasterParID(sub_mf));
	    MF_Set_GlobalID(new_mf,MF_GlobalID(sub_mf));
	    
	    MF_to_list_id[h*max_nf+MF_ID(sub_mf)-1] = MF_ID(new_mf)-1;
	    List_Add(faces,sub_mf);
	    List_Add(added_faces,new_mf);
	    MEnt_Mark(sub_mf,mkfid);
	    
	    mfedges = MF_Edges(sub_mf,1,0);
	    nfe = List_Num_Entries(mfedges);
	    for(j = 0; j < nfe; j++) {
	      sub_me = List_Entry(mfedges,j);
	      fedirs[j] = MF_EdgeDir_i(sub_mf,j) == 1 ? 1 : 0;
	      new_me = NULL;
	      if(MEnt_IsMarked(sub_me,mkeid)) /* first check if it is already in mesh */
		new_me = MESH_Edge(mesh,ME_to_list_id[h*max_ne+ME_ID(sub_me)-1]); 
	      else 
		new_me = (MEdge_ptr)entity_on_list(sub_me,&added_edges); /* check if it is already added */
	      if(new_me)
		if(MV_GlobalID(ME_Vertex(new_me,0)) != MV_GlobalID(ME_Vertex(sub_me,0)))
		  fedirs[j] = 1 - fedirs[j];  /* if the edge dir is not the same, reverse the edge dir */
	      
	      if(!new_me)  {                 /* if this is really a new edge */
		new_me = ME_New(mesh);      /* add new edge and copy information */
		ME_Set_GEntDim(new_me,ME_GEntDim(sub_me));
		ME_Set_GEntID(new_me,ME_GEntID(sub_me));
		ME_Set_PType(new_me,PGHOST);
		ME_Set_MasterParID(new_me,ME_MasterParID(sub_me));
		ME_Set_GlobalID(new_me,ME_GlobalID(sub_me));
		
		ME_to_list_id[h*max_ne+ME_ID(sub_me)-1] = ME_ID(new_me)-1;
		List_Add(edges,sub_me);
		List_Add(added_edges,new_me);
		MEnt_Mark(sub_me,mkeid);
		for(k = 0; k < 2; k++) {
		  sub_mv = ME_Vertex(sub_me,k);
		  new_mv = NULL;
		  if(MEnt_IsMarked(sub_mv,mkvid)) 
		    new_mv = MESH_Vertex(mesh,MV_to_list_id[h*max_nv+MV_ID(sub_mv)-1]);
		  else
		    new_mv = (MVertex_ptr)entity_on_list(sub_mv,&added_verts);
		  
		  if(!new_mv) {
		    new_mv = MV_New(mesh);  /* add new vertex and copy information */
		    MV_Set_GEntDim(new_mv,MV_GEntDim(sub_mv));
		    MV_Set_GEntID(new_mv,MV_GEntID(sub_mv));
		    MV_Set_PType(new_mv,PGHOST);
		    MV_Set_MasterParID(new_mv,MV_MasterParID(sub_mv));
		    MV_Set_GlobalID(new_mv,MV_GlobalID(sub_mv));
		    MV_Coords(sub_mv,coor);
		    MV_Set_Coords(new_mv,coor);
		    
		    MV_to_list_id[h*max_nv+MV_ID(sub_mv)-1] = MV_ID(new_mv)-1;
		    List_Add(verts,sub_mv);
		    List_Add(added_verts,new_mv);
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
	  rfaces[i] = new_mf;
	}
	MR_Set_Faces(new_mr,nrf,rfaces,rfdirs); /* set region-face */
	List_Delete(mrfaces);
      }
  }      
      
  List_Unmark(faces,mkfid);
  List_Unmark(edges,mkeid);
  List_Unmark(verts,mkvid);
  List_Unmark(faces,mkeid2);
  List_Unmark(edges,mkeid2);
  List_Unmark(verts,mkvid2);
  List_Delete(boundary_faces);
  List_Delete(boundary_edges);
  List_Delete(boundary_verts);
  List_Delete(added_faces);
  List_Delete(added_edges);
  List_Delete(added_verts);

  List_Delete(faces);
  List_Delete(edges);
  List_Delete(verts);
  MSTK_free(MV_to_list_id);
  MSTK_free(ME_to_list_id);
  MSTK_free(MF_to_list_id);
  MSTK_free(MV_global_id);
  MSTK_free(ME_global_id);
  MSTK_free(MF_global_id);
  MSTK_free(fedges);
  MSTK_free(fedirs);
  MSTK_free(rfaces);
  MSTK_free(rfdirs);

  return 1;
}

#ifdef __cplusplus
}
#endif

