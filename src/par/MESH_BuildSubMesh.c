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


int MESH_BuildSubMesh(Mesh_ptr mesh, Mesh_ptr submesh) {
  int nf, nr;
  RepType rtype;
  /* basic mesh information */
  rtype = MESH_RepType(mesh);
  nf = MESH_Num_Faces(mesh);
  nr = MESH_Num_Regions(mesh);
  if (nr)
    MESH_BuildSubMesh_Region(mesh, submesh);
  else if(nf) 
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
  RepType rtype;
  int nv, ne, nf, nfe, i, j, k;
  MVertex_ptr mv, new_mv;
  MEdge_ptr me, new_me, *fedges;
  MFace_ptr mf, new_mf;
  List_ptr mfedges;
  List_ptr faces, edges, verts;
  int adj, idx, *fedirs;
  int mkvid, mkeid, mkfid;
  double coor[3];
  nv = MESH_Num_Vertices(mesh);
  ne = MESH_Num_Edges(mesh);
  nf = MESH_Num_Faces(mesh);
  mkvid = MSTK_GetMarker();
  mkeid = MSTK_GetMarker();
  mkfid = MSTK_GetMarker();
  faces = List_New(10);
  edges = List_New(10);
  verts = List_New(10);

  int *MV_to_list_id = (int *)MSTK_malloc((nv+1)*sizeof(int));
  int *ME_to_list_id = (int *)MSTK_malloc((ne+1)*sizeof(int));
  //  int *MF_to_list_id = (int *)MSTK_malloc((nf+1)*sizeof(int));
  /* build the 1-ring outside layer send mesh */
  idx = 0;

  fedges = (MEdge_ptr *) malloc(MAXPV2*sizeof(MEdge_ptr));
  fedirs = (int *) malloc(MAXPV2*sizeof(int));
  while(mf = MESH_Next_Face(mesh, &idx)) {
    if(MF_PType(mf) != POVERLAP) continue;
    if(MEnt_IsMarked(mf,mkfid)) continue;
    List_Add(faces,mf);
    MEnt_Mark(mf,mkfid);
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
      if(MEnt_IsMarked(me,mkeid)) new_me = MESH_Edge(submesh,ME_to_list_id[ME_ID(me)]-1);
      else {
	new_me = ME_New(submesh);    /* add new edge and copy information */
	ME_Set_GEntDim(new_me,ME_GEntDim(me));
	ME_Set_GEntID(new_me,ME_GEntID(me));
	ME_Set_PType(new_me,ME_PType(me));
	ME_Set_MasterParID(new_me,ME_MasterParID(me));
	ME_Set_GlobalID(new_me,ME_GlobalID(me));

	ME_to_list_id[ME_ID(me)] = ME_ID(new_me);
	List_Add(edges,me);
	MEnt_Mark(me,mkeid);
 	for(k = 0; k < 2; k++) {
	  mv = ME_Vertex(me,k);
	  if(MEnt_IsMarked(mv,mkvid)) new_mv = MESH_Vertex(submesh,MV_to_list_id[MV_ID(mv)]-1);
	  else {
	    new_mv = MV_New(submesh);  /* add new vertex and copy information */
	    MV_Set_GEntDim(new_mv,MV_GEntDim(mv));
	    MV_Set_GEntID(new_mv,MV_GEntID(mv));
	    MV_Set_PType(new_mv,MV_PType(mv));
	    MV_Set_MasterParID(new_mv,MV_MasterParID(mv));
	    MV_Set_GlobalID(new_mv,MV_GlobalID(mv));
	    MV_Coords(mv,coor);
	    MV_Set_Coords(new_mv,coor);

	    MV_to_list_id[MV_ID(mv)] = MV_ID(new_mv);
	    List_Add(verts,mv);
	    MEnt_Mark(mv,mkvid);
	  }
	  ME_Set_Vertex(new_me,k,new_mv);
	}
      }
      
      fedges[j] = new_me;
      fedirs[j] = MF_EdgeDir_i(mf,j) == 1 ? 1 : 0;
    }
    MF_Set_Edges(new_mf,nfe,fedges,fedirs);
    printf("new face %d with %d edges\n",MF_ID(new_mf),nfe);
  }

  List_Unmark(faces,mkfid);
  List_Unmark(edges,mkeid);
  List_Unmark(verts,mkvid);
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

int MESH_BuildSubMesh_Region(Mesh_ptr mesh, Mesh_ptr submesh) {
  return 1;
}

#ifdef __cplusplus
}
#endif

