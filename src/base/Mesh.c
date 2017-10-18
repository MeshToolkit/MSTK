#define _H_Mesh_Private

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Mesh.h"
#include "MSTK.h"
#include "MSTK_private.h"


#ifdef __cplusplus
extern "C" {
#endif

  RepType MESH_rtype[5] = {F1, F4, R1, R2, R4};
  char MESH_rtype_str[5][3] = {"F1\0","F4\0","R1\0","R2\0","R4\0"};

Mesh_ptr MESH_New(RepType type) {
  Mesh_ptr newmesh;

  newmesh = (Mesh_ptr) malloc(sizeof(Mesh));

  newmesh->reptype = type;
  newmesh->nv = newmesh->ne = newmesh->nf = newmesh->nr = 0;
  newmesh->mvertex = (List_ptr) NULL;
  newmesh->medge = (List_ptr) NULL;
  newmesh->mface = (List_ptr) NULL;
  newmesh->mregion = (List_ptr) NULL;

#ifdef MSTK_HAVE_MPI
  newmesh->mypartn  = 0;
  newmesh->numpartns = 0;
  newmesh->ghvertex = (List_ptr) NULL;
  newmesh->ghedge = (List_ptr) NULL;
  newmesh->ghface = (List_ptr) NULL;
  newmesh->ghregion = (List_ptr) NULL;
  newmesh->ovvertex = (List_ptr) NULL;
  newmesh->ovedge = (List_ptr) NULL;
  newmesh->ovface = (List_ptr) NULL;
  newmesh->ovregion = (List_ptr) NULL;
  newmesh->par_adj_current = 0;
  newmesh->par_adj_flags = NULL;
  newmesh->par_recv_info = NULL;

  newmesh->max_ghvid = newmesh->max_gheid = newmesh->max_ghfid = newmesh->max_ghrid = 0;
  newmesh->gid_sorted_mvlist = newmesh->gid_sorted_melist =
    newmesh->gid_sorted_mflist = newmesh->gid_sorted_mrlist = NULL;
#endif

  newmesh->geom = (GModel_ptr) NULL;
  newmesh->AttribList = (List_ptr) NULL;
  newmesh->MSetList = (List_ptr) NULL;
  newmesh->autolock = 0;

  if (type >= R1 && type <= R4) newmesh->hedge = Hash_New(0, 1);
  else newmesh->hedge = (Hash_ptr) NULL;
 
  if (type >= R1 && type <= R2) newmesh->hface = Hash_New(0, 1);
  else newmesh->hface = (Hash_ptr) NULL;

  newmesh->max_vid = newmesh->max_eid = newmesh->max_fid = newmesh->max_rid = 0;


  return newmesh;
}

void MESH_Delete(Mesh_ptr mesh) {
  int i, nv, ne, nf, nr;
  MVertex_ptr mv, ghv;
  MEdge_ptr me, ghe;
  MFace_ptr mf, ghf;
  MRegion_ptr mr, ghr;
  MAttrib_ptr attrib;
  MSet_ptr mset;

#ifdef DEBUGHIGH
  if (mesh->hedge) {
    Hash_Print(mesh->hedge);
  }
  if (mesh->hface) {
    Hash_Print(mesh->hface);
  }
#endif 

#ifdef MSTK_HAVE_MPI

  if(mesh->ovvertex)
    List_Delete(mesh->ovvertex);
  if(mesh->ovedge)
    List_Delete(mesh->ovedge);
  if(mesh->ovface)
    List_Delete(mesh->ovface);
  if(mesh->ovregion)
    List_Delete(mesh->ovregion);
  if (mesh->par_adj_flags)
    free(mesh->par_adj_flags);
  if (mesh->par_recv_info)
    free(mesh->par_recv_info);
  if (mesh->gid_sorted_mvlist)
    List_Delete(mesh->gid_sorted_mvlist);
  if (mesh->gid_sorted_melist)
    List_Delete(mesh->gid_sorted_melist);
  if (mesh->gid_sorted_mflist)
    List_Delete(mesh->gid_sorted_mflist);
  if (mesh->gid_sorted_mrlist)
    List_Delete(mesh->gid_sorted_mrlist);

#endif

  if (mesh->mregion) {
    nr = mesh->nr;
    i = 0;
    while ((mr = List_Next_Entry(mesh->mregion,&i))) {
      MR_Destroy_For_MESH_Delete(mr);
    }
    List_Delete(mesh->mregion);
  }
  if (mesh->hface) {
    Hash_Delete(mesh->hface);
    if (mesh->mface) {
      /* Hash_Delete already cleaned all implicit face entities */
      List_Delete(mesh->mface);
    }
  } else {
    if (mesh->mface) {
      nf = mesh->nf;
      i = 0;
      while ((mf = List_Next_Entry(mesh->mface,&i))) {
	MF_Destroy_For_MESH_Delete(mf);
      }
      List_Delete(mesh->mface);
    }
  }

  if (mesh->hedge) {
    Hash_Delete(mesh->hedge);
    if (mesh->medge) {
      /* Hash_Delete already cleaned all implicit edge entities */
      List_Delete(mesh->medge);
    }
  } else {
    if (mesh->medge) {
      ne = mesh->ne;
      i = 0;
      while ((me = List_Next_Entry(mesh->medge,&i))) {
	ME_Destroy_For_MESH_Delete(me);
      }
      List_Delete(mesh->medge);
    }
  }
  if (mesh->mvertex) {
    nv = mesh->nv;
    i = 0;
    while ((mv = List_Next_Entry(mesh->mvertex,&i))) {
      MV_Destroy_For_MESH_Delete(mv);
    }
    List_Delete(mesh->mvertex);
  }

#ifdef MSTK_HAVE_MPI
  /*
    FOR NOW IT SEEMS THAT GHOST REGIONS ARE ALSO ENCOUNTERED WHEN
    GOING THROUGH THE REGULAR MESH AND SO WE DON'T NEED TO DELETE THEM
    SEPARATELY. NOT SURE THAT IT IS THE BEHAVIOR WE WANT */

  if (mesh->ghregion) {
    /*
      i = 0;
      while ((ghr = List_Next_Entry(mesh->ghregion,&i))) {
      MR_Destroy_For_MESH_Delete(ghr);
      }
    */
    List_Delete(mesh->ghregion);
  }
  if (mesh->ghface) {
    /*
      i = 0;
      while ((ghf = List_Next_Entry(mesh->ghface,&i))) {
      MF_Destroy_For_MESH_Delete(ghf);
      }
    */
    List_Delete(mesh->ghface);
  }
  if (mesh->ghedge) {
    /*
    i = 0;
    while ((ghe = List_Next_Entry(mesh->ghedge,&i))) {
      ME_Destroy_For_MESH_Delete(ghe);
    }
    */
    List_Delete(mesh->ghedge);
  }
  if (mesh->ghvertex) {
    /*
    i = 0;
    while ((ghv = List_Next_Entry(mesh->ghvertex,&i))) {
      MV_Destroy_For_MESH_Delete(ghv);
    }
    */
    List_Delete(mesh->ghvertex);
  }
#endif
  
  if (mesh->AttribList) {
    i = 0;
    while ((attrib = List_Next_Entry(mesh->AttribList,&i)))
      MAttrib_Destroy_For_MESH_Delete(attrib);
    List_Delete(mesh->AttribList);
  }

  if (mesh->MSetList) {
    i = 0;
    while ((mset = List_Next_Entry(mesh->MSetList,&i)))
      MSet_Delete(mset);
    List_Delete(mesh->MSetList);
  }

  free(mesh);
}

int MESH_SetRepType(Mesh_ptr mesh, RepType type) {

  if (mesh->reptype != type) {
    mesh->reptype = type;

    if (type >= R1 && type <= R4) 
      mesh->hedge = Hash_New(0, 1);
    else {
      if (mesh->hedge)
	Hash_Delete(mesh->hedge);
      mesh->hedge = (Hash_ptr) NULL;
    }
 
    if (type >= R1 && type <= R2) 
      mesh->hface = Hash_New(0, 1);
    else {
      if (mesh->hface)
	Hash_Delete(mesh->hface);
      mesh->hface = (Hash_ptr) NULL;
    }
  }

  return 1;
}

void MESH_SetGModel(Mesh_ptr mesh, GModel_ptr gm) {
  mesh->geom = gm;
}

GModel_ptr MESH_GModel(Mesh_ptr mesh) {
  return mesh->geom;
}

RepType MESH_RepType(Mesh_ptr mesh) {
  return mesh->reptype;
}

char *MESH_RepType_Str(Mesh_ptr mesh) {
  char *rstr = (char *) malloc(3*sizeof(char));

  strcpy(rstr,MESH_rtype_str[mesh->reptype]);
  return rstr;
}


int MESH_Num_Attribs(Mesh_ptr mesh) {
  if (mesh->AttribList)
    return List_Num_Entries(mesh->AttribList);
  else
    return 0;
}

MAttrib_ptr MESH_Attrib(Mesh_ptr mesh, int i) {
  if (mesh->AttribList)
    return List_Entry(mesh->AttribList,i);
  else
    return NULL;
}
  
MAttrib_ptr MESH_Next_Attrib(Mesh_ptr mesh, int *index) {
  if (mesh->AttribList)
    return List_Next_Entry(mesh->AttribList,index);
  else
    return NULL;
}

MAttrib_ptr MESH_AttribByName(Mesh_ptr mesh, const char *name) {
  if (mesh->AttribList) {
    int idx = 0;
    MAttrib_ptr attrib;
    char attname[256];
    
    while ((attrib = List_Next_Entry(mesh->AttribList,&idx))) {
      MAttrib_Get_Name(attrib,attname);
      if (strcmp(name,attname) == 0)
	return attrib;
    }
  }
  
  return NULL;
}

void MESH_Add_Attrib(Mesh_ptr mesh, MAttrib_ptr attrib) {
  if (!mesh->AttribList)
    mesh->AttribList = List_New(3);
  
  List_Add(mesh->AttribList,attrib);
}

void MESH_Rem_Attrib(Mesh_ptr mesh, MAttrib_ptr attrib) {
  int idx;
  MVertex_ptr mv;
  MEdge_ptr me;
  MFace_ptr mf;
  MRegion_ptr mr;

  if (mesh->AttribList) {
    List_Rem(mesh->AttribList,attrib);

    switch (MAttrib_Get_EntDim(attrib)) {
    case MVERTEX:
      idx = 0;
      while ((mv = MESH_Next_Vertex(mesh,&idx)))
	MEnt_Rem_AttVal(mv,attrib);
      break;

    case MEDGE:
      idx = 0;
      while ((me = MESH_Next_Edge(mesh,&idx)))
	MEnt_Rem_AttVal(me,attrib);
      break;

    case MFACE:
      idx = 0;
      while ((mf = MESH_Next_Face(mesh,&idx)))
	MEnt_Rem_AttVal(mf,attrib);
      break;

    case MREGION:
      idx = 0;
      while ((mr = MESH_Next_Region(mesh,&idx)))
	MEnt_Rem_AttVal(mr,attrib);
      break;

    case MALLTYPE:
      idx = 0;
      while ((mv = MESH_Next_Vertex(mesh,&idx)))
	MEnt_Rem_AttVal(mv,attrib);

      idx = 0;
      while ((me = MESH_Next_Edge(mesh,&idx)))
	MEnt_Rem_AttVal(me,attrib);
      
      idx = 0;
      while ((mf = MESH_Next_Face(mesh,&idx)))
	MEnt_Rem_AttVal(mf,attrib);
      
      idx = 0;
      while ((mr = MESH_Next_Region(mesh,&idx)))
	MEnt_Rem_AttVal(mr,attrib);
      break;

    default:
      break;
    }
  }
}


void MESH_Clear_Attrib(Mesh_ptr mesh, MAttrib_ptr attrib) {

  if (mesh->AttribList) {
    int idx = 0;
    MVertex_ptr mv;
    while ((mv = MESH_Next_Vertex(mesh,&idx)))
      MEnt_Clear_AttVal(mv,attrib);

    idx = 0;
    MEdge_ptr me;
    while ((me = MESH_Next_Edge(mesh,&idx)))
      MEnt_Clear_AttVal(me,attrib);

    idx = 0;
    MFace_ptr mf;
    while ((mf = MESH_Next_Face(mesh,&idx)))
      MEnt_Clear_AttVal(mf,attrib);

    idx = 0;
    MRegion_ptr mr;
    while ((mr = MESH_Next_Region(mesh,&idx)))
      MEnt_Clear_AttVal(mr,attrib);
  }
}


int MESH_Num_MSets(Mesh_ptr mesh) {
  if (mesh->MSetList)
    return List_Num_Entries(mesh->MSetList);
  else
    return 0;
}

MSet_ptr MESH_MSet(Mesh_ptr mesh, int i) {
  if (mesh->MSetList)
    return List_Entry(mesh->MSetList,i);
  else
    return NULL;
}
  
MSet_ptr MESH_Next_MSet(Mesh_ptr mesh, int *index) {
  if (mesh->MSetList)
    return List_Next_Entry(mesh->MSetList,index);
  else
    return NULL;
}

MSet_ptr MESH_MSetByName(Mesh_ptr mesh, const char *name) {
  if (mesh->MSetList) {
    int idx = 0;
    MSet_ptr mset;
    char attname[256];
    
    while ((mset = List_Next_Entry(mesh->MSetList,&idx))) {
      MSet_Name(mset,attname);
      if (strcmp(name,attname) == 0)
	return mset;
    }
  }
  
  return NULL;
}

void MESH_Add_MSet(Mesh_ptr mesh, MSet_ptr mset) {
  if (!mesh->MSetList)
    mesh->MSetList = List_New(3);
  
  List_Add(mesh->MSetList,mset);
}

void MESH_Rem_MSet(Mesh_ptr mesh, MSet_ptr mset) {
  if (mesh->MSetList)
    List_Rem(mesh->MSetList,mset);
}


void MESH_FillHash_Edges(Mesh_ptr mesh) {
  MRegion_ptr region;
  MEdge_ptr edge;
  List_ptr redges;
  int i, locks;

  locks = MESH_AutoLock(mesh);
  /* Force AutoLocking */
  MESH_Set_AutoLock(mesh, 1);
#ifdef DEBUG
  MSTK_Report("Mesh_FillHash_Edges","Inefficient to call routines like MESH_Num_Edges, MESH_Edge, MESH_Next_Edge with this representation\n",MSTK_WARN);
#endif
  i = 0;
  while ((region = MESH_Next_Region(mesh, &i))) {
    redges = MR_Edges(region);
    List_Delete(redges);
  }
  if (mesh->medge) free(mesh->medge);
  mesh->medge = Hash_Entries(mesh->hedge);
  mesh->ne = List_Num_Entries(mesh->medge);
  for (i = 0; i < mesh->ne; i++) {
    edge = MESH_Edge(mesh,i);
    ME_Set_ID(edge,i+1);
  }
  MESH_Set_AutoLock(mesh, locks);
}

void MESH_FillHash_Faces(Mesh_ptr mesh) {
  MRegion_ptr region;
  MFace_ptr face;
  List_ptr rfaces;
  int i, locks;

  locks = MESH_AutoLock(mesh);
  /* Force AutoLocking */
  MESH_Set_AutoLock(mesh, 1);
#ifdef DEBUG
  MSTK_Report("Mesh_FillHash_Edges","Inefficient to call routines like MESH_Num_Faces, MESH_Face, MESH_Next_Face with this representation\n",MSTK_WARN);
#endif
  i = 0;
  while ((region = MESH_Next_Region(mesh, &i))) {
    rfaces = MR_Faces(region);
    List_Delete(rfaces);
  }
  if (mesh->mface) free(mesh->mface);
  mesh->mface = Hash_Entries(mesh->hface);
  mesh->nf = List_Num_Entries(mesh->mface);
  for (i = 0; i < mesh->nf; i++) {
    face = MESH_Face(mesh,i);
    ME_Set_ID(face,i+1);
  }
  MESH_Set_AutoLock(mesh, locks);
}

int MESH_Num_Vertices(Mesh_ptr mesh) {
  return mesh->nv;
}

int MESH_Num_Edges(Mesh_ptr mesh) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R4)) {
    if ((mesh->ne == 0) || (mesh->ne != Hash_Num_Entries(mesh->hedge))) {
      MESH_FillHash_Edges(mesh);
    }
  }
  return mesh->ne;
}

int MESH_Num_Faces(Mesh_ptr mesh) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R2)) {
    if (mesh->nr == 0)   /* 2D mesh */
      return mesh->nf;   
    else if ((mesh->nf == 0) || (mesh->nf != Hash_Num_Entries(mesh->hface))) {
      MESH_FillHash_Faces(mesh);
    }
  }
  return mesh->nf;
}

int MESH_Num_Regions(Mesh_ptr mesh) {
  return mesh->nr;
}

/*
int        MESH_Num_Elements(Mesh_ptr mesh) {
  return mesh->nel;
}
*/

MVertex_ptr MESH_Vertex(Mesh_ptr mesh, int i) {
  if (i >= mesh->nv) {
#ifdef DEBUG
    MSTK_Report("MESH_Vertex","Non-existent vertex requested\n",MSTK_ERROR);
#endif
    return (MVertex_ptr) NULL;
  }
  else
    return (MVertex_ptr) List_Entry(mesh->mvertex, i);
}

MEdge_ptr MESH_Edge(Mesh_ptr mesh, int i) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R4)) {
    if ((mesh->ne == 0) || (mesh->ne != Hash_Num_Entries(mesh->hedge))) {
      MESH_FillHash_Edges(mesh);
    }
  }
  if (i >= mesh->ne) {
#ifdef DEBUG
    MSTK_Report("MESH_Edge","Non-existent edge requested\n",MSTK_ERROR);
#endif
    return (MEdge_ptr) NULL;
  }
  else 
    return (MEdge_ptr) List_Entry(mesh->medge, i);
}

MFace_ptr MESH_Face(Mesh_ptr mesh, int i) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R2)) {
    if (mesh->nr) {
      if ((mesh->nf == 0) || (mesh->nf != Hash_Num_Entries(mesh->hface))) {
	MESH_FillHash_Faces(mesh);
      }
    }
  }
  if (i >= mesh->nf) {
#ifdef DEBUG
    MSTK_Report("MESH_Face","Non-existent face requested\n",MSTK_ERROR);
#endif
    return (MFace_ptr) NULL;
  }
  else
    return (MFace_ptr) List_Entry(mesh->mface, i);
}

MRegion_ptr MESH_Region(Mesh_ptr mesh, int i) {

  if (i >= mesh->nr) {
#ifdef DEBUG
    MSTK_Report("MESH_Region","Non-existent region requested\n",MSTK_ERROR);
#endif
    return (MRegion_ptr) NULL;
  }
  else
    return (MRegion_ptr) List_Entry(mesh->mregion, i);
}

MVertex_ptr MESH_Next_Vertex(Mesh_ptr mesh, int *index) {
  if (mesh->mvertex)
    return (MVertex_ptr) List_Next_Entry(mesh->mvertex, index);
  else
    return NULL;
}

MEdge_ptr MESH_Next_Edge(Mesh_ptr mesh, int *index) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R4)) {
    if ((mesh->ne == 0) || (mesh->ne != Hash_Num_Entries(mesh->hedge))) {
      MESH_FillHash_Edges(mesh);
    }
  }
  if (mesh->medge)
    return (MEdge_ptr) List_Next_Entry(mesh->medge, index);
  else
    return NULL;
}

MFace_ptr MESH_Next_Face(Mesh_ptr mesh, int *index) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R2)) {
    if (mesh->nr) {
      if ((mesh->nf == 0) || (mesh->nf != Hash_Num_Entries(mesh->hface))) {
	MESH_FillHash_Faces(mesh);
      }
    }
  }
  if (mesh->mface)
    return (MFace_ptr) List_Next_Entry(mesh->mface, index);
  else
    return NULL;
}

MRegion_ptr MESH_Next_Region(Mesh_ptr mesh, int *index) {
  if (mesh->mregion)
    return (MRegion_ptr) List_Next_Entry(mesh->mregion, index);
  else
    return NULL;
}

  /* The current implementation of following four functions,
     MESH_*FromID, relies on the fact the the mesh entities are stored
     in linear arrays. So go to the (ID-1)'th element in the list. If
     it is a static mesh, this will be the element with the right
     ID. If we get NULL (the ID is higher than the list size due to
     deletions) or the wrong ID, then brute force search the list */


MVertex_ptr MESH_VertexFromID(Mesh_ptr mesh, int id) {
  int istart, j;
  MVertex_ptr mv;
  int list_size = List_Size_Raw(mesh->mvertex); 

  if (id < 1 || id > list_size || id > mesh-> max_vid)
    return NULL;

  /* Typically this will work for static meshes */

  istart = id-1;
  mv = (MVertex_ptr) List_Entry_Raw(mesh->mvertex,istart);
  if (mv && MV_ID(mv) == id)
    return mv;

  /* If entities have been deleted and new ones introduced, we have to
     search. Its a linear search but somewhat smart because it assumes
     that the entity with the desired ID will be close to where we
     expected to find it in a static list */

  for (j = istart; j >= 0; j--) {
    mv = (MVertex_ptr) List_Entry_Raw(mesh->mvertex,j);
    if (mv && MV_ID(mv) == id)
      return mv;
  }

  for (j = istart; j < list_size; j++) {
    mv = (MVertex_ptr) List_Entry_Raw(mesh->mvertex,j);
    if (mv && MV_ID(mv) == id)
      return mv;
  }

  return NULL;
}

MEdge_ptr MESH_EdgeFromID(Mesh_ptr mesh, int id) {
  int istart, j;
  MEdge_ptr me;
  int list_size = List_Size_Raw(mesh->medge); 

  if (id < 1 || id > list_size || id > mesh->max_eid)
    return NULL;

  /* Typically this will work for static meshes */

  istart = id-1;
  me = (MEdge_ptr) List_Entry_Raw(mesh->medge,istart);
  if (me && ME_ID(me) == id)
    return me;
    
  /* If entities have been deleted and new ones introduced, we have to
     search. Its a linear search but somewhat smart because it assumes
     that the entity with the desired ID will be close to where we
     expected to find it in a static list */

  for (j = istart; j >= 0; j--) {
    me = (MEdge_ptr) List_Entry_Raw(mesh->medge,j);
    if (me && ME_ID(me) == id)
      return me;
  }

  for (j = istart; j < list_size; j++) {
    me = (MEdge_ptr) List_Entry_Raw(mesh->medge,j);
    if (me && ME_ID(me) == id)
      return me;
  }

  return NULL;
}

MFace_ptr MESH_FaceFromID(Mesh_ptr mesh, int id) {
  int istart, j;
  MFace_ptr mf;
  int list_size = List_Size_Raw(mesh->mface); 

  if (id < 1 || id > list_size || id > mesh->max_fid)
    return NULL;

  /* Typically this will work for static meshes */

  istart = id-1;
  mf = (MFace_ptr) List_Entry_Raw(mesh->mface,istart);
  if (mf && MF_ID(mf) == id)
    return mf;
    
  /* If entities have been deleted and new ones introduced, we have to
     search. Its a linear search but somewhat smart because it assumes
     that the entity with the desired ID will be close to where we
     expected to find it in a static list */

  for (j = istart; j >= 0; j--) {
    mf = (MFace_ptr) List_Entry_Raw(mesh->mface,j);
    if (mf && MF_ID(mf) == id)
      return mf;
  }

  for (j = istart; j < list_size; j++) {
    mf = (MFace_ptr) List_Entry_Raw(mesh->mface,j);
    if (mf && MF_ID(mf) == id)
      return mf;
  }

  return NULL;
}

MRegion_ptr MESH_RegionFromID(Mesh_ptr mesh, int id) {
  int istart, j;
  MRegion_ptr mr;
  int list_size = List_Size_Raw(mesh->mregion); 

  if (id < 1 || id > list_size || id > mesh->max_rid)
    return NULL;

  /* Typically this will work for static meshes */

  istart = id-1;
  mr = (MRegion_ptr) List_Entry_Raw(mesh->mregion,istart);
  if (mr && MR_ID(mr) == id)
    return mr;
    
  /* If entities have been deleted and new ones introduced, we have to
     search. Its a linear search but somewhat smart because it assumes
     that the entity with the desired ID will be close to where we
     expected to find it in a static list */

  for (j = istart; j >= 0; j--) {
    mr = (MRegion_ptr) List_Entry_Raw(mesh->mregion,j);
    if (MR_ID(mr) == id)
      return mr;
  }

  for (j = istart; j < list_size; j++) {
    mr = (MRegion_ptr) List_Entry_Raw(mesh->mregion,j);
    if (MR_ID(mr) == id)
      return mr;
  }

  return NULL;
}

 
MEntity_ptr MESH_EntityFromID(Mesh_ptr mesh, int mtype, int id) {

  switch (mtype) {
  case 0:
    return MESH_VertexFromID(mesh,id);
  case 1:
    return MESH_EdgeFromID(mesh,id);
  case 2:
    return MESH_FaceFromID(mesh,id);
  case 3:
    return MESH_RegionFromID(mesh,id);
  default:
    MSTK_Report("MESH_EntityFromID","Unrecognized entity type",MSTK_ERROR);
    return NULL;
  }

}

void MESH_Add_Vertex(Mesh_ptr mesh, MVertex_ptr v) {
  if (mesh->mvertex == (List_ptr) NULL)
    mesh->mvertex = List_New(10);

  mesh->mvertex = List_Add(mesh->mvertex, (void *) v);
  mesh->nv = List_Num_Entries(mesh->mvertex);

  if (MV_ID(v) == 0) { /* New Vertex */
    (mesh->max_vid)++;
    MV_Set_ID(v,mesh->max_vid);
  }
}

void MESH_Add_Edge(Mesh_ptr mesh, MEdge_ptr e){
  /* Have to check if edges exist in this type of representation */
  if (mesh->reptype >= R1 && mesh->reptype <= MSTK_MAXREP)
    return;

  if (mesh->medge == (List_ptr) NULL)
    mesh->medge = List_New(10);

  mesh->medge = List_Add(mesh->medge, (void *) e);
  mesh->ne = List_Num_Entries(mesh->medge);
  if (ME_ID(e) == 0) { /* New edge */
    (mesh->max_eid)++;
    ME_Set_ID(e,mesh->max_eid);
  }
}    
     
void MESH_Add_Face(Mesh_ptr mesh, MFace_ptr f){
  /* Have to check if faces exist in this type of representation */
  if (mesh->nr && (mesh->reptype == R1 || mesh->reptype == R2))
    return;

  if ((mesh->reptype == R4) && (MF_Region(f,0) || MF_Region(f,1))) {
#ifdef DEBUG
    MSTK_Report("MESH_Add_Face","Can add disconnected faces only",MSTK_ERROR);
#endif
    return;
  }

  if (mesh->mface == (List_ptr) NULL)
    mesh->mface = List_New(10);
  
  mesh->mface = List_Add(mesh->mface, (void *) f);
  mesh->nf = List_Num_Entries(mesh->mface);

  if (MF_ID(f) == 0) { /* New face */
    (mesh->max_fid)++;
    MF_Set_ID(f,mesh->max_fid);

  }
}    
     
void MESH_Add_Region(Mesh_ptr mesh, MRegion_ptr r){
  if (mesh->mregion == (List_ptr) NULL)
    mesh->mregion = List_New(10);

  mesh->mregion = List_Add(mesh->mregion, (void *) r);
  mesh->nr = List_Num_Entries(mesh->mregion);

  if (MR_ID(r) == 0) { /* New region */
    (mesh->max_rid)++;
    MR_Set_ID(r,mesh->max_rid);
  }
}    
     
void MESH_Rem_Vertex(Mesh_ptr mesh, MVertex_ptr v) {
  int fnd=0, i, id;

  if (mesh->mvertex == (List_ptr) NULL) {
#ifdef DEBUG
    MSTK_Report("Mesh_Rem_Vertex","No vertices in mesh to remove", MSTK_ERROR);
#endif
    return;
  }

  /* If the list of vertices has not been compressed or the mesh has not
     been renumbered, the real position of the entry will be ID-1. If it
     is then delete it directly */

  id = MV_ID(v);

  i = id-1;
  if (List_Entry_Raw(mesh->mvertex,i) == v) {
    List_Remi_Raw(mesh->mvertex,i);
    fnd = 1;
  }

  /* If not, look for the entry assuming the list of vertices is
   * sorted according to the passed function (in this case
   * MV_ID). This is true because we never put new entities in holes
   * in the entity list created by deleted entities - they always go
   * at the end. */
  
  if (!fnd)
    fnd = List_RemSorted(mesh->mvertex,v,&(MV_ID));

  if (!fnd)
    MSTK_Report("MESH_Rem_Vertex","Vertex not found in list",MSTK_FATAL);

  mesh->nv = List_Num_Entries(mesh->mvertex);

#ifdef MSTK_HAVE_MPI
  if (MV_PType(v) == POVERLAP) {
    fnd = List_RemSorted(mesh->ovvertex,v,&(MV_GlobalID));
    if (!fnd)
      MSTK_Report("MESH_Rem_Vertex","Vertex not found in overlap list",MSTK_FATAL);
    MESH_Mark_ParallelAdj_Stale(mesh);
  }
#endif


  return;
}    
     
void MESH_Rem_Edge(Mesh_ptr mesh, MEdge_ptr e) {
  int fnd=0, i, id;

  if (mesh->medge == (List_ptr) NULL) {
#ifdef DEBUG
    MSTK_Report("Mesh_Rem_Edge","No Edges in mesh to remove",MSTK_ERROR);
#endif
    return;
  }

  /* If the list of edges has not been compressed or the mesh has not
     been renumbered, the real position of the entry will be ID-1. If it
     is then delete it directly */

  id = ME_ID(e);

  i = id-1;
  if (List_Entry_Raw(mesh->medge,i) == e) {
    List_Remi_Raw(mesh->medge,i);
    fnd = 1;
  }
  
  /* If not, look for the entry assuming the list of vertices is
   * sorted according to the passed function (in this case
   * MV_ID). This is true because we never put new entities in holes
   * in the entity list created by deleted entities - they always go
   * at the end. */
  
  if (!fnd)
    fnd = List_RemSorted(mesh->medge,e,&(ME_ID));

  if (!fnd)
    MSTK_Report("MESH_Rem_Edge","Edge not found in list",MSTK_FATAL);

  mesh->ne = List_Num_Entries(mesh->medge);

#ifdef MSTK_HAVE_MPI
  if (ME_PType(e) == POVERLAP) {
    fnd = List_RemSorted(mesh->ovedge,e,&(ME_GlobalID));
    if (!fnd)
      MSTK_Report("MESH_Rem_Edge","Edge not found in overlap list",MSTK_FATAL);
    MESH_Mark_ParallelAdj_Stale(mesh);
  }
#endif

  return;
}    
     
void MESH_Rem_Face(Mesh_ptr mesh, MFace_ptr f){
  int fnd=0, i, id;

  if (mesh->mface == (List_ptr) NULL) {
#ifdef DEBUG
    MSTK_Report("Mesh_Rem_Face","No Faces in mesh to remove",MSTK_ERROR);
#endif
    return;
  }

  /* If the list of faces has not been compressed or the mesh has not
     been renumbered, the real position of the entry will be ID-1. If it
     is then delete it directly */
  
  id = MF_ID(f);

  i = id-1;
  if (List_Entry_Raw(mesh->mface,i) == f) {
    List_Remi_Raw(mesh->mface,i);
    fnd = 1;
  }

  /* If not, look for the entry assuming the list of vertices is
   * sorted according to the passed function (in this case
   * MV_ID). This is true because we never put new entities in holes
   * in the entity list created by deleted entities - they always go
   * at the end. */
  
  if (!fnd)
    fnd = List_RemSorted(mesh->mface,f,&(MF_ID));

  if (!fnd)
    MSTK_Report("MESH_Rem_Face","Face not found in list",MSTK_FATAL);

  mesh->nf = List_Num_Entries(mesh->mface);

#ifdef MSTK_HAVE_MPI
  if (MF_PType(f) == POVERLAP) {
    fnd = List_RemSorted(mesh->ovface,f,&(MF_GlobalID));
    if (!fnd)
      MSTK_Report("MESH_Rem_Face","Face not found in overlap list",MSTK_FATAL);
    MESH_Mark_ParallelAdj_Stale(mesh);
  }
#endif

  return;
}    
     
void MESH_Rem_Region(Mesh_ptr mesh, MRegion_ptr r){
  int fnd=0, i, id;

  if (mesh->mregion == (List_ptr) NULL) {
    MSTK_Report("Mesh_Rem_Region","No regions in mesh to remove",MSTK_ERROR);
    return;
  }

  /* If the list of regions has not been compressed or the mesh has not
     been renumbered, the real position of the entry will be ID-1. If it
     is then delete it directly */

  
  id = MR_ID(r);

  i = id-1;
  if (List_Entry_Raw(mesh->mregion,i) == r) {
    List_Remi_Raw(mesh->mregion,i);
    fnd = 1;
  }

  /* If not, look for the entry assuming the list of vertices is
   * sorted according to the passed function (in this case
   * MV_ID). This is true because we never put new entities in holes
   * in the entity list created by deleted entities - they always go
   * at the end. */
  
  if (!fnd)
    fnd = List_RemSorted(mesh->mregion,r,&(MR_ID));

  if (!fnd)
    MSTK_Report("MESH_Rem_Region","Region not found in list",MSTK_FATAL);

  mesh->nr = List_Num_Entries(mesh->mregion);

#ifdef MSTK_HAVE_MPI
  if (MR_PType(r) == POVERLAP) {
    fnd = List_RemSorted(mesh->ovregion,r,&(MR_GlobalID));
    if (!fnd)
      MSTK_Report("MESH_Rem_Region","Region not found in overlap list",MSTK_FATAL);
    MESH_Mark_ParallelAdj_Stale(mesh);
  }
#endif
}    

List_ptr   MESH_Vertex_List(Mesh_ptr mesh) {
  return mesh->mvertex;
}
List_ptr   MESH_Edge_List(Mesh_ptr mesh) {
  return mesh->medge;
}
List_ptr   MESH_Face_List(Mesh_ptr mesh) {
  return mesh->mface;
}
List_ptr   MESH_Region_List(Mesh_ptr mesh) {
  return mesh->mregion;
}



#ifdef MSTK_HAVE_MPI

  /* Even though some of these routine names are uncharacteristically
     long, they describe exactly what they do, so use them */

  /*
  void MESH_Set_Num_Partitions(Mesh_ptr mesh, unsigned int numpartitions) {
    mesh->numpartns = numpartitions;
    if (!mesh->par_adj_flags)
      mesh->par_adj_flags = (unsigned int *) calloc(numpartitions,sizeof(unsigned int));
  }
  */

  /* Also, for an explanation of the fields of par_adj_flags and
     par_adj_info see the note at the in Mesh.h file */

  void MESH_Set_Prtn(Mesh_ptr mesh, unsigned int partition, unsigned int numpartitions) {
    mesh->mypartn = partition;
    mesh->numpartns = numpartitions;

    if (!mesh->par_adj_flags)
      mesh->par_adj_flags = (unsigned int *)
        calloc(numpartitions,sizeof(unsigned int));
  }

  unsigned int MESH_Prtn(Mesh_ptr mesh) {
    return mesh->mypartn;
  }


  int MESH_ParallelAdj_Current(Mesh_ptr mesh) {
    return mesh->par_adj_current;
  }

  void MESH_Mark_ParallelAdj_Current(Mesh_ptr mesh) {
    mesh->par_adj_current = 1;
  }

  void MESH_Mark_ParallelAdj_Stale(Mesh_ptr mesh) {
    mesh->par_adj_current = 0;
  }


  void MESH_Flag_Has_Ghosts_From_Prtn(Mesh_ptr mesh, unsigned int prtn, MType mtype) {
#ifdef DEBUG
    if (!mesh->par_adj_flags)
      MSTK_Report("MESH_Flag_Has_Overlaps_On_Prtn",
                  "Call MESH_Set_Prtn before calling this routine for the first time",
                  MSTK_FATAL);
#endif
#ifdef MSTK_HAVE_MPI

    if (prtn == mesh->mypartn) return;

    if (mtype == MUNKNOWNTYPE || mtype == MANYTYPE) 
      return;
    else if (mtype == MALLTYPE) {
      mesh->par_adj_flags[prtn] |= 1;
      mesh->par_adj_flags[prtn] |= 1<<2;
      mesh->par_adj_flags[prtn] |= 1<<4;
      mesh->par_adj_flags[prtn] |= 1<<6;
    }
    else
      mesh->par_adj_flags[prtn] |= 1<<(2*mtype);
#endif
  }
  
  unsigned int MESH_Has_Ghosts_From_Prtn(Mesh_ptr mesh, unsigned int prtn, MType mtype) {
#ifdef DEBUG
    if (!mesh->par_adj_flags)
      MSTK_Report("MESH_Flag_Has_Overlaps_On_Prtn",
                  "Call MESH_Set_Prtn before calling this routine for the first time",
                  MSTK_FATAL);
#endif
#ifdef MSTK_HAVE_MPI
    if (mtype == MUNKNOWNTYPE || mesh->par_adj_flags == NULL)
      return 0;
    else if (mtype == MANYTYPE)
      return ((mesh->par_adj_flags[prtn] & 1) |
              (mesh->par_adj_flags[prtn]>>2 & 1) |
              (mesh->par_adj_flags[prtn]>>4 & 1) |
              (mesh->par_adj_flags[prtn]>>6 & 1));

    else if (mtype == MALLTYPE)
      return ((mesh->par_adj_flags[prtn] & 1) &
              (mesh->par_adj_flags[prtn]>>2 & 1) &
              (mesh->par_adj_flags[prtn]>>4 & 1) &
              (mesh->par_adj_flags[prtn]>>6 & 1));
    else      
      return ((mesh->par_adj_flags[prtn])>>(2*mtype) & 1);
#else
   return 0;
#endif
  }
  
  void MESH_Flag_Has_Overlaps_On_Prtn(Mesh_ptr mesh, unsigned int prtn, MType mtype) {
#ifdef DEBUG
    if (!mesh->par_adj_flags)
      MSTK_Report("MESH_Flag_Has_Overlaps_On_Prtn",
                  "Call MESH_Set_Prtn before calling this routine for the first time",
                  MSTK_FATAL);
#endif

#ifdef MSTK_HAVE_MPI
    if (mtype == MUNKNOWNTYPE || mtype == MANYTYPE) 
      return;
    else if (mtype == MALLTYPE) {
      mesh->par_adj_flags[prtn] |= 1<<1;
      mesh->par_adj_flags[prtn] |= 1<<3;
      mesh->par_adj_flags[prtn] |= 1<<5;
      mesh->par_adj_flags[prtn] |= 1<<7;
    }
    else
      mesh->par_adj_flags[prtn] |= 1<<(2*mtype+1);
#endif
  }
  
  unsigned int MESH_Has_Overlaps_On_Prtn(Mesh_ptr mesh, unsigned int prtn, MType mtype) {
#ifdef DEBUG
    if (!mesh->par_adj_flags)
      MSTK_Report("MESH_Flag_Has_Overlaps_On_Prtn",
                  "Call MESH_Set_Prtn before calling this routine for the first time",
                  MSTK_FATAL);
#endif
#ifdef MSTK_HAVE_MPI
    if (mtype == MUNKNOWNTYPE || mesh->par_adj_flags == NULL) 
      return 0;
    else if (mtype == MANYTYPE)
      return ((mesh->par_adj_flags[prtn]>>1 & 1) |
              (mesh->par_adj_flags[prtn]>>3 & 1) |
              (mesh->par_adj_flags[prtn]>>5 & 1) |
              (mesh->par_adj_flags[prtn]>>7 & 1));

    else if (mtype == MALLTYPE)
      return ((mesh->par_adj_flags[prtn]>>1 & 1) &
              (mesh->par_adj_flags[prtn]>>3 & 1) &
              (mesh->par_adj_flags[prtn]>>5 & 1) &
              (mesh->par_adj_flags[prtn]>>7 & 1));
    else      
      return ((mesh->par_adj_flags[prtn])>>(2*mtype+1) & 1);
#else
    return 0;
#endif
  }


  unsigned int MESH_Num_GhostPrtns(Mesh_ptr mesh) {
    if (!mesh->par_recv_info) MESH_Init_Par_Recv_Info(mesh);
    return (unsigned int) mesh->par_recv_info[0];
  }


  void MESH_GhostPrtns(Mesh_ptr mesh, unsigned int *pnums) {
    unsigned int ghnum = mesh->par_recv_info[0];
    memcpy(pnums,mesh->par_recv_info+1,ghnum*sizeof(unsigned int));
  }


  void MESH_Init_Par_Recv_Info(Mesh_ptr mesh) {
    int i;

    if (mesh->par_recv_info) { /* delete old information */
      free(mesh->par_recv_info);
    }
    int ngp = 0;
      
    /* count the number of ghost processor on this processor */
    
    for (i = 0; i < mesh->numpartns; i++)
      if (MESH_Has_Ghosts_From_Prtn(mesh, i, MANYTYPE))
	ngp++;
    
    mesh->par_recv_info = (unsigned int *) calloc((1+5*ngp),sizeof(unsigned int));
    
    mesh->par_recv_info[0] = ngp;
    
    ngp = 0;
    for (i = 0; i < mesh->numpartns; i++)
      if (MESH_Has_Ghosts_From_Prtn(mesh, i, MANYTYPE)) {
	mesh->par_recv_info[1+ngp] = i;
	ngp++;
      }
    
  }
  
  void MESH_Set_Num_Recv_From_Prtn(Mesh_ptr mesh, unsigned int prtn, MType mtype, unsigned int numrecv) {
    int found, i, ghnum;

    if (!mesh->par_recv_info) 
      MESH_Init_Par_Recv_Info(mesh);
 
    if (mtype < MVERTEX || mtype > MREGION) {
      MSTK_Report("MESH_Set_Num_RecvEnts_On_Prtn","Invalid entity type",MSTK_ERROR);
      return;
    }

    found = 0;
    ghnum = mesh->par_recv_info[0]; /* Number of ghost procs */
    /*    printf("set rank %d has %d processors to receive type %d\n",prtn,ghnum,mtype); */
    for (i = 0; i < ghnum; i++) {
      /* printf("   par_recv_info[%d]=%d\n",i,mesh->par_recv_info[i]); */
      if (mesh->par_recv_info[1+i] == prtn) {
        found = 1;
        break;
      }
    }
    if (!found) {
      char mesg[256];
      sprintf(mesg,"This partition (%-d) does not have ghost entities %d from partition %-d",mesh->mypartn,mtype,prtn);
      MSTK_Report("MESH_Set_Num_Recv_From_Prtn",mesg,MSTK_ERROR);
      return;
    }

    mesh->par_recv_info[1+ghnum+4*i+mtype] = numrecv;
  }

  unsigned int MESH_Num_Recv_From_Prtn(Mesh_ptr mesh, unsigned int prtn, MType mtype) {
    int found, i, ghnum;

    if (mtype < MVERTEX || mtype > MREGION) {
      MSTK_Report("MESH_Num_Recv_From_Prtn","Invalid entity type",MSTK_ERROR);
      return 0;
    }

    if (!mesh->par_recv_info) {
      MSTK_Report("MESH_Num_Recv_From_Prtn","Coding error. Data not initialized",MSTK_FATAL);
    }

    if (prtn == mesh->mypartn) return 0;
    

    found = 0;
    ghnum = mesh->par_recv_info[0]; /* Number of ghost procs */
    for (i = 0; i < ghnum; i++) 
      if (mesh->par_recv_info[1+i] == prtn) {
        found = 1;
        break;
      }

    if (!found) {
      char mesg[256];
      sprintf(mesg,"This partition (%-d) does not have ghost entities from partition %-d",mesh->mypartn,prtn);
      MSTK_Report("MESH_Num_Recv_From_Prtn",mesg,MSTK_ERROR);
      return 0;
    }

    return mesh->par_recv_info[1+ghnum+4*i+mtype];
  }



void MESH_Enable_GlobalIDSearch(Mesh_ptr mesh) {

  if (mesh->nv && !mesh->gid_sorted_mvlist) {
    mesh->gid_sorted_mvlist = List_Copy(mesh->mvertex);
    List_Compress(mesh->gid_sorted_mvlist);  // can't sort list with gaps
    List_Sort(mesh->gid_sorted_mvlist,mesh->nv,sizeof(MVertex_ptr),compareGlobalID);
  }

  if (mesh->ne && !mesh->gid_sorted_melist) {
    mesh->gid_sorted_melist = List_Copy(mesh->medge);
    List_Compress(mesh->gid_sorted_melist);
    List_Sort(mesh->gid_sorted_melist,mesh->ne,sizeof(MEdge_ptr),compareGlobalID);
  }

  if (mesh->nf && !mesh->gid_sorted_mflist) {
    mesh->gid_sorted_mflist = List_Copy(mesh->mface);
    List_Compress(mesh->gid_sorted_mflist);
    List_Sort(mesh->gid_sorted_mflist,mesh->nf,sizeof(MFace_ptr),compareGlobalID);
  }

  if (mesh->nr && !mesh->gid_sorted_mrlist) {
    mesh->gid_sorted_mrlist = List_Copy(mesh->mregion);
    List_Compress(mesh->gid_sorted_mrlist);
    List_Sort(mesh->gid_sorted_mrlist,mesh->nr,sizeof(MRegion_ptr),compareGlobalID);
  }

}

void MESH_Sort_GlobalIDSearch_Lists(Mesh_ptr mesh) {

  if (mesh->nv && mesh->gid_sorted_mvlist) {
    List_Compress(mesh->gid_sorted_mvlist);  // can't sort list with gaps
    List_Sort(mesh->gid_sorted_mvlist,mesh->nv,sizeof(MVertex_ptr),compareGlobalID);
  }

  if (mesh->ne && mesh->gid_sorted_melist) {
    List_Compress(mesh->gid_sorted_melist);
    List_Sort(mesh->gid_sorted_melist,mesh->ne,sizeof(MEdge_ptr),compareGlobalID);
  }

  if (mesh->nf && mesh->gid_sorted_mflist) {
    List_Compress(mesh->gid_sorted_mflist);
    List_Sort(mesh->gid_sorted_mflist,mesh->nf,sizeof(MFace_ptr),compareGlobalID);
  }

  if (mesh->nr && mesh->gid_sorted_mrlist) {
    List_Compress(mesh->gid_sorted_mrlist);
    List_Sort(mesh->gid_sorted_mrlist,mesh->nr,sizeof(MRegion_ptr),compareGlobalID);
  }

}


void MESH_Disable_GlobalIDSearch(Mesh_ptr mesh) {

  if (mesh->gid_sorted_mvlist) {
    List_Delete(mesh->gid_sorted_mvlist);
    mesh->gid_sorted_mvlist = NULL;
  }

  if (mesh->gid_sorted_melist) {
    List_Delete(mesh->gid_sorted_melist);
    mesh->gid_sorted_melist = NULL;
  }
  if (mesh->gid_sorted_mflist) {
    List_Delete(mesh->gid_sorted_mflist);
    mesh->gid_sorted_mflist = NULL;
  }
  if (mesh->gid_sorted_mrlist) {
    List_Delete(mesh->gid_sorted_mrlist);
    mesh->gid_sorted_mrlist = NULL;
  }

}
                  
  /* right now using linear search, a hash table should be more efficient */
MVertex_ptr MESH_VertexFromGlobalID(Mesh_ptr mesh, int global_id) {
  int idx;
  MVertex_ptr mv;

  if (!mesh->gid_sorted_mvlist)
    MESH_Enable_GlobalIDSearch(mesh);

  MVertex_ptr mv_tmp = MV_New(NULL);
  MEnt_Set_GlobalID(mv_tmp,global_id);

  mv = (MVertex_ptr) List_Search(mesh->gid_sorted_mvlist,(void *)&mv_tmp,
                                 mesh->nv,sizeof(MVertex_ptr),compareGlobalID);

  MV_Delete(mv_tmp,0);

  return mv;
}

MEdge_ptr MESH_EdgeFromGlobalID(Mesh_ptr mesh, int global_id) {
  int idx;
  MEdge_ptr me;

  if (!mesh->gid_sorted_melist)
    MESH_Enable_GlobalIDSearch(mesh);

  MEdge_ptr me_tmp = ME_New(NULL);
  MEnt_Set_GlobalID(me_tmp,global_id);

  me = (MEdge_ptr) List_Search(mesh->gid_sorted_melist,(void *)&me_tmp,
                               mesh->ne,sizeof(MEdge_ptr),compareGlobalID);

  ME_Delete(me_tmp,0);
  
  return me;
}

MFace_ptr MESH_FaceFromGlobalID(Mesh_ptr mesh, int global_id) {
  int idx;
  MFace_ptr mf;

  if (!mesh->gid_sorted_mflist)
    MESH_Enable_GlobalIDSearch(mesh);

  MFace_ptr mf_tmp = MF_New(NULL);
  MEnt_Set_GlobalID(mf_tmp,global_id);

  mf = (MFace_ptr) List_Search(mesh->gid_sorted_mflist,(void *)&mf_tmp,
                               mesh->nf,sizeof(MFace_ptr),compareGlobalID);

  MF_Delete(mf_tmp,0);
  
  return mf;
}

MRegion_ptr MESH_RegionFromGlobalID(Mesh_ptr mesh, int global_id) {
  int idx;
  MRegion_ptr mr;

  if (!mesh->gid_sorted_mrlist)
    MESH_Enable_GlobalIDSearch(mesh);

  MRegion_ptr mr_tmp = MR_New(NULL);
  MEnt_Set_GlobalID(mr_tmp,global_id);

  mr = (MRegion_ptr) List_Search(mesh->gid_sorted_mrlist,(void *)&mr_tmp,
                                 mesh->nr,sizeof(MRegion_ptr),compareGlobalID);

  MR_Delete(mr_tmp,0);

  return mr;
}

MEntity_ptr MESH_EntityFromGlobalID(Mesh_ptr mesh, int mtype, int id) {

  switch (mtype) {
  case 0:
    return MESH_VertexFromGlobalID(mesh,id);
  case 1:
    return MESH_EdgeFromGlobalID(mesh,id);
  case 2:
    return MESH_FaceFromGlobalID(mesh,id);
  case 3:
    return MESH_RegionFromGlobalID(mesh,id);
  default:
    MSTK_Report("MESH_EntityFromGlobalID","Unrecognized entity type",MSTK_ERROR);
    return NULL;
  }

}


int MESH_Sort_GhostLists(Mesh_ptr mesh, 
                         int (*compfunc)(const void*, const void*)) {

  if (mesh->ghvertex) {
    List_Compress(mesh->ghvertex);  /* cannot sort list with gaps - no cost if no gaps */
    List_Sort(mesh->ghvertex,List_Num_Entries(mesh->ghvertex),sizeof(MVertex_ptr),compfunc);
  }
  if (mesh->ghedge) {
    List_Compress(mesh->ghedge);
    List_Sort(mesh->ghedge,List_Num_Entries(mesh->ghedge),sizeof(MEdge_ptr),compfunc);
  }
  if (mesh->ghface) {
    List_Compress(mesh->ghface);
    List_Sort(mesh->ghface,List_Num_Entries(mesh->ghface),sizeof(MFace_ptr),compfunc);
  }
  if (mesh->ghregion) {
    List_Compress(mesh->ghregion);
    List_Sort(mesh->ghregion,List_Num_Entries(mesh->ghregion),sizeof(MRegion_ptr),compfunc);
  }

  if (mesh->ovvertex) {
    List_Compress(mesh->ovvertex);
    List_Sort(mesh->ovvertex,List_Num_Entries(mesh->ovvertex),sizeof(MVertex_ptr),compfunc);
  }
  if (mesh->ovedge) {
    List_Compress(mesh->ovedge);
    List_Sort(mesh->ovedge,List_Num_Entries(mesh->ovedge),sizeof(MEdge_ptr),compfunc);
  }
  if (mesh->ovface) {
    List_Compress(mesh->ovface);
    List_Sort(mesh->ovface,List_Num_Entries(mesh->ovface),sizeof(MFace_ptr),compfunc);
  }
  if (mesh->ovregion) {
    List_Compress(mesh->ovregion);
    List_Sort(mesh->ovregion,List_Num_Entries(mesh->ovregion),sizeof(MRegion_ptr),compfunc);
  }

  return 1;
}

int MESH_Num_GhostVertices(Mesh_ptr mesh) {
  if(mesh->ghvertex)
    return List_Num_Entries(mesh->ghvertex);
  else
    return 0;
}
int MESH_Num_OverlapVertices(Mesh_ptr mesh) {
  if(mesh->ovvertex)
    return List_Num_Entries(mesh->ovvertex);
  else
    return 0;
}

int MESH_Num_InteriorVertices(Mesh_ptr mesh) {
  if(mesh->mvertex)
    return List_Num_Entries(mesh->mvertex) -\
      List_Num_Entries(mesh->ovvertex)- \
      List_Num_Entries(mesh->ghvertex);
  else
    return 0;
}

List_ptr MESH_GhostVertex_List(Mesh_ptr mesh) {
  return (mesh->ghvertex);
}
List_ptr MESH_OverlapVertex_List(Mesh_ptr mesh) {
  return (mesh->ovvertex);
}

int MESH_Num_GhostEdges(Mesh_ptr mesh) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R4)) {
    MSTK_Report("MESH_Num_GhostEdges",
		"No ghost edges in reduced representations",MSTK_WARN);
    return 0;
  }
  if (mesh->ghedge)
    return List_Num_Entries(mesh->ghedge);
  else
    return 0;
}
int MESH_Num_OverlapEdges(Mesh_ptr mesh) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R4)) {
    MSTK_Report("MESH_Num_OverlapEdges",
		"No ghost edges in reduced representations",MSTK_WARN);
    return 0;
  }
  if (mesh->ovedge)
    return List_Num_Entries(mesh->ovedge);
  else
    return 0;
}
int MESH_Num_InteriorEdges(Mesh_ptr mesh) {
  RepType rtype;
  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R4)) {
    MSTK_Report("MESH_Num_InteriorEdges",
		"No interior edges in reduced representations",MSTK_WARN);
    return 0;
  }
  else 
    return List_Num_Entries(mesh->medge) -\
      List_Num_Entries(mesh->ovedge)- \
      List_Num_Entries(mesh->ghedge);
}

List_ptr MESH_GhostEdge_List(Mesh_ptr mesh) {
  return (mesh->ghedge);
}
List_ptr MESH_OverlapEdge_List(Mesh_ptr mesh) {
  return (mesh->ovedge);
}


int MESH_Num_GhostFaces(Mesh_ptr mesh) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R2)) {
    if ((mesh->nr == 0) && (mesh->ghface))  /* 2D mesh */
      return List_Num_Entries(mesh->ghface);   
    else {
      MSTK_Report("MESH_Num_GhostFaces",
		  "No ghost faces in reduced representation",MSTK_WARN);
      return 0;
    }
  }
  if (mesh->ghface)
    return List_Num_Entries(mesh->ghface);
  else
    return 0;
}
int MESH_Num_OverlapFaces(Mesh_ptr mesh) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R2)) {
    if ((mesh->nr == 0) && (mesh->ovface))    /* 2D mesh */
      return List_Num_Entries(mesh->ovface);   
    else {
      MSTK_Report("MESH_Num_OverlapFaces",
		  "No ghost faces in reduced representation",MSTK_WARN);
      return 0;
    }
  }
  if (mesh->ovface)
    return List_Num_Entries(mesh->ovface);
  else
    return 0;
}

int MESH_Num_InteriorFaces(Mesh_ptr mesh) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R2)) {
    if ((mesh->nr == 0) && (mesh->mface))    /* 2D mesh */
      return List_Num_Entries(mesh->mface)\
	-List_Num_Entries(mesh->ovface)\
	-List_Num_Entries(mesh->ghface);   
    else {
      MSTK_Report("MESH_Num_InteriorFaces",
		  "No interior faces in reduced representation",MSTK_WARN);
      return 0;
    }
  }
  if (mesh->mface)
    return List_Num_Entries(mesh->mface)\
      -List_Num_Entries(mesh->ovface)\
      -List_Num_Entries(mesh->ghface);   
  else
    return 0;
}

List_ptr MESH_GhostFace_List(Mesh_ptr mesh) {
  return (mesh->ghface);
}
List_ptr MESH_OverlapFace_List(Mesh_ptr mesh) {
  return (mesh->ovface);
}

int MESH_Num_GhostRegions(Mesh_ptr mesh) {
  if (mesh->ghregion)
    return List_Num_Entries(mesh->ghregion);
  else
    return 0;
}
int MESH_Num_OverlapRegions(Mesh_ptr mesh) {
  if(mesh->ovregion)
    return List_Num_Entries(mesh->ovregion);
  else
    return 0;
}
int MESH_Num_InteriorRegions(Mesh_ptr mesh) {
  if(mesh->mregion)
    return List_Num_Entries(mesh->mregion)\
      -List_Num_Entries(mesh->ovregion)
      -List_Num_Entries(mesh->ghregion);
  else
    return 0;
}
List_ptr MESH_GhostRegion_List(Mesh_ptr mesh) {
  return (mesh->ghregion);
}
List_ptr MESH_OverlapRegion_List(Mesh_ptr mesh) {
  return (mesh->ovregion);
}


MVertex_ptr MESH_GhostVertex(Mesh_ptr mesh, int i) {
  return (MVertex_ptr) List_Entry(mesh->ghvertex, i);
}
MVertex_ptr MESH_OverlapVertex(Mesh_ptr mesh, int i) {
  return (MVertex_ptr) List_Entry(mesh->ovvertex, i);
}
MEdge_ptr MESH_GhostEdge(Mesh_ptr mesh, int i) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R4)) {
    MSTK_Report("MESH_GhostEdge",
		"No ghost edges in reduced representation",MSTK_FATAL);
  }
  return (MEdge_ptr) List_Entry(mesh->ghedge, i);
}

MEdge_ptr MESH_OverlapEdge(Mesh_ptr mesh, int i) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R4)) {
    MSTK_Report("MESH_OverlapEdge",
		"No ghost edges in reduced representation",MSTK_FATAL);
  }
  return (MEdge_ptr) List_Entry(mesh->ovedge, i);
}

MFace_ptr MESH_GhostFace(Mesh_ptr mesh, int i) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R2)) {
    if (mesh->nr) {
      MSTK_Report("MESH_GhostFace",
		  "No ghost faces in reduced representation",MSTK_FATAL);
    }
  }

  return (MFace_ptr) List_Entry(mesh->ghface, i);
}

MFace_ptr MESH_OverlapFace(Mesh_ptr mesh, int i) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R2)) {
    if (mesh->nr) {
      MSTK_Report("MESH_OverlapFace",
		  "No ghost faces in reduced representation",MSTK_FATAL);
    }
  }

  return (MFace_ptr) List_Entry(mesh->ovface, i);
}

MRegion_ptr MESH_GhostRegion(Mesh_ptr mesh, int i) {
  return (MRegion_ptr) List_Entry(mesh->ghregion, i);
}
MRegion_ptr MESH_OverlapRegion(Mesh_ptr mesh, int i) {
  return (MRegion_ptr) List_Entry(mesh->ovregion, i);
}

MVertex_ptr MESH_Next_GhostVertex(Mesh_ptr mesh, int *index) {
  if (mesh->ghvertex)
    return (MVertex_ptr) List_Next_Entry(mesh->ghvertex, index);
  else
    return NULL;
}


MVertex_ptr MESH_Next_OverlapVertex(Mesh_ptr mesh, int *index) {
  if (mesh->ovvertex)
    return (MVertex_ptr) List_Next_Entry(mesh->ovvertex, index);
  else
    return NULL;
}

MVertex_ptr MESH_Next_InteriorVertex(Mesh_ptr mesh, int *index) {
  MVertex_ptr mv;
  if (mesh->mvertex) {
    mv = List_Next_Entry(mesh->mvertex,index);
    while( MV_PType(mv) != PINTERIOR )
      mv = List_Next_Entry(mesh->mvertex,index);
    return mv;
  }
  else
    return NULL;
}

MEdge_ptr MESH_Next_GhostEdge(Mesh_ptr mesh, int *index) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R4)) {
    MSTK_Report("MESH_Next_GhostEdge",
		"No ghost edges in reduced representation",MSTK_WARN);
    return NULL;
  }
  if (mesh->ghedge)
    return (MEdge_ptr) List_Next_Entry(mesh->ghedge, index);
  else
    return NULL;
}
MEdge_ptr MESH_Next_OverlapEdge(Mesh_ptr mesh, int *index) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R4)) {
    MSTK_Report("MESH_Next_OverlapEdge",
		"No ghost edges in reduced representation",MSTK_WARN);
    return NULL;
  }
  if (mesh->ovedge)
    return (MEdge_ptr) List_Next_Entry(mesh->ovedge, index);
  else
    return NULL;
}

MEdge_ptr MESH_Next_InteriorEdge(Mesh_ptr mesh, int *index) {
  RepType rtype;
  MEdge_ptr me;
  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R4)) {
    MSTK_Report("MESH_Next_InteriorEdge",
		"No interior edges in reduced representation",MSTK_WARN);
    return NULL;
  }
  if (mesh->medge) {
    me = List_Next_Entry(mesh->medge,index);
    while( ME_PType(me) != PINTERIOR )
      me = List_Next_Entry(mesh->medge,index);
    return me;
  }
  else
    return NULL;
}

MFace_ptr MESH_Next_GhostFace(Mesh_ptr mesh, int *index) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R2)) {
    if (mesh->nr) {
      MSTK_Report("MESH_Next_GhostFace",
		  "No ghost faces in reduced representation",MSTK_WARN);
      return NULL;
    }
  }
  if (mesh->ghface)
    return (MFace_ptr) List_Next_Entry(mesh->ghface, index);
  else
    return NULL;
}
MFace_ptr MESH_Next_OverlapFace(Mesh_ptr mesh, int *index) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R2)) {
    if (mesh->nr) {
      MSTK_Report("MESH_Next_OverlapFace",
		  "No ghost faces in reduced representation",MSTK_WARN);
      return NULL;
    }
  }
  if (mesh->ovface)
    return (MFace_ptr) List_Next_Entry(mesh->ovface, index);
  else
    return NULL;
}

MFace_ptr MESH_Next_InteriorFace(Mesh_ptr mesh, int *index) {
  RepType rtype;
  MFace_ptr mf;
  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R2)) {
    if (mesh->nr) {
      MSTK_Report("MESH_Next_InteriorFace",
		  "No interior faces in reduced representation",MSTK_WARN);
      return NULL;
    }
  }
  if (mesh->mface) {
    mf = List_Next_Entry(mesh->mface,index);
    while( MF_PType(mf) != PINTERIOR )
      mf = List_Next_Entry(mesh->mface,index);
    return mf;
  }
  else
    return NULL;
}

MRegion_ptr MESH_Next_GhostRegion(Mesh_ptr mesh, int *index) {
  if (mesh->ghregion)
    return (MRegion_ptr) List_Next_Entry(mesh->ghregion, index);
  else
    return NULL;
}
MRegion_ptr MESH_Next_OverlapRegion(Mesh_ptr mesh, int *index) {
  if (mesh->ovregion)
    return (MRegion_ptr) List_Next_Entry(mesh->ovregion, index);
  else
    return NULL;
}
MRegion_ptr MESH_Next_InteriorRegion(Mesh_ptr mesh, int *index) {
  MRegion_ptr mr;

  if (mesh->mregion) {
    mr = List_Next_Entry(mesh->mregion,index);
    while( MR_PType(mr) != PINTERIOR )
      mr = List_Next_Entry(mesh->mregion,index);
    return mr;
  }
  else
    return NULL;
}

void MESH_Add_GhostVertex(Mesh_ptr mesh, MVertex_ptr v) {
  if (mesh->ghvertex == (List_ptr) NULL)
    mesh->ghvertex = List_New(10);

  mesh->ghvertex = List_Add(mesh->ghvertex, (void *) v);
  
  if (MV_ID(v) == 0) {
    (mesh->max_ghvid)++;
    MV_Set_ID(v,mesh->max_ghvid);
  }
  MESH_Mark_ParallelAdj_Stale(mesh);
}
void MESH_Add_OverlapVertex(Mesh_ptr mesh, MVertex_ptr v) {
  if (mesh->ovvertex == (List_ptr) NULL)
    mesh->ovvertex = List_New(10);
  mesh->ovvertex = List_Add(mesh->ovvertex, (void *) v);
  MESH_Mark_ParallelAdj_Stale(mesh);
}

void MESH_Add_GhostEdge(Mesh_ptr mesh, MEdge_ptr e){
  /* Have to check if edges exist in this type of representation */
  if (mesh->reptype >= R1 && mesh->reptype <= MSTK_MAXREP)
    return;

  if (mesh->ghedge == (List_ptr) NULL)
    mesh->ghedge = List_New(10);

  mesh->ghedge = List_Add(mesh->ghedge, (void *) e);

  if (ME_ID(e) == 0) {
    (mesh->max_gheid)++;
    ME_Set_ID(e,mesh->max_gheid);
  }
  MESH_Mark_ParallelAdj_Stale(mesh);
}    
void MESH_Add_OverlapEdge(Mesh_ptr mesh, MEdge_ptr e){
  /* Have to check if edges exist in this type of representation */
  if (mesh->reptype >= R1 && mesh->reptype <= MSTK_MAXREP)
    return;

  if (mesh->ovedge == (List_ptr) NULL)
    mesh->ovedge = List_New(10);

  mesh->ovedge = List_Add(mesh->ovedge, (void *) e);
  MESH_Mark_ParallelAdj_Stale(mesh);
}    
     
void MESH_Add_GhostFace(Mesh_ptr mesh, MFace_ptr f){
  /* Have to check if faces exist in this type of representation */
  if (mesh->nr && (mesh->reptype == R1 || mesh->reptype == R2))
    return;

  if ((mesh->reptype == R4) && (MF_Region(f,0) || MF_Region(f,1))) {
#ifdef DEBUG
    MSTK_Report("MESH_Add_Face","Can add disconnected faces only",MSTK_ERROR);
#endif
    return;
  }

  if (mesh->ghface == (List_ptr) NULL)
    mesh->ghface = List_New(10);
  
  mesh->ghface = List_Add(mesh->ghface, (void *) f);

  if (MF_ID(f) == 0) {
    (mesh->max_ghfid)++;
    MF_Set_ID(f,mesh->max_ghfid);
  }
  MESH_Mark_ParallelAdj_Stale(mesh);
}    

void MESH_Add_OverlapFace(Mesh_ptr mesh, MFace_ptr f){
  /* Have to check if faces exist in this type of representation */
  if (mesh->nr && (mesh->reptype == R1 || mesh->reptype == R2))
    return;

  if ((mesh->reptype == R4) && (MF_Region(f,0) || MF_Region(f,1))) {
#ifdef DEBUG
    MSTK_Report("MESH_Add_Face","Can add disconnected faces only",MSTK_ERROR);
#endif
    return;
  }

  if (mesh->ovface == (List_ptr) NULL)
    mesh->ovface = List_New(10);
  
  mesh->ovface = List_Add(mesh->ovface, (void *) f);
  MESH_Mark_ParallelAdj_Stale(mesh);
}    
     
void MESH_Add_GhostRegion(Mesh_ptr mesh, MRegion_ptr r){
  if (mesh->ghregion == (List_ptr) NULL)
    mesh->ghregion = List_New(10);

  mesh->ghregion = List_Add(mesh->ghregion, (void *) r);

  if (MR_ID(r) == 0) {
    (mesh->max_ghrid)++;
    MR_Set_ID(r,mesh->max_ghrid);
  }
  MESH_Mark_ParallelAdj_Stale(mesh);
}    
void MESH_Add_OverlapRegion(Mesh_ptr mesh, MRegion_ptr r){
  if (mesh->ovregion == (List_ptr) NULL)
    mesh->ovregion = List_New(10);
  mesh->ovregion = List_Add(mesh->ovregion, (void *) r);
  MESH_Mark_ParallelAdj_Stale(mesh);
}    
     
void MESH_Rem_GhostVertex(Mesh_ptr mesh, MVertex_ptr v) {
  int fnd=0, i, id;

  if (mesh->ghvertex == (List_ptr) NULL) {
#ifdef DEBUG
    MSTK_Report("Mesh_Rem_Vertex","No vertices in mesh to remove", MSTK_ERROR);
#endif
    return;
  }

  /* If the list of vertices has not been compressed or the mesh has not
     been renumbered, the real position of the entry will be ID-1. If it
     is then delete it directly */

  id = MV_ID(v);

  i = id-1;
  if (List_Entry_Raw(mesh->mvertex,i) == v) {
    List_Remi_Raw(mesh->mvertex,i);
    fnd = 1;
  }

  if (!fnd)
    fnd = List_RemSorted(mesh->mvertex,v,&(MV_ID));

  if (!fnd)
    MSTK_Report("MESH_Rem_GhostVertex","Vertex not found in mesh vertex list",MSTK_FATAL);

  mesh->nv = List_Num_Entries(mesh->mvertex);

  fnd = List_RemSorted(mesh->ghvertex,v,&(MV_GlobalID));
  if (!fnd)
    MSTK_Report("MESH_Rem_GhostVertex","Vertex not found in ghost vertex list",MSTK_FATAL);

  MESH_Mark_ParallelAdj_Stale(mesh);
}    
     
void MESH_Rem_GhostEdge(Mesh_ptr mesh, MEdge_ptr e) {
  int fnd=0, i, id;

  if (mesh->ghedge == (List_ptr) NULL) {
#ifdef DEBUG
    MSTK_Report("Mesh_Rem_Edge","No Edges in mesh to remove",MSTK_ERROR);
#endif
    return;
  }

  /* If the list of edges has not been compressed or the mesh has not
     been renumbered, the real position of the entry will be ID-1. If it
     is then delete it directly */

  id = ME_ID(e);

  i = id-1;
  if (List_Entry_Raw(mesh->medge,i) == e) {
    List_Remi_Raw(mesh->medge,i);
    fnd = 1;
  }
  
  if (!fnd)
    fnd = List_RemSorted(mesh->medge,e,&(ME_ID));

  if (!fnd)
    MSTK_Report("MESH_Rem_GhostEdge","Edge not found in mesh edge list",MSTK_FATAL);

  mesh->ne = List_Num_Entries(mesh->medge);

  fnd = List_RemSorted(mesh->ghedge,e,&(ME_GlobalID));
  if (!fnd)
    MSTK_Report("MESH_Rem_GhostEdge","Edge not found in ghost edge list",MSTK_FATAL);

  MESH_Mark_ParallelAdj_Stale(mesh);
}    
     
void MESH_Rem_GhostFace(Mesh_ptr mesh, MFace_ptr f){
  int fnd=0, i, id;

  if (mesh->ghface == (List_ptr) NULL) {
#ifdef DEBUG
    MSTK_Report("Mesh_Rem_GhostFace","No Faces in mesh to remove",MSTK_ERROR);
#endif
    return;
  }

  /* If the list of faces has not been compressed or the mesh has not
     been renumbered, the real position of the entry will be ID-1. If it
     is then delete it directly */
  
  id = MF_ID(f);

  i = id-1;
  if (List_Entry_Raw(mesh->mface,i) == f) {
    List_Remi_Raw(mesh->mface,i);
    fnd = 1;
  }

  if (!fnd)
    fnd = List_RemSorted(mesh->mface,f,&(MF_ID));

  if (!fnd)
    MSTK_Report("MESH_Rem_Face","Face not found in mesh face list",MSTK_FATAL);

  mesh->nf = List_Num_Entries(mesh->mface);

  fnd = List_RemSorted(mesh->ghface,f,&(MF_GlobalID));
  if (!fnd)
    MSTK_Report("MESH_Rem_Face","Face not found in ghost list",MSTK_FATAL);

  MESH_Mark_ParallelAdj_Stale(mesh);
}    
     
void MESH_Rem_GhostRegion(Mesh_ptr mesh, MRegion_ptr r){
  int fnd=0, i, id;

  if (mesh->ghregion == (List_ptr) NULL) {
    MSTK_Report("Mesh_Rem_GhostRegion","No regions in mesh to remove",MSTK_ERROR);
    return;
  }

  /* If the list of regions has not been compressed or the mesh has not
     been renumbered, the real position of the entry will be ID-1. If it
     is then delete it directly */

  
  id = MR_ID(r);

  i = id-1;
  if (List_Entry_Raw(mesh->mregion,i) == r) {
    List_Remi_Raw(mesh->mregion,i);
    fnd = 1;
  }

  if (!fnd)
    fnd = List_RemSorted(mesh->mregion,r,&(MR_ID));

  if (!fnd)
    MSTK_Report("MESH_Rem_GhostRegion","Region not found in mesh region list",MSTK_FATAL);

  mesh->nr = List_Num_Entries(mesh->mregion);

  fnd = List_RemSorted(mesh->ghregion,r,&(MR_GlobalID));
  if (!fnd)
    MSTK_Report("MESH_Rem_GhostRegion","Region not found in ghost region list",MSTK_FATAL);

  MESH_Mark_ParallelAdj_Stale(mesh);
}    


/* using binary search to find vertex on ghost list based on global_id */
  MVertex_ptr MESH_GhostVertexFromGlobalID(Mesh_ptr mesh, int global_id) {
    MVertex_ptr mv;
    int idx = 0;
    while(  (mv = MESH_Next_GhostVertex(mesh,&idx)) ){
	if (MV_GlobalID(mv) == global_id)
	  return mv;
    }
    idx = 0;
    while(  (mv = MESH_Next_OverlapVertex(mesh,&idx)) ){
	if (MV_GlobalID(mv) == global_id)
	  return mv;
    }

    /*
    MVertex_ptr *loc_mv, mv = MV_New(NULL);
    int iloc;
    MV_Set_GlobalID(mv,global_id);
    //    printf("global id %d, num of gh vertex %d\n",global_id, MESH_Num_GhostVertices(mesh));
    loc_mv = (MVertex_ptr *)bsearch(&mv,
				    List_Entries(mesh->ghvertex),
				    MESH_Num_GhostVertices(mesh),
				    sizeof(MVertex_ptr),
				    compareGlobalID);
    if(loc_mv) {
      // printf("found global id %d\n",MV_GlobalID(*loc_mv));
    return *loc_mv;
    }
    //      MV_Delete(mv,0);
    */
    return NULL;
  }


#endif /* MSTK_HAVE_MPI */


     
void MESH_Set_GModel(Mesh_ptr mesh, GModel_ptr geom){
  mesh->geom = geom;
}    
     
void MESH_SetRepTypeIni(Mesh_ptr mesh, RepType reptype){
  if (!(mesh->reptype >= R1 && mesh->reptype <= R4)) {
    if (reptype >= R1 && reptype <= R4) {
      mesh->hedge = Hash_New(0, 1);
    }
  }
  if (!(mesh->reptype >= R1 && mesh->reptype <= R2)) {
    if (reptype >= R1 && reptype <= R2) {
      mesh->hface = Hash_New(0, 1);
    }
  }
  mesh->reptype = reptype;
}    
     
int  MESH_ChangeRepType(Mesh_ptr mesh, RepType nureptype){
  /* Only certain types are allowed */
  return 1;
}


  /* Extra functionality for hash-tables */

  Hash_ptr MESH_Hash_Edges(Mesh_ptr mesh) {
    return mesh->hedge;
  }
  
  Hash_ptr MESH_Hash_Faces(Mesh_ptr mesh) {
    return mesh->hface;
  }

  int MESH_AutoLock(Mesh_ptr mesh) {
     return mesh->autolock;
  } 

  void MESH_Set_AutoLock(Mesh_ptr mesh, int autolock) {
     mesh->autolock = autolock;
  } 

#ifdef __cplusplus
}
#endif
