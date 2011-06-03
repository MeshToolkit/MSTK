#define _H_Mesh_Private

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Mesh.h"
#include "MSTK.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"


#ifdef __cplusplus
extern "C" {
#endif

  RepType MESH_rtype[5] = {F1, F4, R1, R2, R4};
  char MESH_rtype_str[5][3] = {"F1\0","F4\0","R1\0","R2\0","R4\0"};

Mesh_ptr MESH_New(RepType type) {
  Mesh_ptr newmesh;

  newmesh = (Mesh_ptr) MSTK_malloc(sizeof(Mesh));

  newmesh->reptype = type;
  newmesh->nv = newmesh->ne = newmesh->nf = newmesh->nr = 0;
  newmesh->mvertex = (List_ptr) NULL;
  newmesh->medge = (List_ptr) NULL;
  newmesh->mface = (List_ptr) NULL;
  newmesh->mregion = (List_ptr) NULL;

#ifdef MSTK_HAVE_MPI
  newmesh->ghvertex = (List_ptr) NULL;
  newmesh->ghedge = (List_ptr) NULL;
  newmesh->ghface = (List_ptr) NULL;
  newmesh->ghregion = (List_ptr) NULL;
  newmesh->ovvertex = (List_ptr) NULL;
  newmesh->ovedge = (List_ptr) NULL;
  newmesh->ovface = (List_ptr) NULL;
  newmesh->ovregion = (List_ptr) NULL;
  newmesh->global_info = NULL;
  newmesh->local_info = NULL;
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

#ifdef MSTK_HAVE_MPI
  newmesh->max_ghvid = newmesh->max_gheid = newmesh->max_ghfid = newmesh->max_ghrid = 0;
#endif

  return newmesh;
}

#ifdef MSTK_HAVE_MPI
int* MESH_GlobalInfo(Mesh_ptr mesh) {
  return mesh->global_info;
}

int* MESH_LocalInfo(Mesh_ptr mesh) {
  return mesh->local_info;
}

void MESH_Set_GlobalInfo(Mesh_ptr mesh, int *global_info) {
  mesh->global_info = global_info;
}

void MESH_Set_LocalInfo(Mesh_ptr mesh, int *local_info) {
  mesh->local_info = local_info;
}
#endif

void MESH_Delete(Mesh_ptr mesh) {
  int i, nv, ne, nf, nr;
  MVertex_ptr mv, ghv;
  MEdge_ptr me, ghe;
  MFace_ptr mf, ghf;
  MRegion_ptr mr, ghr;
  MAttrib_ptr attrib;

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

  /*
  if (mesh->ghregion) {
    i = 0;
    while ((ghr = List_Next_Entry(mesh->ghregion,&i))) {
      MR_Destroy_For_MESH_Delete(ghr);
    }
    List_Delete(mesh->ghregion);
  }
  if (mesh->ghface) {
    i = 0;
    while ((ghf = List_Next_Entry(mesh->ghface,&i))) {
      MF_Destroy_For_MESH_Delete(ghf);
    }
    List_Delete(mesh->ghface);
  }
  if (mesh->ghedge) {
    i = 0;
    while ((ghe = List_Next_Entry(mesh->ghedge,&i))) {
      ME_Destroy_For_MESH_Delete(ghe);
    }
    List_Delete(mesh->ghedge);
  }
  if (mesh->ghvertex) {
    i = 0;
    while ((ghv = List_Next_Entry(mesh->ghvertex,&i))) {
      MV_Destroy_For_MESH_Delete(ghv);
    }
    List_Delete(mesh->ghvertex);
  }
  */
#endif /* MSTK_HAVE_MPI */
  
  if (mesh->AttribList) {
    i = 0;
    while ((attrib = List_Next_Entry(mesh->AttribList,&i)))
      MAttrib_Delete(attrib);
    List_Delete(mesh->AttribList);
  }

  MSTK_free(mesh);
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
  char *rstr = (char *) MSTK_malloc(3*sizeof(char));

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
  if (mesh->AttribList)
    List_Rem(mesh->AttribList,attrib);
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
  MSTK_Report("Mesh_FillHash_Edges","Inefficient to call routines like MESH_Num_Edges, MESH_Edge, MESH_Next_Edge with this representation\n",WARN);
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
  MSTK_Report("Mesh_FillHash_Edges","Inefficient to call routines like MESH_Num_Faces, MESH_Face, MESH_Next_Face with this representation\n",WARN);
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
    MSTK_Report("Mesh_Vertex","Non-existent vertex requested\n",ERROR);
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
    MSTK_Report("Mesh_Edge","Non-existent edge requested\n",ERROR);
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
    MSTK_Report("Mesh_Face","Non-existent face requested\n",ERROR);
#endif
    return (MFace_ptr) NULL;
  }
  else
    return (MFace_ptr) List_Entry(mesh->mface, i);
}

MRegion_ptr MESH_Region(Mesh_ptr mesh, int i) {

  if (i >= mesh->nr) {
#ifdef DEBUG
    MSTK_Report("Mesh_Region","Non-existent region requested\n",ERROR);
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
     ID. If not, search before this (ID-1)'th upto the beginning of
     the list. Chances are some elements got deleted and the list got
     compressed. If still not found, search after the (ID-1)'th entry
     upto the end of the list. This should be quite efficient for
     static meshes and modestly efficient for dynamic meshes. However,
     if we use a different data structure to store mesh entities (like
     a tree), the efficiency may decrease. So, use with care!!!! */


MVertex_ptr MESH_VertexFromID(Mesh_ptr mesh, int id) {
  int istart, j;
  MVertex_ptr mv;

  if (id < 1)
    return NULL;

  istart = id-1;
  if (istart < mesh->nv) {
    mv = (MVertex_ptr) List_Entry(mesh->mvertex,istart);
    if (MV_ID(mv) == id)
      return mv;
  }
  else
    istart = mesh->nv-1;
  
    
  for (j = istart; j >= 0; j--) {
    mv = (MVertex_ptr) List_Entry(mesh->mvertex,j);
    if (MV_ID(mv) == id)
      return mv;
  }

  for (j = istart; j < mesh->nv; j++) {
    mv = (MVertex_ptr) List_Entry(mesh->mvertex,j);
    if (MV_ID(mv) == id)
      return mv;
  }

  return NULL;
}

MEdge_ptr MESH_EdgeFromID(Mesh_ptr mesh, int id) {
  int istart, j;
  MEdge_ptr me;

  if (id < 1)
    return NULL;

  istart = id-1;
  if (istart < mesh->ne) {
    me = (MEdge_ptr) List_Entry(mesh->medge,istart);
    if (ME_ID(me) == id)
      return me;
  }
  else
    istart = mesh->ne-1;
  
    
  for (j = istart; j >= 0; j--) {
    me = (MEdge_ptr) List_Entry(mesh->medge,j);
    if (ME_ID(me) == id)
      return me;
  }

  for (j = istart; j < mesh->ne; j++) {
    me = (MEdge_ptr) List_Entry(mesh->medge,j);
    if (ME_ID(me) == id)
      return me;
  }

  return NULL;
}

MFace_ptr MESH_FaceFromID(Mesh_ptr mesh, int id) {
  int istart, j;
  MFace_ptr mf;

  if (id < 1)
    return NULL;

  istart = id-1;
  if (istart < mesh->nf) {
    mf = (MFace_ptr) List_Entry(mesh->mface,istart);
    if (MF_ID(mf) == id)
      return mf;
  }
  else
    istart = mesh->nf-1;
  
    
  for (j = istart; j >= 0; j--) {
    mf = (MFace_ptr) List_Entry(mesh->mface,j);
    if (MF_ID(mf) == id)
      return mf;
  }

  for (j = istart; j < mesh->nf; j++) {
    mf = (MFace_ptr) List_Entry(mesh->mface,j);
    if (MF_ID(mf) == id)
      return mf;
  }

  return NULL;
}

MRegion_ptr MESH_RegionFromID(Mesh_ptr mesh, int id) {
  int istart, j;
  MRegion_ptr mr;

  if (id < 1)
    return NULL;

  istart = id-1;
  if (istart < mesh->nr) {
    mr = (MRegion_ptr) List_Entry(mesh->mregion,istart);
    if (MR_ID(mr) == id)
      return mr;
  }
  else
    istart = mesh->nr-1;
  
    
  for (j = istart; j >= 0; j--) {
    mr = (MRegion_ptr) List_Entry(mesh->mregion,j);
    if (MR_ID(mr) == id)
      return mr;
  }

  for (j = istart; j < mesh->nr; j++) {
    mr = (MRegion_ptr) List_Entry(mesh->mregion,j);
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
    MSTK_Report("MESH_EntityFromID","Unrecognized entity type",ERROR);
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
    MSTK_Report("MESH_Add_Face","Can add disconnected faces only",ERROR);
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
    MSTK_Report("Mesh_Rem_Vertex","No vertices in mesh to remove", ERROR);
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
    MSTK_Report("MESH_Rem_Vertex","Vertex not found in list",FATAL);

  mesh->nv = List_Num_Entries(mesh->mvertex);

  return;
}    
     
void MESH_Rem_Edge(Mesh_ptr mesh, MEdge_ptr e) {
  int fnd=0, i, id;

  if (mesh->medge == (List_ptr) NULL) {
#ifdef DEBUG
    MSTK_Report("Mesh_Rem_Edge","No Edges in mesh to remove",ERROR);
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
    MSTK_Report("MESH_Rem_Edge","Edge not found in list",FATAL);

  mesh->ne = List_Num_Entries(mesh->medge);

  return;
}    
     
void MESH_Rem_Face(Mesh_ptr mesh, MFace_ptr f){
  int fnd=0, i, id;

  if (mesh->mface == (List_ptr) NULL) {
#ifdef DEBUG
    MSTK_Report("Mesh_Rem_Face","No Faces in mesh to remove",ERROR);
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
    MSTK_Report("MESH_Rem_Face","Face not found in list",FATAL);

  mesh->nf = List_Num_Entries(mesh->mface);

  return;
}    
     
void MESH_Rem_Region(Mesh_ptr mesh, MRegion_ptr r){
  int fnd=0, i, id;

  if (mesh->mregion == (List_ptr) NULL) {
    MSTK_Report("Mesh_Rem_Region","No regions in mesh to remove",ERROR);
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
    MSTK_Report("MESH_Rem_Region","Region not found in list",FATAL);

  mesh->nr = List_Num_Entries(mesh->mregion);
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

  int MESH_Sort_GhostLists(Mesh_ptr mesh, 
	       int (*compfunc)(MEntity_ptr, MEntity_ptr)) {

    if (mesh->ghvertex)
      List_Sort(mesh->ghvertex,List_Num_Entries(mesh->ghvertex),sizeof(MVertex_ptr),compfunc);
    if (mesh->ghedge)
      List_Sort(mesh->ghedge,List_Num_Entries(mesh->ghedge),sizeof(MEdge_ptr),compfunc);
    if (mesh->ghface)
      List_Sort(mesh->ghface,List_Num_Entries(mesh->ghface),sizeof(MFace_ptr),compfunc);
    if (mesh->ghregion)
      List_Sort(mesh->ghregion,List_Num_Entries(mesh->ghregion),sizeof(MRegion_ptr),compfunc);

    if (mesh->ovvertex)
      List_Sort(mesh->ovvertex,List_Num_Entries(mesh->ovvertex),sizeof(MVertex_ptr),compfunc);
    if (mesh->ovedge)
      List_Sort(mesh->ovedge,List_Num_Entries(mesh->ovedge),sizeof(MEdge_ptr),compfunc);
    if (mesh->ovface)
      List_Sort(mesh->ovface,List_Num_Entries(mesh->ovface),sizeof(MFace_ptr),compfunc);
    if (mesh->ovregion)
      List_Sort(mesh->ovregion,List_Num_Entries(mesh->ovregion),sizeof(MRegion_ptr),compfunc);
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

int MESH_Num_GhostEdges(Mesh_ptr mesh) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R4)) {
    MSTK_Report("MESH_Num_GhostEdges",
		"No ghost edges in reduced representations",WARN);
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
		"No ghost edges in reduced representations",WARN);
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
		"No interior edges in reduced representations",WARN);
    return 0;
  }
  else 
    return List_Num_Entries(mesh->medge) -\
      List_Num_Entries(mesh->ovedge)- \
      List_Num_Entries(mesh->ghedge);
}

int MESH_Num_GhostFaces(Mesh_ptr mesh) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R2)) {
    if ((mesh->nr == 0) && (mesh->ghface))  /* 2D mesh */
      return List_Num_Entries(mesh->ghface);   
    else {
      MSTK_Report("MESH_Num_GhostFaces",
		  "No ghost faces in reduced representation",WARN);
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
		  "No ghost faces in reduced representation",WARN);
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
		  "No interior faces in reduced representation",WARN);
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
		"No ghost edges in reduced representation",FATAL);
  }
  return (MEdge_ptr) List_Entry(mesh->ghedge, i);
}

MEdge_ptr MESH_OverlapEdge(Mesh_ptr mesh, int i) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R4)) {
    MSTK_Report("MESH_OverlapEdge",
		"No ghost edges in reduced representation",FATAL);
  }
  return (MEdge_ptr) List_Entry(mesh->ovedge, i);
}

MFace_ptr MESH_GhostFace(Mesh_ptr mesh, int i) {
  RepType rtype;

  rtype = MESH_RepType(mesh);
  if ((rtype >= R1) && (rtype <= R2)) {
    if (mesh->nr) {
      MSTK_Report("MESH_GhostFace",
		  "No ghost faces in reduced representation",FATAL);
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
		  "No ghost faces in reduced representation",FATAL);
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
		"No ghost edges in reduced representation",WARN);
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
		"No ghost edges in reduced representation",WARN);
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
		"No interior edges in reduced representation",WARN);
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
		  "No ghost faces in reduced representation",WARN);
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
		  "No ghost faces in reduced representation",WARN);
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
		  "No interior faces in reduced representation",WARN);
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
}
void MESH_Add_OverlapVertex(Mesh_ptr mesh, MVertex_ptr v) {
  if (mesh->ovvertex == (List_ptr) NULL)
    mesh->ovvertex = List_New(10);
  mesh->ovvertex = List_Add(mesh->ovvertex, (void *) v);
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
}    
void MESH_Add_OverlapEdge(Mesh_ptr mesh, MEdge_ptr e){
  /* Have to check if edges exist in this type of representation */
  if (mesh->reptype >= R1 && mesh->reptype <= MSTK_MAXREP)
    return;

  if (mesh->ovedge == (List_ptr) NULL)
    mesh->ovedge = List_New(10);

  mesh->ovedge = List_Add(mesh->ovedge, (void *) e);
}    
     
void MESH_Add_GhostFace(Mesh_ptr mesh, MFace_ptr f){
  /* Have to check if faces exist in this type of representation */
  if (mesh->nr && (mesh->reptype == R1 || mesh->reptype == R2))
    return;

  if ((mesh->reptype == R4) && (MF_Region(f,0) || MF_Region(f,1))) {
#ifdef DEBUG
    MSTK_Report("MESH_Add_Face","Can add disconnected faces only",ERROR);
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
}    

void MESH_Add_OverlapFace(Mesh_ptr mesh, MFace_ptr f){
  /* Have to check if faces exist in this type of representation */
  if (mesh->nr && (mesh->reptype == R1 || mesh->reptype == R2))
    return;

  if ((mesh->reptype == R4) && (MF_Region(f,0) || MF_Region(f,1))) {
#ifdef DEBUG
    MSTK_Report("MESH_Add_Face","Can add disconnected faces only",ERROR);
#endif
    return;
  }

  if (mesh->ovface == (List_ptr) NULL)
    mesh->ovface = List_New(10);
  
  mesh->ovface = List_Add(mesh->ovface, (void *) f);
}    
     
void MESH_Add_GhostRegion(Mesh_ptr mesh, MRegion_ptr r){
  if (mesh->ghregion == (List_ptr) NULL)
    mesh->ghregion = List_New(10);

  mesh->ghregion = List_Add(mesh->ghregion, (void *) r);

  if (MR_ID(r) == 0) {
    (mesh->max_ghrid)++;
    MR_Set_ID(r,mesh->max_ghrid);
  }
}    
void MESH_Add_OverlapRegion(Mesh_ptr mesh, MRegion_ptr r){
  if (mesh->ovregion == (List_ptr) NULL)
    mesh->ovregion = List_New(10);

  mesh->ovregion = List_Add(mesh->ovregion, (void *) r);
}    
     
void MESH_Rem_GhostVertex(Mesh_ptr mesh, MVertex_ptr v) {
  int fnd=0, i, id;

  if (mesh->ghvertex == (List_ptr) NULL) {
#ifdef DEBUG
    MSTK_Report("Mesh_Rem_Vertex","No vertices in mesh to remove", ERROR);
#endif
    return;
  }

  /* If the list of vertices has not been compressed or the mesh has not
     been renumbered, the real position of the entry will be ID-1. If it
     is then delete it directly */

  id = MV_ID(v);

  i = id-1;
  if (List_Entry_Raw(mesh->ghvertex,i) == v) {
    List_Remi_Raw(mesh->ghvertex,i);
    fnd = 1;
  }

  if (!fnd)
    fnd = List_RemSorted(mesh->ghvertex,v,&(MV_GlobalID));

  if (!fnd)
    MSTK_Report("MESH_Rem_GhostVertex","Vertex not found in list",FATAL);

  return;
}    
     
void MESH_Rem_GhostEdge(Mesh_ptr mesh, MEdge_ptr e) {
  int fnd=0, i, id;

  if (mesh->ghedge == (List_ptr) NULL) {
#ifdef DEBUG
    MSTK_Report("Mesh_Rem_Edge","No Edges in mesh to remove",ERROR);
#endif
    return;
  }

  /* If the list of edges has not been compressed or the mesh has not
     been renumbered, the real position of the entry will be ID-1. If it
     is then delete it directly */

  id = ME_ID(e);

  i = id-1;
  if (List_Entry_Raw(mesh->ghedge,i) == e) {
    List_Remi_Raw(mesh->ghedge,i);
    fnd = 1;
  }
  
  if (!fnd)
    fnd = List_RemSorted(mesh->ghedge,e,&(ME_ID));

  if (!fnd)
    MSTK_Report("MESH_Rem_GhostEdge","Edge not found in list",FATAL);

  return;
}    
     
void MESH_Rem_GhostFace(Mesh_ptr mesh, MFace_ptr f){
  int fnd=0, i, id;

  if (mesh->ghface == (List_ptr) NULL) {
#ifdef DEBUG
    MSTK_Report("Mesh_Rem_GhostFace","No Faces in mesh to remove",ERROR);
#endif
    return;
  }

  /* If the list of faces has not been compressed or the mesh has not
     been renumbered, the real position of the entry will be ID-1. If it
     is then delete it directly */
  
  id = MF_ID(f);

  i = id-1;
  if (List_Entry_Raw(mesh->ghface,i) == f) {
    List_Remi_Raw(mesh->ghface,i);
    fnd = 1;
  }

  if (!fnd)
    fnd = List_RemSorted(mesh->ghface,f,&(MF_ID));

  if (!fnd)
    MSTK_Report("MESH_Rem_Face","Face not found in list",FATAL);

  return;
}    
     
void MESH_Rem_GhostRegion(Mesh_ptr mesh, MRegion_ptr r){
  int fnd=0, i, id;

  if (mesh->ghregion == (List_ptr) NULL) {
    MSTK_Report("Mesh_Rem_GhostRegion","No regions in mesh to remove",ERROR);
    return;
  }

  /* If the list of regions has not been compressed or the mesh has not
     been renumbered, the real position of the entry will be ID-1. If it
     is then delete it directly */

  
  id = MR_ID(r);

  i = id-1;
  if (List_Entry_Raw(mesh->ghregion,i) == r) {
    List_Remi_Raw(mesh->ghregion,i);
    fnd = 1;
  }

  if (!fnd)
    fnd = List_RemSorted(mesh->ghregion,r,&(MR_ID));

  if (!fnd)
    MSTK_Report("MESH_Rem_GhostRegion","Region not found in list",FATAL);

  return;
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

/* File name should have the .mstk extension - will not check here */
int MESH_InitFromFile(Mesh_ptr mesh, const char *filename) {
  FILE *fp;
  char inp_rtype[16], temp_str[256], fltype_str[16], rltype_str[16];
  char attname[256], atttype_str[256], attent_str[256];
  int i, j, found, NV=0, NE=0, NF=0, NR=0, nav, nar, gdim, gid;
  int id, dim, vid1, vid2, eid, fid, rid, adjvid, adjrid, adjv_flag;
  int nfv, max_nfv=0, nfe, max_nfe=0, nrv, max_nrv=0, nrf, max_nrf=0;
  int *fedirs, *rfdirs, ival, ncomp, nent, done, status;
  int processed_vertices=0, processed_adjv=0, processed_edges=0;
  int processed_faces=0, processed_regions=0, processed_adjr=0;
  double ver, xyz[3], rval;
  double *rval_arr;
  MVertex_ptr mv, ev1, ev2, adjv, *fverts, *rverts;
  MEdge_ptr me, *fedges;
  MFace_ptr mf, *rfaces;
  MRegion_ptr mr, adjr;
  MEntity_ptr ent;
  MAttrib_ptr attrib;
  MType attent;
  MAttType atttype;
  RepType file_reptype;


  if (!(fp = fopen(filename,"r"))) {
    MSTK_Report("MESH_InitFromFile","Cannot open file",ERROR);
    return 0;
  }

  status = fscanf(fp,"%s %lf",temp_str,&ver);
  if (strcmp(temp_str,"MSTK") != 0) {
    MSTK_Report("MESH_InitFromFile","Not a MSTK file",ERROR);
    fclose(fp);
    return 0;
  }
  if (status == EOF)
    MSTK_Report("MESH_InitFromFile",
		"Premature end of file before any mesh data is read",FATAL);

  if (ver != MSTK_VER) {
    MSTK_Report("MESH_InitFromFile","Version mismatch",WARN);
  }

  status = fscanf(fp,"%s %d %d %d %d\n",inp_rtype,
		  &(mesh->nv),&(mesh->ne),&(mesh->nf),&(mesh->nr));
  if (status == EOF)
    MSTK_Report("MESH_InitFromFile",
		"Premature end of file before any mesh data is read",FATAL);


  found = 0;
  for (i = 0; i < MSTK_MAXREP; i++) {
    if (strncmp(inp_rtype,MESH_rtype_str[i],2) == 0) {
      file_reptype = MESH_rtype[i];
      found = 1;

      if (mesh->reptype == UNKNOWN_REP)
	MESH_SetRepType(mesh,file_reptype); /* function does additional things
					       that a simple assignement to
					       mesh->reptype does not */
      break;
    }
  }

  if (!found) {
    MSTK_Report("MESH_InitFromFile","Unrecognized representation type",ERROR);
    fclose(fp);
    return 0;
  }


  /* For now, the reduced representations only allow region and
     face elements - no edge elements are allowed */

  if (file_reptype >= R1 && file_reptype <= R4) {
    if (mesh->ne)
      MSTK_Report("Mesh_InitFromFile",
		  "Representation does not allow edges to be specified",WARN);
  }


  status = fscanf(fp,"%s",temp_str);
  if (status == EOF)
    MSTK_Report("MESH_InitFromFile",
		"Premature end of file while looking for vertex data",FATAL);
  else if (status == 0)
    MSTK_Report("MESH_InitFromFile",
		"Error in reading vertex data",FATAL);

  if (strncmp(temp_str,"vertices",8) == 0) {

    NV = mesh->nv;
    mesh->mvertex = List_New(NV);

    for (i = 0; i < NV; i++) {
      status = fscanf(fp,"%lf %lf %lf %d %d",&xyz[0],&xyz[1],&xyz[2],&gdim,&gid);
      if (status == EOF)
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file while reading vertices",FATAL);
      else if (status == 0)
	MSTK_Report("MESH_InitFromFile",
		    "Error in reading vertex data",FATAL);

      mv = MV_New(mesh);
      MV_Set_Coords(mv,xyz);
      
      MV_Set_GEntID(mv,gid);
      MV_Set_GEntDim(mv,gdim);
      
      MV_Set_ID(mv,(i+1));
    }

    processed_vertices = 1;
  }
  else {
    MSTK_Report("MESH_InitFromFile",
		"Vertex information should be listed first",ERROR);
    fclose(fp);
    return 0;
  }


  status = fscanf(fp,"%s",temp_str);
  if (status == EOF) {
    if (mesh->ne || mesh->nf || mesh->nr)
      MSTK_Report("MESH_InitFromFile",
		  "Premature end of file after vertex data",FATAL);
    else
      return 0;
  }



  /* ADJACENT VERTEX DATA */

  if (strncmp(temp_str,"adjvertices",11) == 0) {
    adjv_flag = 1;
    if (mesh->reptype == R2 || mesh->reptype == R4) {

      /* Mesh representation supports storing of adjacent vertex info */

      for (i = 0; i < NV; i++) {
	mv = List_Entry(mesh->mvertex, i);
	status = fscanf(fp,"%d",&nav);
	if (status == EOF)
	  MSTK_Report("MESH_InitFromFile",
		      "Premature end of file while reading adjacent vertices",
		      FATAL);
	else if (status == 0)
	  MSTK_Report("MESH_InitFromFile",
		      "Error in reading adjacent vertex data",FATAL);

	for (j = 0; j < nav; j++) {
	  status = fscanf(fp,"%d",&adjvid);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading adjacent vertices"
			,FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading adjacent vertex data",FATAL);

	  adjv = List_Entry(mesh->mvertex,adjvid-1);
#ifdef DEBUG
	  if (MV_ID(adjv) != adjvid)
	    MSTK_Report("MESH_InitFromFile","Adjacent vertex ID mismatch",
			ERROR);
#endif
	  
	  MV_Add_AdjVertex(mv,adjv);
	}
      }
    }
    else {

      /* Mesh representation does not support storing of adjacent
	 vertex information. Just read and discard */

      for (i = 0; i < NV; i++) {
	status = fscanf(fp,"%d",&nav);
	if (status == EOF)
	  MSTK_Report("MESH_InitFromFile",
		      "Premature end of file while reading adjacent vertices",
		      FATAL);
	else if (status == 0)
	  MSTK_Report("MESH_InitFromFile",
		      "Error in reading adjacent vertex data",FATAL);
	
	for (j = 0; j < nav; j++)
	  fscanf(fp,"%d",&adjvid);
      }
    }

    processed_adjv = 1;
  }
  else {
    adjv_flag = 0;
    if (file_reptype == R2 || file_reptype == R4) {
      MSTK_Report("MESH_InitFromFile",
		  "Expected adjacent vertex information",ERROR);
    }
  }
    

  if (processed_adjv) {
    status = fscanf(fp,"%s",temp_str);
    if (status == EOF) {
      if (mesh->ne || mesh->nf || mesh->nr)
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file after adjacent vertex data",FATAL);
      else {
        fclose(fp);
	return 1;
      }
    }
  }



  /* EDGE DATA */

  if (strncmp(temp_str,"edges",5) == 0) {

    NE = mesh->ne;
    mesh->medge = List_New(NE);

    if (mesh->reptype >= F1 && mesh->reptype <= F4) {

      /* Mesh representation supports edges */

      for (i = 0; i < NE; i++) {
	status = fscanf(fp,"%d %d %d %d",&vid1,&vid2,&gdim,&gid);
	if (status == EOF)
	  MSTK_Report("MESH_InitFromFile",
		      "Premature end of file while reading edges",FATAL);
	else if (status == 0)
	  MSTK_Report("MESH_InitFromFile",
		      "Error in reading edge data",FATAL);

	ev1 = List_Entry(mesh->mvertex,vid1-1);
	ev2 = List_Entry(mesh->mvertex,vid2-1);
#ifdef DEBUG
	if (MV_ID(ev1) != vid1 || MV_ID(ev2) != vid2)
	  MSTK_Report("MESH_InitFromFile","Mesh vertex ID mismatch",ERROR);
#endif

	me = ME_New(mesh);

	ME_Set_Vertex(me,0,ev1);
	ME_Set_Vertex(me,1,ev2);

	ME_Set_GEntID(me,gid);
	ME_Set_GEntDim(me,gdim);

	ME_Set_ID(me, i+1);
      }
    }
    else {

      /* Mesh representation does not support edges. Read and discard */

      for (i = 0; i < NE; i++) {
	status = fscanf(fp,"%d %d %d %d",&vid1,&vid2,&gdim,&gid);
	if (status == EOF)
	  MSTK_Report("MESH_InitFromFile",
		      "Premature end of file while reading edges",FATAL);
	else if (status == 0)
	  MSTK_Report("MESH_InitFromFile",
		      "Error in reading edge data",FATAL);
      }
    }

    processed_edges = 1;
  }
  else {
    if (file_reptype >= F1 && file_reptype <= F4) {
      MSTK_Report("MESH_InitFromFile","Expected edge information",ERROR);
      return 0;
    }
  }
  
  if (processed_edges) {
    status = fscanf(fp,"%s",temp_str);
    if (status == EOF) {
      if (mesh->nf || mesh->nr) {
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file after edge data",FATAL);
      }
      else {
        fclose(fp);
	return 1;
      }
    }
  }



  /* FACE DATA */

  /* It is possible for someone to specify faces for a solid mesh specified 
   in the R1 format - do we have to idiot-proof? */

  if (strncmp(temp_str,"face",4) == 0) {
    if (strncmp(temp_str,"faces",5) != 0) 
      MSTK_Report("MESH_InitFromFile","Expected keyword \"faces\"",ERROR);
    
    NF = mesh->nf;
    mesh->mface = List_New(NF);

    status = fscanf(fp,"%s",fltype_str);
    if (status == EOF)
      MSTK_Report("MESH_InitFromFile",
		  "Premature end of file while reading faces",FATAL);

    if (file_reptype >= R1 && file_reptype <= R4) {

      if (strncmp(fltype_str,"vertex",6) == 0) {      

	/* Mesh representation supports faces for surface meshes */

	fverts = NULL;
	for (i = 0; i < NF; i++) {
	  mf = MF_New(mesh);

	  status = fscanf(fp,"%d",&nfv);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading faces",FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading face data",FATAL);

	  if (fverts) {
	    if (nfv > max_nfv) {
	      max_nfv = nfv;
	      fverts = MSTK_realloc(fverts,max_nfv*sizeof(MVertex_ptr));
	    }
	  }
	  else {
	    max_nfv = nfv;
	    fverts = MSTK_malloc(nfv*sizeof(MVertex_ptr));
	  }
	  
	  for (j = 0; j < nfv; j++) {
	    status = fscanf(fp,"%d",&vid1);
	    if (status == EOF)
	      MSTK_Report("MESH_InitFromFile",
			  "Premature end of file while reading face data",
			  FATAL);
	    else if (status == 0)
	      MSTK_Report("MESH_InitFromFile",
			  "Error in reading edge data",FATAL);

	    fverts[j] = List_Entry(mesh->mvertex,vid1-1);
#ifdef DEBUG
	    if (MV_ID(fverts[j]) != vid1)
	      MSTK_Report("MESH_InitFromFile","Mesh vertex ID mismatch",ERROR);
#endif
	  }
	  MF_Set_Vertices(mf,nfv,fverts);  

	  status = fscanf(fp,"%d %d",&gdim,&gid);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading face data",FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading face data",FATAL);
	  
	  MF_Set_GEntID(mf,gid);
	  MF_Set_GEntDim(mf,gdim);
	  
	  MF_Set_ID(mf,i+1);
	}
	if (fverts)
	  MSTK_free(fverts);

	processed_faces = 1;
      }
      else {
	MSTK_Report("MESH_InitFromFile",
		    "Expected face description in terms of vertices",ERROR);
        fclose(fp);
	return 0;
      }

    }
    else if (file_reptype >= F1 && file_reptype <= F4) { 

      if (strncmp(fltype_str,"edge",4) == 0) {

	fedges = NULL;	fedirs = NULL;
	for (i = 0; i < NF; i++) {
	  mf = MF_New(mesh);

	  status = fscanf(fp,"%d",&nfe);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading face edge data",
			FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading face data",FATAL);

	  if (fedges) {
	    if (nfe > max_nfe) {
	      max_nfe = nfe;
	      fedges = MSTK_realloc(fedges,max_nfe*sizeof(MVertex_ptr));
	      fedirs = MSTK_realloc(fedirs,max_nfe*sizeof(int));
	    }
	  }
	  else {
	    max_nfe = nfe;
	    fedges = MSTK_malloc(max_nfe*sizeof(MVertex_ptr));
	    fedirs = MSTK_malloc(max_nfe*sizeof(MVertex_ptr));
	  }
	  
	  for (j = 0; j < nfe; j++) {
	    status = fscanf(fp,"%d",&eid);
	    if (status == EOF)
	      MSTK_Report("MESH_InitFromFile",
			  "Premature end of file while reading face edge data",
			  FATAL);
	    else if (status == 0)
	      MSTK_Report("MESH_InitFromFile",
			  "Error in reading face data",FATAL);

	    fedirs[j] = eid > 0 ? 1 : 0;
	    fedges[j] = List_Entry(mesh->medge,abs(eid)-1);
#ifdef DEBUG
	    if (ME_ID(fedges[j]) != abs(eid))
	      MSTK_Report("MESH_InitFromFile","Mesh edge ID mismatch",ERROR);
#endif
	  }
	  
	  if (mesh->reptype >= F1 && mesh->reptype <= F4) {

	    MF_Set_Edges(mf,nfe,fedges,fedirs);

	  }
	  else {

	    MSTK_Report("MESH_InitFromFile","Need to convert edge list to vertex list",ERROR);
	    fclose(fp);
	    return 0;

	  }
	  
	  status = fscanf(fp,"%d %d",&gdim,&gid);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading faces",FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading face data",FATAL);
	  
	  MF_Set_GEntDim(mf,gdim);
	  MF_Set_GEntID(mf,gid);
	  
	  MF_Set_ID(mf,i+1);
	}
	if (fedges) {
	  MSTK_free(fedges);
	  MSTK_free(fedirs);
	}
	
	processed_faces = 1;
      }
      else {
	MSTK_Report("MESH_InitFromFile",
		    "Expect face description in terms of edges",ERROR);
	return 0;
      }
    }
  }
  else {
    if (file_reptype >= F1 && file_reptype <= F4) {
      MSTK_Report("MESH_InitFromFile","Expected face information",ERROR);
      fclose(fp);
      return 0;
    }
  }
  
  
  if (processed_faces) {
    status = fscanf(fp,"%s",temp_str);
    if (status == EOF) {
      if (mesh->nr != 0)
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file after face data",FATAL);
      else {
        fclose(fp);
	return 1;
      }
    }
    else if (status == 0)
      MSTK_Report("MESH_InitFromFile",
		  "Error in reading region data",FATAL);
  }


  /* REGION DATA */

  if (strncmp(temp_str,"region",6) == 0) {
    if (strncmp(temp_str,"regions",7) != 0) 
      MSTK_Report("MESH_InitFromFile","Expected keyword \"regions\"",ERROR);

    NR = mesh->nr;
    mesh->mregion = List_New(NR);

    status = fscanf(fp,"%s",rltype_str);
    if (status == EOF)
      MSTK_Report("MESH_InitFromFile",
		  "Premature end of file while reading regions",FATAL);
    else if (status == 0)
      MSTK_Report("MESH_InitFromFile",
		  "Error in reading region data",FATAL);

    if (file_reptype == R1 || file_reptype == R2) {
      if (strncmp(rltype_str,"vertex",6) == 0) {
	rverts = NULL;
	for (i = 0; i < NR; i++) {
	  mr = MR_New(mesh);
	  
	  status = fscanf(fp,"%d",&nrv);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading region data"
			,FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading region data",FATAL);

	  if (rverts) {
	    if (nrv > max_nrv) {
	      max_nrv = nrv;
	      rverts = MSTK_realloc(rverts,max_nrv*sizeof(MVertex_ptr));
	    }
	  }
	  else {
	    max_nrv = nrv;
	    rverts = MSTK_malloc(nrv*sizeof(MVertex_ptr));
	  }
	  
	  for (j = 0; j < nrv; j++) {
	    status = fscanf(fp,"%d",&vid1);
	    if (status == EOF)
	      MSTK_Report("MESH_InitFromFile",
		      "Premature end of file while reading region data",
			  FATAL);
	    else if (status == 0)
	      MSTK_Report("MESH_InitFromFile",
			  "Error in reading region data",FATAL);

	    rverts[j] = List_Entry(mesh->mvertex,vid1-1);
#ifdef DEBUG
	    if (MV_ID(rverts[j]) != vid1)
	      MSTK_Report("MESH_InitFromFile","Mesh vertex ID mismatch",ERROR);
#endif
	  }
	  MR_Set_Vertices(mr,nrv,rverts,0,NULL);  
	  
	  status = fscanf(fp,"%d %d",&gdim,&gid);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading regions",FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading region data",FATAL);
	  
	  MR_Set_GEntDim(mr,gdim);

	  MR_Set_GEntID(mr,gid);
	  
	  MR_Set_ID(mr,i+1);
	}
	if (rverts)
	  MSTK_free(rverts);

	processed_regions = 1;
      }
      else {
	MSTK_Report("MESH_InitFromFile",
		    "Expected region description in terms of faces",FATAL);
        fclose(fp);
	return 0;
      }
    }
    else {
      if (strncmp(rltype_str,"face",4) == 0) {
	rfaces = NULL;
	rfdirs = NULL;
	for (i = 0; i < NR; i++) {
	  mr = MR_New(mesh);
	  
	  status = fscanf(fp,"%d",&nrf);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading region data",
			FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading region data",FATAL);
	  
	  if (rfaces) {
	    if (nrf > max_nrf) {
	      max_nrf = nrf;
	      rfaces = MSTK_realloc(rfaces,max_nrf*sizeof(MFace_ptr));
	      rfdirs = MSTK_realloc(rfdirs,max_nrf*sizeof(int));
	    }
	  }
	  else {
	    max_nrf = nrf;
	    rfaces = MSTK_malloc(nrf*sizeof(MFace_ptr));
	    rfdirs = MSTK_malloc(nrf*sizeof(int));
	  }
	  
	  for (j = 0; j < nrf; j++) {
	    status = fscanf(fp,"%d",&fid);
	    if (status == EOF)
	      MSTK_Report("MESH_InitFromFile",
			  "Premature end of file while reading region data",
			  FATAL);
	    else if (status == 0)
	      MSTK_Report("MESH_InitFromFile",
			  "Error in reading region data",FATAL);
	    
	    rfdirs[j] = fid > 0 ? 1 : 0;
	    rfaces[j] = List_Entry(mesh->mface,abs(fid)-1);
#ifdef DEBUG
	    if (MF_ID(rfaces[j]) != abs(fid))
	      MSTK_Report("MESH_InitFromFile","Mesh face ID mismatch",ERROR);
#endif
	  }

	  if (mesh->reptype != R1 && mesh->reptype != R2) {

	    MR_Set_Faces(mr,nrf,rfaces,rfdirs);

	  }
	  else {

	    MSTK_Report("MESH_InitFromFile","Need to convert face list to vertex list", ERROR);
	    fclose(fp);
	    return 0;

	  }
	  
	  status = fscanf(fp,"%d %d",&gdim,&gid);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			  "Premature end of file while reading faces",FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading region data",FATAL);
	  
	  MR_Set_GEntDim(mr,gdim);
	  
	  MR_Set_GEntID(mr,gid);
	  
	  MR_Set_ID(mr,i+1);
	}
	if (rfaces) {
	  MSTK_free(rfaces);
	  MSTK_free(rfdirs);
	}
	
	processed_regions = 1;
      }
      else {
	MSTK_Report("MESH_InitFromFile",
		    "Expected region description in terms of vertices",ERROR);
	fclose(fp);
	return 0;
      }
    }
  }

  if (processed_regions) {
    status = fscanf(fp,"%s",temp_str);
    if (status == EOF) {
      if (mesh->reptype == R2)
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file after reading regions",FATAL);
      else {
        fclose(fp);
	return 1;
      }
    }
  }


  /* ADJACENT REGION DATA */

  if (strncmp(temp_str,"adjregions",10) == 0) {
    if (mesh->reptype == R2) {
      for (i = 0; i < NR; i++) {
	mr = List_Entry(mesh->mregion,i);

	status = fscanf(fp,"%d",&nar);
	if (status == EOF)
	  MSTK_Report("MESH_InitFromFile",
		      "Premature end of file while reading adjacent regions",
		      FATAL);
	else if (status == 0)
	  MSTK_Report("MESH_InitFromFile",
		      "Error in reading adjacent region data",FATAL);

	for (j = 0; j < nar; j++) {
	  status = fscanf(fp,"%d",&adjrid);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading ajdacent regions",
			FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading adjacent region data",FATAL);
	  
	  if (adjrid == 0)
	    continue;
	  adjr = List_Entry(mesh->mregion,adjrid-1);
#ifdef DEBUG
	  if (MR_ID(adjr) != adjrid)
	    MSTK_Report("MESH_InitFromFile",
			"Adjacent region ID mismatch",ERROR);
#endif

	  MR_Add_AdjRegion(mr,j,adjr);
	}
      }

      processed_adjr = 1;
    }
    else {
      for (i = 0; i < NR; i++) {
	status = fscanf(fp,"%d",&nar);
	if (status == EOF)
	  MSTK_Report("MESH_InitFromFile",
		      "Premature end of file while reading adjacent regions",
		      FATAL);
	else if (status == 0)
	  MSTK_Report("MESH_InitFromFile",
		      "Error in reading adjacent region data",FATAL);
	
	for (j = 0; j < nar; j++) {
	  fscanf(fp,"%d",&rid);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading adjacent regions",
			FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading adjacent region data",FATAL);
	}
      }

      processed_adjr = 1;
    }
  }

  if (processed_adjr) {
    status = fscanf(fp,"%s",temp_str);
    if (status == EOF) {
      if (mesh->reptype == R2)
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file after reading regions",FATAL);
      else {
        fclose(fp);
	return 1;
      }
    }
  }

  /* ATTRIBUTE DATA */

  if (strncmp(temp_str,"attributes",10) == 0) {

    done = 0;
    while (!done) {
      status = fscanf(fp,"%s",attname);
      if (status == EOF) {
	done = 1;
	continue;
      }
      else if (status == 0)
	MSTK_Report("MESH_InitFromFile",
		    "Error in reading attribute data",FATAL);

      status = fscanf(fp,"%s",atttype_str);
      if (status == EOF)
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file while reading attributes",FATAL);
      else if (status == 0)
	MSTK_Report("MESH_InitFromFile",
		    "Error in reading attribute data",FATAL);

      if (strncmp(atttype_str,"INT",3) == 0)
	atttype = INT;
      else if (strncmp(atttype_str,"DOUBLE",6) == 0)
	atttype = DOUBLE;
      else if (strncmp(atttype_str,"POINTER",7) == 0)
	MSTK_Report("MESH_InitFromFile",
		    "Cannot specify POINTER attributes in file",FATAL);
      else if (strncmp(atttype_str,"VECTOR",6) == 0)
	atttype = VECTOR;
      else if (strncmp(atttype_str,"TENSOR",6) == 0)
	atttype = TENSOR;
      else {
	sprintf(temp_str,"%-s not a recognized attribute type",atttype_str);
	MSTK_Report("MESH_InitFromFile",temp_str,FATAL);
      }

      status = fscanf(fp,"%d",&ncomp);
      if (status == EOF)
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file while reading attributes",FATAL);
      else if (status == 0)
	MSTK_Report("MESH_InitFromFile",
		    "Error in reading attribute data",FATAL);

      if ((atttype == INT || atttype == DOUBLE) && (ncomp != 1)) 
	MSTK_Report("MESH_InitFromFile","Number of components should be 1 for attributes of type INT or DOUBLE",WARN);
      else if ((atttype == VECTOR || atttype == TENSOR) && (ncomp == 0))
	MSTK_Report("MESH_InitFromFile","Number of components should be non-zero for attributes of type VECTOR or TENSOR",FATAL);



      status = fscanf(fp,"%s",attent_str);
      if (status == EOF)
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file while reading attributes",FATAL);
      else if (status == 0)
	MSTK_Report("MESH_InitFromFile",
		    "Error in reading attribute data",FATAL);
      if (strncmp(attent_str,"MVERTEX",7) == 0) 
	attent = MVERTEX;
      else if (strncmp(attent_str,"MEDGE",5) == 0)
	attent = MEDGE;
      else if (strncmp(attent_str,"MFACE",5) == 0)
	attent = MFACE;
      else if (strncmp(attent_str,"MREGION",7) == 0)
	attent = MREGION;
      else if (strncmp(attent_str,"MALLTYPE",8) == 0)
	attent = MALLTYPE;
      else {
	sprintf(temp_str,"%s not a recognized entity type for attributes",
		attent_str);
	MSTK_Report("MESH_InitFromFile",temp_str,FATAL);
      }

      status = fscanf(fp,"%d",&nent);
      if (status == EOF)
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file while reading attributes",FATAL);
      else if (status == 0)
	MSTK_Report("MESH_InitFromFile",
		    "Error in reading attribute data",FATAL);

      if (nent < 1) 
	MSTK_Report("MESH_InitFromFile",
		    "Attribute applied on no entities?",FATAL);

      if (atttype == INT || atttype == DOUBLE || atttype == POINTER)
	attrib = MAttrib_New(mesh,attname,atttype,attent);
      else
	attrib = MAttrib_New(mesh,attname,atttype,attent,ncomp);
      
      for (i = 0; i < nent; i++) {
	fscanf(fp,"%d %d",&dim,&id);
	if (status == EOF)
	  MSTK_Report("MESH_InitFromFile",
		      "Premature end of file while reading attributes",FATAL);
	else if (status == 0)
	  MSTK_Report("MESH_InitFromFile",
		      "Error in reading attribute data",FATAL);
	
	if (attent != dim && attent != MALLTYPE) {
	  MSTK_Report("MESH_InitFromFile",
		      "Attribute not applicable to this type of entity",WARN);
	  if (atttype == INT)
	    fscanf(fp,"%d",&ival);
	  else
	    for (j = 0; j < ncomp; j++)
	      fscanf(fp,"%lf",&rval);
	}
	else {

	  if (atttype == INT)
	    fscanf(fp,"%d",&ival);
	  else if (atttype == DOUBLE)
	    fscanf(fp,"%lf",&rval);
	  else if (atttype == VECTOR || atttype == TENSOR) {
	    ival = ncomp;
	    rval_arr = (double *) MSTK_malloc(ncomp*sizeof(double));
	    for (j = 0; j < ncomp; j++) {
	      fscanf(fp,"%lf",&(rval_arr[j]));
	    }
	  }

	  switch (dim) {
	  case MVERTEX:
	    ent = MESH_VertexFromID(mesh,id);
	    break;
	  case MEDGE:
	    ent = MESH_EdgeFromID(mesh,id);
	    break;
	  case MFACE:
	    ent = MESH_FaceFromID(mesh,id);
	    break;
	  case MREGION:
	    ent = MESH_RegionFromID(mesh,id);
	    break;
	  default:
	    ent = NULL;
	    MSTK_Report("MESH_InitFromFile","Invalid entity type",FATAL);
	  }
	  
	  MEnt_Set_AttVal(ent,attrib,ival,rval,(void *)rval_arr);
	  
	}

      } /* for (i = 0; i < nent; i++) */

    } /* while (!done) */

  }

  fclose(fp);

  return 1;
}

int MESH_WriteToFile(Mesh_ptr mesh, const char *filename, RepType rtype) {
  FILE *fp;
  char mesg[80], attname[256];
  int i, j, k, idx;
  int gdim, gid;
  int mvid, mvid0, mvid1, mvid2, mrid2, meid, mfid, mrid;
  int nav, nar, nfe, nfv, nrf, nrv, dir=0;
  int nv, ne, nf, nr;
  int natt, ncomp, ival, nent;
  double xyz[3], rval, rdummy, *rval_arr;
  void *pval, *pdummy;
  MVertex_ptr mv, mv0, mv1, mv2;
  MEdge_ptr me;
  MFace_ptr mf;
  MRegion_ptr mr, mr2;
  List_ptr adjverts, mfedges, mfverts, mrfaces, mrverts, adjregs;
  RepType reptype;
  MAttrib_ptr attrib, vidatt, eidatt, fidatt, ridatt;
  MType attentdim;
  MAttType atttype;
  

  if (!(fp = fopen(filename,"w"))) {
    sprintf(mesg,"Cannot open file %-s for writing",filename);
    MSTK_Report("MESH_WriteToFile",mesg,ERROR);
    return 0;
  }

  if (rtype != UNKNOWN_REP) {
    reptype = rtype;
  }
  else {
    reptype = MESH_RepType(mesh);
  }

  nv = MESH_Num_Vertices(mesh);
  ne = MESH_Num_Edges(mesh);
  nf = MESH_Num_Faces(mesh);
  nr = MESH_Num_Regions(mesh);

  fprintf(fp,"MSTK %-2.1lf\n",MSTK_VER);
  fprintf(fp,"%s %d %d %d %d\n",
	  MESH_rtype_str[reptype], 
	  nv, 
	  (reptype >= R1 && reptype <= R4)?0:ne, 
	  (reptype >= R1 && reptype <= R2 && nr)?0:nf, 
	  nr);

  vidatt = MAttrib_New(mesh,"vidatt",INT,MVERTEX);
  eidatt = MAttrib_New(mesh,"eidatt",INT,MEDGE);
  fidatt = MAttrib_New(mesh,"fidatt",INT,MFACE);
  ridatt = MAttrib_New(mesh,"ridatt",INT,MREGION);

  idx = 0; i = 0;
  while ((mv = MESH_Next_Vertex(mesh,&idx)))
    MEnt_Set_AttVal(mv,vidatt,++i,0.0,NULL);

  idx = 0; i = 0;
  while ((me = MESH_Next_Edge(mesh,&idx)))
    MEnt_Set_AttVal(me,eidatt,++i,0.0,NULL);

  idx = 0; i = 0;
  while ((mf = MESH_Next_Face(mesh,&idx)))
    MEnt_Set_AttVal(mf,fidatt,++i,0.0,NULL);

  idx = 0; i = 0;
  while ((mr = MESH_Next_Region(mesh,&idx)))
    MEnt_Set_AttVal(mr,ridatt,++i,0.0,NULL);
  
  
  fprintf(fp,"vertices\n");
  idx = 0;
  while ((mv = MESH_Next_Vertex(mesh,&idx))) {

    MV_Coords(mv,xyz);

    gdim = MV_GEntDim(mv);
    gid = MV_GEntID(mv);

    fprintf(fp,"%24.16lf %24.16lf %24.16lf   %d %d\n",
	    xyz[0],xyz[1],xyz[2],gdim,gid);
    
  }

  if (reptype == R2 || reptype == R4) {
    fprintf(fp,"adjvertices\n");

    idx = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx))) {

      nav = MV_Num_AdjVertices(mv);
      fprintf(fp,"%d ",nav);
      
      adjverts = MV_AdjVertices(mv);
      for (j = 0; j < nav; j++) {
	mv2 = List_Entry(adjverts,j);
	MEnt_Get_AttVal(mv2,vidatt,&mvid2,&rval,&pval);
	fprintf(fp,"%d ",mvid2);
      }
      fprintf(fp,"\n");
      List_Delete(adjverts);
    }
  }



  if (reptype <= F4 && ne) {
    fprintf(fp,"edges\n");

    idx = 0;
    while ((me = MESH_Next_Edge(mesh,&idx))) {

      mv0 = ME_Vertex(me,0);
      MEnt_Get_AttVal(mv0,vidatt,&mvid0,&rval,&pval);
      mv1 = ME_Vertex(me,1);
      MEnt_Get_AttVal(mv1,vidatt,&mvid1,&rval,&pval);

      gdim = ME_GEntDim(me);
      gid = ME_GEntID(me);

      fprintf(fp,"%d %d \t%d %d\n",mvid0,mvid1,gdim,gid);
    }
  }



  if (reptype <= F4) {

    /* For full representations, always write out faces in terms of edges */

    fprintf(fp,"faces edge\n");
    
    idx = 0;
    while ((mf = MESH_Next_Face(mesh,&idx))) {
      
      nfe = MF_Num_Edges(mf);
      fprintf(fp,"%d ",nfe);
      
      mfedges = MF_Edges(mf,1,0);
      for (j = 0; j < nfe; j++) {
	me = List_Entry(mfedges,j);
	dir = MF_EdgeDir_i(mf,j);
	MEnt_Get_AttVal(me,eidatt,&meid,&rval,&pval);
	if (dir != 1) meid = -meid;
	fprintf(fp,"%d ",meid);
      }
      List_Delete(mfedges);
      
      gdim = MF_GEntDim(mf);
      /*
	gent = MF_GEntity(mf);
	gid = gent ? -99 : 0;
      */
      gid = MF_GEntID(mf);
      
      fprintf(fp,"\t%d %d\n",gdim,gid);
    }
  }
  else {

    /* For reduced representations, R3 and R4 always write out faces
       in terms of vertices. For reduced representations, R1 and R2
       write out faces in terms of vertices only when there are no
       regions (i.e. faces are the highest level mesh entities) */

    if ((reptype > R2) || (nr == 0)) {

      fprintf(fp,"faces vertex\n");

      idx = 0;
      while ((mf = MESH_Next_Face(mesh,&idx))) {
	
	nfv = MF_Num_Edges(mf);
	fprintf(fp,"%d ",nfv);
	
	mfverts = MF_Vertices(mf,1,0);
	for (j = 0; j < nfv; j++) {
	  mv = List_Entry(mfverts,j);
	  MEnt_Get_AttVal(mv,vidatt,&mvid,&rval,&pval);
	  fprintf(fp,"%d ",mvid);
	}
	List_Delete(mfverts);

	gdim = MF_GEntDim(mf);
	gid = MF_GEntID(mf);
	
	fprintf(fp,"\t%d %d\n",gdim,gid);
      }
    }
	
  }


  if (nr) {
    if (reptype <= F4 || reptype >= R2) {
      fprintf(fp,"regions face\n");

      idx = 0;
      while ((mr = MESH_Next_Region(mesh,&idx))) {

	nrf = MR_Num_Faces(mr);
	fprintf(fp,"%d ",nrf);

	mrfaces = MR_Faces(mr);
	for (j = 0; j < nrf; j++) {
	  mf = List_Entry(mrfaces,j);
	  dir = MR_FaceDir_i(mr,j);
	  MEnt_Get_AttVal(mf,fidatt,&mfid,&rval,&pval);
	  if (dir != 1) mfid = -mfid;
	  fprintf(fp,"%d ",mfid);
	}
	List_Delete(mrfaces);
	
	gdim = MF_GEntDim(mr);
	gid = MR_GEntID(mr);

	fprintf(fp,"\t%d %d\n",gdim,gid);
      }
    }
    else {
      fprintf(fp,"regions vertex\n");

      idx = 0;
      while ((mr = MESH_Next_Region(mesh,&idx))) {

	nrv = MR_Num_Vertices(mr);
	fprintf(fp,"%d ",nrv);

	mrverts = MR_Vertices(mr);
	for (j = 0; j < nrv; j++) {
	  mv = List_Entry(mrverts,j);
	  MEnt_Get_AttVal(mv,vidatt,&mvid,&rval,&pval);
	  fprintf(fp,"%d ",mvid);
	}
	List_Delete(mrverts);
	
	gdim = MR_GEntDim(mr);
	gid = MR_GEntID(mr);

	fprintf(fp,"\t%d %d\n",gdim,gid);
      }
    }

    if (reptype == R2 || reptype == R4) {
      fprintf(fp,"adjregions\n");
      
      idx = 0;
      while ((mr = MESH_Next_Region(mesh,&idx))) {

	nar = MR_Num_Faces(mr);
	fprintf(fp,"%d ",nar);

	adjregs = MR_AdjRegions(mr);

	for (j = 0; j < nar; j++) {
	  mr2 = List_Entry(adjregs,j);
	  if ((long) mr2 == -1) 
	    fprintf(fp,"%d ",0);
	  else {
	    MEnt_Get_AttVal(mr2,ridatt,&mrid2,&rval,&pval);
	    fprintf(fp,"%d ",mrid2);
	  }
	}
	fprintf(fp,"\n");
	List_Delete(adjregs);
      }
    }
  }


  /* Write out attributes if there are more than the 4 that we created 
    in this routine */


  if ((natt = MESH_Num_Attribs(mesh)) > 4) {

    fprintf(fp,"attributes\n");

    for (i = 0; i < natt; i++) {
      
      attrib = MESH_Attrib(mesh,i);

      /* Don't write out attribs we created for the internal use of 
	 this routine */
      if (attrib == vidatt || attrib == eidatt || attrib == fidatt || 
	  attrib == ridatt) continue;
      
      MAttrib_Get_Name(attrib,attname);

      atttype = MAttrib_Get_Type(attrib);
      if (atttype == POINTER) continue;  /* cannot write it out */

      ncomp = MAttrib_Get_NumComps(attrib);

      attentdim = MAttrib_Get_EntDim(attrib);


      /* First count how many entities actually have the attribute assigned */

      nent = 0;
      switch(attentdim) {
      case MVERTEX:
	idx = 0;
	while ((mv = MESH_Next_Vertex(mesh,&idx)))
	  if (MEnt_Get_AttVal(mv,attrib,&ival,&rval,&pval)) nent++;
	break;
      case MEDGE:
	idx = 0;
	while ((me = MESH_Next_Edge(mesh,&idx)))
	  if (MEnt_Get_AttVal(me,attrib,&ival,&rval,&pval)) nent++;
	break;
      case MFACE:
	idx = 0;
	while ((mf = MESH_Next_Face(mesh,&idx)))
	  if (MEnt_Get_AttVal(mf,attrib,&ival,&rval,&pval)) nent++;	    
	break;
      case MREGION: 
	idx = 0;
	while ((mr = MESH_Next_Region(mesh,&idx)))
	  if (MEnt_Get_AttVal(mr,attrib,&ival,&rval,&pval)) nent++;
	break;
      case MALLTYPE:
	idx = 0;
	while ((mv = MESH_Next_Vertex(mesh,&idx)))
	  if (MEnt_Get_AttVal(mv,attrib,&ival,&rval,&pval)) nent++;
	idx = 0;
	while ((me = MESH_Next_Edge(mesh,&idx)))
	  if (MEnt_Get_AttVal(me,attrib,&ival,&rval,&pval)) nent++;
	idx = 0;
	while ((mf = MESH_Next_Face(mesh,&idx)))
	  if (MEnt_Get_AttVal(mf,attrib,&ival,&rval,&pval)) nent++;	    
	idx = 0;
	while ((mr = MESH_Next_Region(mesh,&idx)))
	  if (MEnt_Get_AttVal(mr,attrib,&ival,&rval,&pval)) nent++;
	break;	
      default:
	break;
      } /* switch (attentdim) */


      /* No point in writing out attribute if no entity uses it! Or is there? */

      if (!nent) continue;



      fprintf(fp,"%-s\n",attname);

      switch(atttype) {
      case INT:
	fprintf(fp,"INT\n");
	break;
      case DOUBLE:
	fprintf(fp,"DOUBLE\n");
	break;
      case VECTOR:
	fprintf(fp,"VECTOR\n");
	break;
      case TENSOR:
	fprintf(fp,"TENSOR\n");
	break;
      default:
	MSTK_Report("MESH_WriteToFile",
		    "Unrecognizable or unprintable attribute type\n",WARN);
	continue;	
      }

      fprintf(fp,"%-d\n",ncomp);

      switch(attentdim) {
      case MVERTEX:
	fprintf(fp,"MVERTEX\n");
	break;
      case MEDGE:
	fprintf(fp,"MEDGE\n");
	break;
      case MFACE:
	fprintf(fp,"MFACE\n");
	break;
      case MREGION:
	fprintf(fp,"MREGION\n");
	break;
      case MALLTYPE:
	fprintf(fp,"MALLTYPE\n");
	break;
      default:
	MSTK_Report("Mesh_WriteToFile","Unrecognized entity type",WARN);
	break;
      }

      fprintf(fp,"%-d\n",nent);


      switch(attentdim) {
      case MVERTEX:
	idx = 0;
	while ((mv = MESH_Next_Vertex(mesh,&idx))) {
	  if (MEnt_Get_AttVal(mv,attrib,&ival,&rval,&pval)) {
	    MEnt_Get_AttVal(mv,vidatt,&mvid,&rdummy,&pdummy);
	    fprintf(fp,"0 %-d ",mvid);
	    switch (atttype) {
	    case INT:
	      fprintf(fp," %-d",ival);
	      break;
	    case DOUBLE: 
	      fprintf(fp," %-lf ",rval);
	      break;
	    case VECTOR: case TENSOR:
	      rval_arr = (double *) pval;
	      for (k = 0; k < ncomp; k++)
		fprintf(fp," %-lf ",rval_arr[k]);
	      break;
	    default:
	      break;
	    }
	    fprintf(fp,"\n");
	  }
	}
	break;
      case MEDGE:
	idx = 0;
	while ((me = MESH_Next_Edge(mesh,&idx))) {
	  if (MEnt_Get_AttVal(me,attrib,&ival,&rval,&pval)) {
	    MEnt_Get_AttVal(me,eidatt,&meid,&rdummy,&pdummy);
	    fprintf(fp,"1 %-d ",meid);
	    switch (atttype) {
	    case INT:
	      fprintf(fp," %-d",ival);
	      break;
	    case DOUBLE: 
	      fprintf(fp," %-lf ",rval);
	      break;
	    case VECTOR: case TENSOR:
	      rval_arr = (double *) pval;
	      for (k = 0; k < ncomp; k++)
		fprintf(fp," %-lf ",rval_arr[k]);
	      break;
	    default:
	      break;
	    }
	    fprintf(fp,"\n");
	  }
	}
	break;
      case MFACE:
	idx = 0;
	while ((mf = MESH_Next_Face(mesh,&idx))) {
	  if (MEnt_Get_AttVal(mf,attrib,&ival,&rval,&pval)) {
	    MEnt_Get_AttVal(mf,fidatt,&mfid,&rdummy,&pdummy);
	    fprintf(fp,"2 %-d ",mfid);
	    switch (atttype) {
	    case INT:
	      fprintf(fp," %-d",ival);
	      break;
	    case DOUBLE: 
	      fprintf(fp," %-lf ",rval);
	      break;
	    case VECTOR: case TENSOR:
	      rval_arr = (double *) pval;
	      for (k = 0; k < ncomp; k++)
		fprintf(fp," %-lf ",rval_arr[k]);
	      break;
	    default:
	      break;
	    }
	    fprintf(fp,"\n");
	  }
	}
	break;

      case MREGION: 
	idx = 0;
	while ((mr = MESH_Next_Region(mesh,&idx))) {
	  if (MEnt_Get_AttVal(mr,attrib,&ival,&rval,&pval)) {
	    MEnt_Get_AttVal(mr,ridatt,&mrid,&rdummy,&pdummy);
	    fprintf(fp,"3 %-d ",mrid);
	    switch (atttype) {
	    case INT:
	      fprintf(fp," %-d",ival);
	      break;
	    case DOUBLE: 
	      fprintf(fp," %-lf ",rval);
	      break;
	    case VECTOR: case TENSOR:
	      rval_arr = (double *) pval;
	      for (k = 0; k < ncomp; k++)
		fprintf(fp," %-lf ",rval_arr[k]);
	      break;
	    default:
	      break;
	    }
	    fprintf(fp,"\n");
	  }
	}
	break;

      case MALLTYPE:
	idx = 0;
	while ((mv = MESH_Next_Vertex(mesh,&idx))) {
	  if (MEnt_Get_AttVal(mv,attrib,&ival,&rval,&pval)) {
	    MEnt_Get_AttVal(mv,vidatt,&mvid,&rdummy,&pdummy);
	    fprintf(fp,"0 %-d ",mvid);
	    switch (atttype) {
	    case INT:
	      fprintf(fp," %-d",ival);
	      break;
	    case DOUBLE: 
	      fprintf(fp," %-lf ",rval);
	      break;
	    case VECTOR: case TENSOR:
	      rval_arr = (double *) pval;
	      for (k = 0; k < ncomp; k++)
		fprintf(fp," %-lf ",rval_arr[k]);
	      break;
	    default:
	      break;
	    }
	    fprintf(fp,"\n");
	  }
	}
	idx = 0;
	while ((me = MESH_Next_Edge(mesh,&idx))) {
	  if (MEnt_Get_AttVal(me,attrib,&ival,&rval,&pval)) {
	    MEnt_Get_AttVal(me,eidatt,&meid,&rdummy,&pdummy);
	    fprintf(fp,"1 %-d ",meid);
	    switch (atttype) {
	    case INT:
	      fprintf(fp," %-d",ival);
	      break;
	    case DOUBLE: 
	      fprintf(fp," %-lf ",rval);
	      break;
	    case VECTOR: case TENSOR:
	      rval_arr = (double *) pval;
	      for (k = 0; k < ncomp; k++)
		fprintf(fp," %-lf ",rval_arr[k]);
	      break;
	    default:
	      break;
	    }
	    fprintf(fp,"\n");
	  }
	}
	idx = 0;
	while ((mf = MESH_Next_Face(mesh,&idx))) {
	  if (MEnt_Get_AttVal(mf,attrib,&ival,&rval,&pval)) {
	    MEnt_Get_AttVal(mf,fidatt,&mfid,&rdummy,&pdummy);
	    fprintf(fp,"2 %-d ",mfid);
	    switch (atttype) {
	    case INT:
	      fprintf(fp," %-d",ival);
	      break;
	    case DOUBLE: 
	      fprintf(fp," %-lf ",rval);
	      break;
	    case VECTOR: case TENSOR:
	      rval_arr = (double *) pval;
	      for (k = 0; k < ncomp; k++)
		fprintf(fp," %-lf ",rval_arr[k]);
	      break;
	    default:
	      break;
	    }
	    fprintf(fp,"\n");
	  }
	}
	idx = 0;
	while ((mr = MESH_Next_Region(mesh,&idx))) {
	  if (MEnt_Get_AttVal(mr,attrib,&ival,&rval,&pval)) {
	    MEnt_Get_AttVal(mr,ridatt,&mrid,&rdummy,&pdummy);
	    fprintf(fp,"3 %-d ",mrid);
	    switch (atttype) {
	    case INT:
	      fprintf(fp," %-d",ival);
	      break;
	    case DOUBLE: 
	      fprintf(fp," %-lf ",rval);
	      break;
	    case VECTOR: case TENSOR:
	      rval_arr = (double *) pval;
	      for (k = 0; k < ncomp; k++)
		fprintf(fp," %-lf ",rval_arr[k]);
	      break;
	    default:
	      break;
	    }
	    fprintf(fp,"\n");
	  }
	}
	break;	
      default:
	break;
      } /* switch (attentdim) */

    } /* for (i = 0; i < natt) */
    
  } /* if (Mesh_Num_Attribs(mesh)) */
  

  idx = 0; i = 0;
  while ((mv = MESH_Next_Vertex(mesh,&idx)))
    MEnt_Rem_AttVal(mv,vidatt);

  idx = 0; i = 0;
  while ((me = MESH_Next_Edge(mesh,&idx)))
    MEnt_Rem_AttVal(me,eidatt);

  idx = 0; i = 0;
  while ((mf = MESH_Next_Face(mesh,&idx)))
    MEnt_Rem_AttVal(mf,fidatt);

  idx = 0; i = 0;
  while ((mr = MESH_Next_Region(mesh,&idx)))
    MEnt_Rem_AttVal(mr,ridatt);
  
  MAttrib_Delete(vidatt);
  MAttrib_Delete(eidatt);
  MAttrib_Delete(fidatt);
  MAttrib_Delete(ridatt);




  fclose(fp);

  return 1;
}


  /* Enforce continuous numbering for mesh entities */

void MESH_Renumber(Mesh_ptr mesh) {
  MVertex_ptr mv;
  MEdge_ptr me;
  MFace_ptr mf;
  MRegion_ptr mr;
  int idx, n;

  idx = 0; n = 0;
  while ((mv = MESH_Next_Vertex(mesh,&idx)))
    MV_Set_ID(mv,++n);

  idx = 0; n = 0;
  while ((me = MESH_Next_Edge(mesh,&idx)))
    ME_Set_ID(me,++n);

  idx = 0; n = 0;
  while ((mf = MESH_Next_Face(mesh,&idx)))
    MF_Set_ID(mf,++n);

  idx = 0; n = 0;
  while ((mr = MESH_Next_Region(mesh,&idx)))
    MR_Set_ID(mr,++n);

  return;
}

  /*

Mesh_ptr MESH_Copy(Mesh_ptr oldmesh) {
  MVertex_ptr mv;
  MEdge_ptr me;
  MFace_ptr mf;
  MRegion_ptr mr;
  int idx, n;

  newmesh = MESH_New(oldmesh->repType);

  newmesh->nv = oldmesh->nv;

  return;
}

  */


  /* Time for this function to go */


  int MESH_Init_ParAtts(Mesh_ptr mesh) {
    
    plnumatt = MAttrib_New(mesh,"paratts",POINTER,MALLTYPE);

    if (plnumatt)
      return 1;
    else
      return 0;
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
