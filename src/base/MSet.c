#define _H_MSet_Private

#include <string.h>
#include <stdarg.h>
#include "MSet.h"
#include "Mesh.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif


  MSet_ptr MSet_New(Mesh_ptr mesh, const char *set_name, MType entdim) {
    MSet_ptr set;
    int i, nset;

#ifdef DEBUG
    nset = MESH_Num_MSets(mesh);
    for (i = 0; i < nset; i++) {
      set = MESH_MSet(mesh,i);
      if (strcmp(set_name,set->name) == 0) {
	MSTK_Report("MSet_New","Set with given name exists",WARN);
	return set;
      }
    }
#endif

    set = (MSet_ptr) MSTK_malloc(sizeof(MSet));
    set->name = (char *) MSTK_malloc((strlen(set_name)+1)*sizeof(char));
    strcpy(set->name,set_name);
    set->entdim = entdim; /* what type of entity can this be attached to? */

    set->entlist = List_New(0);

    set->mesh = mesh;
    MESH_Add_MSet(mesh,set);

    return set;
  }


  char *MSet_Name(MSet_ptr set, char *set_name) {
    strcpy(set_name,set->name);
    return set_name;
  }

  MType MSet_EntDim(MSet_ptr set) {
    return set->entdim;
  }

  Mesh_ptr MSet_Mesh(MSet_ptr set) {
    return set->mesh;
  }

  void MSet_Delete(MSet_ptr set) {
    MESH_Rem_MSet(set->mesh,set);
    MSTK_free(set->name);
    List_Delete(set->entlist);
    MSTK_free(set);
  }

  MSet_ptr    MSet_Add(MSet_ptr set, void *entry) {
    return List_Add(set->entlist,entry);
  }

  MSet_ptr    MSet_ChknAdd(MSet_ptr set, void *entry) {
    return List_ChknAdd(set->entlist,entry);
  }

  MSet_ptr    MSet_Insert(MSet_ptr set, void *nuentry, void *b4entry) {
    return List_Insert(set->entlist,nuentry,b4entry);
  }


  MSet_ptr    MSet_Inserti(MSet_ptr set, void *nuentry, int i) {
    return List_Inserti(set->entlist,nuentry,i);
  }

  int         MSet_Rem(MSet_ptr set, void *entry) {
    return List_Rem(set->entlist,entry);
  }

  int         MSet_Remi(MSet_ptr set, int i) {
    return List_Remi(set->entlist,i);
  }

  int         MSet_Replace(MSet_ptr set, void *entry, void *nuentry) {
    return List_Replace(set->entlist,entry,nuentry);
  }

  int         MSet_Replacei(MSet_ptr set, int i, void *nuentry) {
    return List_Replacei(set->entlist,i,nuentry);
  }

  int         MSet_Contains(MSet_ptr set, void *entry) {
    return List_Contains(set->entlist, entry);
  }

  int         MSet_Locate(MSet_ptr set, void *entry) {
    return List_Contains(set->entlist, entry);
  }

  void       *MSet_Entry(MSet_ptr set, int i) {
    return List_Entry(set->entlist, i);
  }

  void       *MSet_Next_Entry(MSet_ptr set, int *i) {
    return List_Next_Entry(set->entlist, i);
  }

  int         MSet_Num_Entries(MSet_ptr set) {
    return List_Num_Entries(set->entlist);
  }


  MSet_ptr    MSet_Copy(MSet_ptr src) {
    MSet_ptr newset;
    char name[256];

    newset = MSet_New(MSet_Mesh(src),MSet_Name(src,name),MSet_EntDim(src));

    newset->entlist = List_Copy(src->entlist);

    return newset;
  }


  MSet_ptr    MSets_Cat(MSet_ptr dest, MSet_ptr src) {    
    if (MSet_Mesh(dest) != MSet_Mesh(src) ||
	MSet_EntDim(dest) != MSet_EntDim(dest)) {
      MSTK_Report("MSets_Cat","MSets must belong to same mesh and have same type of entities",FATAL);
    }

    dest->entlist = List_Cat(dest->entlist, src->entlist);
    
    return dest;
  }

  /* The entities in the two sets must be valid entities (not just
     some pointers) because we are going to be marking them. Also they
     should belong to the same mesh */

  MSet_ptr    MSets_Union(MSet_ptr s1, MSet_ptr s2) {
    MSet_ptr newset;
    char *newname, s1name[256], s2name[256];
    int idx, mkid;
    MEntity_ptr ent;

    if (MSet_Mesh(s1) != MSet_Mesh(s2)) {

#ifdef DEBUG
      MSTK_Report("MSets_Union","MSets must belong to same mesh",WARN); 
#endif

      return NULL;
    }

    MSet_Name(s1,s1name);
    MSet_Name(s2,s2name);

    newname = (char *) MSTK_malloc((strlen(s1name)+strlen(s2name)+7)*sizeof(char));
    strcpy(newname,"in_");
    strcat(newname,s1name);
    strcat(newname,"_or_");
    strcat(newname,s2name);

    if (MSet_EntDim(s1) != MSet_EntDim(s2)) {

      /* If entities in the two sets are of different types, the only way 
	 we can accommodate them in the new set is to make it accommodate
	 all entity types */

      newset = MSet_New(MSet_Mesh(s1),newname,MALLTYPE);
    }
    else
      newset = MSet_New(MSet_Mesh(s1),newname,MSet_EntDim(s1));

    
    newset->entlist = List_Copy(s1->entlist);
    
    mkid = MSTK_GetMarker();

    MSet_Mark(newset,mkid);

    idx = 0;
    while ((ent = (MEntity_ptr) MSet_Next_Entry(s2,&idx))) {
      if (!MEnt_IsMarked(ent,mkid)) {
	MSet_Add(newset,ent);
	MEnt_Mark(ent,mkid);
      }
    }

    MSet_Unmark(newset,mkid);
    MSTK_FreeMarker(mkid);

    return newset;
  }


  /* Returns NULL pointer if the intersection is NULL */

  MSet_ptr    MSets_Intersect(MSet_ptr s1, MSet_ptr s2) {
    MSet_ptr newset;
    char *newname, s1name[256], s2name[256];
    int mkid1, mkid2, idx;
    MEntity_ptr ent;
    MSet_ptr sa, sb;
    
    if (MSet_Mesh(s1) != MSet_Mesh(s2)) {

#ifdef DEBUG
      MSTK_Report("MSets_Intersect","MSets must belong to same mesh",WARN); 
#endif

      return NULL;
    }

    MSet_Name(s1,s1name);
    MSet_Name(s2,s2name);

    newname = (char *) MSTK_malloc((strlen(s1name)+strlen(s2name)+8)*sizeof(char));
    strcpy(newname,"in_");
    strcat(newname,s1name);
    strcat(newname,"_and_");
    strcat(newname,s2name);

    if (MSet_EntDim(s1) != MSet_EntDim(s2)) {

      /* If entities in the two sets are of different types, the only way 
	 we can accommodate them in the new set is to make it accommodate
	 all entity types */

      newset = MSet_New(MSet_Mesh(s1),newname,MALLTYPE);
    }
    else
      newset = MSet_New(MSet_Mesh(s1),newname,MSet_EntDim(s1));


    mkid1 = MSTK_GetMarker();
    MSet_Mark(s1,mkid1);

    mkid2 = MSTK_GetMarker();
    MSet_Mark(s2,mkid2);

    if (MSet_Num_Entries(s1) < MSet_Num_Entries(s2)) {
      sa = s1;
      sb = s2;
    }
    else {
      sa = s2;
      sb = s1;
    }

    idx = 0;
    while ((ent = MSet_Next_Entry(sa,&idx))) {
      if (MEnt_IsMarked(ent,mkid1) && MEnt_IsMarked(ent,mkid2))
	MSet_Add(newset,ent);
    }

    MSet_Unmark(s1,mkid1);
    MSTK_FreeMarker(mkid1);
    MSet_Unmark(s2,mkid2);
    MSTK_FreeMarker(mkid2);
    

    if (!MSet_Num_Entries(newset)) {
      MSet_Delete(newset);
      newset = NULL;
    }

    return newset;
  }


  /* Will return NULL if resulting set is of zero length */

  MSet_ptr    MSets_Subtract(MSet_ptr s1, MSet_ptr s2) {

    MSet_ptr newset;
    char *newname, s1name[256], s2name[256];
    int mkid1, mkid2, idx;
    MEntity_ptr ent;
    
    if (MSet_Mesh(s1) != MSet_Mesh(s2)) {

#ifdef DEBUG
      MSTK_Report("MSets_Intersect","MSets must belong to same mesh",WARN); 
#endif

      return NULL;
    }

    MSet_Name(s1,s1name);
    MSet_Name(s2,s2name);

    newname = (char *) MSTK_malloc((strlen(s1name)+strlen(s2name)+10)*sizeof(char));
    strcpy(newname,"in_");
    strcat(newname,s1name);
    strcat(newname,"_minus_");
    strcat(newname,s2name);


    newset = MSet_New(MSet_Mesh(s1),newname,MSet_EntDim(s1));


    mkid1 = MSTK_GetMarker();
    MSet_Mark(s1,mkid1);

    mkid2 = MSTK_GetMarker();
    MSet_Mark(s2,mkid2);

    idx = 0;
    while ((ent = MSet_Next_Entry(s1,&idx))) {
      if (MEnt_IsMarked(ent,mkid1) && !MEnt_IsMarked(ent,mkid2))
	MSet_Add(newset,ent);
    }

    MSet_Unmark(s1,mkid1);
    MSTK_FreeMarker(mkid1);
    MSet_Unmark(s2,mkid2);
    MSTK_FreeMarker(mkid2);
    

    if (!MSet_Num_Entries(newset)) {
      MSet_Delete(newset);
      newset = NULL;
    }

    return newset;

  }


  void MSet_Mark(MSet_ptr list, int markerID) {
    MEntity_ptr ent;
    int idx = 0;
    
    while ((ent = (MEntity_ptr) List_Next_Entry(list,&idx)))
      MEnt_Mark(ent,markerID);
  }

  void MSet_Unmark(MSet_ptr list, int markerID) {
    MEntity_ptr ent;
    int idx = 0;
    
    while ((ent = (MEntity_ptr) List_Next_Entry(list,&idx)))
      MEnt_Unmark(ent,markerID);
  }


#ifdef __cplusplus
}
#endif
