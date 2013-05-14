#ifndef _H_MSet
#define _H_MSet


#include "MSTK_defines.h"
#include "MSTK_types.h"


#ifdef __cplusplus
extern "C" {
#endif

#ifdef _H_MSet_Private
  typedef struct MSet { /* Set Definition */
    char *name;
    Mesh_ptr mesh;
    MType entdim;
    List_ptr entlist;
  } MSet, *MSet_ptr;
#else
  typedef void *MSet_ptr;
#endif



  MSet_ptr    MSet_New(Mesh_ptr mesh, const char *set_name, MType entdim);
  char       *MSet_Name(MSet_ptr set, char *set_name);
  Mesh_ptr    MSet_Mesh(MSet_ptr);
  MType       MSet_EntDim(MSet_ptr set);
  void        MSet_Delete(MSet_ptr set);

  void        MSet_Rename(MSet_ptr set, char *newname);
  MSet_ptr    MSet_Add(MSet_ptr set, void *entry);
  MSet_ptr    MSet_ChknAdd(MSet_ptr set, void *entry);
  MSet_ptr    MSet_Insert(MSet_ptr set, void *nuentry, void *b4entry);
  MSet_ptr    MSet_Inserti(MSet_ptr set, void *nuentry, int i);
  int         MSet_Rem(MSet_ptr set, void *entry);
  int         MSet_Remi(MSet_ptr set, int i);
  int         MSet_Replace(MSet_ptr set, void *entry, void *nuentry);
  int         MSet_Replacei(MSet_ptr set, int i, void *nuentry);
  int         MSet_Contains(MSet_ptr set, void *entry);
  int         MSet_Locate(MSet_ptr set, void *entry);
  void       *MSet_Entry(MSet_ptr set, int i);
  void       *MSet_Next_Entry(MSet_ptr set, int *i);
  int         MSet_Num_Entries(MSet_ptr set);
  MSet_ptr    MSet_Cat(MSet_ptr dest, MSet_ptr src);
  MSet_ptr    MSet_Copy(MSet_ptr oldset);

  MSet_ptr    MSets_Union(MSet_ptr s1, MSet_ptr s2);
  MSet_ptr    MSets_Intersect(MSet_ptr s1, MSet_ptr s2);
  MSet_ptr    MSets_Subtract(MSet_ptr s1, MSet_ptr s2);

#ifdef DEBUG
  void        MSet_Print(MSet_ptr set);
#endif



#ifdef __cplusplus
}
#endif

#endif
