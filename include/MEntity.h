#ifndef _H_MEntity
#define _H_MEntity

#ifdef __cplusplus
extern "C" {
#endif

#include "MSTK_defines.h"
#include "MSTK_types.h"
#include "List.h"

#ifdef _H_MEntity_Private
  typedef struct MEntity {
    int id;
    int marker;
    Mesh_ptr mesh;
    char dim;
    char gdim;
    int gid;
    GEntity_ptr gent;
  } MEntity, *MEntity_ptr;
#else
  typedef void *MEntity_ptr;
#endif

  void MEnt_Set_GEntity(MEntity_ptr mentity, GEntity_ptr gent);
  void MEnt_Set_GEntDim(MEntity_ptr mentity, int gdim);
  void MEnt_Set_ID(MEntity_ptr mentity, int id);

  int MEnt_ID(MEntity_ptr mentity);
  int MEnt_Dim(MEntity_ptr mentity);
  int MEnt_OrigDim(MEntity_ptr mentity);
  Mesh_ptr MEnt_Mesh(MEntity_ptr mentity);
  int MEnt_GEntDim(MEntity_ptr mentity);
  GEntity_ptr MEnt_GEntity(MEntity_ptr mentity);

  void MEnt_Delete(MEntity_ptr mentity, int keep);

  void MEnt_Mark(MEntity_ptr ent, int mkr);
  int  MEnt_IsMarked(MEntity_ptr ent, int mkr);
  void MEnt_UnMark(MEntity_ptr ent, int mkr);
  void List_Mark(List_ptr list, int mkr);
  void List_Unmark(List_ptr list, int mkr);

  
#ifdef __cplusplus
}
#endif

#endif
