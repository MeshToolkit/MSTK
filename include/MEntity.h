#ifndef _H_MEntity
#define _H_MEntity

#ifdef __cplusplus
extern "C" {
#endif

#include "MSTK_defines.h"
#include "MSTK_types.h"

#ifdef _H_MEntity_Private

  typedef struct MEntity_Data {
    Mesh_ptr mesh;
    unsigned int dim_id;
    unsigned int rtype_gdim_gid;
    unsigned int marker;
    List_ptr AttInsList;
  } MEntity_Data, *MEntity_Data_ptr;

  typedef struct MEntity {
    MEntity_Data entdat;
  } MEntity, *MEntity_ptr;

#else
  typedef void *MEntity_ptr;
#endif

  void MEnt_Init_CmnData(MEntity_ptr ent);
  void MEnt_Free_CmnData(MEntity_ptr ent);

  void  MEnt_Set_Dim(MEntity_ptr ent, MType dim);
  MType MEnt_Dim(MEntity_ptr mentity);
  MType MEnt_OrigDim(MEntity_ptr mentity);

  void MEnt_Set_Volatile(MEntity_ptr mentity);
  int  MEnt_IsVolatile(MEntity_ptr mentity);

  void MEnt_Set_ID(MEntity_ptr mentity, int id);
  int  MEnt_ID(MEntity_ptr mentity);


  void MEnt_Set_GEntity(MEntity_ptr mentity, GEntity_ptr gent);
  void MEnt_Set_GEntID(MEntity_ptr mentity, int gid);
  void MEnt_Set_GEntDim(MEntity_ptr mentity, int gdim);

  GEntity_ptr MEnt_GEntity(MEntity_ptr mentity);
  int         MEnt_GEntID(MEntity_ptr mentity);
  int         MEnt_GEntDim(MEntity_ptr mentity);


  void     MEnt_Set_Mesh(MEntity_ptr mentity, Mesh_ptr mesh);
  Mesh_ptr MEnt_Mesh(MEntity_ptr mentity);

  void    MEnt_Set_RepType_Data(MEntity_ptr mentity, RepType rtype);
  void    MEnt_Set_RepType(MEntity_ptr mentity, RepType rtype);
  RepType MEnt_RepType(MEntity_ptr mentity);


  void MEnt_Set_DelFlag(MEntity_ptr ent);
  void MEnt_Rem_DelFlag(MEntity_ptr ent);
  void MEnt_Delete(MEntity_ptr mentity, int keep);


  /* Entity Marking */

  void MEnt_Mark(MEntity_ptr ent, int mkr);
  int  MEnt_IsMarked(MEntity_ptr ent, int mkr);
  void MEnt_UnMark(MEntity_ptr ent, int mkr);
  void List_Mark(List_ptr list, int mkr);
  void List_Unmark(List_ptr list, int mkr);


  /* Attributes */

  void  MEnt_Set_AttVal(MEntity_ptr ent, MAttrib_ptr attrib, int ival, 
			double lval, void *pval);
  void  MEnt_Rem_AttVal(MEntity_ptr ent, MAttrib_ptr attrib);
  int   MEnt_Get_AttVal(MEntity_ptr ent, MAttrib_ptr attrib, int *ival, 
			double *lval, void **pval);

  /* Adhoc routines for parallel output of MSTK files */

  int MEnt_NumProcs(MEntity_ptr ent);
  void MEnt_Set_ProcIDs(MEntity_ptr ent, int np, int *procids);
  int MEnt_ProcIDs(MEntity_ptr ent, int *np, int *procids);
  void MEnt_Set_LocalID(MEntity_ptr ent, int procid, int lnum);
  int MEnt_LocalID(MEntity_ptr ent, int procid);

#ifdef __cplusplus
}
#endif

#endif
