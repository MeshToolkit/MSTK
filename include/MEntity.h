#ifndef _H_MEntity
#define _H_MEntity

#include "MSTK_defines.h"
#include "MSTK_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _H_MEntity_Private

  typedef struct MEntity_Data {

    Mesh_ptr mesh;
    List_ptr AttInsList;

    /* The first bit (from the right) in ent->dim_id contains flag
       indicating if the entity is alive or deleted. Bit 2 indicates
       if the entity is temporary/volatile or permanent (edges in all
       reduced representations are temporary).Bits 3,4,5 combined
       contain the dimension of the entity. If this number is greater
       than or equal to 4, it means the dimension of the entity is
       undefined. The rest of the bits encode the entity ID. Since
       this is an unsigned int (32 bits), the max ID can be 2^(32-5)-1
       = 134217727 (approx. 134.2 million entities of each type) */

    unsigned int deleted:1;   /* flag indicating if entity is deleted */
    unsigned int temporary:1; /* flag indicating if entity is temporary */
    unsigned int dim:3;       /* dimension of entity */
    unsigned int id:27;       /* ID of entity */

    unsigned int rtype:3;     /* Representation type */
    unsigned int gdim:3;      /* Geometric entity dimension */
    unsigned int gid:26;      /* Geometric entity ID */

    unsigned int marker;

#ifdef MSTK_HAVE_MPI
    unsigned int ptype:2;    /* parallel entity type */
    unsigned int parbdry:1;   /* is or is not on partition boundary */
    unsigned int masterparid:29; /* This is 2^29-1 partitions - enough, no? */
    unsigned int globalid;  
#endif

  } MEntity_Data, *MEntity_Data_ptr;

  typedef struct MEntity {
    MEntity_Data entdat;
  } MEntity, *MEntity_ptr;

#else
  typedef void *MEntity_ptr;
#endif

  PType MEnt_PType(MEntity_ptr ent);

  /* Is the entity on the partition boundary? */
  int MEnt_OnParBoundary(MEntity_ptr ent);

  int   MEnt_MasterParID(MEntity_ptr ent);

  int   MEnt_GlobalID(MEntity_ptr ent);

#ifdef MSTK_HAVE_MPI
  void  MEnt_Set_PType(MEntity_ptr ent, PType ptype);

  /* Mark/Unmark the entity as being on the partition boundary */
  void MEnt_Flag_OnParBoundary(MEntity_ptr ent);
  void MEnt_Unflag_OnParBoundary(MEntity_ptr ent);
  void  MEnt_Set_MasterParID(MEntity_ptr ent, int masterparid);
  void  MEnt_Set_GlobalID(MEntity_ptr ent, int globalid);

  /* local id of owning entity on master processor */
  /*
  int   MEnt_MasterLocalID(MEntity_ptr ent);
  void   MEnt_Set_MasterLocalID(MEntity_ptr ent, int localid);
  */
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
  void  MEnt_Clear_AttVal(MEntity_ptr ent, MAttrib_ptr attrib);
  int   MEnt_Get_AttVal(MEntity_ptr ent, MAttrib_ptr attrib, int *ival, 
			double *lval, void **pval);

#ifdef __cplusplus
}
#endif

#endif
