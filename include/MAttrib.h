#ifndef _H_MATTRIB
#define _H_MATTRIB


#ifdef __cplusplus
extern "C" {
#endif

#include "MSTK_defines.h"
#include "MSTK_types.h"
#include "List.h"


#ifdef _H_MAttrib_Private
  typedef struct MAttrib { /* Attribute Definition */
    char *name;
    MAttType type;
    Mesh_ptr mesh;
  } MAttrib, *MAttrib_ptr;

  typedef struct MAttIns { /* Instance of an attribute */
    MAttrib_ptr attrib;
    union {
      int ival;
      double lval;
      void *pval;
    } att_val;
  } MAttIns, *MAttIns_ptr;
#endif


  MAttIns_ptr MAttIns_New(MAttrib_ptr attrib);
  MAttrib_ptr MAttIns_Attrib(MAttIns_ptr att);
  void        MAttIns_Set_Value(MAttIns_ptr att, int ival, double rval, void *pval);
  int         MAttIns_Get_Value(MAttIns_ptr att, int *ival, double *rval, void **pval);
  void        MAttIns_Delete(MAttIns_ptr att);


  MAttrib_ptr MAttrib_New(Mesh_ptr mesh, char *att_name, MAttType att_type);
  char       *MAttrib_Get_Name(MAttrib_ptr attrib, char *att_name);
  MAttType    MAttrib_Get_Type(MAttrib_ptr attrib);
  void        MAttrib_Delete(MAttrib_ptr attrib);




#ifdef __cplusplus
}
#endif

#endif
