/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#ifndef _H_MATTRIB
#define _H_MATTRIB


#include "MSTK_defines.h"
#include "MSTK_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _H_MAttrib_Private
  typedef struct MAttrib { /* Attribute Definition */
    char *name;
    MAttType type;
    int ncomp;
    MType entdim;
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


  MAttrib_ptr MAttrib_New(Mesh_ptr mesh, const char *att_name, MAttType att_type, MType entdim, ...);
  char       *MAttrib_Get_Name(MAttrib_ptr attrib, char *att_name);
  MAttType    MAttrib_Get_Type(MAttrib_ptr attrib);
  MType       MAttrib_Get_EntDim(MAttrib_ptr attrib);
  int         MAttrib_Get_NumComps(MAttrib_ptr attrib);
  void        MAttrib_Delete(MAttrib_ptr attrib);
  void        MAttrib_Clear(MAttrib_ptr attrib);




#ifdef __cplusplus
}
#endif

#endif
