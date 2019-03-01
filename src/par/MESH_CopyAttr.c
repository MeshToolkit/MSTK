/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "MSTK.h"

#ifdef __cplusplus
extern "C" {
#endif

 /* 
    *** WARNING WARNING:***

    THIS IS AN INTERNAL MSTK ROUTINE AND MUST BE USED UNDER SPECIFIC
    CIRCUMSTANCES AND IS APPROPRIATE TO BE CALLED BY
    'MSTK_Mesh_Distribute' or 'MSTK_Mesh_Partition' ONLY. IT RELIES ON
    THE ASSUMPTION THAT THERE ARE LINKS FROM THE ENTITIES OF THE
    GLOBAL MESH TO CORRESPONDING ENTITIES OF THE SUBMESHES THROUGH AN
    ATTRIBUTE CALLED 'g2latt'
    

    this function copy attribute information from mesh to submesh
    assume mesh and submesh have the same attributes list. 

    attr_name: the attribute name

 */

int MESH_CopyAttr(Mesh_ptr mesh, int num, Mesh_ptr *submeshes, const char *attr_name) {
  int i, ncomp, ival, idx;
  double rval, *pval_arr=NULL;
  void *pval;
  MType mtype;
  MAttrib_ptr *local_attrib, global_attrib, g2latt;
  MAttType atttype;


  global_attrib = MESH_AttribByName(mesh,attr_name);
  if(!global_attrib) {
    MSTK_Report("MESH_CopyAttr","Root mesh has no attribute of given name",MSTK_WARN);
    return 0;
  }

  g2latt = MESH_AttribByName(mesh,"Global2Local");

  local_attrib = (MAttrib_ptr *) malloc(num*sizeof(MAttrib_ptr));
  for (i = 0; i < num; i++) {
    local_attrib[i] = MESH_AttribByName(submeshes[i],attr_name);
    /* if there is no such attribute */
    if(!local_attrib[i]) {
      MSTK_Report("MESH_CopyAttr","Submesh has no attribute of given name",MSTK_WARN);
      return 0;
    }
  }
    

  /* get attribute properties */
  ncomp = MAttrib_Get_NumComps(global_attrib);
  mtype = MAttrib_Get_EntDim(global_attrib);
  atttype = MAttrib_Get_Type(global_attrib);

  if (atttype == POINTER) {
    free(local_attrib);
    return 0; /* Don't see why one would want to
                 transmit pointer info to the submesh */
  }
  
  /* attribute entity type */
  
  if (mtype == MVERTEX || mtype == MALLTYPE) {
    MVertex_ptr gmv, lmv;
    List_ptr lmvlist;

    idx = 0;
    while ((gmv = MESH_Next_Vertex(mesh,&idx))) {

      MEnt_Get_AttVal(gmv,global_attrib,&ival,&rval,&pval);

      /* Don't copy null values */

      if ((atttype == INT && ival == 0) || 
          (atttype == DOUBLE  && rval == 0.0) ||
          (atttype == POINTER && pval == NULL))
        continue;

      int dummy_ival;
      double dummy_rval;
      MEnt_Get_AttVal(gmv,g2latt,&dummy_ival,&dummy_rval,&lmvlist);
      if (!lmvlist) continue;

      int idx2 = 0;
      while ((lmv = List_Next_Entry(lmvlist,&idx2))) {
        if (ncomp > 1) {
          pval_arr = (void *)malloc(ncomp*sizeof(double));
          memcpy(pval_arr,pval,ncomp*sizeof(double));
        }
        else
          pval_arr = NULL;

        Mesh_ptr entmesh = MEnt_Mesh(lmv);
        for (i = 0; i < num; i++)
          if (entmesh == submeshes[i]) {
            MEnt_Set_AttVal(lmv,local_attrib[i],ival,rval,(void *)pval_arr);
            break;
          }
      }
    }
  }

  if (mtype == MEDGE || mtype == MALLTYPE) {
    MVertex_ptr gme, lme;
    List_ptr lmelist;

    idx = 0;
    while ((gme = MESH_Next_Edge(mesh,&idx))) {

      MEnt_Get_AttVal(gme,global_attrib,&ival,&rval,&pval);

      /* Don't copy null values */

      if ((atttype == INT && ival == 0) || 
          (atttype == DOUBLE  && rval == 0.0) ||
          (atttype == POINTER && pval == NULL))
        continue;

      int dummy_ival;
      double dummy_rval;
      MEnt_Get_AttVal(gme,g2latt,&dummy_ival,&dummy_rval,&lmelist);
      if (!lmelist) continue;

      int idx2 = 0;
      while ((lme = List_Next_Entry(lmelist,&idx2))) {
        if (ncomp > 1) {
          pval_arr = (void *)malloc(ncomp*sizeof(double));
          memcpy(pval_arr,pval,ncomp*sizeof(double));
        }
        else
          pval_arr = NULL;

        Mesh_ptr entmesh = MEnt_Mesh(lme);
        for (i = 0; i < num; i++)
          if (entmesh == submeshes[i]) {
            MEnt_Set_AttVal(lme,local_attrib[i],ival,rval,(void *)pval_arr);
            break;
          }
      }
    }
  }

  if (mtype == MFACE || mtype == MALLTYPE) {
    MFace_ptr gmf, lmf;
    List_ptr lmflist;

    idx = 0;
    while ((gmf = MESH_Next_Face(mesh,&idx))) {

      MEnt_Get_AttVal(gmf,global_attrib,&ival,&rval,&pval);

      /* Don't copy null values */

      if ((atttype == INT && ival == 0) || 
          (atttype == DOUBLE  && rval == 0.0) ||
          (atttype == POINTER && pval == NULL))
        continue;

      int dummy_ival;
      double dummy_rval;
      MEnt_Get_AttVal(gmf,g2latt,&dummy_ival,&dummy_rval,&lmflist);
      if (!lmflist) continue;

      int idx2 = 0;
      while ((lmf = List_Next_Entry(lmflist,&idx2))) {
        if (ncomp > 1) {
          pval_arr = (void *)malloc(ncomp*sizeof(double));
          memcpy(pval_arr,pval,ncomp*sizeof(double));
        }
        else
          pval_arr = NULL;

        Mesh_ptr entmesh = MEnt_Mesh(lmf);
        for (i = 0; i < num; i++)
          if (entmesh == submeshes[i]) {
            MEnt_Set_AttVal(lmf,local_attrib[i],ival,rval,(void *)pval_arr);
            break;
          }
      }
    }
  }


  if (mtype == MREGION || mtype == MALLTYPE) {
    MRegion_ptr gmr, lmr;
    List_ptr lmrlist;

    idx = 0;
    while ((gmr = MESH_Next_Region(mesh,&idx))) {

      MEnt_Get_AttVal(gmr,global_attrib,&ival,&rval,&pval);

      /* Don't copy null values */

      if ((atttype == INT && ival == 0) || 
          (atttype == DOUBLE  && rval == 0.0) ||
          (atttype == POINTER && pval == NULL))
        continue;

      int dummy_ival;
      double dummy_rval;
      MEnt_Get_AttVal(gmr,g2latt,&dummy_ival,&dummy_rval,&lmrlist);
      if (!lmrlist) continue;

      int idx2 = 0;
      while ((lmr = List_Next_Entry(lmrlist,&idx2))) {
        if (ncomp > 1) {
          pval_arr = (void *)malloc(ncomp*sizeof(double));
          memcpy(pval_arr,pval,ncomp*sizeof(double));
        }
        else
          pval_arr = NULL;

        Mesh_ptr entmesh = MEnt_Mesh(lmr);
        for (i = 0; i < num; i++)
          if (entmesh == submeshes[i]) {
            MEnt_Set_AttVal(lmr,local_attrib[i],ival,rval,(void *)pval_arr);
            break;
          }
      }
    }
  }

  free(local_attrib);
  
  return 1;
}
  
#ifdef __cplusplus
}
#endif

