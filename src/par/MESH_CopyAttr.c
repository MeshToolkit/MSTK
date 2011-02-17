#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "MSTK.h"

#ifdef __cplusplus
extern "C" {
#endif

 /* 
    this function copy attribute information from mesh to submesh
    assume mesh and submesh have the same attributes list

    attr_name: the attribute name

 */

int MESH_CopyAttr(Mesh_ptr mesh, Mesh_ptr submesh, const char *attr_name) {
  int j, num, ncomp, ival;
  double rval, *pval_arr;
  void *pval;
  MType mtype;

  MEntity_ptr local_ment, global_ment;
  MAttrib_ptr local_attrib, global_attrib;

  global_attrib = MESH_AttribByName(mesh,attr_name);
  local_attrib = MESH_AttribByName(submesh,attr_name);
  /* if there is no such attribute */
  if(!global_attrib || !local_attrib) {
    MSTK_Report("MESH_CopyAttr","There is no such attribute",WARN);
    return 0;
  }

  /* get attribute properties */
  ncomp = MAttrib_Get_NumComps(global_attrib);
  mtype = MAttrib_Get_EntDim(global_attrib);


  /* attribute entity type */
  switch (mtype) {
  case MVERTEX:
    num = MESH_Num_Vertices(submesh);
    break;
  case MEDGE:
    num = MESH_Num_Edges(submesh);
    break;
  case MFACE:
    num = MESH_Num_Faces(submesh);
    break;
  case MREGION:
    num = MESH_Num_Regions(submesh);
    break;
  default:
    num = 0;
    MSTK_Report("MESH_CopyAttr()","Invalid entity type",FATAL);
  }
  
  /* collect data */
  for(j = 0; j < num; j++) {
    switch (mtype) {
    case MVERTEX:
      local_ment = MESH_Vertex(submesh,j);
      global_ment = MESH_VertexFromID(mesh,MEnt_GlobalID(local_ment));
      break;
    case MEDGE:
      local_ment = MESH_Edge(submesh,j);
      global_ment = MESH_EdgeFromID(mesh,MEnt_GlobalID(local_ment));
      break;
    case MFACE:
      local_ment = MESH_Face(submesh,j);
      global_ment = MESH_FaceFromID(mesh,MEnt_GlobalID(local_ment));
      break;
    case MREGION:
      local_ment = MESH_Region(submesh,j);
      global_ment = MESH_RegionFromID(mesh,MEnt_GlobalID(local_ment));
      break;
    default:
      MSTK_Report("MESH_CopyAttr()","Invalid entity type",FATAL);
    }
    
    MEnt_Get_AttVal(global_ment,global_attrib,&ival,&rval,&pval);
    
    if(ncomp > 1) {
      pval_arr = (double *)MSTK_malloc(ncomp*sizeof(double));
      pval_arr = (double *)pval;
    }	  
    MEnt_Set_AttVal(local_ment,local_attrib,ival,rval,(void*) pval_arr);
  }
  
  return 1;
}
  
#ifdef __cplusplus
}
#endif

