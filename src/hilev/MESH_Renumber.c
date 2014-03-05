#include <stdio.h>
#include <stdlib.h>
#include "MSTK_private.h"
#include "MSTK.h"

/* Enforce continuous and possibly optimized numbering for mesh entities */
/* renum_type = 0   --- Sequential renumbering 
              = 1   --- Reverse Cuthill Mckee 
*/



void MESH_Renumber(Mesh_ptr mesh, int renum_type, MType mtype) {
  MVertex_ptr mv, v0;
  MEdge_ptr me;
  MFace_ptr mf;
  MRegion_ptr mr;
  int idx, idx2;
  int nv, ne, nf, nr;
  int done;
  MAttrib_ptr vidatt;
  List_ptr vlist;
  double xyz[3];

  if (renum_type == 0) {
    if (mtype == MVERTEX || mtype == MALLTYPE) {
      idx = 0; nv = 0;
      while ((mv = MESH_Next_Vertex(mesh,&idx)))
	MV_Set_ID(mv,++nv);
    }
    
    if (mtype == MEDGE || mtype == MALLTYPE) {
      idx = 0; ne = 0;
      while ((me = MESH_Next_Edge(mesh,&idx)))
	ME_Set_ID(me,++ne);
    }
    
    if (mtype == MFACE || mtype == MALLTYPE) {
      idx = 0; nf = 0;
      while ((mf = MESH_Next_Face(mesh,&idx)))
	MF_Set_ID(mf,++nf);
    }
    
    if (mtype == MREGION || mtype == MALLTYPE) {
      idx = 0; nr = 0;
      while ((mr = MESH_Next_Region(mesh,&idx)))
	MR_Set_ID(mr,++nr);
    }
  }
  else if (renum_type == 1) {
  }
  

  vidatt = MAttrib_New(mesh,"vidrcm",INT,MVERTEX);
  idx = 0;
  while ((mv = MESH_Next_Vertex(mesh,&idx))) {
    fprintf(stderr,"MV_ID = %d\n",MV_ID(mv));
    MEnt_Set_AttVal(mv,vidatt,MV_ID(mv),0.0,NULL);
  }
    
  return;
}


