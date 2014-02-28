#include <stdio.h>
#include <stdlib.h>
#include "MSTK.h"

/* Enforce continuous and possibly optimized numbering for mesh entities */
/* type = 0   --- Just remove holes
        = 1   --- Reverse Cuthill Mckee
*/



void MESH_Renumber(Mesh_ptr mesh, int type) {
  MVertex_ptr mv;
  MEdge_ptr me;
  MFace_ptr mf;
  MRegion_ptr mr;
  int idx, idx2;
  int nv, ne, nf, nr, nve, nve2;
  MAttrib_ptr vidatt;

  if (type == 0) {
    idx = 0; nv = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx)))
      MV_Set_ID(mv,++nv);
    
    idx = 0; ne = 0;
    while ((me = MESH_Next_Edge(mesh,&idx)))
      ME_Set_ID(me,++ne);
    
    idx = 0; nf = 0;
    while ((mf = MESH_Next_Face(mesh,&idx)))
      MF_Set_ID(mf,++nf);
    
    idx = 0; nr = 0;
    while ((mr = MESH_Next_Region(mesh,&idx)))
      MR_Set_ID(mr,++nr);
  }
  else if (type == 1) {

    mkid = MSTK_GetMarker();

    /* Renumber vertices using the Reverse Cuthill Mckee algorithm */

    /* First renumber the vertices */
    /* Start with the vertex in the lower leftmost corner */

    minx = miny = minz = 1.0e+12;
    v0 = NULL;
    idx = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx))) {
      MV_Coords(mv,&xyz);
      if (xyz[0] <= minx && xyz[1] <= miny && xyz[2] <= minz) {
        minx = xyz[0];
        miny = xyz[1];
        minz = xyz[2];
        v0 = mv;
      }
    }

    
    vlist = List_New(MV_Num_Vertices(mesh));
    List_Add(vlist,v0);

    nv = 0;
    done = 0;
    while ((mv = List_Next_Entry(vlist,&idx))) {

      MV_Set_ID(mv,++nv);

      /* Add the neighbors of the vertex to the list in increasing
         order of their valence */

      vedges = MV_Edges(mv);
      nve = List_Num_Entries(vedges);
      oppvlist = List_New(nve);

      idx2 = 0;
      nve2 = 0;
      while ((ve = List_Next_Entry(oppvlist,&idx2))) {
        oppv = MV_OppVertex(mv);
        if (MV_PType(oppv) != PGHOST &&
            !MEnt_IsMarked(oppv,mkid)) {
          List_Add(oppvlist,oppv);
          MEnt_Mark(oppv,mkid);
          nve2++;
        }
      }
      List_Delete(vedges);

      List_Sort(oppvlist,nve2,sizeof(void *),compareValence);

      vlist = List_Cat(vlist,oppvlist);
      List_Delete(oppvlist);
    }

    /* Number the vertices in reverse to get the Reverse Cuthill McKee
       ordering */
    
    idx = 0;
    while ((mv = List_Next_Entry(vlist,&idx))) {
      MV_Set_ID(mv,nv--);
    }
  }
  

  vidatt = MAttrib_New(mesh,"vidrcm",INT,MVERTEX);
  idx = 0;
  while ((mv = MESH_Next_Vertex(mesh,&idx)))
    MEnt_Set_AttVal(mv,vidatt,MV_ID(mv),0.0,NULL);
    
  return;
}


