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

  if (type == 0) {
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
  else if (type == 1) {
    double minx, miny, minz;
    MVertex_ptr v0;
    int mkid = MSTK_GetMarker();


    /* Renumber using the Reverse Cuthill Mckee algorithm */

    if (mtype == MVERTEX || mtype == MALLTYPE) { 
      /* First renumber the vertices */
      /* Start with the vertex in the lower leftmost corner */
      
      minx = miny = minz = 1.0e+12;
      v0 = NULL;
      idx = 0;
      while ((mv = MESH_Next_Vertex(mesh,&idx))) {
	MV_Coords(mv,xyz);
	if (xyz[0] <= minx && xyz[1] <= miny && xyz[2] <= minz) {
	  minx = xyz[0];
	  miny = xyz[1];
	  minz = xyz[2];
	  v0 = mv;
	}
      }

    
      vlist = List_New(MESH_Num_Vertices(mesh));
      List_Add(vlist,v0);

      nv = 0;
      done = 0;
      idx = 0;
      while ((mv = List_Next_Entry(vlist,&idx))) {
	List_ptr oppvlist, vedges;
	MVertex_ptr oppv;
	MEdge_ptr ve;
	int nve, nve2;

	MV_Set_ID(mv,++nv);

	/* Add the neighbors of the vertex to the list in increasing
	   order of their valence */

	vedges = MV_Edges(mv);
	nve = List_Num_Entries(vedges);
	oppvlist = List_New(nve);

	idx2 = 0;
	nve2 = 0;
	while ((ve = List_Next_Entry(vedges,&idx2))) {
	  oppv = ME_OppVertex(ve,mv);
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
      while ((mv = List_Next_Entry(vlist,&idx)))
	MV_Set_ID(mv,nv--);
    }
 
    /* Reorder edges according to a breadth first algorithm applied to
       edges (differs from RCM in that it does not add adjacent nodes
       in ascending order of their valence) */

    if (mtype == MEDGE || mtype == MALLTYPE) {

      if (mtype == MALLTYPE) { 
	/* RCM algorithm already applied on the vertices. Use an edge
	   connected to the starting vertex as the first edge */

	vedges = MV_Edges(v0);
	e0 = List_Entry(vedges,0);
	List_Delete(vedges);
      }
      else {
	/* Find the edge whose mid point is a minimum point */
	minx = miny = minz = 1.0e+12;
	e0 = NULL;
	idx = 0;
	while ((me = MESH_Next_Edge(mesh,&idx))) {
	  MV_Coords(ME_Vertex(me,0),exyz[0]);
	  MV_Coords(ME_Vertex(me,1),exyz[1]);
	  xyz[0] = (exyz[0][0]+exyz[1][0])/2.0;
	  xyz[1] = (exyz[0][1]+exyz[1][1])/2.0;
	  xyz[2] = (exyz[0][2]+exyz[1][2])/2.0;
	  if (xyz[0] < minx && xyz[1] < miny && xyz[2] < minz) {
	    minx = xyz[0];
	    miny = xyz[1];
	    minz = xyz[2];
	    e0 = me;
	  }
	}
      }

      

    }

    MSTK_FreeMarker(mkid);
  }
  

  vidatt = MAttrib_New(mesh,"vidrcm",INT,MVERTEX);
  idx = 0;
  while ((mv = MESH_Next_Vertex(mesh,&idx))) {
    fprintf(stderr,"MV_ID = %d\n",MV_ID(mv));
    MEnt_Set_AttVal(mv,vidatt,MV_ID(mv),0.0,NULL);
  }
    
  return;
}


