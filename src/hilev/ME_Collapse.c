#include "MSTK_private.h"


/*                                                                      */
/* Collapse an edge                                                     */
/*   All elements connected to this edge will be collapsed/deleted and  */
/* entities that become overlapping will be merged                      */
/*                                                                      */
/* THIS ROUTINE DOES NOT CHECK IF THE RESULTING ELEMENTS ARE            */
/* GEOMETRICALLY VALID!                                                 */
/*                                                                      */
/* e is the edge to be collapsed                                        */
/*                                                                      */
/* vkeep_in is the vertex of the edge that must be preserved            */
/* If vkeep_in is NULL, then either vertex might be collapsed to the    */
/* other.                                                               */
/*                                                                      */
/* if topoflag is 1, then topological constraints are respected during  */
/* the collapse. That means that the collapse will take place only if   */
/* it does not violate topological conformity with the geometric        */
/* model. For example, the routine will not collapse the edge, if its   */
/* vertices are classified on two different model entities of the same  */
/* dimension (say two different model faces) or if the vertex to be     */
/* retained is classified on a higher dimensional entity than the       */
/* vertex to be deleted (vdel). If topoflag is 0, then the collapse is  */
/* performed without consideration to the classification of entities    */
/*                                                                      */
/* The routine returns the retained vertex if collapse was successful   */


MVertex_ptr ME_Collapse(MEdge_ptr e, MVertex_ptr vkeep_in, int topoflag) {
  MVertex_ptr vdel, vkeep, ev00, ev01, ev10, ev11, vert;
  MEdge_ptr edge, edge2, oldedges[3], nuedges[2];
  MFace_ptr face, face2, rface1, rface2;
  MRegion_ptr reg, reg2;
  List_ptr vedges, efaces, eregs, fedges, rfaces, fverts1, fverts2, vfaces;
  int idx1, idx2, idx3, dir, status, nfe, nrf, allfound, degenerate;
  int i, j, nfe2, nfv1, nfv2;

  status = 1;

  if (vkeep_in == NULL) {
    vdel = ME_Vertex(e,0);
    vkeep = ME_Vertex(e,1);
  }
  else {
    vkeep = vkeep_in;
    vdel = ME_OppVertex(e,vkeep);
  }

  if (topoflag == 1) {  /* keep topological conformity with geometric model */
    int dimkeep, dimdel;
    
    dimkeep = MV_GEntDim(vkeep); /* Model entity dim of vertex to keep */
    dimdel = MV_GEntDim(vdel);   /* Model entity dim  of vertex to delete */
    
    if (dimkeep == dimdel) {
      if (MV_GEntID(vkeep) != MV_GEntID(vdel))
	status = 0;                 /* cannot allow since it will cause 
				       a dimensional reduction in mesh */
    }
    else if (dimdel < dimkeep) {
      status = 0;                  /* boundary of mesh will get messed up */

      if (vkeep_in == NULL) {
	
	/* If no preference was indicated on which vertex to retain,
	   we can collapse in the other direction */

	vdel = ME_Vertex(e,1);
	vkeep = ME_Vertex(e,0);
	status = 1;
      }
    }
  }


  if (status == 0)
    return NULL;   /* Cannot collapse due to constraints of topological
		   conformity with geometric model */


  /* Need to collect this in advance because the info gets messed up later */

  efaces = ME_Faces(e);
  eregs = ME_Regions(e);


  /* Replace vdel with vkeep in all edges connected to vdel */

  vedges = MV_Edges(vdel);
  idx1 = 0;
  while ((edge = List_Next_Entry(vedges,&idx1))) {
    ME_Replace_Vertex(edge,vdel,vkeep);
  }
  List_Delete(vedges);

  /* Remove edge 'e' from all faces connected to e */
  /* This part of the code is using some reliance on the internal
     implementation of MF_Edges. While unlikely, it _might_ break if
     the innards of MF_Edges are changed */

  
  idx1 = 0;
  while ((face = List_Next_Entry(efaces,&idx1))) {

    fedges = MF_Edges(face,1,0);
    nfe = List_Num_Entries(fedges);

    /* Find the edge before and after e in the face */

    oldedges[0] = oldedges[2] = NULL;
    for (i = 0; i < nfe; i++) {
      edge = List_Entry(fedges,i);
      if (edge == e) continue;

      dir = MF_EdgeDir_i(face,i);

      if (ME_Vertex(edge,dir) == vkeep)
	oldedges[0] = edge;
      else if (ME_Vertex(edge,!dir) == vkeep)
	oldedges[2] = edge;
    }
    oldedges[1] = e;

    nuedges[0] = oldedges[0];
    nuedges[1] = oldedges[2];


    /* Replace oldedges[0], oldedges[1] (=e), oldedges[2] with 
       oldedges[0], oldedges[2] since e is degenerate */

    MF_Replace_Edges(face,3,oldedges,2,nuedges);

    List_Delete(fedges);
  }



  /* Delete topologically degenerate regions */
  /* Defined as two faces of the regions having the same vertices */

  idx1 = 0;
  while ((reg = List_Next_Entry(eregs,&idx1))) {

    rfaces = MR_Faces(reg);
    nrf = List_Num_Entries(rfaces);

    if (nrf == 4) {
      
      /* This is a tet - it will become degenerate */

      MR_Delete(reg,0);

    }
    else {

      degenerate = 0;

      for (i = 0; i < nrf; i++) {

	rface1 = List_Entry(rfaces,i);
	
	fverts1 = MF_Vertices(rface1,1,0);
	nfv1 = List_Num_Entries(fverts1);
		
	for (j = i+1; j < nrf; j++) {
	  
	  rface2 = List_Entry(rfaces,j);

	  fverts2 = MF_Vertices(rface2,1,0);
	  nfv2 = List_Num_Entries(fverts2);
	  
	  if (nfv1 != nfv2) {
	    List_Delete(fverts2);
	    continue;             /* can't be exactly coincident */
	  }

	  allfound = 1;
	  idx2 = 0;
	  while ((vert = List_Next_Entry(fverts2,&idx2))) {
	    if (!List_Contains(fverts1,vert)) {
	      allfound = 0;
	      break;
	    }
	  }
	  
	  List_Delete(fverts2);
	  
	  if (allfound) {
	    degenerate = 1;
	    break;
	  }
	  
	} /* for (j = i+1 ... */

	List_Delete(fverts1);

	if (degenerate) break;

      } /* for (i = 0; i < nrf;.... */

      if (degenerate) 
	MR_Delete(reg,0);

    } /* if (nrf == 4) .. else ... */

    List_Delete(rfaces);

  } /* while ((reg = ...)) */



  /* Delete topologically degenerate faces */

  idx1 = 0;
  while ((face = List_Next_Entry(efaces,&idx1))) {

    fedges = MF_Edges(face,1,0);

    if (List_Num_Entries(fedges) == 2) {

      /* Disconnect the regions from the face before deleting */

      List_ptr fregs = MF_Regions(face);

      idx2 = 0;
      while ((reg = List_Next_Entry(fregs,&idx2)))
	MF_Rem_Region(face,reg);

      List_Delete(fregs);


      MF_Delete(face,0);
    }

    List_Delete(fedges);
  }



  /* Now merge edges which have the same end vertices */

  vedges = MV_Edges(vkeep);
  idx1 = 0; 
  while ((edge = List_Next_Entry(vedges,&idx1))) {
    if (edge == e) continue;
    
    ev00 = ME_Vertex(edge,0);
    ev01 = ME_Vertex(edge,1);

    idx2 = 0;
    while ((edge2 = List_Next_Entry(vedges,&idx2))) {
      if (edge == e || edge == edge2) continue;

      ev10 = ME_Vertex(edge2,0);
      ev11 = ME_Vertex(edge2,1);

      if ((ev00 == ev10 && ev01 == ev11) ||
	  (ev00 == ev11 && ev10 == ev01)) {
	    
	MEs_Merge(edge,edge2);
	
	List_Rem(vedges,edge2);
	
	break;
      }
    }
  }
  List_Delete(vedges);



  /* Merge faces with the same set of edges */
      
  vfaces = MV_Faces(vkeep);
  
  idx1 = 0;
  while ((face = List_Next_Entry(vfaces,&idx1))) {
    
    fedges = MF_Edges(face,1,0);
    nfe = List_Num_Entries(fedges);
    
    idx2 = 0;
    while ((face2 = List_Next_Entry(vfaces,&idx2))) {
      List_ptr fedges2;

      if (face2 == face) continue;

      fedges2 = MF_Edges(face2,1,0);
      nfe2 = List_Num_Entries(fedges2);

      if (nfe != nfe2) {
	List_Delete(fedges2);
	continue;
      }

      allfound = 1;

      for (i = 0; i < nfe2; i++) {
	edge = List_Entry(fedges2,i);
	if (!List_Contains(fedges,edge)) {
	  allfound = 0;
	  break;
	}
      }
      List_Delete(fedges2);

      if (allfound) {
	MFs_Merge(face,face2);	
	List_Rem(vfaces,face2);
	break;
      }
	
    } /* while (face2 = List_Next_Entry(vfaces,... */

    List_Delete(fedges);

  } /* while (face = List_Next_Entry(vfaces,... */
  List_Delete(vfaces);


  /* Now actually delete the collapse edge and the to-be-merged vertex */

  ME_Delete(e,0);
  MV_Delete(vdel,0);


  List_Delete(efaces);
  List_Delete(eregs);

  return vkeep;
}
  
