#define _H_MFace_Private

#include "MFace.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MF_Set_Edges_FN(MFace_ptr f, int n, MEdge_ptr *e, int *dir) {
    int i;
    MFace_DownAdj_FN *downadj;

#ifdef DEBUG
    if (n > 32)
      MSTK_Report("MF_Set_Edges_FN","Currently, only 32 edges supported per face",ERROR);
#endif

    downadj = (MFace_DownAdj_FN *) f->downadj

    downadj->ne = n;
    downadj->edirs = 0;
    downadj->fedges = List_New(n);
    
    for (i = 0; i < n; i++) {
      downadj->edirs = downadj->edirs | (dir[i] << i);
      List_Add(downadj->fedges,e[i]);
      ME_Add_Face(e[i],f);
    }
  }



  /* Replace the i'th edge in a face with a set of edges. Assumption
     is that the new edges are consecutive */

/*   void MF_Replace_Edge_i_FN(MFace_ptr f, int i, int nnu, MEdge_ptr *nuedges, int *nudirs) { */
/*     MFace_DownAdj_FN *downadj; */
/*     MEdge_ptr fedge; */
/*     int j, k, dir; */

/*     downadj = (MFace_DownAdj_FN *) f->downadj; */

/*     if (downadj->ne == 0) */
/*       MSTK_Report("MF_Replace_Edge_i","No initial set of edges for face", */
/* 		  ERROR); */

/*     fedge = List_Entry(downadj->fedges,i); */
/*     downadj->edirs = (downadj->edirs & ~(1<<i)); /\* set bit i to 0 *\/ */

/*     /\* First one is easy - Just replace the old edge with the new one *\/ */
/*     List_Replacei(downadj->fedges,i,nuedges[0]); */
/*     downadj->edirs = (downadj->edirs | (nudirs[0]<<i)); /\* set to dir *\/ */

/*     /\* Insert the rest of the nuedges at the right place *\/ */
/*     for (j = 1; j < nnu; j++) */
/*       List_Inserti(downadj->fedges,nuedges[j],i+j); */

/*     /\* Move the direction bits after the i'th bit to the left by */
/*        (nnu-1) spaces and insert the direction bits for the rest of */
/*        the inserted edges *\/ */

/*     for (k = downadj->ne-1; k > i; k--) { */

/*       /\* set bit (k+nnu-1) to 0 *\/ */
/*       downadj->edirs = downadj->edirs & ~(1<<(k+nnu-1)); */

/*       /\* get bit k *\/ */
/*       dir = (downadj->edirs>>k) & 1; */

/*       /\* move bit k to (k+nnu-1)'th position *\/ */
/*       downadj->edirs = downadj->edirs | (dir<<(k+nnu-1)); */

/*       /\* set bit k to 0 *\/ */
/*       downadj->edirs = downadj->edirs & ~(1<<k); */
/*     } */

/*     /\* set bit j according to input *\/ */
/*     for (j = 1; j < nnu; j++) */
/*       downadj->edirs = downadj->edirs | (nudirs[j]<<(i+j)); */

/*     /\* One edge replaced an existing edge, others were added *\/ */
/*     downadj->ne += nnu-1;  */

/*     /\* Update upward adjacencies *\/ */

/*     ME_Rem_Face(fedge,f); */
/*     for (i = 0; i < nnu; i++) */
/*       ME_Add_Face(nuedges[i],f); */
/*   } */



  /* Replace an edge in a face with another set of edges. Assumption
     is that the edges in the new set are consecutive */

/*   void MF_Replace_Edge_FN(MFace_ptr f, MEdge_ptr e, int nnu, MEdge_ptr *nuedges, int *nudirs) { */
/*     int i, found = 0; */
/*     MFace_DownAdj_FN *downadj; */
/*     MEdge_ptr fedge; */

/*     downadj = (MFace_DownAdj_FN *) f->downadj; */


/*     if (downadj->ne == 0) */
/*       MSTK_Report("MF_Replace_Edge_F1","No initial set of edges for face", */
/* 		  ERROR); */
      
/*     for (i = 0; i < downadj->ne; i++) */
/*       if ((fedge = List_Entry(downadj->fedges,i)) == e) { */
/* 	found = 1; */
/* 	break; */
/*       } */
/*     if (!found) { */
/*       MSTK_Report("MF_Replace_Edge","Edge not found in face",ERROR); */
/*       return; */
/*     } */

/*     MF_Replace_Edge_i_FN(f, i, nnu, nuedges, nudirs); */
/*   } */



  /* Replace a set of 'n' edges, starting with edge 'i', in a face
     with another set of edges. Assumption is that the edges are
     consecutive */

  void MF_Replace_Edges_i_FN(MFace_ptr f, int nold, int i, int nnu, 
			     MEdge_ptr *nuedges) {
    MFace_DownAdj_FN *downadj;
    MEdge_ptr fedge, *oldedges, *newedges;
    MVertex_ptr vold_0, vold_1, lastv;
    int j, k, ne, ncom, dir, *olddirs, *newdirs, rev;

    downadj = (MFace_DownAdj_FN *) f->downadj;
    ne = downadj->ne;

    if (ne == 0)
      MSTK_Report("MF_Replace_Edge_i","No initial set of edges for face",
		  ERROR);

    oldedges = (MEdge_ptr *) malloc(nold*sizeof(MEdge_ptr));
    olddirs  = (int *) malloc(nold*sizeof(int));
    for (j = 0; j < nold; j++) {
      k = (i+j)%ne;
      oldedges[j] = List_Entry(downadj->fedges,k);
      olddirs[j] = (downadj->edirs>>k) & 1;
    }

    rev = 0;
    vold_0 = ME_Vertex(oldedges[0],!olddirs[0]);
    vold_1 = ME_Vertex(oldedges[nold-1],olddirs[nold-1]);

    /* Assume that nuedges[0] or nuedges[nnu-1] are connected to vold_0 */
    /* (likewise for vold_1) and that the other edges are in sequence  */

    if (ME_UsesEntity(nuedges[0],vold_0,MVERTEX) && ME_UsesEntity(nuedges[nnu-1],vold_1,MVERTEX))
      rev = 0;
    else if (ME_UsesEntity(nuedges[0],vold_1,MVERTEX) && ME_UsesEntity(nuedges[nnu-1],vold_0,MVERTEX))
      rev = 1;
    else
      MSTK_Report("MF_Replace_Edge_i","Mismatched set of edges",ERROR);


    newedges = (MEdge_ptr *) malloc(nnu*sizeof(MEdge_ptr));
    newdirs  = (int *) malloc(nnu*sizeof(MEdge_ptr));
    lastv = vold_0;
    for (j = 0; j < nnu; j++) {
      k = rev ? nnu-1-j : j;
      newedges[j] = nuedges[k];      
      newdirs[j] = (ME_Vertex(newedges[j],0) == lastv) ? 1 : 0;
      lastv = ME_Vertex(newedges[j],newdirs[j]);	
    }

    /* Number of common edges */

    ncom = nold < nnu ? nold : nnu;


    /* Replace old edges with new ones (up to the common number of edges) */

    for (j = 0; j < ncom; j++) {
      k = (i+j)%ne; 
      List_Replacei(downadj->fedges,k,newedges[j]);		      
      downadj->edirs = (downadj->edirs & ~(1<<k)); /* set bit k to 0 */
      downadj->edirs = (downadj->edirs | (newdirs[j]<<k)); /* set to dir */
    }



    if (nold < nnu) {
      /* Insert the rest of the nuedges at the right place */

      for (j = ncom; j < nnu; j++) {
	List_Inserti(downadj->fedges,newedges[j],(i+j)%ne);

      
	/* Move the direction bits after the (i+ncom)'th bit to the
	 left by (nnu-ncom) spaces and insert the direction bits for
	 the rest of the inserted edges */

	
	for (k = ne-1; k >= (i+j)%ne; k--) {

	  /* set bit k+1 to 0 */
	  downadj->edirs = downadj->edirs & ~(1<<(k+1));

	  /* get bit k */
	  dir = (downadj->edirs>>k) & 1;

	  /* move bit k to (k+1)'th position */
	  downadj->edirs = downadj->edirs | (dir<<(k+1));	
	}

	/* set bit j according to input */
	downadj->edirs = downadj->edirs & ~(1<<(i+j)%ne); /* Set bit to 0 */
	downadj->edirs = downadj->edirs | (newdirs[j]<<((i+j)%ne)); /* set to newdir */

	ne++;
      }

    }
    else {
      /* Delete the remaining old edges */

      for (j = ncom; j < nold; j++)
	List_Rem(downadj->fedges,oldedges[j]);

      /* Delete (nold-ncom) direction bits after the (i+ncom)'th bit and
	 move the remaining bits to the right */
      
      for (j = 0; j < nold-ncom; j++) {
	for (k = (i+ncom+j)%ne; k < ne-1; k++) {
	  /* set bit k to 0 */
	  downadj->edirs = downadj->edirs & ~(1<<k);

	  /* get bit k+1 */
	  dir = (downadj->edirs>>(k+1)) & 1;

	  /* move bit k+1 to k'th position */
	  downadj->edirs = downadj->edirs | (dir<<k);
	}

	ne--;
      }
    }
      
    /* ncom edges were replaced existing edges, others were added */
    downadj->ne = ne;  /* also equal to "downadj->ne - nold + nnu" */ 

    /* Update upward adjacencies */
    for (j = 0; j < nold; j++) {
      ME_Rem_Face(oldedges[j],f);
    }
    for (j = 0; j < nnu; j++) {
      ME_Add_Face(newedges[j],f);
    }

    free(oldedges);
    free(olddirs);
    free(newedges);
    free(newdirs);
  }


  /* Replace a set of edges in a face with another set of
     edges. Assumption is that the edges in each set (old and new) are
     consecutive */

  void MF_Replace_Edges_FN(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges) {
    int i, j, ne, found = 0, imin1, istart, *eindex, nxtidx, ipos1;
    MFace_DownAdj_FN *downadj;
    MEdge_ptr fedge;

    downadj = (MFace_DownAdj_FN *) f->downadj;
    ne = downadj->ne;

    if (ne == 0)
      MSTK_Report("MF_Replace_Edge_F1","No initial set of edges for face",
		  ERROR);

    /* Order the old edges in the direction of the face */

    eindex = (int *) malloc(nold*sizeof(int));

    for (i = 0; i < nold; i++) {           
      for (j = 0, found = 0; j < ne; j++)
	if ((fedge = List_Entry(downadj->fedges,j)) == oldedges[i]) {
	  found = 1;
	  eindex[i] = j;
	  break;
	}
      if (!found) {
	MSTK_Report("MF_Replace_Edge","Edge not found in face",ERROR);
	free(eindex);
	return;
      }
    }

    /* Find the lowest index; set it to be the start index for
       replacing edges */

    imin1 = 10000; ipos1 = -1;
    for (i = 0; i < nold; i++)
      if (eindex[i] < imin1) {
	imin1 = eindex[i];
	ipos1 = i;
      }
    istart = imin1;
    
    
    /* Check if there is a break in the edge indices sorted in
       ascending order. If there is a break, that means the edges are
       spanning the start vertex in the face vertex cyle.*/
	
    nxtidx = imin1+1;
    for (i = 1; i < nold; i++) {
      j = (ipos1+i)%nold;

      if (eindex[j] != nxtidx) {
	/* There is a break in the edge indices and this is the next
	   smallest index. Use this as the start index for replacing
	   the edges */

	istart = eindex[j];
	break;
      }
      else
	nxtidx++;
    }
      
	  
    MF_Replace_Edges_i_FN(f, nold, istart, nnu, nuedges);

    free(eindex);
  }





  int MF_Num_Edges_FN(MFace_ptr f) {
    return ((MFace_DownAdj_FN *) f->downadj)->ne;
  }	

  List_ptr MF_Edges_FN(MFace_ptr f, int dir, MVertex_ptr v0) {
    int i, k=0, n, fnd, edir;
    List_ptr fedges;
    MEdge_ptr e;
    MFace_DownAdj_FN *downadj;

    downadj = (MFace_DownAdj_FN *) f->downadj;


    if (v0 == NULL) {
      if (dir)
	return List_Copy(downadj->fedges);
      else {
	n = downadj->ne;
	fedges = List_New(n);
	for (i = n-1; i >= 0; i--)
	  List_Add(fedges,List_Entry(downadj->fedges,i));
	return fedges;
      }
    }
    else {
      n = downadj->ne;
      fnd = 0;
      for (i = 0; i < n; i++) {
	e = List_Entry(downadj->fedges,i);
	edir = ((downadj->edirs)>>i) & 1;
	if (ME_Vertex(e,edir^dir) == v0) {
	  fnd = 1;
	  k = i;
	  break;
	}
      }

      if (!fnd)
	MSTK_Report("MF_Edges_F1","Cannot find vertex in face!!",FATAL);
	
      fedges = List_New(n);
      for (i = 0; i < n; i++) {
	e = dir ? List_Entry(downadj->fedges,(k+i)%n) :
	  List_Entry(downadj->fedges,(k+n-i)%n);
	List_Add(fedges,e);
      }	
    }

    return fedges;
  }

  int MF_EdgeDir_FN(MFace_ptr f, MEdge_ptr e) {
    int n, i;
    MFace_DownAdj_FN *downadj;

    downadj = (MFace_DownAdj_FN *) f->downadj;
    
    n = downadj->ne;
    for (i = 0; i < n; i++) {
      if (List_Entry(downadj->fedges,i) == e)
	return ((downadj->edirs)>>i) & 1;
    }

    MSTK_Report("MF_Edges_FN","Cannot find edge in face!!",FATAL);

    return 0;
  }

  int MF_EdgeDir_i_FN(MFace_ptr f, int i) {
    MFace_DownAdj_FN *downadj;

#ifdef DEBUG
    if (i > 31)
      MSTK_Report("MF_Set_Edges_FN","Currently, only 32 edges supported per face",ERROR);
#endif

    downadj = (MFace_DownAdj_FN *) f->downadj;

    return ((downadj->edirs)>>i) & 1;
  }
			
  int MF_UsesEdge_FN(MFace_ptr f, MEdge_ptr e) {
    MFace_DownAdj_FN *downadj;

    downadj = (MFace_DownAdj_FN *) f->downadj;

    return List_Contains(downadj->fedges,e);
  }

#ifdef __cplusplus
}
#endif
