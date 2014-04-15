#define _H_MFace_Private

#include "MFace.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  int MF_Set_GInfo_Auto_F1F3(MFace_ptr f) {
    int i, ne, same, fgdim, fgid, egdim, egid, egdim0, egid0;
    MEdge_ptr e;
    MFace_Adj_F1F3 *adj;

    adj = (MFace_Adj_F1F3 *) f->adj;
    ne = List_Num_Entries(adj->fedges);

    same = 1;
    fgdim = -1;
    fgid = -1;
    egid = -1;

    e = List_Entry(adj->fedges,0);    
    egid0 = ME_GEntID(e);
    egdim0 = ME_GEntDim(e);

    for (i = 1; i < ne; i++) {
      e = List_Entry(adj->fedges,i);
      egid = ME_GEntID(e);
      egdim = ME_GEntDim(e);
      if (egdim == egdim0 && egid == egid0)
	continue; /* all edges have same classification so far */
      else {
	same = 0;
	/* I don't think this is right always
	if (egdim > fgdim) {
	  fgdim = egdim;
	  fgid = egid;
	}
	*/
	break;
      }
    }
    if (same) {
      fgdim = egdim0;
      fgid = egid;
    }
     
    if (fgdim == -1 || fgdim < 2) {
      List_ptr fregions;

      /* We are unable to find proper classification info from the
	 vertices. Lets look at the number of regions connected to the
	 face and their classification */

      fregions = MF_Regions(f);
      
      if (fregions == NULL || List_Num_Entries(fregions) == 1) {
	
	/* In a complete mesh, this face must be on a model face */

	fgdim = 2;

      }
      else {
	MRegion_ptr fregion0, fregion1;
	int rgid0, rgid1;

	/* Internal face. Check if it is a mesh region or on an
	   interior interface */

	fregion0 = List_Entry(fregions,0); rgid0 = MEnt_GEntID(fregion0);
	fregion1 = List_Entry(fregions,1); rgid1 = MEnt_GEntID(fregion1);

	if (rgid0 == -1 || rgid1 == -1) {
	  
	  /* One of the regions is not classified properly. Just
	     assume this is an internal face */

	  fgdim = 3;	  

	}
	else {

	  fgdim = (rgid0 == rgid1) ? 3 : 2;

	}
      }

      if (fregions) List_Delete(fregions);
    }

    MEnt_Set_GEntDim((MEntity_ptr) f,fgdim);
    MEnt_Set_GEntID((MEntity_ptr) f,fgid);

    if (fgdim == 4)
      return 0;
    else
      return 1;
  }


  void MF_Set_Edges_F1F3(MFace_ptr f, int n, MEdge_ptr *e, int *dir) {
    int i;
    MFace_Adj_F1F3 *adj;

#ifdef DEBUG
    if (n > 32)
      MSTK_Report("MF_Set_Edges_F1F3","Currently, only 32 edges supported per face",MSTK_ERROR);
#endif

    adj = (MFace_Adj_F1F3 *) f->adj;
    adj->edirs = 0UL;
    if (adj->fedges)
	List_Delete(adj->fedges);
    adj->fedges = List_New(n);
    
    for (i = 0; i < n; i++) {
#ifdef DEBUG
      if (MF_Mesh(f) != ME_Mesh(e[i]))
	MSTK_Report("MF_Set_Edges_F1F3",
		    "Face and edge belong to different meshes",MSTK_FATAL);
#endif

      adj->edirs = adj->edirs | (dir[i] << i);
      List_Add(adj->fedges,e[i]);
      ME_Add_Face(e[i],f);
    }
  }



  /* Remove an edge from a face */

  void MF_Rem_Edge_F1F3(MFace_ptr f, MEdge_ptr remedge) {
    int i, j, ne, eindex, found;
    MFace_Adj_F1F3 *adj;

    adj = (MFace_Adj_F1F3 *) f->adj;
    ne = List_Num_Entries(adj->fedges);

    if (ne == 0)
      MSTK_Report("MF_Remove_Edges_F1F3","No initial set of edges for face",
		  MSTK_ERROR);

    /* Order the old edges in the direction of the face */

    for (j = 0, found = 0; j < ne; j++)
      if (List_Entry(adj->fedges,j) == remedge) {
	found = 1;
	eindex = j;
	break;
      }
    if (!found) {
      MSTK_Report("MF_Remove_Edge_F1F3","Edge not found in face",MSTK_ERROR);
      return;
    }


    List_Remi(adj->fedges,eindex);

    /* Move the remaining bits to the right */
    int k;
    long long dir;
    for (k = eindex; k < ne-1; k++) {
      /* set bit k to 0 */
      adj->edirs = adj->edirs & ~(1UL<<k);

      /* get bit k+1 */
      dir = (adj->edirs>>(k+1)) & 1UL;

      /* move bit k+1 to k'th position */
      adj->edirs = adj->edirs | (dir<<k);
    }

    ne--;

      
    /* Tell the edge that it is not connected to the face anymore */

    ME_Rem_Face(remedge,f);

  }




  /* Replace a set of 'n' edges, starting with edge 'i', in a face
     with another set of edges. Assumption is that the edges are
     consecutive */

  void MF_Replace_Edges_i_F1F3(MFace_ptr f, int nold, int i, int nnu, 
			     MEdge_ptr *nuedges) {
    MFace_Adj_F1F3 *adj;
    MEdge_ptr *oldedges, *newedges;
    MVertex_ptr vold_0, vold_1, lastv;
    int j, k, ne, ncom, dir, *olddirs, *newdirs, rev;

    adj = (MFace_Adj_F1F3 *) f->adj;
    ne = List_Num_Entries(adj->fedges);

    if (ne == 0)
      MSTK_Report("MF_Replace_Edge_i","No initial set of edges for face",
		  MSTK_ERROR);

    oldedges = (MEdge_ptr *) MSTK_malloc(nold*sizeof(MEdge_ptr));
    olddirs  = (int *) MSTK_malloc(nold*sizeof(int));
    for (j = 0; j < nold; j++) {
      k = (i+j)%ne;
      oldedges[j] = List_Entry(adj->fedges,k);
      olddirs[j] = (adj->edirs>>k) & 1UL;
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
      MSTK_Report("MF_Replace_Edge_i","Mismatched set of edges",MSTK_ERROR);


    newedges = (MEdge_ptr *) MSTK_malloc(nnu*sizeof(MEdge_ptr));
    newdirs  = (int *) MSTK_malloc(nnu*sizeof(MEdge_ptr));
    lastv = vold_0;
    for (j = 0; j < nnu; j++) {
      k = rev ? nnu-1-j : j;
      newedges[j] = nuedges[k];      
      newdirs[j] = (ME_Vertex(newedges[j],0) == lastv) ? 1 : 0;
      lastv = ME_Vertex(newedges[j],newdirs[j]);	

#ifdef DEBUG
      if (MF_Mesh(f) != ME_Mesh(nuedges[k]))
	  MSTK_Report("MF_Replace_Edges_i_F1F3",
		      "Face and edge belong to different meshes",MSTK_FATAL);
#endif
    }

    /* Number of common edges */

    ncom = nold < nnu ? nold : nnu;


    /* Replace old edges with new ones (up to the common number of edges) */

    for (j = 0; j < ncom; j++) {
      k = (i+j)%ne; 
      List_Replacei(adj->fedges,k,newedges[j]);		      
      adj->edirs = (adj->edirs & ~(1UL<<k)); /* set bit k to 0 */
      adj->edirs = (adj->edirs | (newdirs[j]<<k)); /* set to dir */
    }



    if (nold < nnu) {
      /* Insert the rest of the nuedges at the right place */

      for (j = ncom; j < nnu; j++) {
	List_Inserti(adj->fedges,newedges[j],i+j);

      
	/* Move the direction bits after the (i+ncom)'th bit to the
	 left by (nnu-ncom) spaces and insert the direction bits for
	 the rest of the inserted edges */

	
	for (k = ne-1; k >= i+j; k--) {

	  /* set bit k+1 to 0 */
	  adj->edirs = adj->edirs & ~(1UL<<(k+1));

	  /* get bit k */
	  dir = (adj->edirs>>k) & 1UL;

	  /* move bit k to (k+1)'th position */
	  adj->edirs = adj->edirs | (dir<<(k+1));	
	}

	/* set bit j according to input */
	adj->edirs = adj->edirs & ~(1UL<<(i+j)); /* Set bit to 0 */
	adj->edirs = adj->edirs | (newdirs[j]<<(i+j)); /* set to newdir */

	ne++;
      }

    }
    else {
      /* Delete the remaining old edges */

      for (j = ncom; j < nold; j++)
	List_Rem(adj->fedges,oldedges[j]);

      /* Delete (nold-ncom) direction bits after the (i+ncom)'th bit and
	 move the remaining bits to the right */
      
      for (j = 0; j < nold-ncom; j++) {
	for (k = (i+ncom+j)%ne; k < ne-1; k++) {
	  /* set bit k to 0 */
	  adj->edirs = adj->edirs & ~(1UL<<k);

	  /* get bit k+1 */
	  dir = (adj->edirs>>(k+1)) & 1UL;

	  /* move bit k+1 to k'th position */
	  adj->edirs = adj->edirs | (dir<<k);
	}

	ne--;
      }
    }
      
    /* Update upward adjacencies */
    for (j = 0; j < nold; j++) {
      ME_Rem_Face(oldedges[j],f);
    }
    for (j = 0; j < nnu; j++) {
      ME_Add_Face(newedges[j],f);
    }

    MSTK_free(oldedges);
    MSTK_free(olddirs);
    MSTK_free(newedges);
    MSTK_free(newdirs);
  }


  /* Replace a set of edges in a face with another set of
     edges. Assumption is that the edges in each set (old and new) are
     consecutive */

  void MF_Replace_Edges_F1F3(MFace_ptr f, int nold, MEdge_ptr *oldedges, int nnu, MEdge_ptr *nuedges) {
    int i, j, ne, found = 0, imin1, istart, *eindex, nxtidx, ipos1;
    MFace_Adj_F1F3 *adj;
    MEdge_ptr fedge;

    adj = (MFace_Adj_F1F3 *) f->adj;
    ne = List_Num_Entries(adj->fedges);

    if (ne == 0)
      MSTK_Report("MF_Replace_Edge_F1","No initial set of edges for face",
		  MSTK_ERROR);

    /* Order the old edges in the direction of the face */

    eindex = (int *) MSTK_malloc(nold*sizeof(int));

    for (i = 0; i < nold; i++) {           
      for (j = 0, found = 0; j < ne; j++)
	if ((fedge = List_Entry(adj->fedges,j)) == oldedges[i]) {
	  found = 1;
	  eindex[i] = j;
	  break;
	}
      if (!found) {
	MSTK_Report("MF_Replace_Edge","Edge not found in face",MSTK_ERROR);
	MSTK_free(eindex);
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
      
	  
    MF_Replace_Edges_i_F1F3(f, nold, istart, nnu, nuedges);

    MSTK_free(eindex);
  }


  int MF_Rev_EdgeDir_i_F1F3(MFace_ptr f, int i) {
    MFace_Adj_F1F3 *adj;

    adj = (MFace_Adj_F1F3 *) f->adj;
    adj->edirs ^= (1UL << i);

    return 1;
  }


  int MF_Rev_EdgeDir_F1F3(MFace_ptr f, MEdge_ptr e) {
    int n, i;
    MFace_Adj_F1F3 *adj;

    adj = (MFace_Adj_F1F3 *) f->adj;

    n = List_Num_Entries(adj->fedges);
    for (i = 0; i < n; i++) {
      if (List_Entry(adj->fedges,i) == e)
        MF_Rev_EdgeDir_i_F1F3(f, i);
    }

    return 1;
  }


  int MF_Num_Edges_F1F3(MFace_ptr f) {
    List_ptr fedges = ((MFace_Adj_F1F3 *) f->adj)->fedges;
    return List_Num_Entries(fedges);
  }	

  List_ptr MF_Edges_F1F3(MFace_ptr f, int dir, MVertex_ptr v0) {
    int i, k=0, n, fnd, edir;
    List_ptr fedges;
    MEdge_ptr e;
    MFace_Adj_F1F3 *adj;

    adj = (MFace_Adj_F1F3 *) f->adj;


    if (v0 == NULL) {
      if (dir)
	return List_Copy(adj->fedges);
      else {
	n = List_Num_Entries(adj->fedges);
	fedges = List_New(n);
	for (i = n-1; i >= 0; i--)
	  List_Add(fedges,List_Entry(adj->fedges,i));
	return fedges;
      }
    }
    else {
      n = List_Num_Entries(adj->fedges);
      fnd = 0;
      for (i = 0; i < n; i++) {
	e = List_Entry(adj->fedges,i);
	edir = ((adj->edirs)>>i) & 1UL;
	if (ME_Vertex(e,edir^dir) == v0) {
	  fnd = 1;
	  k = i;
	  break;
	}
      }

      if (!fnd)
	MSTK_Report("MF_Edges_F1","Cannot find vertex in face!!",MSTK_FATAL);
	
      fedges = List_New(n);
      for (i = 0; i < n; i++) {
	e = dir ? List_Entry(adj->fedges,(k+i)%n) :
	  List_Entry(adj->fedges,(k+n-i)%n);
	List_Add(fedges,e);
      }	
    }

    return fedges;
  }

  void MF_EdgeIDs_F1F3(MFace_ptr f, int dir, int startvertid, int *nfe,
                       int *fedgeids) {
    int i, k=0, fnd, edir;
    List_ptr fedges;
    MEdge_ptr e;
    MFace_Adj_F1F3 *adj;

    adj = (MFace_Adj_F1F3 *) f->adj;

    *nfe = List_Num_Entries(adj->fedges);

    if (startvertid == 0) {
      if (dir) {
        for (i = 0; i < *nfe; i++)
          fedgeids[i] = MEnt_ID(List_Entry(adj->fedges,i));
      }
      else {
	for (i = 0; i < *nfe; i++)
	  fedgeids[i] = MEnt_ID(List_Entry(adj->fedges,*nfe-i-1));
      }
    }
    else {
      fnd = 0;
      for (i = 0; i < *nfe; i++) {
	e = List_Entry(adj->fedges,i);
	edir = ((adj->edirs)>>i) & 1UL;
	if (ME_VertexID(e,edir^dir) == startvertid) {
	  fnd = 1;
	  k = i;
	  break;
	}
      }

      if (!fnd)
	MSTK_Report("MF_Edges_F1","Cannot find vertex in face!!",MSTK_FATAL);
	
      for (i = 0; i < *nfe; i++) {
	e = dir ? List_Entry(adj->fedges,(k+i)%(*nfe)) :
	  List_Entry(adj->fedges,(k+(*nfe)-i)%(*nfe));
	fedgeids[i] = MEnt_ID(e);
      }	
    }
  }


  int MF_EdgeDir_F1F3(MFace_ptr f, MEdge_ptr e) {
    int n, i;
    MFace_Adj_F1F3 *adj;

    adj = (MFace_Adj_F1F3 *) f->adj;
    
    n = List_Num_Entries(adj->fedges);
    for (i = 0; i < n; i++) {
      if (List_Entry(adj->fedges,i) == e)
	return ((adj->edirs)>>i) & 1UL;
    }

    MSTK_Report("MF_Edges_F1F3","Cannot find edge in face!!",MSTK_FATAL);

    return 0;
  }

  int MF_EdgeDir_i_F1F3(MFace_ptr f, int i) {
    MFace_Adj_F1F3 *adj;

#ifdef DEBUG
    if (i > 31)
      MSTK_Report("MF_Set_Edges_F1F3","Currently, only 32 edges supported per face",MSTK_ERROR);
#endif

    adj = (MFace_Adj_F1F3 *) f->adj;

    return ((adj->edirs)>>i) & 1UL;
  }
			
  int MF_UsesEdge_F1F3(MFace_ptr f, MEdge_ptr e) {
    MFace_Adj_F1F3 *adj;

    adj = (MFace_Adj_F1F3 *) f->adj;

    return List_Contains(adj->fedges,e);
  }


#ifdef __cplusplus
}
#endif
