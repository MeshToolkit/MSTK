#define _H_MVertex_Private

#include <stdlib.h>

#include "MVertex.h"
#include "MVertex_jmp.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif


  void MV_Set_RepType_F4(MVertex_ptr v) {
    MVertex_Adj_F1F4 *adj;

    adj = v->adj = (MVertex_Adj_F1F4 *) malloc(sizeof(MVertex_Adj_F1F4));
    adj->vedges = List_New(5);
  }

  void MV_Delete_F4(MVertex_ptr v, int keep) {
    MVertex_Adj_F1F4 *adj;

    if (!keep) {
      adj = (MVertex_Adj_F1F4 *) v->adj;
      if (adj) {
	if (adj->vedges)
	  List_Delete(adj->vedges);
	free(adj);
      }
    }
  }

  void MV_Restore_F4(MVertex_ptr v) {
    MEnt_Set_Dim((MEntity_ptr) v,MVERTEX);
  }

  void MV_Destroy_For_MESH_Delete_F4(MVertex_ptr v) {
    MVertex_Adj_F1F4 *adj;

    adj = (MVertex_Adj_F1F4 *) v->adj;
    if (adj) {
      if (adj->vedges)
	List_Delete(adj->vedges);
      free(adj);
    }
  }

  int MV_Num_AdjVertices_F4(MVertex_ptr v) {
    List_ptr vedges = ((MVertex_Adj_F1F4 *) v->adj)->vedges;
    return List_Num_Entries(vedges);
  }

  int MV_Num_Edges_F4(MVertex_ptr v) {
    List_ptr vedges = ((MVertex_Adj_F1F4 *) v->adj)->vedges;
    return List_Num_Entries(vedges);
  }

  int MV_Num_Faces_F4(MVertex_ptr v) {
    List_ptr vfaces;
    int nf;

#ifdef DEBUG
    MSTK_Report("MV_Num_Faces",
		"Ineficient to call this routine with this representation",
		MSTK_MESG);
#endif
    
    vfaces = MV_Faces_F4(v);
    nf = List_Num_Entries(vfaces);
    List_Delete(vfaces);
    return nf;
  }

  int MV_Num_Regions_F4(MVertex_ptr v) {
    List_ptr vregions;
    int nr;

#ifdef DEBUG
    MSTK_Report("MV_Num_Regions",
		"Inefficient to call this routine with this representation",
		MSTK_MESG);
#endif
	
    vregions = MV_Regions_F4(v);
    nr = List_Num_Entries(vregions);
    List_Delete(vregions);
    return nr;
  }

  List_ptr MV_AdjVertices_F4(MVertex_ptr v) {
    MVertex_Adj_F1F4 *adj;
    List_ptr vedges, adjv;
    int ne, i;
    MEdge_ptr vedge;
    MVertex_ptr ov;

    adj = (MVertex_Adj_F1F4 *) v->adj;
    vedges = adj->vedges;
    if (vedges == 0)
      return 0;

    ne = List_Num_Entries(vedges);
    adjv = List_New(ne);
    for (i = 0; i < ne; i++) {
      vedge = List_Entry(vedges,i);
      ov = ME_OppVertex(vedge,v);
      List_Add(adjv,ov);
    }

    return adjv;
  }


  void MV_AdjVertexIDs_F4(MVertex_ptr v, int *nvadj, int *adjvids) {
    MVertex_Adj_F1F4 *adj;
    List_ptr vedges, adjv;
    int ne, i;
    MEdge_ptr vedge;
    MVertex_ptr ov;

    adj = (MVertex_Adj_F1F4 *) v->adj;

    *nvadj = 0;
    vedges = adj->vedges;
    if (vedges) {
      *nvadj = List_Num_Entries(vedges);
      for (i = 0; i < *nvadj; i++) {
        vedge = List_Entry(vedges,i);
        ov = ME_OppVertex(vedge,v);
        adjvids[i] = MEnt_ID((MEntity_ptr)ov);
      }      
    }
  }  
    

  List_ptr MV_Edges_F4(MVertex_ptr v) {
    MVertex_Adj_F1F4 *adj;
    List_ptr vedges;

    adj = (MVertex_Adj_F1F4 *) v->adj;
    vedges = List_Copy(adj->vedges);
    return vedges;
  }

  void MV_EdgeIDs_F4(MVertex_ptr v, int *nve, int *vedgeids) {
    MVertex_Adj_F1F4 *adj;
    adj = (MVertex_Adj_F1F4 *) v->adj;

    if (adj->vedges) {
      int i;
      *nve = List_Num_Entries(adj->vedges);
      for (i = 0; i < *nve; i++) 
        vedgeids[i] = MEnt_ID(List_Entry(adj->vedges,i));
    }
    else
      *nve = 0;
  }


  List_ptr MV_Faces_F4(MVertex_ptr v) {
    MVertex_Adj_F1F4 *adj;
    int i, j, k, ne, nf, nr, n;
    List_ptr vedges, eregions, rfaces, efaces, vfaces;
    MEdge_ptr edge;
    MFace_ptr face;
    MRegion_ptr region;

    adj = (MVertex_Adj_F1F4 *) v->adj;
    vedges = adj->vedges;
    ne = List_Num_Entries(vedges);

    n = 0;
    vfaces = List_New(ne);

    for (i = 0; i < ne; i++) {
      edge = List_Entry(vedges,i);
      eregions = ME_Regions(edge);
      if (eregions) {
	nr = List_Num_Entries(eregions);
	
	for (j = 0; j < nr; j++) {
	  region = List_Entry(eregions,j);
	  
	  rfaces = MR_Faces(region);
	  nf = List_Num_Entries(rfaces);
	  
	  for (k = 0; k < nf; k++) {
	    face = List_Entry(rfaces,k);
            int fmarked;
            fmarked = List_Contains(vfaces,face);
	    if (!fmarked) {
	      if (MF_UsesEntity(face,(MEntity_ptr) v,0)) {		
		List_Add(vfaces,face);
		n++;
	      }
	    }
	  }
	  List_Delete(rfaces);
	}
	List_Delete(eregions);
      }
      else {
	/* perhaps the edge has boundary faces (not connected to regions) */
	efaces = ME_Faces(edge);
	if (efaces) {
	  nf = List_Num_Entries(efaces);
	  
	  for (k = 0; k < nf; k++) {
	    face = List_Entry(efaces,k);
            int fmarked;
            fmarked = List_Contains(vfaces,face);
            if (!fmarked) {
	      if (MF_UsesEntity(face,(MEntity_ptr) v,0)) {
		List_Add(vfaces,face);
		n++;
	      }
	    }
	  }
	  List_Delete(efaces);
	}
      }
    }
    if (n > 0)
      return vfaces;
    else {
      List_Delete(vfaces);
      return 0;
    }
  }

  void MV_FaceIDs_F4(MVertex_ptr v, int *nvf, int *vfaceids) {
    List_ptr vfaces = MV_Faces_F4(v);
    if (vfaces) {
      int i;
      *nvf = List_Num_Entries(vfaces);
      for (i = 0; i < *nvf; i++)
        vfaceids[i] = MEnt_ID(List_Entry(vfaces,i));
      List_Delete(vfaces);
    }
    else
      *nvf = 0;
  }


  List_ptr MV_Regions_F4(MVertex_ptr v) {
    MVertex_Adj_F1F4 *adj;
    int i, j, ne, nr, n;
    List_ptr vedges, eregions, vregions;
    MEdge_ptr edge;
    MRegion_ptr region;

    adj = (MVertex_Adj_F1F4 *) v->adj;
    vedges = adj->vedges;
    ne = List_Num_Entries(vedges);

    n = 0;
    vregions = List_New(ne);

    for (i = 0; i < ne; i++) {
      edge = List_Entry(vedges,i);

      eregions = ME_Regions(edge);
      if (eregions) {
	nr = List_Num_Entries(eregions);
	
	for (j = 0; j < nr; j++) {
	  region = List_Entry(eregions,j);

          int rmarked;
          rmarked = List_Contains(vregions,region);
          if (!rmarked) {
	    List_Add(vregions,region);
	    n++;
	  }
	}
	
	List_Delete(eregions);
      }
    }

    if (n > 0)
      return vregions;
    else {
      List_Delete(vregions);
      return 0;
    }
  }

  void MV_RegionIDs_F4(MVertex_ptr v, int *nvr, int *vregionids) {
    List_ptr vregions = MV_Regions_F4(v);
    if (vregions) {
      int i;
      *nvr = List_Num_Entries(vregions);
      for (i = 0; i < *nvr; i++)
        vregionids[i] = MEnt_ID(List_Entry(vregions,i));
      List_Delete(vregions);
    }
    else
      *nvr = 0;
  }

  void MV_Add_AdjVertex_F4(MVertex_ptr v, MVertex_ptr av) {
#ifdef DEBUG
    MSTK_Report("MV_Add_AdjVertex","Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void MV_Rem_AdjVertex_F4(MVertex_ptr v, MVertex_ptr av) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_AdjVertex","Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void MV_Add_Edge_F4(MVertex_ptr v, MEdge_ptr e) {
    MVertex_Adj_F1F4 *adj;

    adj = (MVertex_Adj_F1F4 *) v->adj;

    if (adj->vedges == NULL)
      adj->vedges = List_New(10);
    List_ChknAdd(adj->vedges,e);
  }

  void MV_Rem_Edge_F4(MVertex_ptr v, MEdge_ptr e) {
    MVertex_Adj_F1F4 *adj;

    adj = (MVertex_Adj_F1F4 *) v->adj;

    if (adj->vedges == NULL)
      return;
    List_Rem(adj->vedges,e);
  }

  void MV_Add_Face_F4(MVertex_ptr v, MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("MV_Add_Face","Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void MV_Rem_Face_F4(MVertex_ptr v, MFace_ptr f) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_Face","Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void MV_Add_Region_F4(MVertex_ptr v, MRegion_ptr f) {
#ifdef DEBUG
    MSTK_Report("MV_Add_Region","Function call not suitable for this representation",MSTK_WARN);
#endif
  }

  void MV_Rem_Region_F4(MVertex_ptr v, MRegion_ptr f) {
#ifdef DEBUG
    MSTK_Report("MV_Rem_Region","Function call not suitable for this representation",MSTK_WARN);
#endif
  }

#ifdef __cplusplus
}
#endif
