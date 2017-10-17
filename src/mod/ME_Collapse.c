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
/*                                                                      */
/* 'deleted_entities' contains a list of all entities deleted due to    */
/* this operation. 'merged_entity_pairs' contains info about merged     */
/* entities by storing the deleted and kept entity in each merged pair  */
/* (in that order)                                                      */
/*                                                                      */
/* NOTE: THIS ROUTINE ONLY MARKS THE DELETED ENTITY AS MDELETED BUT     */
/* DOES NOT ACTUALLY DELETE THEM. IF YOU WANT THEM REALLY DELETED, YOU  */
/* CALL THE MEnt_Delete(entity, 0) ROUTINE ON THEM                      */

MVertex_ptr ME_Collapse(MEdge_ptr e, MVertex_ptr vkeep_in, int topoflag,
                        List_ptr *deleted_entities,
                        List_ptr *merged_entity_pairs) {
  MVertex_ptr vdel, vkeep, ev00, ev01, ev10, ev11, vert;
  MEdge_ptr edge, edge2, oldedges[3], nuedges[2];
  MFace_ptr face, face2, rface1, rface2;
  MRegion_ptr reg, reg2;
  List_ptr vedges, efaces, eregs, fedges, rfaces, fverts1, fverts2, vfaces;
  int idx1, idx2, idx3, dir, status, nfe, nrf, nrv, allfound, degenerate;
  int i, j, nfe2, nfv1, nfv2;
  MEntity_ptr delent;

  status = 1;

  if (vkeep_in == NULL) {
    vdel = ME_Vertex(e,0);
    vkeep = ME_Vertex(e,1);
  }
  else {
    vkeep = vkeep_in;
    vdel = ME_OppVertex(e,vkeep);
  }

  int dimkeep, dimdel;
    
  dimkeep = MV_GEntDim(vkeep); /* Model entity dim of vertex to keep */
  dimdel = MV_GEntDim(vdel);   /* Model entity dim  of vertex to delete */
  
  if (topoflag == 1) {
    if (dimkeep == dimdel) {
      
      if (MV_GEntID(vkeep) != MV_GEntID(vdel))
        status = 0;                 /* cannot allow since it will cause 
                                       a dimensional reduction in mesh */
      else {
        /* No restriction because of topology - pick the one with the
         * lower global ID to delete - important for parallel
         * operation consistency */
        int vdel_gid = MV_GlobalID(vdel);
        int vkeep_gid = MV_GlobalID(vkeep);
        if (vdel_gid < vkeep_gid) {
          MVertex_ptr vtemp = vdel;
          vdel = vkeep;
          vkeep = vtemp;
        }
      }
    }
    else if (dimdel < dimkeep) {
      
      if (vkeep_in == NULL) {
        /* If no preference was indicated on which vertex to retain,
	   we can collapse in the other direction */
        MVertex_ptr vtemp = vdel;
	vdel = vkeep;
	vkeep = vtemp;
      }
      else
        status = 0; /* can't reverse order or vertices and boundary of
                       mesh will get messed up if we go through as is */
    }
  }
  else if (vkeep_in == NULL) { 

    /* If no preference was indicated for the kept vertex and
       topological conformity with the underlying geometric model was
       not requested, we prefer to keep an external boundary vertex
       over an interior vertex or interior boundary vertex. This is
       because it is more likely that the external boundary vertex
       would have a boundary condition applied to it. If a preference
       was indicated, we just have to respect that. */

    int vdel_external = 0;

    /* Check if any edges connected to vdel have only one connected face */

    vedges = MV_Edges(vdel);    

    idx1 = 0;
    while ((edge = (MEdge_ptr) List_Next_Entry(vedges,&idx1))) {
      List_ptr efaces = ME_Faces(edge);
      int nef = List_Num_Entries(efaces);
      List_Delete(efaces);
      if (nef < 2) {
        vdel_external = 1;
        break;
      }
    }                                          
    
    List_Delete(vedges);

    /* check if any face connected to vdel has only one region
       connected to it */

    if (!vdel_external) {
      vfaces = MV_Faces(vdel);

      idx1 = 0;
      while ((face = (MFace_ptr) List_Next_Entry(vfaces,&idx1))) {
        List_ptr fregs = MF_Regions(face);
        int nfr = fregs ? List_Num_Entries(fregs) : 0;
        if (fregs) List_Delete(fregs);

        if (nfr == 1) {
          vdel_external = 0;
          break;
        }
      }

      List_Delete(vfaces);
    }

    if (vdel_external) {
      /* swap the vertices in the hope that vkeep is not also on an
         external boundary. Since we have to go through with the
         collapse anyway, there is no use of doing a detailed check
         for whether vkeep is also on an external boundary */

      MVertex_ptr vtemp = vdel;
      vdel = vkeep;
      vkeep = vtemp;
    }
    
  }


  if (status == 0)
    return NULL;   /* Cannot collapse due to constraints of topological
		   conformity with geometric model */

  int num_merged_pairs = 0;
  *merged_entity_pairs = List_New(10);

  List_Add(*merged_entity_pairs, vdel);
  List_Add(*merged_entity_pairs, vkeep);

  int num_deleted = 0;
  *deleted_entities = List_New(10);

  /* Need to collect this in advance because the info gets messed up later */

  efaces = ME_Faces(e);
  eregs = ME_Regions(e);

  /* Delete topologically degenerate tetrahedrons */
  
  if (eregs) {
    idx1 = 0;
    while ((reg = List_Next_Entry(eregs,&idx1))) {
      
      List_ptr rverts = MR_Vertices(reg);
      nrv = List_Num_Entries(rverts);
      
      if (nrv == 4) {
        List_Add(*deleted_entities,reg);
        MR_Delete(reg,1);    /* This is a tet that degenerated */
      }
      
      List_Delete(rverts);
    }
  }
  
  /* Delete topologically degenerate faces */
  
  if (efaces) {
    idx1 = 0;
    while ((face = List_Next_Entry(efaces,&idx1))) {
      
      fedges = MF_Edges(face,1,0);
      
      if (List_Num_Entries(fedges) == 3) {
        
        /* Disconnect the regions from the face before deleting */
        
        List_ptr fregs = MF_Regions(face);
        if (fregs) {
          idx2 = 0;
          while ((reg = List_Next_Entry(fregs,&idx2)))
            MR_Rem_Face(reg,face);
          
          List_Delete(fregs);
        }
        
        List_Add(*deleted_entities,face);
        MF_Delete(face,1);
      }
      
      List_Delete(fedges);
    }
  }
  
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

  
  if (efaces) {
    idx1 = 0;
    while ((face = List_Next_Entry(efaces,&idx1))) {
      
      if (MEnt_Dim(face) == MDELETED)
        continue;
      
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
    List_Delete(efaces);
  }
  

  /* Delete topologically degenerate regions */
  /* Defined as two faces of the regions having the same vertices */

  if (eregs) {
    idx1 = 0;
    while ((reg = List_Next_Entry(eregs,&idx1))) {

      if (MEnt_Dim(reg) == MDELETED)
        continue;
      
      rfaces = MR_Faces(reg);
      nrf = List_Num_Entries(rfaces);

	degenerate = 0;

	for (i = 0; i < nrf; i++) {

	  rface1 = List_Entry(rfaces,i);
          if (MEnt_Dim(rface1) == MDELETED)
            continue;
          
	  fverts1 = MF_Vertices(rface1,1,0);
	  nfv1 = List_Num_Entries(fverts1);
		
	  for (j = i+1; j < nrf; j++) {
	  
	    rface2 = List_Entry(rfaces,j);
            if (MEnt_Dim(rface2) == MDELETED)
              continue;
            
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

	if (degenerate) {
          List_Add(*deleted_entities,reg);
          MR_Delete(reg,1);
        }

      List_Delete(rfaces);
    } /* while ((reg = ...)) */
  }

  /* Now merge edges which have the same end vertices */
  /* Prefer to preserve edges on external boundaries over internal edges */

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

        int edge_on_boundary, edge2_on_boundary;
        int edim = 4;
	
        edge_on_boundary = 0;
        edim = ME_GEntDim(edge);
        if (edim == 1 || edim == 2)
          edge_on_boundary = 1;
        else if (edim == 4 && ME_PType(edge) != PGHOST) {
          /* Don't know classification of edge so check if its a
             boundary edge based on topology. However, do this check
             only if the edge is not a ghost, because edges on the
             outer boundary of a ghost layer can mimic an external
             edge */
          
          efaces = ME_Faces(edge);
          int nef = List_Num_Entries(efaces);
          if (nef == 1) {
            edge_on_boundary = 1;
          } else {
            idx3 = 0;
            while ((face = List_Next_Entry(efaces,&idx3))) {
              List_ptr fregs = MF_Regions(face);
              if (fregs) {
                int nfr = List_Num_Entries(fregs);
                if (nfr == 1 || 
                    (nfr == 2 && 
                     MR_GEntID(List_Entry(fregs,0)) != MR_GEntID(List_Entry(fregs,1))))
                  edge_on_boundary = 1;
                List_Delete(fregs);
              }
              if (edge_on_boundary) break;
            }
          }
          List_Delete(efaces);
        }
    
        edge2_on_boundary = 0;
        edim = ME_GEntDim(edge2);
        if (edim == 1 || edim == 2)
          edge_on_boundary = 1;
        else if (edim == 4 && ME_PType(edge2) != PGHOST) {
          /* Don't know classification of edge so check if its a 
             boundary edge based on topology. However, do this check 
             only if the edge is not a ghost, because edges on the outer 
             boundary of a ghost layer can mimic an external edge */
          
          efaces = ME_Faces(edge2);
          int nef = List_Num_Entries(efaces);
          if (nef == 1) {
            edge2_on_boundary = 1;
          } else {
            idx3 = 0;
            while ((face = List_Next_Entry(efaces,&idx3))) {
              List_ptr fregs = MF_Regions(face);
              if (fregs) {
                int nfr = List_Num_Entries(fregs);
                if (nfr == 1 || 
                    (nfr == 2 && 
                     MR_GEntID(List_Entry(fregs,0)) != MR_GEntID(List_Entry(fregs,1))))
                  edge_on_boundary = 1;
                List_Delete(fregs);
              }
              if (edge2_on_boundary) break;
            }
          }
          List_Delete(efaces);
        }
    
        if (edge_on_boundary == edge2_on_boundary) {
          /* if both edges are boundary edges or both are interior, then we
             don't have to worry about topological consistency of the
             mesh with the model which way we merge. In that case,
             merge 'edge2' onto 'edge' only if 'edge' has the lower
             with the lower global ID - this preserves parallel
             consistency; if not, continue on and 'edge' will get
             merged onto 'edge2' when we are processing 'edge2' in the
             outer loop */

          int egid = ME_GlobalID(edge);
          int egid2 = ME_GlobalID(edge2);
          if (egid < egid2) {
            MEs_Merge(edge,edge2,topoflag);	
            List_Rem(vedges,edge2);	
            List_Add(*deleted_entities,edge2);
            List_Add(*merged_entity_pairs, edge2);
            List_Add(*merged_entity_pairs, edge);
            break;
          }
        } else if (edge_on_boundary) {
          /* 'edge' is external while 'edge2' is not; merge 'edge2'
           * onto 'edge' regardless of global IDs because boundary
           * entities may be in boundary sets for boundary
           * conditions. Any other processor encountering these two
           * edges will make the same decision */

          MEs_Merge(edge,edge2,topoflag);	
          List_Rem(vedges,edge2);	
          List_Add(*deleted_entities,edge2);
          List_Add(*merged_entity_pairs, edge2);
          List_Add(*merged_entity_pairs, edge);
          break;
        }
      }
    }
  }
  List_Delete(vedges);


  /* Merge faces with the same set of edges */
      
  vfaces = MV_Faces(vkeep);

  if (vfaces) {
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

          /* Proceed with merge (which will delete face2) only if
             face2 is not an external face or both face and face2 are
             external - HOWEVER, SINCE WE DELETED THE DEGENERATE
             REGION IN BETWEEN THE TWO FACES, AN EXTERNAL FACE WILL
             HAVE 0 REGIONS CONNECTED TO IT */

          List_ptr fregs = MF_Regions(face);
          int face_on_boundary = ((MF_PType(face) != PGHOST) && !fregs) ? 1 : 0;
          if (fregs) List_Delete(fregs);

          List_ptr fregs2 = MF_Regions(face2);
          int face2_on_boundary = ((MF_PType(face2) != PGHOST) && !fregs2) ? 1 : 0;
          if (fregs2) List_Delete(fregs2);

          if (face_on_boundary == face2_on_boundary) {
          /* if both faces are boundary faces or both are interior, then we
             don't have to worry about topological consistency of the
             mesh with the model which way we merge. In that case,
             merge 'face2' onto 'face' only if 'face' has the lower
             with the lower global ID - this preserves parallel
             consistency; if not, continue on and 'face' will get
             merged onto 'face2' when we are processing 'face2' in the
             outer loop */

            int fgid = MF_GlobalID(face);
            int fgid2 = MF_GlobalID(face2);
            if (fgid < fgid2) {
              MFs_Merge(face,face2,topoflag);	
              List_Rem(vfaces,face2);
              List_Add(*deleted_entities,face2);
              List_Add(*merged_entity_pairs, face2);
              List_Add(*merged_entity_pairs, face);
              break;
            }
          } else if (face_on_boundary) {
            MFs_Merge(face,face2,topoflag);	
            List_Rem(vfaces,face2);
            List_Add(*deleted_entities,face2);
            List_Add(*merged_entity_pairs, face2);
            List_Add(*merged_entity_pairs, face);
            break;
          }
        }
	
      } /* while (face2 = List_Next_Entry(vfaces,... */

      List_Delete(fedges);

    } /* while (face = List_Next_Entry(vfaces,... */
    List_Delete(vfaces);
  }

  /* Now actually delete the collapse edge and the to-be-merged vertex */

  ME_Delete(e,1);
  List_Add(*deleted_entities,e);

  MV_Delete(vdel,1);
  List_Add(*deleted_entities,vdel);

  if (eregs) {
    idx1 = 0;
    while ((reg = List_Next_Entry(eregs,&idx1))) {
      if (MEnt_Dim(reg) == MDELETED)
          continue;
      MR_Update_ElementType(reg);
    }
    List_Delete(eregs);
  }

  return vkeep;
}
  
