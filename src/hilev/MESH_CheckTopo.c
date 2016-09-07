#include <stdio.h>
#include <stdlib.h>

#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif

  /* Routine to check that the topological structure of the mesh is valid */
  /* Rao Garimella - 12/16/2010                                           */

  int MESH_CheckTopo(Mesh_ptr mesh) {
    int valid = 1;
    char mesg[256], funcname[32] = "MESH_CheckTopo";
    int idx1, idx2, idx3, idx4;
    MVertex_ptr mv;
    MEdge_ptr me, ve, fe, re;
    MFace_ptr mf, vf, ef, rf;
    MRegion_ptr mr, vr, er, fr;
    int found, done;
    int dir;
    int i, j, k;
    int nfe;
    int vid, eid, fid, rid;
    int gvid, geid, gfid, grid;
    int gvdim, gedim, gfdim, grdim;
    int maxiter = 1000;
    List_ptr vedges, vfaces, vregions;
    List_ptr efaces;
    List_ptr fverts, fedges, fregs, fregs1;
    List_ptr rverts, redges, rfaces;



    /*****************************************************************/
    /* Vertices                                                      */
    /*****************************************************************/    

    
    /* Check that edges connected to vertices reference the vertices */
    /* Also check that the classification of the vertex is consistent
       with respect to the edge */

    int first_unknown_classfn = 1;
    idx1 = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx1))) {

#ifdef MSTK_HAVE_MPI
      if (MV_PType(mv) == PGHOST) continue;
#endif      

      vid = MV_ID(mv);
      gvdim = MV_GEntDim(mv);
      gvid  = MV_GEntID(mv);

      if (gvdim == 4 && first_unknown_classfn) {
        sprintf(mesg, "Vertex %-d - classification unknown\n", vid);
        MSTK_Report(funcname, mesg, MSTK_WARN);
        first_unknown_classfn = 0;
      }
      vedges = MV_Edges(mv);
      if (!vedges) {
        sprintf(mesg,"Vertex %-d does not have any connected edges\n",vid);
        MSTK_Report(funcname,mesg,MSTK_WARN);
        continue;
      }

      idx2 = 0;
      while ((ve = List_Next_Entry(vedges,&idx2))) {
	
	eid = ME_ID(ve);
	     
	if (ME_Vertex(ve,0) != mv && ME_Vertex(ve,1) != mv) {
	  sprintf(mesg,"Vertex %-d connected to edge %-d but edge does not use vertex",vid,eid);
	  MSTK_Report(funcname,mesg,MSTK_ERROR);
	  valid = 0;
	}

      }


      if (gvdim == 1) {

	/* If vertex is classified on a model edge, then it should be
	   connected to two and only two edges that are classified on
	   the same model edge */

	int ne = 0;
	idx2 = 0;
	while ((ve = List_Next_Entry(vedges,&idx2))) {	  
	  gedim = ME_GEntDim(ve);
	  geid  = ME_GEntID(ve);	  
	  if (gedim == 1 && geid == gvid) ne++;
	}
	  
	if (ne != 2) {
	  sprintf(mesg,"Vertex %-d classified on model edge %-d but it is not \n connected to two edges classified on this model edge",vid,gvid);
	  MSTK_Report(funcname,mesg,MSTK_WARN);
	}
      }

      List_Delete(vedges);



      if (gvdim == 2) {
	MEdge_ptr e0, ecur, enxt;
	MFace_ptr fcur;
        int flipped = 0;

	/* If vertex is classified on a model face, then we should be
	   able to find a ring of faces classified on that model
	   face */

	vfaces = MV_Faces(mv);

	found = 0;
	idx2 = 0;
	while ((vf = List_Next_Entry(vfaces,&idx2))) {
	  if (MF_GEntDim(vf) == 2) {
	    found = 1;
	    break;
	  }
	}
        List_Delete(vfaces);

	if (!found) {
	  sprintf(mesg,"Vertex %-d classified on model face %-d but could not \n find connected face classified on this model face",vid,gvid);
	  MSTK_Report(funcname,mesg,MSTK_WARN);
	  valid = 0;
	}

	fcur = vf;

	fedges = MF_Edges(fcur,1,mv);
	nfe = List_Num_Entries(fedges);
	e0 = List_Entry(fedges,0);
	ecur = e0;
	enxt = List_Entry(fedges,nfe-1);
	List_Delete(fedges);
	
	done = 0; i = 0;
	while (!done) {
	  ecur = enxt;
	  efaces = ME_Faces(ecur);
	  found = 0;
	  idx3 = 0;
	  while ((ef = List_Next_Entry(efaces,&idx3))) {
	    if (ef != fcur && MF_GEntDim(ef) == 2 && MF_GEntID(ef) == gvid) {
	      fcur = ef;
	      found = 1;
	      break;
	    }
	  }
	  List_Delete(efaces);

	  if (!found) {
	    sprintf(mesg,"Could not find next boundary face connected to vertex %-d",vid);
	    MSTK_Report(funcname,mesg,MSTK_WARN);
	    valid = 0;
	    break;
	  }
	   
	  fedges = MF_Edges(fcur,1,mv);
	  nfe = List_Num_Entries(fedges);

          if (List_Entry(fedges,0) == ecur)
            enxt = List_Entry(fedges,nfe-1);
          else if (List_Entry(fedges,nfe-1) == ecur) {
            enxt = List_Entry(fedges,0);
            flipped = 1;
          }
          else {
            sprintf(mesg,"Could not find next edge while traversing around vertex %-d on model face %-d",vid,gvid);
            MSTK_Report(funcname,mesg,MSTK_ERROR);
          }

	  List_Delete(fedges);

	  if (enxt == e0)
	    done = 1;

	  if (++i > maxiter)
	    break;
	}

	if (!done) {
	  sprintf(mesg,"Vertex %-d classified on model face %-d but could not  find ring of faces classified on this model face",vid,gvid);
	  MSTK_Report(funcname,mesg,MSTK_WARN);
	}

        if (done && flipped) {
          sprintf(mesg,"Inconsistent orientations of boundary faces around vertex %-d",vid);
          MSTK_Report(funcname,mesg,MSTK_WARN);
        }
      }


    } /* while ((mv = MESH_Next_Vertex(mesh,&idx1))) */




    /*****************************************************************/
    /* Edges                                                      */
    /*****************************************************************/    

    first_unknown_classfn = 1;
    idx1 = 0;
    while ((me = MESH_Next_Edge(mesh,&idx1))) {

#ifdef MSTK_HAVE_MPI
      if (ME_PType(me) == PGHOST) continue;
#endif

      eid = ME_ID(me);
      gedim = ME_GEntDim(me);
      geid = ME_GEntID(me);

      if (gedim == 4 && first_unknown_classfn) {
        sprintf(mesg, "Edge %-d - unknown classification", eid);
        MSTK_Report(funcname, mesg, MSTK_WARN);
        first_unknown_classfn = 0;
      }

      if (ME_Vertex(me,0) == ME_Vertex(me,1)) {
        sprintf(mesg,"Edge %-d has repeated vertices",eid);
        MSTK_Report(funcname,mesg,MSTK_ERROR);
      }

      for (i = 0; i < 2; i++) {
	MVertex_ptr ev = ME_Vertex(me,i);

	vid = MV_ID(ev);
	gvid = MV_GEntID(ev);
	gvdim = MV_GEntDim(ev);

        if (gvdim != 4 && gvdim != 4) {  /* vertex and edge classifn is known */
          if (gedim < gvdim) {
            sprintf(mesg,"Edge %-d classified on lower dimensional entity than  connected vertex %-d",eid,vid);
            MSTK_Report(funcname,mesg,MSTK_WARN);
            valid = 0;
          }
          else if (gedim == gvdim && geid != gvid) {
            sprintf(mesg,"Edge %-d and its vertex %-d classified on different entities of the same dimension",eid,vid);
            MSTK_Report(funcname,mesg,MSTK_WARN);
            valid = 0;
          }
        }

	vedges = MV_Edges(ev);
	if (!List_Contains(vedges,me)) {
	  sprintf(mesg,"Edge %-d sees vertex %-d but not vice versa",eid,vid);
	  MSTK_Report(funcname,mesg,MSTK_ERROR);
	  valid = 0;
	}
	List_Delete(vedges);

	if (gedim == 2) {
	  MFace_ptr ebf[2], fcur, fnxt;
	  MRegion_ptr rcur;
	  int nf, nfr;
	  List_ptr eregs;

	  /* Edge is classified on model face - it should be connected
	     to two and only two faces also classified on this model
	     face */

	  ebf[0] = ebf[1] = NULL;
	  nf = 0;
	  efaces = ME_Faces(me);
	  idx2 = 0;
	  while ((ef = List_Next_Entry(efaces,&idx2))) {
	    fid = MF_ID(ef);

	    if (MF_GEntDim(ef) == 2) {
	      nf++;
	      if (gedim == 2 && MF_GEntID(ef) != geid) {
		sprintf(mesg,"Face %-d connected to edge %-d classified on different model face",fid,eid);
		MSTK_Report(funcname,mesg,MSTK_WARN);
		valid = 0;
	      }

	      if (ebf[0] == NULL)
		ebf[0] = ef;
	      else
		ebf[1] = ef;
	    }
	  }
	  List_Delete(efaces);

	  if (nf != 2) {
	    sprintf(mesg,"Boundary edge %-d is not connected to exactly two\n  faces classified on the boundary",eid);
	    MSTK_Report(funcname,mesg,MSTK_ERROR);
	    valid = 0;
	  }


	  eregs = ME_Regions(me);
	  if (!eregs) 
	    continue;
	  else
	    List_Delete(eregs);

	  /* Can we go from f0 to f1 in one or two dirs? */

	  fcur = ebf[0];
          fnxt = NULL;
	  fregs = MF_Regions(fcur);
	  if (!fregs) {
	    fid = MF_ID(fcur);
	    sprintf(mesg,"Edge %-d connected to regions but face %-d is not",eid,fid);
	    MSTK_Report(funcname,mesg,MSTK_ERROR);
	    valid = 0;
	  }

	  nfr = List_Num_Entries(fregs);
	  for (i = 0; i < nfr; i++) {
	    rcur = List_Entry(fregs,i);
	  
	    rfaces = MR_Faces(rcur);
	    idx3 = 0;
	    found = 0;
	    while ((rf = List_Next_Entry(rfaces,&idx3))) {
	      if (rf != fcur && MF_UsesEntity(rf,me,1)) {
		found = 1;
		fnxt = rf;
		break;
	      }
	    }
	    List_Delete(rfaces);
	    
	    if (!found) {
	      rid = MR_ID(rcur);
	      sprintf(mesg,"Could not find second face in region %-d using edge %-d",rid,eid);
	    }
	    
	    
	    done = 0; j = 0;
	    while (!done) {
	      fcur = fnxt;
	      fid = MF_ID(fcur);
	      
	      if (fnxt == ebf[1]) {
		done = 1;
		break;
	      }
	      
	      fregs1 = MF_Regions(fcur);
	      idx3 = 0;
	      while ((fr = List_Next_Entry(fregs1,&idx3))) {
		if (fr != rcur) {
		  rcur = fr;
		  found = 1;
		  break;
		}
	      }
	      List_Delete(fregs1);

	      if (!found) {
		sprintf(mesg,"Could not find next region around edge %-d",eid);
		MSTK_Report(funcname,mesg,MSTK_ERROR);
		valid = 0;
		break;
	      }
	      
	      
	      rfaces = MR_Faces(rcur);
	      idx3 = 0;
	      found = 0;
	      while ((rf = List_Next_Entry(rfaces,&idx3))) {
		if (rf != fcur && MF_UsesEntity(rf,me,1)) {
		  found = 1;
		  fnxt = rf;
		  break;
		}
	      }
	      List_Delete(rfaces);
	      
	      if (!found) {
		rid = MR_ID(rcur);
		sprintf(mesg,"Could not find second face in region %-d using edge %-d",rid,eid);
	      }
	      
	      if (++j > maxiter)
		break;
	    } /* while (!done) */

	    if (!done) {
	      sprintf(mesg,"Could not traverse around edge %-d from face %-d to face %-d",eid,MF_ID(ebf[0]),MF_ID(ebf[1]));
	      MSTK_Report(funcname,mesg,MSTK_ERROR);
	      valid = 0;
	    }
	  } /* for (i = 0; i < nfr; i++) */
          List_Delete(fregs);

	} /* if (geid == 2) */

      } /* for (i = 0; i < 2; i++) */

    } /* while ((me = MESH_Next_Edge(mesh,&idx1))) */



    /*****************************************************************/
    /* Faces                                                      */
    /*****************************************************************/    

    first_unknown_classfn = 1;
    idx1 = 0;
    while ((mf = MESH_Next_Face(mesh,&idx1))) {

#ifdef MSTK_HAVE_MPI
      if (MF_PType(mf) == PGHOST) continue;
#endif

      fid = MF_ID(mf);
      gfid = MF_GEntID(mf);
      gfdim = MF_GEntDim(mf);

      if (gfdim == 4 && first_unknown_classfn) {
        sprintf(mesg, "Face %-d - unknown classification", fid);
        MSTK_Report(funcname, mesg, MSTK_WARN);
        first_unknown_classfn = 0;
      }

      fedges = MF_Edges(mf,1,0);

      if (List_Num_Entries(fedges) < 3) {
        sprintf(mesg,"Face %-d has less than 3 edges",fid);
        MSTK_Report(funcname,mesg,MSTK_ERROR);
      }

      idx2 = 0;
      while ((fe = List_Next_Entry(fedges,&idx2))) {
	eid = ME_ID(fe);
	geid = ME_GEntID(fe);
	gedim = ME_GEntDim(fe);

        if (gedim != 4 && gfdim != 4) {  /* Edge, Face classfn is known */
          if (gfdim < gedim) {
            sprintf(mesg,"Face %-d classified on lower order entity than edge %-d",fid,ME_ID(fe));
            MSTK_Report(funcname,mesg,MSTK_WARN);
            valid = 0;
          }
          else if (gedim == gfdim && geid != gfid) {
            sprintf(mesg,"Face %-d and edge %-d classified on different\n entities of the same dimension",fid,eid);
            MSTK_Report(funcname,mesg,MSTK_WARN);
          }
        }

	efaces = ME_Faces(fe);
	if (!List_Contains(efaces,mf)) {
	  sprintf(mesg,"Face %-d refers to edge %-d but not vice versa",fid,ME_ID(fe));
	  MSTK_Report(funcname,mesg,MSTK_ERROR);
	  valid = 0;
	}
	List_Delete(efaces);
      }
      List_Delete(fedges);
      

      fregs = MF_Regions(mf);

      if (gfdim == 3) {
	if (!fregs || List_Num_Entries(fregs) != 2) {
	  sprintf(mesg,"Interior face %-d does not have two connected regions",fid);
	  MSTK_Report(funcname,mesg,MSTK_ERROR);
	  valid = 0;
	}
      }


      if (fregs) {
	if (List_Num_Entries(fregs) == 2) {
	  if (MR_FaceDir(List_Entry(fregs,0),mf) == MR_FaceDir(List_Entry(fregs,1),mf)) {
	    sprintf(mesg,"Both regions using face %-d in the same sense",fid);
	    MSTK_Report(funcname,mesg,MSTK_ERROR);
	    valid = 0;
	  }
	}
	List_Delete(fregs);      
      }


    } /* while ((mf = MESH_Next_Face(mesh,&idx1))) */



    /*****************************************************************/
    /* Regions                                                      */
    /*****************************************************************/    

    idx1 = 0;
    while ((mr = MESH_Next_Region(mesh,&idx1))) {
      
#ifdef MSTK_HAVE_MPI
      if (MR_PType(mr) == PGHOST) continue;
#endif

      rid = MR_ID(mr);
      grid = MR_GEntID(mr);

      rfaces = MR_Faces(mr);
      int nrf = List_Num_Entries(rfaces);

      if (nrf < 4) {
        sprintf(mesg,"Region %-d has less than 4 faces",rid);
        MSTK_Report(funcname,mesg,MSTK_ERROR);
      }

      /* Check that face to region and region to face links are consistent
         with each other */
      int *rfdirs = (int *) malloc(nrf*sizeof(int));

      i = 0;
      idx2 = 0;
      while ((rf = List_Next_Entry(rfaces,&idx2))) {
	rfdirs[i] = MR_FaceDir_i(mr,i);
	if (mr != MF_Region(rf,!rfdirs[i])) {
	  sprintf(mesg,"Region %-d to face %-d dir inconsistent with \n face to region dir",rid,MF_ID(rf));
	  MSTK_Report(funcname,mesg,MSTK_ERROR);
	  valid = 0;
	}

        i++;
      }


      /* Check that faces of a region have consistent orientation in
         the region with respect to each other */
      
      for (i = 0; i < nrf; i++) {
        MFace_ptr rf, rf2;

        rf = List_Entry(rfaces,i);

        fedges = MF_Edges(rf,1,0);
        nfe = List_Num_Entries(fedges);
        
        for (j = 0; j < nfe; j++) {
          fe = List_Entry(fedges,j);
          int fedir = MF_EdgeDir_i(rf,j);
          
          /* Find adjacent face in the region */
          found = 0;
          for (k = 0; k < nrf; k++) {
            rf2 = List_Entry(rfaces,k);
            if (rf != rf2 && MF_UsesEntity(rf2,fe,MEDGE)) {
              found = 1;
              break;
            }
          }

          if (!found) {
            sprintf(mesg,"Cannot find another face in region %-d sharing edge %-d (ID = %-d) of face with ID = %-d",MR_ID(mr),j,ME_ID(fe),MF_ID(rf));
            MSTK_Report(funcname,mesg,MSTK_ERROR);
            valid = 0;
          }
          
          int fedir_adj = MF_EdgeDir(rf2,fe);
          
          /* If the two faces use the edge in opposite directions then
             the region should use the faces in the same direction and
             vice versa */

          if (((fedir_adj == fedir) && (rfdirs[i] == rfdirs[k])) ||
              ((fedir_adj != fedir) && (rfdirs[i] != rfdirs[k]))) {
            sprintf(mesg,"Region %-d faces are inconsistently oriented",MR_ID(mr));
            MSTK_Report(funcname,mesg,MSTK_ERROR);
            valid = 0;
          }
        }
        List_Delete(fedges);
           
      }
      List_Delete(rfaces);
      free(rfdirs);
    }

    return valid;

  } /* int MESH_CheckTopo */


#ifdef __cplusplus
}
#endif
