#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif


int MESH_ExportToGMV(Mesh_ptr mesh, const char *filename, const int natt, 
		     const char **attnames) {
  int			gentid, *gentities;
  MFType		ftype;
  MRType                rtype;
  List_ptr		rverts, fverts, rfaces, fregs, efaces;
  MVertex_ptr           vertex;
  MEdge_ptr             edge;
  MFace_ptr	        face;
  MRegion_ptr           region;
  MAttrib_ptr           attrib, *outattribs;
  MAttType              atttype;
  char                  attname[256], matname[256], date_str[256];
  int                   jv, jr, jf, nf, nr, gmodel, nrf, nrv, nfv, dir;
  int			i, found, k, len, suff, nmeshatt, noutatt, ival;
  int                   nalloc, ngent, fnum, fid, vid, nv, icr;
  int                   attentdim, j, ncells, polygons=0, cellmk, idx;
  int                   ncells1, ncells2;
  double		vxyz[3], rval;
  void                 *pval;
  FILE		        *fp;
  time_t                ctime;

  int			ftmpl[2][4] =   {{0,1,2,-1},{0,1,2,3}};

  int			rtmpl[5][8] =   {{0,2,1,3,-1,-1,-1,-1},
					 {1,2,3,4,0,-1,-1,-1},
					 {3,4,5,0,1,2,-1,-1},
					 {-1,-1,-1,-1,-1,-1,-1,-1},
					 {4,5,6,7,0,1,2,3}};

  gmodel = 1;
  
  if (!(fp = fopen(filename,"w"))) {
    fprintf(stderr,"mstk2gmv: Couldn't open output file %s\n",filename);
    exit(2);
  }


  fprintf(fp,"gmvinput ascii\n");
  nv = MESH_Num_Vertices(mesh);

  fprintf(fp,"codename MSTK_V_1.3\n");
  ctime = time(&ctime);
  strftime(date_str,sizeof(date_str),"%m/%d/%Y",localtime(&ctime));
  fprintf(fp,"simdate %s\n",date_str);

  fprintf(fp,"nodev %d\n",nv);
  for (jv = 0; jv < nv; jv++) {
    vertex = MESH_Vertex(mesh,jv);
    MV_Coords(vertex,vxyz);
    fprintf(fp,"% 5.9f % 5.9f % 5.9f\n",vxyz[0],vxyz[1],vxyz[2]);
    MV_Set_ID(vertex,(jv+1));
  }

  
  cellmk = MSTK_GetMarker();

  /* number of regions */

  ncells = MESH_Num_Regions(mesh); 

  /* May have problems if we want to write polygons as 'faces' */
  /* number of faces not connected to a region */
  ncells2 = 0;
  idx = 0;
  while ((face = MESH_Next_Face(mesh,&idx))) {
    if ((fregs = MF_Regions(face)))
      List_Delete(fregs);
    else {
      MEnt_Mark(face,cellmk);
      ncells++;
      ncells2++;
    }
  }
  
  /* number of edges not connected to a face */
  ncells1 = 0;
  idx = 0;
  while ((edge = MESH_Next_Edge(mesh,&idx))) {
    if ((efaces = ME_Faces(edge)))
      List_Delete(efaces);
    else {
      MEnt_Mark(edge,cellmk);
      ncells++;
      ncells1++;
    }
  }
  

  /* write region connectivity for entire mesh */

  fprintf(fp,"cells %d\n",ncells);

  ngent = 0;
  nalloc = 10;
  gentities = (int *) MSTK_malloc(nalloc*sizeof(int));
  
  idx = 0;
  while ((region = MESH_Next_Region(mesh,&idx))) {
    rverts = MR_Vertices(region);
    nrv = List_Num_Entries(rverts);
    nrf = MR_Num_Faces(region);
    switch (nrv) {
    case 4:
      rtype = TET;
      fprintf(fp,"tet 4 ");
      break;
    case 5:
      if (nrf == 5) {
	rtype = PYRAMID;
	fprintf(fp,"pyramid 5 ");
      }
      else 
	rtype = POLYHED;
      break;
    case 6:
      if (nrf == 5) { /* need additional check to fully verify */
	rtype = PRISM;
	fprintf(fp,"prism 6 ");
      }
      else
	rtype = POLYHED;
      break;
    case 8:
      if (nrf == 6) {
	rtype = HEX;
	fprintf(fp,"hex 8 ");
      }
      else
	rtype = POLYHED;
      break;
    default:
      rtype = POLYHED;
      break;
    }
    if (rtype != POLYHED) {
      for (jv = 0; jv < nrv; jv++) {
	vertex = List_Entry(rverts,rtmpl[nrv-4][jv]);
	fprintf(fp," % 8d",MV_ID(vertex));
      }
      fprintf(fp,"\n");
      List_Delete(rverts);
    }
    else {
      rfaces = MR_Faces(region);
      nrf = MR_Num_Faces(region);
      fprintf(fp,"general %d\n",nrf);
      for (jf = 0; jf < nrf; jf++) {
	face = List_Entry(rfaces,jf);
	fprintf(fp,"%d ",MF_Num_Edges(face)); /* assuming linear elements */
      }
      fprintf(fp,"\n");
      for (jf = 0; jf < nrf; jf++) {
	face = List_Entry(rfaces,jf);
	dir = MR_FaceDir_i(region,jf);
	
	fverts = MF_Vertices(face,dir,0);
	nfv = List_Num_Entries(fverts);
	
	for (jv = 0; jv < nfv; jv++) {
	  vertex = List_Entry(fverts,jv);
	  fprintf(fp," %8d",MV_ID(vertex));
	}
	fprintf(fp,"\n");
	List_Delete(fverts);
      }
      List_Delete(rfaces);
    }
    
    gentid = MR_GEntID(region);

    found = 0;
    for (i = 0; i < ngent; i++)
      if (gentities[i] == gentid) {
	found = 1;
	break;
      }
    
    if (!found) {
      if (ngent+1 >= nalloc) {
	nalloc *= 2;
	gentities = (int *) MSTK_realloc(gentities,nalloc*sizeof(int));
      }
      gentities[ngent] = gentid;
      ngent++;
    }
  }

  if (ncells2) {
    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {
      if (!MEnt_IsMarked(face,cellmk))
	continue;
      
      if (MF_Num_Vertices(face) > 4) {
	polygons = 1;
	break;
      }
    }
      
    /* Turn off writing polygons as faces for now */
    polygons = 0;      
    if (polygons) 
      fprintf(fp,"faces %d 0\n",ncells2);
      
    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {
      if (!MEnt_IsMarked(face,cellmk))
	continue;

      fnum++;
      fid = MF_ID(face);
      
      fverts = MF_Vertices(face,1,0);
      nfv = List_Num_Entries(fverts);
      if (polygons) {
	fprintf(fp,"%-d ",nfv);
	
	for (jv = 0; jv < nfv; jv++) {
	  vertex = List_Entry(fverts,jv);
	  vid = MV_ID(vertex);
	  fprintf(fp," % 8d",vid);
	}
	
	fprintf(fp," 0  0 ");
	
	fprintf(fp,"\n");
	List_Delete(fverts);
      }
      else {
	switch (nfv) {
	case 3:
	  ftype = TRI;
	  fprintf(fp,"tri 3 ");
	  break;
	case 4:
	  ftype = QUAD;
	  fprintf(fp,"quad 4 ");
	  break;
	default:
	  ftype = POLYGON;
	  fprintf(fp,"general 1\n");
	  fprintf(fp,"%-d ",nfv);
	  break;
	}
	
	for (jv = 0; jv < nfv; jv++) {
	  vertex = List_Entry(fverts,jv);
	  vid = MV_ID(vertex);
	  fprintf(fp," % 8d",vid);
	}
	
	fprintf(fp,"\n");
	List_Delete(fverts);
      }
	
      gentid = MF_GEntID(face);
	
      found = 0;
      for (i = 0; i < ngent; i++)
	if (gentities[i] == gentid) {
	  found = 1;
	  break;
	}
      
      if (!found) {
	if (ngent+1 >= nalloc) {
	  nalloc *= 2;
	  gentities = (int *) MSTK_realloc(gentities,nalloc*sizeof(int));
	}
	gentities[ngent] = gentid;
	ngent++;
      }
    }
  }

  
  /* Write out stand-alone edges (not connected to any face */

  if (ncells1) {
    idx = 0;
    while ((edge = MESH_Next_Edge(mesh,&idx))) {
      if (!MEnt_IsMarked(edge,cellmk))
	continue;
      
      fprintf(fp,"line 2 ");

      for (jv = 0; jv < 2; jv++) {
	vertex = ME_Vertex(edge,jv);
	vid = MV_ID(vertex);
	fprintf(fp," % 8d",vid);
      }
	  
      fprintf(fp,"\n");
	
      gentid = ME_GEntID(edge);
	  
      found = 0;
      for (i = 0; i < ngent; i++)
	if (gentities[i] == gentid) {
	  found = 1;
	  break;
	}
      
      if (!found) {
	if (ngent+1 >= nalloc) {
	  nalloc *= 2;
	  gentities = (int *) MSTK_realloc(gentities,nalloc*sizeof(int));
	}
	gentities[ngent] = gentid;
	ngent++;
      }
    }
  }


  if (gmodel) {
    fprintf(fp,"material %d 0\n",ngent);
    for (i = 0; i < ngent; i++) {
      sprintf(matname,"mat%-d\n",(i+1));
      fprintf(fp,"%s",matname);
    }
    
    k = 0;
    idx = 0;
    while ((region = MESH_Next_Region(mesh,&idx))) {
      gentid = MR_GEntID(region);
      
      for (i = 0, found = 0; i < ngent; i++) {
	if (gentities[i] == gentid) {
	  found = 1;
	  break;
	}
      }
      fprintf(fp,"%d ",(i+1));
      if ((++k)%10 == 0)
	fprintf(fp,"\n");
    }

    if (ncells2) {
      while ((face = MESH_Next_Face(mesh,&idx))) {
	if (!MEnt_IsMarked(face,cellmk))
	  continue;

	gentid = MF_GEntID(face);
	for (i = 0, found = 0; i < ngent; i++) {
	  if (gentities[i] == gentid) {
	    found = 1;
	    break;
	  }
	}
	fprintf(fp,"%d ",(i+1));
	if ((++k)%10 == 0)
	  fprintf(fp,"\n");
      }
    }

    if (ncells1) {
      while ((edge = MESH_Next_Edge(mesh,&idx))) {
	if (!MEnt_IsMarked(edge,cellmk))
	  continue;

	gentid = ME_GEntID(edge);
	for (i = 0, found = 0; i < ngent; i++) {
	  if (gentities[i] == gentid) {
	    found = 1;
	    break;
	  }
	}
	fprintf(fp,"%d ",(i+1));
	if ((++k)%10 == 0)
	  fprintf(fp,"\n");
      }
    }
    fprintf(fp,"\n");


    /* Other variables related to geometric classification */
    
    fprintf(fp,"variable \n");
    fprintf(fp,"itetclr  0\n");
    
    k = 0;
    idx = 0;
    while ((region = MESH_Next_Region(mesh,&idx))) {
      gentid = MR_GEntID(region);
      fprintf(fp,"%d ",gentid);      
      if ((++k)%10 == 0)
	fprintf(fp,"\n");
    }

    if (ncells2) {
      while ((face = MESH_Next_Face(mesh,&idx))) {
	if (!MEnt_IsMarked(face,cellmk))
	  continue;
	else {
	  gentid = MF_GEntID(face);
	  fprintf(fp,"%d ",gentid);
	  if ((++k)%10 == 0)
	    fprintf(fp,"\n");
	}
      }
    }

    if (ncells1) {
      while ((edge = MESH_Next_Edge(mesh,&idx))) {
	if (!MEnt_IsMarked(edge,cellmk))
	  continue;
	else {
	  gentid = ME_GEntID(edge);
	  fprintf(fp,"%d ",gentid);
	  if ((++k)%10 == 0)
	    fprintf(fp,"\n");
	}
      }
    }

    fprintf(fp,"\n");

    
    fprintf(fp,"icr1   1 \n"); 
    k = 0;
    idx = 0;
    while ((vertex = MESH_Next_Vertex(mesh,&idx))) {
      icr = 3.0-MV_GEntDim(vertex);
      fprintf(fp,"%f ",(float)icr);
      if ((++k)%10 == 0)
	fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
  }



  nmeshatt = MESH_Num_Attribs(mesh);
      
  /* Must write out other attributes that exist on mesh entities */ 
  /* If natt is 0, all attributes are to be written out */
  
  if (natt >= 0 && nmeshatt) {
    
    /* Collect the attributes */

    outattribs = (MAttrib_ptr *) MSTK_malloc(nmeshatt*sizeof(int));
    
    for (i = 0, noutatt = 0; i < nmeshatt; i++) {
      attrib = MESH_Attrib(mesh,i);
      
      attentdim = MAttrib_Get_EntDim(attrib);
      if ((attentdim == MVERTEX)  || (attentdim == MREGION) || 
	  (ncells2 && attentdim == MFACE) || 
	  (ncells1 && attentdim == MEDGE) || (attentdim == MALLTYPE)) {
	
	atttype = MAttrib_Get_Type(attrib);
	if (atttype == INT || atttype == DOUBLE) {
	  
	  MAttrib_Get_Name(attrib,attname);
	  
	  if (natt) {
	    found = 0;
	    for (j = 0; j < natt; j++) {
	      if (strcmp(attname,attnames[j]) == 0) {
		found = 1;
	      }
	    }
	    if (!found)
	      continue;
	  }
	  
	  outattribs[noutatt] = attrib;
	  noutatt++;
	}
      }
    }


    /* write them out */

    if (noutatt && !gmodel)
      fprintf(fp,"variable \n");

    for (i = 0; i < noutatt; i++) {
      attrib = outattribs[i];
      
      MAttrib_Get_Name(attrib,attname);
      atttype = MAttrib_Get_Type(attrib);
      attentdim = MAttrib_Get_EntDim(attrib);

      k = 0;
      if (attentdim == MVERTEX || attentdim == MALLTYPE) {
	fprintf(fp,"%s ",attname);
	fprintf(fp," 1 \n");
	
	idx = 0;
	while ((vertex = MESH_Next_Vertex(mesh,&idx))) {
	  
	  MEnt_Get_AttVal(vertex,attrib,&ival,&rval,&pval);
	  
	  if (atttype == INT) {
	    fprintf(fp,"%d ",ival);
	    if ((k+1)%10 == 0 && k != nv) fprintf(fp,"\n");
	    k++;
	  }
	  else if (atttype == DOUBLE) {
	    fprintf(fp,"%14.7lf ", rval);
	    if ((k+1)%5 == 0 && k != nv) fprintf(fp,"\n");
	    k++;
	  }
	}
      }
      
      if (attentdim != MVERTEX && attentdim != MUNKNOWNTYPE) {
	fprintf(fp,"%s ",attname);
	fprintf(fp," 0 \n");
	
	idx = 0;
	while ((region = MESH_Next_Region(mesh,&idx))) {
	  
	  MEnt_Get_AttVal(region,attrib,&ival,&rval,&pval);
	  
	  if (atttype == INT) {
	    fprintf(fp,"%d ",ival);
	    if ((k+1)%10 == 0) fprintf(fp,"\n");
	    k++;
	  }
	  else if (atttype == DOUBLE) {
	    fprintf(fp,"%lf ", rval);
	    if ((k+1)%5 == 0) fprintf(fp,"\n");
	    k++;
	  }
	}

	if (ncells2) {
	  idx = 0;
	  while ((face = MESH_Next_Face(mesh,&idx))) {
	    if (!MEnt_IsMarked(face,cellmk))
	      continue;

	    MEnt_Get_AttVal(face,attrib,&ival,&rval,&pval);
	    
	    if (atttype == INT) {
	      fprintf(fp,"%d ",ival);
	      if ((k+1)%10 == 0) fprintf(fp,"\n");
	      k++;
	    }
	    else if (atttype == DOUBLE) {
	      fprintf(fp,"%lf ", rval);
	      if ((k+1)%5 == 0) fprintf(fp,"\n");
	      k++;
	    }
	  }
	}
	
	if (ncells1) {
	  idx = 0;
	  while ((edge = MESH_Next_Edge(mesh,&idx))) {
	    if (!MEnt_IsMarked(edge,cellmk))
	      continue;

	    MEnt_Get_AttVal(edge,attrib,&ival,&rval,&pval);
	    
	    if (atttype == INT) {
	      fprintf(fp,"%d ",ival);
	      if ((k+1)%10 == 0) fprintf(fp,"\n");
	      k++;
	    }
	    else if (atttype == DOUBLE) {
	      fprintf(fp,"%lf ", rval);
	      if ((k+1)%5 == 0) fprintf(fp,"\n");
	      k++;
	    }
	  }
	}
      }
    }

    if (noutatt)
      fprintf(fp,"\n");
	  
    MSTK_free(outattribs);
  }
  
  fprintf(fp,"endvars \n");
  
  
  /* Finish Export */

  fprintf(fp,"endgmv\n");

  fclose(fp);



  /* Clean up */

  MSTK_free(gentities);

  if (ncells2) {
    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) 
      MEnt_Unmark(face,cellmk);
  }

  if (ncells1) {
    idx = 0;
    while ((edge = MESH_Next_Edge(mesh,&idx)))
      MEnt_Unmark(edge,cellmk);
  }

  MSTK_FreeMarker(cellmk);


  return 1;
}

#ifdef __cplusplus
  }
#endif

