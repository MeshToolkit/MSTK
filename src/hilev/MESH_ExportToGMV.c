#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif



  /* Function to export MSTK mesh to GMV format */

  /* if natt = 0, all attributes are written out 
     if natt = -1, no attributes are written out 
     if natt > 0, only attributes specified in attnames are written out */

  /* opts is an array of flags that controls how mesh is exported to GMV 
     
     if opts[0] = 0, nodes are written out with the 'nodev' keyword
     and as triplets of (x,y,z) coordinates. If it is 1, nodes are
     written out with the 'nodes' keyword and all the x coordinates
     are written out, followed by all the y coordinates and then the z
     coordinates 

     if opts[1] = 0, cells are written out with the 'cells' keyword
     and are listed normally as tets, hexes, tris etc and point to
     node ids if opts[1] = 1, cells are written out as 'vface3d' or
     'vface2d' which point to 'vfaces' which in turn point to node ids
  */

  /* If you want default options, call the functions as

     MESH_ExportToGMV(mymesh,mygmvfile,0,NULL,NULL);
  */



int MESH_ExportToGMV(Mesh_ptr mesh, const char *filename, const int natt, 
		     const char **attnames, int *opts) {
  int			gentid, *gentities;
  MFType		ftype;
  MRType                rtype;
  List_ptr		rverts, fverts, rfaces, fregs, efaces;
  MVertex_ptr           vertex;
  MEdge_ptr             edge;
  MFace_ptr	        face;
  MRegion_ptr           region, region0, region1, oppreg;
  MAttrib_ptr           attrib, *outattribs, oppatt;
  MAttType              atttype;
  char                  attname[256], matname[256], date_str[256];
  int                   jv, jf, gmodel, nrf, nrv, nfv, dir;
  int			i, found, k, nmeshatt, noutatt, ival;
  int                   nalloc, ngent, fnum, fid, vid, nv, icr;
  int                   attentdim, j, ncells, polygons=0, cellmk, idx;
  int                   ncells1, ncells2, ndup, NFACES, rid0, rid1, rid;
  int                   oppfid;
  double		vxyz[3], rval, *xcoord, *ycoord, *zcoord;
  void                 *pval;
  FILE		        *fp;
  time_t                ctime;

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

  fprintf(fp,"codename MSTK_V_1.4\n");
  ctime = time(&ctime);
  strftime(date_str,sizeof(date_str),"%m/%d/%Y",localtime(&ctime));
  fprintf(fp,"simdate %s\n",date_str);

  if (!opts || opts[0] == 0) {
    fprintf(fp,"nodev %d\n",nv);
    for (jv = 0; jv < nv; jv++) {
      vertex = MESH_Vertex(mesh,jv);
      MV_Coords(vertex,vxyz);
      fprintf(fp,"% 20.12f % 20.12lf % 20.12lf\n",vxyz[0],vxyz[1],vxyz[2]);
      MV_Set_ID(vertex,(jv+1));
    }
  }
  else {
    fprintf(fp,"nodes %d\n",nv);
    xcoord = (double *) MSTK_malloc(nv*sizeof(double));
    ycoord = (double *) MSTK_malloc(nv*sizeof(double));
    zcoord = (double *) MSTK_malloc(nv*sizeof(double));
    for (jv = 0; jv < nv; jv++) {
      vertex = MESH_Vertex(mesh,jv);
      MV_Coords(vertex,vxyz);
      xcoord[jv] = vxyz[0];
      ycoord[jv] = vxyz[1];
      zcoord[jv] = vxyz[2];
      MV_Set_ID(vertex,(jv+1));
    }

    for (jv = 0; jv < nv; jv++) {
      fprintf(fp,"% 20.12lf",xcoord[jv]);
      if (jv%5 == 0 || jv == nv-1)
	fprintf(fp,"\n");
    }
    for (jv = 0; jv < nv; jv++) {
      fprintf(fp,"% 20.12lf",ycoord[jv]);
      if (jv%5 == 0 || jv == nv-1)
	fprintf(fp,"\n");
    }
    for (jv = 0; jv < nv; jv++) {
      fprintf(fp,"% 20.12lf",zcoord[jv]);
      if (jv%5 == 0 || jv == nv-1)
	fprintf(fp,"\n");
    }
    MSTK_free(xcoord);
    MSTK_free(ycoord);
    MSTK_free(zcoord);
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

  if (!opts || opts[1] == 0) { /* write out cells as cells*/
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
  }
  else { /* Write out cells as vface3D and faces as vfaces */
    NFACES = MESH_Num_Faces(mesh);

    oppatt = MAttrib_New(mesh,"oppfaceID",INT,MFACE); 

    fprintf(fp,"cells %d\n",ncells);
    
    ngent = 0;
    nalloc = 10;
    gentities = (int *) MSTK_malloc(nalloc*sizeof(int));
  
    ndup = 0;
    idx = 0;
    while ((region = MESH_Next_Region(mesh,&idx))) {      
      rid = MR_ID(region);

      rfaces = MR_Faces(region);
      nrf = List_Num_Entries(rfaces);
      
      fprintf(fp,"vface3d % d ",nrf);
      
      for (jf = 0; jf < nrf; jf++) {
	face = List_Entry(rfaces,jf);

	/* if face has a region on the other side, we will be writing
	   out a duplicate face for the opposite region. We use the
	   convention that if current region ID from which we are
	   looking at the face is smaller, then the face ID will be
	   written out as is, whereas for the opposite region a new
	   face with ID = face_id+num_mesh_faces will be written
	   out. Vice versa if this region ID is higher */

	fregs = MF_Regions(face);
	if (fregs && List_Num_Entries(fregs) == 2) {
	  
	  oppreg = List_Entry(fregs,0);
	  if (oppreg == region)
	    oppreg = List_Entry(fregs,1);
	  
	  if (rid < MR_ID(oppreg))
	    fprintf(fp," % d ",MF_ID(face));
	  else {
	    ndup++;
	    fprintf(fp," % d ",NFACES+ndup);
	    MEnt_Set_AttVal(face,oppatt,NFACES+ndup,0,NULL);
	  }
	}
	else
	  fprintf(fp,"% d ",MF_ID(face));
      }
      fprintf(fp,"\n");

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
    

    /* Now to write out vfaces */

    fprintf(fp,"vfaces % d\n",NFACES+ndup);

    /* Write out all the faces of the mesh first */
    ndup = 0;
    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {      

      fregs = MF_Regions(face);
      if (!fregs)	
	continue;

      if (List_Num_Entries(fregs) == 2) {
	ndup++;

	region0 = List_Entry(fregs,0);
	region1 = List_Entry(fregs,1);
	  
	/* This is the face written out with its regular ID. The
	   connected region with the lower ID will use this face such
	   that its normal points out of the region */

	rid0 = MR_ID(region0); rid1 = MR_ID(region1);
	if (rid0 < rid1) {
	  region = region0;
	  rid = rid0;
	}
	else {
	  region = region1;
	  rid = rid1;
	}

	dir = MR_FaceDir(region,face);	  
	fverts = MF_Vertices(face,dir,0);

	nfv = List_Num_Entries(fverts);
	
	fid = MF_ID(face);
	MEnt_Get_AttVal(face,oppatt,&oppfid,&rval,&pval);
	
	fprintf(fp,"% d  1  % d 1 % d ",nfv,oppfid,rid);
	for (i = 0; i < nfv; i++)
	  fprintf(fp," % d ",MV_ID(List_Entry(fverts,i)));
	fprintf(fp,"\n");
	
	List_Delete(fverts);
      }
      else {
	region = List_Entry(fregs,0);
	rid = MR_ID(region);
	dir = MR_FaceDir(region,face);	  
	fverts = MF_Vertices(face,dir,0);

	nfv = List_Num_Entries(fverts);
	
	fid = MF_ID(face);
	oppfid = 0;
	
	fprintf(fp,"% d  1  % d 1 % d ",nfv,oppfid,rid);
	for (i = 0; i < nfv; i++)
	  fprintf(fp," % d ",MV_ID(List_Entry(fverts,i)));
	fprintf(fp,"\n");
	
	List_Delete(fverts);
      }
      List_Delete(fregs);
    }

    /* Write out all the faces that are shared by two regions again.
       Since we didn't store the order which the faces were assigned
       the "duplicate face IDs" (and probably we don't want to store
       this information) we need to loop over regions (cells) again,
       to duplicate the constructive process so that the faces are
       written out in the right order.
    */
    
    idx = 0; 
    ndup = 0;

    while ((region = MESH_Next_Region(mesh,&idx))) {      

      rid = MR_ID(region);

      rfaces = MR_Faces(region);
      nrf = List_Num_Entries(rfaces);
      
      for (jf = 0; jf < nrf; jf++) {

	face = List_Entry(rfaces,jf);

	fregs = MF_Regions(face);
 
	if (fregs && List_Num_Entries(fregs) == 2) {

	  oppreg = List_Entry(fregs,0);
	  if (oppreg == region){
	    oppreg = List_Entry(fregs,1);
	  }
	  
	  if ( rid > MR_ID(oppreg)) {
	    fid = MF_ID(face);
	    dir = MR_FaceDir(region,face);	  
	    fverts = MF_Vertices(face,dir,0);
	    nfv = List_Num_Entries(fverts);
	    fprintf(fp,"% d  1  % d 1 % d ",nfv,fid,rid);

	    for (i = 0; i < nfv; i++) {
	      fprintf(fp," % d ",MV_ID(List_Entry(fverts,i)));
	    }
	    fprintf(fp,"\n");
	    List_Delete(fverts);      
	  }  
	}
	List_Delete(fregs);
      }

    }

    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx)))
      MEnt_Rem_AttVal(face,oppatt);
    MAttrib_Delete(oppatt);
  }

  
  if (ncells2) {
    if (opts && opts[1] == 1) {
      fprintf(stderr,"vface2d not yet implemented\n");
      return 0;
    }

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

