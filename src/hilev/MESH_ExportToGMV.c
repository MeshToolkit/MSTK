#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif


int MESH_ExportToGMV(Mesh_ptr mesh, const char *filename, const int natt, 
		     const char **attnames) {
  int			gregid, *gregions, gfaceid, *gfaces;
  MFType		ftype;
  MRType                rtype;
  List_ptr		rverts, fverts, rfaces;
  MVertex_ptr           vertex;
  MFace_ptr	        face;
  MRegion_ptr           region;
  MAttrib_ptr           attrib, *outattribs;
  MAttType              atttype;
  char                  attname[256];
  int                   jv, jr, jf, nf, nr, gmodel, nrf, nrv, nfv, dir;
  int			i, found, k, len, suff, nmeshatt, noutatt, ival;
  int                   nalloc, ngregions, ngfaces,fnum,fid,vid, nv, icr;
  int                   attentdim, j;
  double		vxyz[3], rval;
  void                 *pval;
  char                  matname[256];
  FILE		        *fp;

  /* Although the GMV template for tets is 0,2,1,3, LAGRIT thinks it is
     0,1,2,3 and GMV does not seem to care */

  int			ftmpl[2][8] =   {{0,1,2,-1,-1,-1,-1,-1},
					 {0,1,2,3,-1,-1,-1,-1}};

  int			rtmpl[4][8] =   {{0,1,2,3,-1,-1,-1,-1},
					 {4,0,1,2,3,-1,-1,-1},
					 {3,4,5,0,1,2,-1,-1},
					 {4,5,6,7,0,1,2,3}};

  gmodel = 0;
  
  if (!(fp = fopen(filename,"w"))) {
    fprintf(stderr,"mstk2gmv: Couldn't open output file %s\n",filename);
    exit(2);
  }


  fprintf(fp,"gmvinput ascii\n");
  nv = MESH_Num_Vertices(mesh);

  fprintf(fp,"nodev %d\n",nv);
  vid = 0;
  for (jv = 0; jv < nv; jv++) {
    vertex = MESH_Vertex(mesh,jv);
    MV_Coords(vertex,vxyz);
    fprintf(fp,"% 5.9f % 5.9f % 5.9f\n",vxyz[0],vxyz[1],vxyz[2]);
    MV_Set_ID(vertex,vid++);
  }


  /* write region connectivity for entire mesh */

  if ((nr = MESH_Num_Regions(mesh))) {
    fprintf(fp,"cells %d\n",nr);

    ngregions = 0;
    nalloc = 10;
    gregions = (int *) malloc(nalloc*sizeof(int));

    for (jr = 0; jr < nr; jr++) {
      region = MESH_Region(mesh, jr);
      rverts = MR_Vertices(region);
      nrv = List_Num_Entries(rverts);
      switch (nrv) {
      case 4:
	rtype = TET;
	fprintf(fp,"tet 4 ");
	break;
      case 5:
        rtype = PYRAMID;
	fprintf(fp,"pyramid 5 ");
	break;
      case 6:
	rtype = PRISM;
	fprintf(fp,"prism 6 ");
	break;
      case 8:
	rtype = HEX;
	fprintf(fp,"hex 8 ");
        break;
      default:
	rtype = POLYHED;
	break;
      }
      if (rtype != POLYHED) {
	for (jv = 0; jv < nrv; jv++) {
	  vertex = List_Entry(rverts,rtmpl[rtype-1][jv]);
	  fprintf(fp," % 8d",MV_ID(vertex)+1);
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
	    fprintf(fp," %8d",MV_ID(vertex)+1);
	  }
	  fprintf(fp,"\n");
	  List_Delete(fverts);
	}
	List_Delete(rfaces);
      }
	
      gregid = MR_GEntID(region);
      if (gregid) { 
	/* Assumes that if one region has classification info, all do */
	gmodel = 1;

	found = 0;
	for (i = 0; i < ngregions; i++)
	  if (gregions[i] == gregid) {
	    found = 1;
	    break;
	  }
	
	if (!found) {
	  if (ngregions+1 >= nalloc) {
	    nalloc *= 2;
	    gregions = (int *) realloc(gregions,nalloc*sizeof(int));
	  }
	  gregions[ngregions] = gregid;
	  ngregions++;
	}
      }
    }

    if (gmodel) {
      fprintf(fp,"material %d 0\n",ngregions);
      for (i = 0; i < ngregions; i++) {
	sprintf(matname,"mat%-d\n",(i+1));
	fprintf(fp,"%s",matname);
      }
      
      for (jr = 0, k = 0; jr < nr; jr++) {
	region = MESH_Region(mesh,jr);
	gregid = MR_GEntID(region);

	for (i = 0, found = 0; i < ngregions; i++) {
	  if (gregions[i] == gregid) {
	    found = 1;
	    break;
	  }
	}
	fprintf(fp,"%d ",(i+1));

	k++;
	if (k == 10) {
	  fprintf(fp,"\n");
	  k = 0;
	}
      }
      fprintf(fp,"\n");
      
      fprintf(fp,"variable \n");
      fprintf(fp,"itetclr  0\n");

      for (jr = 0, k = 0; jr < nr; jr++) {
	region = MESH_Region(mesh,jr);
	gregid = MR_GEntID(region);

	for (i = 0, found = 0; i < ngregions; i++) {
	  if (gregions[i] == gregid) {
	    found = 1;
	    break;
	  }
	}
	fprintf(fp,"%d ",(i+1));

	k++;
	if (k == 10) {
	  fprintf(fp,"\n");
	  k = 0;
	}
      }
      fprintf(fp,"\n");

      fprintf(fp,"icr1   1 \n"); k = 0;
      for (jv = 0; jv < nv; jv++) {
	vertex = MESH_Vertex(mesh,jv);
        icr = 3.0-MV_GEntDim(vertex);
	fprintf(fp,"%f ",(float)icr);
	k++;
	if (k%10 == 0)
	  fprintf(fp,"\n");
      }
      fprintf(fp,"\n");

      /* Must write out other attributes that exist on mesh regions 
	 and mesh nodes */

      if (natt >= 0) {
	/* Must write out attributes related to mesh faces and mesh nodes */
	/* If natt is 0, all attributes are to be written out */
	
	nmeshatt = MESH_Num_Attribs(mesh);

	if (nmeshatt) {
	  
	  outattribs = (MAttrib_ptr *) malloc(nmeshatt*sizeof(int));
	  
	  for (i = 0, k = 0, noutatt = 0; i < nmeshatt; i++) {
	    attrib = MESH_Attrib(mesh,i);
	    
	    attentdim = MAttrib_Get_EntDim(attrib);
	    if (attentdim != MVERTEX  && attentdim != MREGION)
	      continue;
	    
	    atttype = MAttrib_Get_Type(attrib);
	    if (atttype != INT && atttype != DOUBLE)
	      continue; /* can't write other types out */
	    
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
	  
	  for (i = 0; i < noutatt; i++) {
	    attrib = outattribs[i];
	    
	    MAttrib_Get_Name(attrib,attname);
	    atttype = MAttrib_Get_Type(attrib);
	    attentdim = MAttrib_Get_EntDim(attrib);
 
	    fprintf(fp,"%s ",attname);
	    if (attentdim == MVERTEX) {
	      fprintf(fp," 1 \n");
	      
	      for (jv = 0, k = 0; jv < nv; jv++) {
		vertex = MESH_Vertex(mesh,jv);
		
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
	    else {
	      fprintf(fp," 0 \n");
	      
	      for (jr = 0, k = 0; jr < nr; jr++) {
		region = MESH_Region(mesh,jr);
		
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
	    }
	  }
	  
	  free(outattribs);
	}
      }

      fprintf(fp,"endvars \n");
    }

    free(gregions);
  }
  else if ((nf = MESH_Num_Faces(mesh))) {
    int polygons = 0;

    for (jf = 0, fnum = 0; jf < nf; jf++) {
      face = MESH_Face(mesh,jf);
      if (MF_Num_Vertices(face) > 4) {
	polygons = 1;
	break;
      }
    }

    /* Turn off writing polygons as faces for now */
    polygons = 0;

    if (polygons) 
      fprintf(fp,"faces %d 0\n",nf);
    else
      fprintf(fp,"cells %d\n",nf);

    ngfaces = 0;
    nalloc = 10;
    gfaces = (int *) malloc(nalloc*sizeof(int));

    for (jf = 0, fnum = 0; jf < nf; jf++) {
      face = MESH_Face(mesh,jf);
      fnum++;
      fid = MF_ID(face);

      fverts = MF_Vertices(face,1,0);
      nfv = List_Num_Entries(fverts);
      if (polygons) {
	fprintf(fp,"%-d ",nfv);

	for (jv = 0; jv < nfv; jv++) {
	  vertex = List_Entry(fverts,jv);
	  vid = MV_ID(vertex);
	  fprintf(fp," % 8d",vid+1);
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
	  fprintf(fp," % 8d",vid+1);
	}

	fprintf(fp,"\n");
	List_Delete(fverts);
      }

      gfaceid = MF_GEntID(face);
      if (gfaceid) {
	/* Assumes that if one face has classification, all do */
	gmodel = 1;

	found = 0;
	for (i = 0; i < ngfaces; i++)
	  if (gfaces[i] == gfaceid) {
	    found = 1;
	    break;
	  }
	
	if (!found) {
	  if (ngfaces+1 >= nalloc) {
	    nalloc *= 2;
	    gfaces = (int *) realloc(gfaces,nalloc*sizeof(int));
	  }
	  gfaces[ngfaces] = gfaceid;
	  ngfaces++;
	}
      }
    }

    if (gmodel) {
      fprintf(fp,"material %d 0\n",ngfaces);
      for (i = 0; i < ngfaces; i++) {
	sprintf(matname,"face%-d\n",(i+1));
	fprintf(fp,"%s",matname);
      }
      fprintf(fp,"\n");
      
      k = 0;
      for (jf = 0; jf < nf; jf++) {
	face = MESH_Face(mesh,jf);
	gfaceid = MF_GEntID(face);
	found = 0;
	for (i = 0; i < ngfaces; i++) {
	  if (gfaces[i] == gfaceid) {
	    found = 1;
	    break;
	  }
	}
	fprintf(fp,"%d ",(i+1));
	k++;
	if (k > 10) {
	  fprintf(fp,"\n");
	  k = 0;
	}
      }
      fprintf(fp,"\n");


      /* FIELD VARIABLES OR ENTITY ATTRIBUTES */
      k = 0;
      fprintf(fp,"variable \n");

      /* SPECIAL MESH DATA */

      fprintf(fp,"itetclr  0\n");
      for (jf = 0; jf < nf; jf++) {
	face = MESH_Face(mesh,jf);
	gfaceid = MF_GEntID(face);
	found = 0;
	for (i = 0; i < ngfaces; i++) {
	  if (gfaces[i] == gfaceid) {
	    found = 1;
	    break;
	  }
	}
	fprintf(fp,"%d ",(i+1));
	k++;
	if (k > 10) {
	  fprintf(fp,"\n");
	  k = 0;
	}
      }
      fprintf(fp,"\n");

      fprintf(fp,"icr1   1 \n"); k = 0;
      for (jv = 0; jv < nv; jv++) {
	vertex = MESH_Vertex(mesh,jv);
        icr = 3-MV_GEntDim(vertex);
	fprintf(fp,"%f ",(float)icr);
	k++;
	if (k%10 == 0)
	  fprintf(fp,"\n");
      }
      fprintf(fp,"\n");


      /* OTHER ATTRIBUTES, MOST PROBABLY REPRESENTING FIELD DATA */

      if (natt >= 0) {
	/* Must write out attributes related to mesh faces and mesh nodes */
	
	nmeshatt = MESH_Num_Attribs(mesh);

	if (nmeshatt) {
	
	  outattribs = (MAttrib_ptr *) malloc(nmeshatt*sizeof(int));
	  
	  for (i = 0, k = 0, noutatt = 0; i < nmeshatt; i++) {
	    attrib = MESH_Attrib(mesh,i);
	    
	    attentdim = MAttrib_Get_EntDim(attrib);
	    if (attentdim != MVERTEX  && attentdim != MFACE)
	      continue;
	    
	    atttype = MAttrib_Get_Type(attrib);
	    if (atttype != INT && atttype != DOUBLE)
	      continue; /* can't write other types out */
	    
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
	  
	  for (i = 0; i < noutatt; i++) {
	    attrib = outattribs[i];
	    
	    MAttrib_Get_Name(attrib,attname);
	    atttype = MAttrib_Get_Type(attrib);
	    attentdim = MAttrib_Get_EntDim(attrib);
	    
	    fprintf(fp,"%s ",attname);
	    if (attentdim == MVERTEX) {
	      fprintf(fp," 1 \n");
	      
	      for (jv = 0, k = 0; jv < nv; jv++) {
		vertex = MESH_Vertex(mesh,jv);
		
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
	    else {
	      fprintf(fp," 0 \n");
	      
	      for (jf = 0, k = 0; jf < nf; jf++) {
		face = MESH_Face(mesh,jf);
		
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
	  }
	  
	  free(outattribs);
	}
      }

      fprintf(fp,"endvars \n");
    }

    free(gfaces);
  }
  fprintf(fp,"endgmv\n");

  fclose(fp);

  return 1;
}

#ifdef __cplusplus
  }
#endif

