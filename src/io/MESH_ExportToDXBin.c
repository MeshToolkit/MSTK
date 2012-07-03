#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MSTK.h"

void MESH_ExportToDXBin(Mesh_ptr mesh, const char *fname) {
  int i, k, k1, bt, bq, bgf, bpy, bpr, bh, nreg, np, offset;
  int nvf, nrf, objnum, fed0, fed1, fed2, rfd0, flip;
  int nbtris=0, nbquads=0, ngenf=0, ntets=0, npyrs=0, nprisms=0, nhexes=0;
  int *idbtri=0, *idbquad=0, *idbtet=0, *idbprm=0, *idbpyr=0, *idbhex=0;
  int *edgarr, *loopstartarr, *facearr, *polcenidx;
  int nfaces, nedges=0, index, maxvf, maxtri=0;
  int genpoly2, id, nfv;
  FILE *fout, *fdat;
  MVertex_ptr vertex, evert0, evert1, lvert[8];
  List_ptr verts, elist, rfaces, fedges;
  MRegion_ptr reg;
  MFace_ptr face, rface, rface0, rface1, rface2;
  MEdge_ptr edge, ecmn;
  float  flxyz[3];
  float *dattri=0, *datquad=0, *datgenf=0, *dattet=0, *datprm=0, *datpyr=0, *dathex=0;
  double xyz[3], fxyz[MAXPV2][3];
  double (*polcen2)[3];
  char   basename[256], datfname[256];
  char   xsb[4];
  int rvtmpl[4][8]={{0,1,2,3,-1,-1,-1,-1},
		    {0,1,2,3,4,-1,-1,-1},
		    {0,1,2,3,4,5,-1,-1},
		    {0,1,3,2,4,5,7,6}};
  int haveModel = 1;
 
  strcpy(basename,fname);
  basename[strlen(fname)-3] = '\0'; 

  strcpy(datfname,basename);
  strcat(datfname,".dat");
  fdat = fopen(datfname,"wb");


  /* Write out vertex coordinates */
  
  index = 0;
  np = 0;
  while( (vertex = MESH_Next_Vertex(mesh,&index)) ){
    MV_Coords(vertex,xyz);
    for (i = 0; i < 3; i++) flxyz[i] = (float) xyz[i];
    fwrite((void *)flxyz,sizeof(float),3,fdat);
    np++;
  }



  nfaces = MESH_Num_Faces(mesh);

  

  /* Check if we have a polygonal mesh and if so, calculate center points */
  /* Calculate center points for the polygons based on least square
     fitting of a sphere */

  /* We are still not checking for the center point being inside the
     polygon which may be non-convex */
  
  polcenidx = (int *) malloc((nfaces+1)*sizeof(int));
  polcen2 = (double (*)[3]) malloc((nfaces+1)*sizeof(double [3]));
  
  index = 0; bgf = 0; maxvf = 0;
  while ((face = MESH_Next_Face(mesh,&index))) {
    if (!MF_Region(face,0) && !MF_Region(face,1)) {
      
      nvf = MF_Num_Edges(face);
      if (nvf > maxvf)
	maxvf = nvf;
      
      if (nvf > 3) {
	genpoly2 = 1;
	
	id = MF_ID(face);
	
	MF_Coords(face,&nfv,fxyz);
	for (k = 0; k < 3; k++) polcen2[id][k] = 0.0;
	for (i = 0; i < nfv; i++) {
	  for (k = 0; k < 3; k++)
	    polcen2[id][k] += fxyz[i][k];
	}
	for (k = 0; k < 3; k++) polcen2[id][k] /= nfv;

	for (k = 0; k < 3; k++)
	  flxyz[k] = polcen2[id][k];
	
	fwrite((void *)flxyz,sizeof(float),3,fdat);
	polcenidx[id] = np;
	np++;
      }
    }
  }
  

  /* Write out triangles, quads and triangular subdivision of polygons */
  
  nbtris = nbquads = 0; bt = bq = 0;
  index = 0;
  while ((face = MESH_Next_Face(mesh,&index)))
    if (!MF_Region(face,0) && !MF_Region(face,1)){
      verts = MF_Vertices(face,1,0);
      nvf = List_Num_Entries(verts);
      if (nvf == 3) {
	if (!idbtri) {
	  maxtri = nfaces;
	  idbtri = (int *)malloc(3*maxtri*sizeof(int));
	  dattri = (float *)malloc(maxtri*sizeof(float));
	}
	nbtris++;
	for (k = 0; k < nvf; k++)
	  idbtri[bt++] = MV_ID(List_Entry(verts,k))-1;	
	dattri[nbtris-1] = MF_GEntID(face);
      }
      /* We have to disallow direct writing of quads since 
	 non-convex quads are not being rendered correctly in DX.
	 They will have to be written as 4 triangles */
      /*
      else if (nvf == 4) {
	if (!idbquad) {
	  idbquad = (int *)malloc(4*nfaces*sizeof(int));
	  datquad = (float *)malloc(nfaces*sizeof(float));
	}
	nbquads++;
	for (k = 0; k < 2; k++)
	  idbquad[bq++] = MV_ID(List_Entry(verts,k))-1;
	idbquad[bq++] = MV_ID(List_Entry(verts,3))-1;
	idbquad[bq++] = MV_ID(List_Entry(verts,2))-1;	
	datquad[nbquads-1] = MF_GEntID(face);
      }
      */
      else {
	if (!idbtri) {
	  maxtri = nfaces;
	  idbtri = (int *)malloc(3*maxtri*sizeof(int));
	  dattri = (float *)malloc(maxtri*sizeof(float));
	}
	
	for (k = 0; k < nvf; k++) {
	  nbtris++;
	  if (nbtris > maxtri) {
	    maxtri = 2*nbtris;
	    idbtri = (int *) realloc(idbtri,3*maxtri*sizeof(int));
	    dattri = (float *) realloc(dattri,maxtri*sizeof(float));
	  }
	  idbtri[bt++] = MV_ID(List_Entry(verts,k))-1;
	  idbtri[bt++] = MV_ID(List_Entry(verts,(k+1)%nvf))-1;
	  idbtri[bt++] = polcenidx[MF_ID(face)];
	  dattri[nbtris-1] = MF_GEntID(face);
	}
      }
      List_Delete(verts);
    }
  
  if (nbtris) {
    fwrite((void *)idbtri,sizeof(int),3*nbtris,fdat);
    if (haveModel)
      fwrite((void *)dattri,sizeof(int),nbtris,fdat);
  }
  
  if (nbquads) {
    fwrite((void *)idbquad,sizeof(int),4*nbquads,fdat);
    if (haveModel)
      fwrite((void *)datquad,sizeof(int),nbquads,fdat);
  }
  
  free(idbtri);
  free(dattri);
  free(idbquad);
  free(datquad);
  
  

  /* If the mesh has more general elements that triangles and quads, then
     write out every mesh faces in terms of edges and loops */

  if (maxvf > 3) {
    datgenf = (float *)malloc(nfaces*sizeof(float));

    edgarr = (int *) malloc(maxvf*nfaces*sizeof(int));
    loopstartarr = (int *) malloc(nfaces*sizeof(int));
    facearr = (int *) malloc(nfaces*sizeof(int));

    ngenf = 0; bgf = 0;
    index = 0;
    while ((face = MESH_Next_Face(mesh,&index)) )
      if (!MF_Region(face,0) && !MF_Region(face,1)) {              
	loopstartarr[ngenf] = bgf;
	datgenf[ngenf] = MF_GEntID(face);
	facearr[ngenf] = ngenf; /* face has no holes */
	ngenf++;

	verts = MF_Vertices(face,1,0);
	nvf = List_Num_Entries(verts);

	for (i = 0; i < nvf; i++)
	  edgarr[bgf++] = MV_ID(List_Entry(verts,i))-1;

	List_Delete(verts);
      }
    nedges = bgf;
    
    if (ngenf) {
      fwrite((void *)edgarr,sizeof(int),nedges,fdat);
      fwrite((void *)loopstartarr,sizeof(int),nfaces,fdat);
      fwrite((void *)facearr,sizeof(int),nfaces,fdat);
      /**/
      if (haveModel)
	fwrite((void *)datgenf,sizeof(float),nfaces,fdat);
      /**/
    } 
  }


  /* Write out mesh regions as necessary */

  if ((nreg = MESH_Num_Regions(mesh))) {

    bt = bpr = bpy = bh = 0;
    index = 0; i = 0; ntets = nprisms = npyrs = nhexes = 0;
    while ((reg = MESH_Next_Region(mesh,&index))) {
      nrf = MR_Num_Faces(reg);
      switch(nrf) {
      case 4:
	if (!idbtet) {
	  idbtet = (int *)malloc(4*nreg*sizeof(int));
	  dattet = (float *)malloc(nreg*sizeof(float));
	}
	verts = MR_Vertices(reg);
        for (k = 0; k < 4; k++) {
	  k1 = rvtmpl[0][k];
	  idbtet[bt++] = MV_ID(List_Entry(verts,k))-1;
	}
	List_Delete(verts);
	ntets++;

	dattet[ntets-1] = MR_GEntID(reg);

	break;

      case 5:
	rfaces = MR_Faces(reg);
	rface0 = List_Entry(rfaces,0);
	rfd0 = MR_FaceDir_i(reg,0);

	switch (MF_Num_Edges(rface0)) {

	case 3: /****** TRIANGULAR PRISM *******/
	  if (!idbprm) {
	    idbprm = (int *)malloc(6*nreg*sizeof(int));
	    datprm = (float *)malloc(nreg*sizeof(float));
	  }

	  verts =  MF_Vertices(rface0,!rfd0,0);
	  for (k = 0; k < 3; k++) 
	    lvert[k] = List_Entry(verts,k);
	  List_Delete(verts);

	  elist = MF_Edges(rface0,!rfd0,lvert[0]);
	  ecmn = List_Entry(elist,0);
	  List_Delete(elist);
	  fed0 = MF_EdgeDir(rface0,ecmn);

	  rface1 = NULL;
	  for (k = 1; k < 6; k++) {
	    rface = List_Entry(rfaces,k);
	    if (rface != rface0 && MF_UsesEntity(rface,ecmn,1)) {
	      rface1 = rface;
	      break;
	    }
	  }
	  fed1 = MF_EdgeDir(rface1,ecmn);

	  fedges = MF_Edges(rface1,1,0);
	  evert0 = evert1 = NULL;
	  for (k = 0; k < 4; k++) {
	    edge = List_Entry(fedges,k);
	    if (edge == ecmn)
	      continue;
	    evert0 = ME_Vertex(edge,0); evert1 = ME_Vertex(edge,1);
	    if (evert0 != lvert[0] && evert1 != lvert[0] &&
		evert0 != lvert[1] && evert1 != lvert[1]) {
	      ecmn = edge;
	      break;
	    }
	  }
	  List_Delete(fedges);
	  fed2 = MF_EdgeDir(rface1,ecmn);

	  flip = 0;
	  if (fed0^(!rfd0)) /* lvert0 and lvert1 are opp to evert00,evert01 */
	    flip = !flip;
	  if (fed1 == fed2) /* edges pointing in different dirs w.r.t face 1 */
	    flip = !flip;
	  lvert[3] = flip==0 ? evert0 : evert1;
	  lvert[4] = flip==0 ? evert1 : evert0;

	  rface2 = 0;
	  for (k = 1; k < 5; k++) {
	    rface = List_Entry(rfaces,k);
	    if (rface != rface0 && rface != rface1 && 
		MF_UsesEntity(rface,ecmn,1)){
	      rface2 = rface;
	      break;
	    }
	  }

	  /* Vertex opposite to edge, ecmn, in face, rface2 */
	  verts = MF_Vertices(rface2,1,0);
	  for (k = 0; k < 3; k++) {
	    vertex = List_Entry(verts,k);
	    if (!ME_UsesEntity(ecmn,vertex,0)) {
	      lvert[5] = vertex;
	    }
	  }
	  List_Delete(verts);
	    
	  idbprm[bpr]   = MV_ID(lvert[0])-1; 
	  idbprm[bpr+1] = MV_ID(lvert[1])-1;
	  idbprm[bpr+2] = MV_ID(lvert[2])-1; 
	  idbprm[bpr+3] = MV_ID(lvert[2])-1;
	  idbprm[bpr+4] = MV_ID(lvert[3])-1; 
	  idbprm[bpr+5] = MV_ID(lvert[4])-1;
	  idbprm[bpr+6] = MV_ID(lvert[5])-1; 
	  idbprm[bpr+7] = MV_ID(lvert[5])-1;
	  bpr += 8;
	  nprisms++;

	  datprm[nprisms-1] = MR_GEntID(reg);
	  break;

	case 4: /***** PYRAMID *****/
	  if (!idbpyr) {
	    idbpyr = (int *)malloc(5*nreg*sizeof(int));
	    datpyr = (float *)malloc(nreg*sizeof(float));
	  }

	  verts =  MF_Vertices(rface0,!MR_FaceDir_i(reg,0),0);
	  for (k = 0; k < 4; k++)
	    lvert[k] = List_Entry(verts,k);
	  List_Delete(verts);
	  rface1 = List_Entry(rfaces,1); /* has to be connected to rface0 */
	  fedges = MF_Edges(rface1,1,0);
	  evert0 = evert1 = NULL;
	  for (k = 0; k < 3; k++) {
	    edge = List_Entry(fedges,k);
	    if (!MF_UsesEntity(rface0,edge,1)) {
	      evert0 = ME_Vertex(edge,0);
	      if (evert0 == lvert[0] || evert0 == lvert[1] || 
		  evert0 == lvert[2] || evert0 == lvert[3])
		lvert[4] = ME_Vertex(edge,1);
	      else
		lvert[4] = evert0;
	      break;
	    }
	  }
	  List_Delete(fedges);
	  idbpyr[bpy]   = MV_ID(lvert[0])-1; 
	  idbpyr[bpy+1] = MV_ID(lvert[1])-1;
	  idbpyr[bpy+2] = MV_ID(lvert[3])-1; 
	  idbpyr[bpy+3] = MV_ID(lvert[2])-1;
	  idbpyr[bpy+4] = idbpyr[bpy+5] = idbpyr[bpy+6] = idbpyr[bpy+7] = 
	    MV_ID(lvert[4])-1;
	  bpy += 8;
	  npyrs++;

	  datpyr[npyrs-1] = MR_GEntID(reg);

	  break;
	}
	break;
      case 6: /******* HEXAHEDRON *********/
	if (!idbhex) {
	  idbhex = (int *)malloc(8*nreg*sizeof(int));
	  dathex = (float *)malloc(nreg*sizeof(float));
	}

	verts = MR_Vertices(reg);
        for (k = 0; k < 8; k++) {
	  k1 = rvtmpl[3][k];
	  idbhex[bh++] = MV_ID(List_Entry(verts,k1))-1;
	}
	List_Delete(verts);
	nhexes++;

	dathex[nhexes-1] = MR_GEntID(reg);

	break;
      default: /* General - have to do it some other way */
	fprintf(stderr,"Not handled\n");
      }
    }

    if (ntets) {
      fwrite((void *)idbtet,sizeof(int),4*ntets,fdat);
      if (haveModel)
	fwrite((void *)dattet,sizeof(float),ntets,fdat);
    }
    if (npyrs) {
      fwrite((void *)idbpyr,sizeof(int),8*npyrs,fdat);
      if (haveModel)
	fwrite((void *)datpyr,sizeof(float),npyrs,fdat);
    }
    if (nprisms) {
      fwrite((void *)idbprm,sizeof(int),8*nprisms,fdat);
      if (haveModel)
	fwrite((void *)datprm,sizeof(float),nprisms,fdat);
    }
    if (nhexes) {
      fwrite((void *)idbhex,sizeof(int),8*nhexes,fdat);
      if (haveModel)
	fwrite((void *)dathex,sizeof(float),nhexes,fdat);
    }
  }
  fclose(fdat);
    
  free(idbtet);

  /* Now write into the .dx file */
  fout = fopen(fname,"w");
  
#ifdef i686_linux
  strcpy(xsb,"lsb");
#else 
  strcpy(xsb,"msb");
#endif

  objnum = 1;
  fprintf(fout,"object %d class array type float rank 1 shape 3 items %d %3s binary\n",objnum++,np,xsb);
  fprintf(fout,"data file %s,0 \n",datfname);
  fprintf(fout," attribute \"dep\" string \"positions\" \n\n");
  offset = 3*np*sizeof(float);


  if (nbtris) {
    fprintf(fout,"object %d class array type int rank 1 shape 3 items %d %3s binary\n",objnum++,nbtris,xsb);
    fprintf(fout,"data file %s,%d\n",datfname,offset);
    fprintf(fout,"attribute \"element type\" string \"triangles\"\n");
    fprintf(fout,"attribute \"ref\" string \"positions\"\n\n");
    offset += 3*nbtris*sizeof(int);
    
    if (haveModel) {
      fprintf(fout,"object %d class array type float rank 0 items %d %3s binary\n",objnum++,nbtris,xsb);
      fprintf(fout,"data file %s,%d\n",datfname,offset);
      fprintf(fout,"attribute \"dep\" string \"connections\"\n");
      offset += nbtris*sizeof(float);
    }
  }
  
  if (nbquads) {
    fprintf(fout,"object %d class array type int rank 1 shape 4 items %d %3s binary\n",objnum++,nbquads,xsb);
    fprintf(fout,"data file %s,%d\n",datfname,offset);
    fprintf(fout,"attribute \"element type\" string \"quads\"\n");
    fprintf(fout,"attribute \"ref\" string \"positions\"\n\n");
    offset += 4*nbquads*sizeof(int);
    
    if (haveModel) {
      fprintf(fout,"object %d class array type float rank 0 items %d %3s binary\n",objnum++,nbquads,xsb);
      fprintf(fout,"data file %s,%d\n",datfname,offset);
      fprintf(fout,"attribute \"dep\" string \"connections\"\n");
      offset += nbquads*sizeof(float);
    }
  }
  
  
  if (ngenf) {
    fprintf(fout,"object %d class array type int rank 0 items %d %3s binary\n",objnum++,nedges,xsb);
    fprintf(fout,"data file %s,%d\n",datfname,offset);
    fprintf(fout,"attribute \"ref\" string \"positions\"\n\n");
    offset += nedges*sizeof(int);
    
    fprintf(fout,"object %d class array type int rank 0 items %d %3s binary\n",objnum++,ngenf,xsb);
    fprintf(fout,"data file %s,%d\n",datfname,offset);
    fprintf(fout,"attribute \"ref\" string \"edges\"\n\n");
    offset += ngenf*sizeof(int);
    
    fprintf(fout,"object %d class array type int rank 0 items %d %3s binary\n",objnum++,ngenf,xsb);
    fprintf(fout,"data file %s,%d\n",datfname,offset);
    fprintf(fout,"attribute \"ref\" string \"loops\"\n\n");
    offset += ngenf*sizeof(int);

    if (haveModel) {
      fprintf(fout,"object %d class array type float rank 0 items %d %3s binary\n",objnum++,ngenf,xsb);
      fprintf(fout,"data file %s,%d\n",datfname,offset);
      fprintf(fout,"attribute \"dep\" string \"faces\"\n\n");
      offset += ngenf*sizeof(float);
    }
  }
  
  
  if (nreg) {
    if (ntets) {
      fprintf(fout,"object %d class array type int rank 1 shape 4 items %d %3s binary\n",objnum++,ntets,xsb);
      fprintf(fout,"data file %s,%d\n",datfname,offset);
      fprintf(fout,"attribute \"element type\" string \"tetrahedra\"\n");
      fprintf(fout,"attribute \"ref\" string \"positions\"\n\n");
      offset += 4*ntets*sizeof(int);

      if (haveModel) {
	fprintf(fout,"object %d class array type float rank 0 items %d %3s binary\n",objnum++,ntets,xsb);
	fprintf(fout,"data file %s,%d\n",datfname,offset);
	fprintf(fout,"attribute \"dep\" string \"connections\"\n");
	offset += ntets*sizeof(float);
      }
    }

    if (npyrs) {
      fprintf(fout,"object %d class array type int rank 1 shape 8 items %d %3s binary\n",objnum++,npyrs,xsb);
      fprintf(fout,"data file %s,%d\n",datfname,offset);
      fprintf(fout,"attribute \"element type\" string \"cubes\"\n");
      fprintf(fout,"attribute \"ref\" string \"positions\"\n\n");
      offset += 8*npyrs*sizeof(int);

      if (npyrs && haveModel) {
	fprintf(fout,"object %d class array type float rank 0 items %d %3s binary\n",objnum++,npyrs,xsb);
	fprintf(fout,"data file %s,%d\n",datfname,offset);
	fprintf(fout,"attribute \"dep\" string \"connections\"\n");
	offset += npyrs*sizeof(float);
      }
    }

    if (nprisms) {
      fprintf(fout,"object %d class array type int rank 1 shape 8 items %d %3s binary\n",objnum++,nprisms,xsb);
      fprintf(fout,"data file %s,%d\n",datfname,offset);
      fprintf(fout,"attribute \"element type\" string \"cubes\"\n");
      fprintf(fout,"attribute \"ref\" string \"positions\"\n\n");
      offset += 8*nprisms*sizeof(int);

      if (nprisms && haveModel) {
	fprintf(fout,"object %d class array type float rank 0 items %d %3s binary\n",objnum++,nprisms,xsb);
	fprintf(fout,"data file %s,%d\n",datfname,offset);
	fprintf(fout,"attribute \"dep\" string \"connections\"\n");
	offset += nprisms*sizeof(float);
      }
    }

    if (nhexes) {
      fprintf(fout,"object %d class array type int rank 1 shape 8 items %d %3s binary\n",objnum++,nhexes,xsb);
      fprintf(fout,"data file %s,%d\n",datfname,offset);
      fprintf(fout,"attribute \"element type\" string \"cubes\"\n");
      fprintf(fout,"attribute \"ref\" string \"positions\"\n\n");
      offset += 8*nhexes*sizeof(int);

      if (nhexes && haveModel) {
	fprintf(fout,"object %d class array type float rank 0 items %d %3s binary\n",objnum++,nhexes,xsb);
	fprintf(fout,"data file %s,%d\n",datfname,offset);
	fprintf(fout,"attribute \"dep\" string \"connections\"\n");
	offset += nhexes*sizeof(float);
      }
    }
  }

  fprintf(fout,"\n");

  objnum = 2;
  if (nbtris) {
    fprintf(fout,"object \"tris\" class field\n");
    fprintf(fout,"component \"positions\" value 1\n");
    fprintf(fout,"component \"connections\" value %d\n",objnum++);
    if (haveModel)
      fprintf(fout,"component \"data\" value %d\n",objnum++);
  }
  
  if (nbquads) {
    fprintf(fout,"object \"quads\" class field\n");
    fprintf(fout,"component \"positions\" value 1\n");
    fprintf(fout,"component \"connections\" value %d\n",objnum++);
    if (haveModel)
      fprintf(fout,"component \"data\" value %d\n",objnum++);
  }
  
  if (ngenf) {
    fprintf(fout,"object \"polygons\" class field\n");
    fprintf(fout,"component \"positions\" value 1\n");
    fprintf(fout,"component \"edges\" value %d\n",objnum++);
    fprintf(fout,"component \"loops\" value %d\n",objnum++);
    fprintf(fout,"component \"faces\" value %d\n",objnum++);
    /**/ if (haveModel) 
      fprintf(fout,"component \"data\" value %d\n",objnum++);
    /**/
  }
  
  if (nreg) {
    if (ntets) {
      fprintf(fout,"object \"tets\" class field\n");
      fprintf(fout,"component \"positions\" value 1\n");
      fprintf(fout,"component \"connections\" value %d\n",objnum++);
      if (haveModel)
	fprintf(fout,"component \"data\" value %d\n",objnum++);
    }

    if (npyrs) {
      fprintf(fout,"object \"pyramids\" class field\n");
      fprintf(fout,"component \"positions\" value 1\n");
      fprintf(fout,"component \"connections\" value %d\n",objnum++);
      if (haveModel)
	fprintf(fout,"component \"data\" value %d\n",objnum++);
    }

    if (nprisms) {
      fprintf(fout,"object \"prisms\" class field\n");
      fprintf(fout,"component \"positions\" value 1\n");
      fprintf(fout,"component \"connections\" value %d\n",objnum++);
      if (haveModel)
	fprintf(fout,"component \"data\" value %d\n",objnum++);
    }

    if (nhexes) {
      fprintf(fout,"object \"hexes\" class field\n");
      fprintf(fout,"component \"positions\" value 1\n");
      fprintf(fout,"component \"connections\" value %d\n",objnum++);
      if (haveModel)
	fprintf(fout,"component \"data\" value %d\n",objnum++);
    }    
  }

  fprintf(fout,"end\n");
  fclose(fout);
}

