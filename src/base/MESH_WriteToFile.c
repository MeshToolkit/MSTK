/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#define _H_Mesh_Private

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Mesh.h"
#include "MSTK.h"
#include "MSTK_private.h"


#ifdef __cplusplus
extern "C" {
#endif

  extern RepType MESH_rtype[5];
  extern char MESH_rtype_str[5][3];

  /* File name should have the .mstk extension - will not check here */
  /* MSTK_Comm argument can be NULL only if MPI is not enabled        */


  int MESH_WriteToFile(Mesh_ptr mesh, const char *filename, RepType rtype, MSTK_Comm comm) {
  FILE *fp;
  char mesg[80], attname[256];
  int i, j, k, idx;
  int gdim, gid;
  int mvid, mvid0, mvid1, mvid2, mrid2, meid, mfid, mrid;
  int nav, nar, nfe, nfv, nrf, nrv, dir=0;
  int nv, ne, nf, nr;
  int natt, ncomp, ival, nent;
  double xyz[3], rval, rdummy, *rval_arr;
  void *pval, *pdummy;
  MVertex_ptr mv, mv0, mv1, mv2;
  MEdge_ptr me;
  MFace_ptr mf;
  MRegion_ptr mr, mr2;
  List_ptr adjverts, mfedges, mfverts, mrfaces, mrverts, adjregs;
  RepType reptype;
  MAttrib_ptr attrib, vidatt, eidatt, fidatt, ridatt;
  MType attentdim;
  MAttType atttype;

  char modfilename[256];
  strcpy(modfilename, filename);
  
  int rank = 0, numprocs = 1;
#ifdef MSTK_HAVE_MPI
  if (comm) {
    MPI_Comm_size((MPI_Comm)comm, &numprocs);
    MPI_Comm_rank((MPI_Comm)comm, &rank);
  }
  if (numprocs > 1) {
    int ndigits = 0;
    int div = 1;
    while (numprocs/div) {div *= 10; ndigits++;}
    sprintf(modfilename,"%s.%d.%0*d",filename,numprocs,ndigits,rank);
  }
#endif
  
  if (!(fp = fopen(modfilename,"w"))) {
    sprintf(mesg,"Cannot open file %-s for writing",modfilename);
    MSTK_Report("MESH_WriteToFile",mesg,MSTK_ERROR);
    return 0;
  }

  if (rtype != UNKNOWN_REP) {
    reptype = rtype;
  }
  else {
    reptype = MESH_RepType(mesh);
  }

  nv = MESH_Num_Vertices(mesh);
  ne = MESH_Num_Edges(mesh);
  nf = MESH_Num_Faces(mesh);
  nr = MESH_Num_Regions(mesh);

  fprintf(fp,"MSTK %-2.1lf\n",MSTK_FILE_VER);
  fprintf(fp,"%s %d %d %d %d\n",
	  MESH_rtype_str[reptype], 
	  nv, 
	  (reptype >= R1 && reptype <= R4)?0:ne, 
	  (reptype >= R1 && reptype <= R2 && nr)?0:nf, 
	  nr);

  vidatt = MAttrib_New(mesh,"vidatt",INT,MVERTEX);
  eidatt = MAttrib_New(mesh,"eidatt",INT,MEDGE);
  fidatt = MAttrib_New(mesh,"fidatt",INT,MFACE);
  ridatt = MAttrib_New(mesh,"ridatt",INT,MREGION);

  idx = 0; i = 0;
  while ((mv = MESH_Next_Vertex(mesh,&idx)))
    MEnt_Set_AttVal(mv,vidatt,++i,0.0,NULL);

  idx = 0; i = 0;
  while ((me = MESH_Next_Edge(mesh,&idx)))
    MEnt_Set_AttVal(me,eidatt,++i,0.0,NULL);

  idx = 0; i = 0;
  while ((mf = MESH_Next_Face(mesh,&idx)))
    MEnt_Set_AttVal(mf,fidatt,++i,0.0,NULL);

  idx = 0; i = 0;
  while ((mr = MESH_Next_Region(mesh,&idx)))
    MEnt_Set_AttVal(mr,ridatt,++i,0.0,NULL);
  
  
  fprintf(fp,"vertices\n");
  idx = 0;
  while ((mv = MESH_Next_Vertex(mesh,&idx))) {

    MV_Coords(mv,xyz);

    gdim = MV_GEntDim(mv);
    gid = MV_GEntID(mv);

    fprintf(fp,"%24.16lf %24.16lf %24.16lf   %d %d\n",
	    xyz[0],xyz[1],xyz[2],gdim,gid);
    
  }

  if (reptype == R2 || reptype == R4) {
    fprintf(fp,"adjvertices\n");

    idx = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx))) {

      nav = MV_Num_AdjVertices(mv);
      fprintf(fp,"%d ",nav);
      
      adjverts = MV_AdjVertices(mv);
      for (j = 0; j < nav; j++) {
	mv2 = List_Entry(adjverts,j);
	MEnt_Get_AttVal(mv2,vidatt,&mvid2,&rval,&pval);
	fprintf(fp,"%d ",mvid2);
      }
      fprintf(fp,"\n");
      List_Delete(adjverts);
    }
  }



  if (reptype <= F4 && ne) {
    fprintf(fp,"edges\n");

    idx = 0;
    while ((me = MESH_Next_Edge(mesh,&idx))) {

      mv0 = ME_Vertex(me,0);
      MEnt_Get_AttVal(mv0,vidatt,&mvid0,&rval,&pval);
      mv1 = ME_Vertex(me,1);
      MEnt_Get_AttVal(mv1,vidatt,&mvid1,&rval,&pval);

      gdim = ME_GEntDim(me);
      gid = ME_GEntID(me);

      fprintf(fp,"%d %d \t%d %d\n",mvid0,mvid1,gdim,gid);
    }
  }



  if (reptype <= F4) {

    /* For full representations, always write out faces in terms of edges */

    fprintf(fp,"faces edge\n");
    
    idx = 0;
    while ((mf = MESH_Next_Face(mesh,&idx))) {
      
      nfe = MF_Num_Edges(mf);
      fprintf(fp,"%d ",nfe);
      
      mfedges = MF_Edges(mf,1,0);
      for (j = 0; j < nfe; j++) {
	me = List_Entry(mfedges,j);
	dir = MF_EdgeDir_i(mf,j);
	MEnt_Get_AttVal(me,eidatt,&meid,&rval,&pval);
	if (dir != 1) meid = -meid;
	fprintf(fp,"%d ",meid);
      }
      List_Delete(mfedges);
      
      gdim = MF_GEntDim(mf);
      /*
	gent = MF_GEntity(mf);
	gid = gent ? -99 : 0;
      */
      gid = MF_GEntID(mf);
      
      fprintf(fp,"\t%d %d\n",gdim,gid);
    }
  }
  else {

    /* For reduced representations, R3 and R4 always write out faces
       in terms of vertices. For reduced representations, R1 and R2
       write out faces in terms of vertices only when there are no
       regions (i.e. faces are the highest level mesh entities) */

    if ((reptype > R2) || (nr == 0)) {

      fprintf(fp,"faces vertex\n");

      idx = 0;
      while ((mf = MESH_Next_Face(mesh,&idx))) {
	
	nfv = MF_Num_Edges(mf);
	fprintf(fp,"%d ",nfv);
	
	mfverts = MF_Vertices(mf,1,0);
	for (j = 0; j < nfv; j++) {
	  mv = List_Entry(mfverts,j);
	  MEnt_Get_AttVal(mv,vidatt,&mvid,&rval,&pval);
	  fprintf(fp,"%d ",mvid);
	}
	List_Delete(mfverts);

	gdim = MF_GEntDim(mf);
	gid = MF_GEntID(mf);
	
	fprintf(fp,"\t%d %d\n",gdim,gid);
      }
    }
	
  }


  if (nr) {
    if (reptype <= F4 || reptype >= R2) {
      fprintf(fp,"regions face\n");

      idx = 0;
      while ((mr = MESH_Next_Region(mesh,&idx))) {

	nrf = MR_Num_Faces(mr);
	fprintf(fp,"%d ",nrf);

	mrfaces = MR_Faces(mr);
	for (j = 0; j < nrf; j++) {
	  mf = List_Entry(mrfaces,j);
	  dir = MR_FaceDir_i(mr,j);
	  MEnt_Get_AttVal(mf,fidatt,&mfid,&rval,&pval);
	  if (dir != 1) mfid = -mfid;
	  fprintf(fp,"%d ",mfid);
	}
	List_Delete(mrfaces);
	
	gdim = MF_GEntDim(mr);
	gid = MR_GEntID(mr);

	fprintf(fp,"\t%d %d\n",gdim,gid);
      }
    }
    else {
      fprintf(fp,"regions vertex\n");

      idx = 0;
      while ((mr = MESH_Next_Region(mesh,&idx))) {

	nrv = MR_Num_Vertices(mr);
	fprintf(fp,"%d ",nrv);

	mrverts = MR_Vertices(mr);
	for (j = 0; j < nrv; j++) {
	  mv = List_Entry(mrverts,j);
	  MEnt_Get_AttVal(mv,vidatt,&mvid,&rval,&pval);
	  fprintf(fp,"%d ",mvid);
	}
	List_Delete(mrverts);
	
	gdim = MR_GEntDim(mr);
	gid = MR_GEntID(mr);

	fprintf(fp,"\t%d %d\n",gdim,gid);
      }
    }

    if (reptype == R2 || reptype == R4) {
      fprintf(fp,"adjregions\n");
      
      idx = 0;
      while ((mr = MESH_Next_Region(mesh,&idx))) {

	nar = MR_Num_Faces(mr);
	fprintf(fp,"%d ",nar);

	adjregs = MR_AdjRegions(mr);

	for (j = 0; j < nar; j++) {
	  mr2 = List_Entry(adjregs,j);
	  if ((long) mr2 == -1) 
	    fprintf(fp,"%d ",0);
	  else {
	    MEnt_Get_AttVal(mr2,ridatt,&mrid2,&rval,&pval);
	    fprintf(fp,"%d ",mrid2);
	  }
	}
	fprintf(fp,"\n");
	List_Delete(adjregs);
      }
    }
  }


  /* Write out attributes if there are more than the 4 that we created 
    in this routine */


  if ((natt = MESH_Num_Attribs(mesh)) > 4) {

    fprintf(fp,"attributes\n");

    for (i = 0; i < natt; i++) {
      
      attrib = MESH_Attrib(mesh,i);

      /* Don't write out attribs we created for the internal use of 
	 this routine */
      if (attrib == vidatt || attrib == eidatt || attrib == fidatt || 
	  attrib == ridatt) continue;
      
      MAttrib_Get_Name(attrib,attname);

      atttype = MAttrib_Get_Type(attrib);
      if (atttype == POINTER) continue;  /* cannot write it out */

      ncomp = MAttrib_Get_NumComps(attrib);

      attentdim = MAttrib_Get_EntDim(attrib);


      /* First count how many entities actually have the attribute assigned */

      nent = 0;
      switch(attentdim) {
      case MVERTEX:
	idx = 0;
	while ((mv = MESH_Next_Vertex(mesh,&idx)))
	  if (MEnt_Get_AttVal(mv,attrib,&ival,&rval,&pval)) nent++;
	break;
      case MEDGE:
	idx = 0;
	while ((me = MESH_Next_Edge(mesh,&idx)))
	  if (MEnt_Get_AttVal(me,attrib,&ival,&rval,&pval)) nent++;
	break;
      case MFACE:
	idx = 0;
	while ((mf = MESH_Next_Face(mesh,&idx)))
	  if (MEnt_Get_AttVal(mf,attrib,&ival,&rval,&pval)) nent++;	    
	break;
      case MREGION: 
	idx = 0;
	while ((mr = MESH_Next_Region(mesh,&idx)))
	  if (MEnt_Get_AttVal(mr,attrib,&ival,&rval,&pval)) nent++;
	break;
      case MALLTYPE:
	idx = 0;
	while ((mv = MESH_Next_Vertex(mesh,&idx)))
	  if (MEnt_Get_AttVal(mv,attrib,&ival,&rval,&pval)) nent++;
	idx = 0;
	while ((me = MESH_Next_Edge(mesh,&idx)))
	  if (MEnt_Get_AttVal(me,attrib,&ival,&rval,&pval)) nent++;
	idx = 0;
	while ((mf = MESH_Next_Face(mesh,&idx)))
	  if (MEnt_Get_AttVal(mf,attrib,&ival,&rval,&pval)) nent++;	    
	idx = 0;
	while ((mr = MESH_Next_Region(mesh,&idx)))
	  if (MEnt_Get_AttVal(mr,attrib,&ival,&rval,&pval)) nent++;
	break;	
      default:
	break;
      } /* switch (attentdim) */


      /* No point in writing out attribute if no entity uses it! Or is there? */

      if (!nent) continue;



      fprintf(fp,"%-s\n",attname);

      switch(atttype) {
      case INT:
	fprintf(fp,"INT\n");
	break;
      case DOUBLE:
	fprintf(fp,"DOUBLE\n");
	break;
      case VECTOR:
	fprintf(fp,"VECTOR\n");
	break;
      case TENSOR:
	fprintf(fp,"TENSOR\n");
	break;
      default:
	MSTK_Report("MESH_WriteToFile",
		    "Unrecognizable or unprintable attribute type\n",MSTK_WARN);
	continue;	
      }

      fprintf(fp,"%-d\n",ncomp);

      switch(attentdim) {
      case MVERTEX:
	fprintf(fp,"MVERTEX\n");
	break;
      case MEDGE:
	fprintf(fp,"MEDGE\n");
	break;
      case MFACE:
	fprintf(fp,"MFACE\n");
	break;
      case MREGION:
	fprintf(fp,"MREGION\n");
	break;
      case MALLTYPE:
	fprintf(fp,"MALLTYPE\n");
	break;
      default:
	MSTK_Report("Mesh_WriteToFile","Unrecognized entity type",MSTK_WARN);
	break;
      }

      fprintf(fp,"%-d\n",nent);


      switch(attentdim) {
      case MVERTEX:
	idx = 0;
	while ((mv = MESH_Next_Vertex(mesh,&idx))) {
	  if (MEnt_Get_AttVal(mv,attrib,&ival,&rval,&pval)) {
	    MEnt_Get_AttVal(mv,vidatt,&mvid,&rdummy,&pdummy);
	    fprintf(fp,"0 %-d ",mvid);
	    switch (atttype) {
	    case INT:
	      fprintf(fp," %-d",ival);
	      break;
	    case DOUBLE: 
	      fprintf(fp," %-lf ",rval);
	      break;
	    case VECTOR: case TENSOR:
	      rval_arr = (double *) pval;
	      for (k = 0; k < ncomp; k++)
		fprintf(fp," %-lf ",rval_arr[k]);
	      break;
	    default:
	      break;
	    }
	    fprintf(fp,"\n");
	  }
	}
	break;
      case MEDGE:
	idx = 0;
	while ((me = MESH_Next_Edge(mesh,&idx))) {
	  if (MEnt_Get_AttVal(me,attrib,&ival,&rval,&pval)) {
	    MEnt_Get_AttVal(me,eidatt,&meid,&rdummy,&pdummy);
	    fprintf(fp,"1 %-d ",meid);
	    switch (atttype) {
	    case INT:
	      fprintf(fp," %-d",ival);
	      break;
	    case DOUBLE: 
	      fprintf(fp," %-lf ",rval);
	      break;
	    case VECTOR: case TENSOR:
	      rval_arr = (double *) pval;
	      for (k = 0; k < ncomp; k++)
		fprintf(fp," %-lf ",rval_arr[k]);
	      break;
	    default:
	      break;
	    }
	    fprintf(fp,"\n");
	  }
	}
	break;
      case MFACE:
	idx = 0;
	while ((mf = MESH_Next_Face(mesh,&idx))) {
	  if (MEnt_Get_AttVal(mf,attrib,&ival,&rval,&pval)) {
	    MEnt_Get_AttVal(mf,fidatt,&mfid,&rdummy,&pdummy);
	    fprintf(fp,"2 %-d ",mfid);
	    switch (atttype) {
	    case INT:
	      fprintf(fp," %-d",ival);
	      break;
	    case DOUBLE: 
	      fprintf(fp," %-lf ",rval);
	      break;
	    case VECTOR: case TENSOR:
	      rval_arr = (double *) pval;
	      for (k = 0; k < ncomp; k++)
		fprintf(fp," %-lf ",rval_arr[k]);
	      break;
	    default:
	      break;
	    }
	    fprintf(fp,"\n");
	  }
	}
	break;

      case MREGION: 
	idx = 0;
	while ((mr = MESH_Next_Region(mesh,&idx))) {
	  if (MEnt_Get_AttVal(mr,attrib,&ival,&rval,&pval)) {
	    MEnt_Get_AttVal(mr,ridatt,&mrid,&rdummy,&pdummy);
	    fprintf(fp,"3 %-d ",mrid);
	    switch (atttype) {
	    case INT:
	      fprintf(fp," %-d",ival);
	      break;
	    case DOUBLE: 
	      fprintf(fp," %-lf ",rval);
	      break;
	    case VECTOR: case TENSOR:
	      rval_arr = (double *) pval;
	      for (k = 0; k < ncomp; k++)
		fprintf(fp," %-lf ",rval_arr[k]);
	      break;
	    default:
	      break;
	    }
	    fprintf(fp,"\n");
	  }
	}
	break;

      case MALLTYPE:
	idx = 0;
	while ((mv = MESH_Next_Vertex(mesh,&idx))) {
	  if (MEnt_Get_AttVal(mv,attrib,&ival,&rval,&pval)) {
	    MEnt_Get_AttVal(mv,vidatt,&mvid,&rdummy,&pdummy);
	    fprintf(fp,"0 %-d ",mvid);
	    switch (atttype) {
	    case INT:
	      fprintf(fp," %-d",ival);
	      break;
	    case DOUBLE: 
	      fprintf(fp," %-lf ",rval);
	      break;
	    case VECTOR: case TENSOR:
	      rval_arr = (double *) pval;
	      for (k = 0; k < ncomp; k++)
		fprintf(fp," %-lf ",rval_arr[k]);
	      break;
	    default:
	      break;
	    }
	    fprintf(fp,"\n");
	  }
	}
	idx = 0;
	while ((me = MESH_Next_Edge(mesh,&idx))) {
	  if (MEnt_Get_AttVal(me,attrib,&ival,&rval,&pval)) {
	    MEnt_Get_AttVal(me,eidatt,&meid,&rdummy,&pdummy);
	    fprintf(fp,"1 %-d ",meid);
	    switch (atttype) {
	    case INT:
	      fprintf(fp," %-d",ival);
	      break;
	    case DOUBLE: 
	      fprintf(fp," %-lf ",rval);
	      break;
	    case VECTOR: case TENSOR:
	      rval_arr = (double *) pval;
	      for (k = 0; k < ncomp; k++)
		fprintf(fp," %-lf ",rval_arr[k]);
	      break;
	    default:
	      break;
	    }
	    fprintf(fp,"\n");
	  }
	}
	idx = 0;
	while ((mf = MESH_Next_Face(mesh,&idx))) {
	  if (MEnt_Get_AttVal(mf,attrib,&ival,&rval,&pval)) {
	    MEnt_Get_AttVal(mf,fidatt,&mfid,&rdummy,&pdummy);
	    fprintf(fp,"2 %-d ",mfid);
	    switch (atttype) {
	    case INT:
	      fprintf(fp," %-d",ival);
	      break;
	    case DOUBLE: 
	      fprintf(fp," %-lf ",rval);
	      break;
	    case VECTOR: case TENSOR:
	      rval_arr = (double *) pval;
	      for (k = 0; k < ncomp; k++)
		fprintf(fp," %-lf ",rval_arr[k]);
	      break;
	    default:
	      break;
	    }
	    fprintf(fp,"\n");
	  }
	}
	idx = 0;
	while ((mr = MESH_Next_Region(mesh,&idx))) {
	  if (MEnt_Get_AttVal(mr,attrib,&ival,&rval,&pval)) {
	    MEnt_Get_AttVal(mr,ridatt,&mrid,&rdummy,&pdummy);
	    fprintf(fp,"3 %-d ",mrid);
	    switch (atttype) {
	    case INT:
	      fprintf(fp," %-d",ival);
	      break;
	    case DOUBLE: 
	      fprintf(fp," %-lf ",rval);
	      break;
	    case VECTOR: case TENSOR:
	      rval_arr = (double *) pval;
	      for (k = 0; k < ncomp; k++)
		fprintf(fp," %-lf ",rval_arr[k]);
	      break;
	    default:
	      break;
	    }
	    fprintf(fp,"\n");
	  }
	}
	break;	
      default:
	break;
      } /* switch (attentdim) */

    } /* for (i = 0; i < natt) */
    
  } /* if (Mesh_Num_Attribs(mesh)) */
  

  idx = 0; i = 0;
  while ((mv = MESH_Next_Vertex(mesh,&idx)))
    MEnt_Rem_AttVal(mv,vidatt);

  idx = 0; i = 0;
  while ((me = MESH_Next_Edge(mesh,&idx)))
    MEnt_Rem_AttVal(me,eidatt);

  idx = 0; i = 0;
  while ((mf = MESH_Next_Face(mesh,&idx)))
    MEnt_Rem_AttVal(mf,fidatt);

  idx = 0; i = 0;
  while ((mr = MESH_Next_Region(mesh,&idx)))
    MEnt_Rem_AttVal(mr,ridatt);
  
  MAttrib_Delete(vidatt);
  MAttrib_Delete(eidatt);
  MAttrib_Delete(fidatt);
  MAttrib_Delete(ridatt);




  fclose(fp);

  return 1;
}

#ifdef __cplusplus
}
#endif
