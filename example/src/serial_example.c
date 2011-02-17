#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MSTK.h"


/* Example program to read a mesh and print out some information about it */


int main(int argc, char *argv[]) {
  int i, k, len, suff, idx, idx2, ok, edir, nv, ne, natt, ival, ncomp;
  double xyz[3], rval, *rval_arr;
  void *pval;
  char mname[256], filename[256], attname[256];
  Mesh_ptr mesh;
  MVertex_ptr v;
  MEdge_ptr e;
  MFace_ptr f;
  MAttrib_ptr attrib;
  GEntity_ptr gent;
  List_ptr fedges;
  MType attentdim;
  MAttType atttype;

  
  fprintf(stderr,"\n");

  if (argc == 1) {
    fprintf(stderr,"Usage: %s testfile.mstk",argv[0]);
    exit(-1);
  }

  /* Initialize MSTK - Always do this even if it does not seem to
     matter in this version of MSTK */

  MSTK_Init();

  /* Load the mesh */

  strcpy(mname,argv[1]);
  len = strlen(mname);
  suff = 0;
  if (len > 5) {
    k = len-5;
    while (!suff && k > 1) {
      if (strncmp(&(mname[k]),".mstk",5) == 0) {
        suff = 1;
      }
      else
        k--;
    }
  }
  if (suff)
    mname[len-5] = '\0';

  strcpy(filename,mname);
  strcat(filename,".mstk");

  mesh = MESH_New(UNKNOWN_REP);
  ok = MESH_InitFromFile(mesh,filename);
  if (!ok) {
    fprintf(stderr,"Cannot open input file %s\n\n\n",filename);
    exit(-1);
  }


  /* Print some info about the mesh */


  nv = MESH_Num_Vertices(mesh);
  for (i = 0; i < nv; i++) {
    v = MESH_Vertex(mesh,i);

    /* Basic info */
    fprintf(stderr,"\n");
    fprintf(stderr,"Vertex: 0x%-x   ID: %-d   ",(unsigned int)v,MV_ID(v));
    
    
    /* Classification w.r.t. geometric model */
      
    if (MV_GEntDim(v) == -1)
      fprintf(stderr,"Unknown Classification\n");
    else {
      fprintf(stderr,"GEntID: %-d    GEntDim: %-d\n", MV_GEntID(v),
	      MV_GEntDim(v));
      if ((gent = MV_GEntity(v)))
	fprintf(stderr,"Model entity pointer: 0x%-x\n",(unsigned int) gent);
    }


    /* Coordinates */
    MV_Coords(v,xyz);
    fprintf(stderr,"Coordinates: %16.8lf %16.8lf %16.8lf\n",xyz[0],xyz[1],xyz[2]);
  }



  idx = 0;
  while ((f = MESH_Next_Face(mesh,&idx))) {

    /* Basic info */
    fprintf(stderr,"\n");
    fprintf(stderr,"Face: 0x%-x   ID: %-d   ",(unsigned int) f,MF_ID(f));
    
    
    /* Classification w.r.t. geometric model */

    if (MF_GEntDim(f) == -1)
      fprintf(stderr,"Unknown Classification\n");
    else {
      fprintf(stderr,"GEntID: %-d    GEntDim: %-d\n", MF_GEntID(f),
	      MF_GEntDim(f));
      if ((gent = MF_GEntity(f)))
	fprintf(stderr,"Model entity pointer: 0x%-x\n",(unsigned int) gent);
    }

    fprintf(stderr,"\n");


    /* Edges of face */
    fedges = MF_Edges(f,1,0);
    ne = List_Num_Entries(fedges);
    fprintf(stderr,"Edges: %-d\n",ne);
    fprintf(stderr,"Object        ID      GEntID   GEntDim    Vertex IDs\n");
    idx2 = 0; i = 0;
    while ((e = List_Next_Entry(fedges,&idx2))) {
      edir = MF_EdgeDir_i(f,i);
      if (edir) 
	fprintf(stderr,"0x%-8x    %-8d %-8d     %-1d       %-d  %-d\n",
		(unsigned int) e,ME_ID(e),ME_GEntID(e),ME_GEntDim(e),
		MV_ID(ME_Vertex(e,0)),MV_ID(ME_Vertex(e,1)));
      else 
	fprintf(stderr,"0x%-8x    %-8d %-8d     %-1d       %-d  %-d\n",
		(unsigned int) e,-ME_ID(e),ME_GEntID(e),ME_GEntDim(e),
		MV_ID(ME_Vertex(e,0)),MV_ID(ME_Vertex(e,1)));
      i++;
    }
    fprintf(stderr,"\n");
    List_Delete(fedges);
  }


  if ((natt = MESH_Num_Attribs(mesh))) {

    fprintf(stderr,"Attributes on the mesh:\n\n");
    
    for (i = 0; i < natt; i++) {

      attrib = MESH_Attrib(mesh,i);
      
      attentdim = MAttrib_Get_EntDim(attrib);

     /* Won't print out edge, face and region based attributes
        but the code should be quite similar */

      if (attentdim != MVERTEX) continue;

      fprintf(stderr,"Attribute Number %-d\n",i+1);

      MAttrib_Get_Name(attrib,attname);
      fprintf(stderr,"Name: %-s\n",attname);

      atttype = MAttrib_Get_Type(attrib);
      switch(atttype) {
      case INT:
	fprintf(stderr,"Type: Integer\n");
	break;
      case DOUBLE:
	fprintf(stderr,"Type: Double\n");
	break;
      case VECTOR:
	fprintf(stderr,"Type: Vector\n");
	break;
      case TENSOR:
	fprintf(stderr,"Type: Tensor\n");
	break;
      default:
	fprintf(stderr,"Unrecognizable or unprintable attribute type\n");
	continue;	
      }

      attentdim = MAttrib_Get_EntDim(attrib);
      switch(attentdim) {
      case MVERTEX:
	fprintf(stderr,"Applicable to vertices only\n");
	break;
      case MEDGE:
      case MFACE:
      case MREGION:
      case MALLTYPE:
      default:
	fprintf(stderr,"Unrecognized entity type\n");
	continue;
      }

      ncomp = MAttrib_Get_NumComps(attrib);
      fprintf(stderr,"Number of components: %-d\n",ncomp);

      switch(attentdim) {
      case MVERTEX:
	idx = 0;
	while ((v = MESH_Next_Vertex(mesh,&idx))) {
	  if (MEnt_Get_AttVal(v,attrib,&ival,&rval,&pval)) {
	    fprintf(stderr,"V %-d: ",MV_ID(v));
	    switch (atttype) {
	    case INT:
	      fprintf(stderr," %-d\n",ival);
	      break;
	    case DOUBLE: 
	      fprintf(stderr," %-lf ",rval);
	      break;
	    case VECTOR: case TENSOR:
	      rval_arr = (double *) pval;
	      for (k = 0; k < ncomp; k++)
		fprintf(stderr," %-lf ",rval_arr[k]);
	      break;
	    default:
	      break;
	    }
	    fprintf(stderr,"\n");
	  }
	}
	break;
      case MEDGE:
      case MFACE:
      case MREGION: 
      case MALLTYPE:
	/* Won't print out edge, face and region based attributes
	   but the code should be quite similar */
	continue;
      default:
	break;
      } /* switch (attentdim) */

    } /* for (i = 0; i < natt) */
    
  } /* if (Mesh_Num_Attribs(mesh)) */
  


  /* Write out a copy of the mesh in the F1 format */

  strcpy(filename,mname);
  strcat(filename,"-copy.mstk");
  MESH_WriteToFile(mesh,filename,F1);

  
  /* Deleting of mesh is not necessary if this the end of the program */

  MESH_Delete(mesh);

  return 1;
}
