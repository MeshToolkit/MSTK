#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MSTK.h"


double Tet_Volume(double (*rxyz)[3]);
double PR_Volume(double (*rxyz)[3], int n, int **rfverts, int *nfv, 
		 int nf, int *star_shaped);


int main(int argc, char **argv) {
  int len, ok, firstwarn1=1,firstwarn2;
  int i, j, idx, dir, nrf, nrv, nfv[MAXPF3];
  int rid, nbad, status, **rfverts=NULL;
  int first=1, star_shaped;
  double vol, rxyz[MAXPV3][3];
  double rcen[3], fcen[3], fxyz[MAXPV2][3], txyz[4][3], tvol;
  int k;
  char infname[256], gmvfname[256], mname[256];
  Mesh_ptr mesh;
  MRegion_ptr mr;
  MFace_ptr rface;
  List_ptr rverts, rfaces, fverts;
  MAttrib_ptr valatt;

  

  switch(argc) {
  case 2:
    strcpy(infname,argv[1]);
    break;
  case 3:
    strcpy(infname,argv[1]);
    break;
  default:
    fprintf(stderr,"Usage: %s infname<.mstk>\n",argv[0]);
    exit(-1);
  }


  len = strlen(infname);
  if (len > 5 && (strncmp(&(infname[len-5]),".mstk",5) == 0)) {
    strcpy(mname,infname);
    mname[len-5] = '\0';
  }
  else
    strcpy(mname,infname);

  strcpy(gmvfname,mname);
  strcat(gmvfname,"-chk.gmv");


  MSTK_Init();

  fprintf(stderr,"\n\nChecking mesh %s\n\n",infname);

  mesh = MESH_New(UNKNOWN_REP);
  ok = MESH_InitFromFile(mesh,infname);
  if (!ok) {
    fprintf(stderr,"Cannot open input file %s\n",infname);
    exit(-1);
  }


  fprintf(stderr,"Checking topological validity of mesh.....");
  status = MESH_CheckTopo(mesh);
  if (status == 1) 
    fprintf(stderr,"OK\n");
  else
    fprintf(stderr,"FAILED\n");


  if (MESH_Num_Regions(mesh) == 0) {
    fprintf(stderr,"Geometric validity checks present only for solid meshes\n");
  }

  fprintf(stderr,"Checking geometric validity of elements....");
  valatt = MAttrib_New(mesh,"valatt",INT,MALLTYPE);

  idx = 0;
  nbad = 0;
  while ((mr = MESH_Next_Region(mesh,&idx))) {
    
    rid = MR_ID(mr);

    rverts = MR_Vertices(mr);
    nrv = List_Num_Entries(rverts);

    for (i = 0; i < nrv; i++) {
      MVertex_ptr rv = List_Entry(rverts,i);
      MV_Coords(rv,rxyz[i]);
    }

    if (nrv == 4 && MR_Num_Faces(mr) == 4) {
      vol = Tet_Volume(rxyz);
      if (vol <= 0.0) {
	if (firstwarn1) {
	  fprintf(stderr,"\n\nElement is invalid - ID %d\n",rid);
	  fprintf(stderr,"Volume = %lf\n",vol);
	  MR_Print(mr,3);
	  fprintf(stderr,"\n\n\n Select/Color elements by cell field \"valatt\" in %s to see all bad elements\n\n\n",gmvfname);

	  firstwarn1 = 0;
          status = 0;
	}
	nbad++;
	MEnt_Set_AttVal(mr,valatt,-1,0.0,NULL);
      }
    }
    else {

      if (first) {
	rfverts = (int **) malloc(MAXPF3*sizeof(int *));
	for (i = 0; i < MAXPF3; i++)
	  rfverts[i] = (int *) malloc(MAXPV2*sizeof(int));
      }

      rfaces = MR_Faces(mr);
      nrf = List_Num_Entries(rfaces);
      
      for (i = 0; i < nrf; i++) {
	rface = List_Entry(rfaces,i);
	dir = MR_FaceDir_i(mr,i);

	fverts = MF_Vertices(rface,!dir,0);
	nfv[i] = List_Num_Entries(fverts);

	for (j = 0; j < nfv[i]; j++) {
	  MVertex_ptr fv = List_Entry(fverts,j);
	  rfverts[i][j] = List_Locate(rverts,fv);
	}

	List_Delete(fverts);
      }
      List_Delete(rfaces);
      

      
      vol = PR_Volume(rxyz, nrv, rfverts, nfv, nrf, &star_shaped);

      if (vol <= 0.0) {
	if (firstwarn1) {
	  fprintf(stderr,"\n\nElement is invalid - ID %d\n",rid);
	  fprintf(stderr,"Volume = %lf\n",vol);	  
	  MR_Print(mr,3);
	  fprintf(stderr,"\n\n\n Select/Color elements by cell field \"valatt\" in %s to see all bad elements\n\n\n",gmvfname);
	  firstwarn1 = 0;
	  status = 0;
	}

	nbad++;
	MEnt_Set_AttVal(mr,valatt,-1,0.0,NULL);
      } /* if (vol <= 0.0) */	    
      else if (!star_shaped) {
	if (firstwarn2) {
	  fprintf(stderr,"\n\nElement is valid but is not star shaped - ID %d\n",rid);
	  fprintf(stderr,"Volume = %lf\n",vol);

	  fprintf(stderr,"Invalid tet in decomposition of polyhedral element.\n");
	  fprintf(stderr,"Perhaps a polyhedral face is severely distorted?\n");
	  MR_Print(mr,3);
	  fprintf(stderr,"\n\n\n Select/Color elements by cell field \"valatt\" in %s to see all bad elements\n\n\n",gmvfname);
	  firstwarn2 = 0;
	  status = 0;
	}

	nbad++;
	MEnt_Set_AttVal(mr,valatt,1,0.0,NULL);
      }
	
    } /* if (MR_NumFaces(mr) == 4) ... else ... */

    List_Delete(rverts);

  } /* while ((mr = MESH_Next_Region)) */

  if (rfverts) {
    for (i = 0; i < MAXPF3; i++)
      free(rfverts[i]);
    free(rfverts);
  }


  if (status)
     fprintf(stderr,"OK\n");

  if (status) {
    fprintf(stderr,"\n\n");
    fprintf(stderr,"***************\n");
    fprintf(stderr,"Mesh is A-OK!!\n");
    fprintf(stderr,"***************\n");
    fprintf(stderr,"\n\n");
    exit(0);
  }
  else {
    if (nbad) {
      fprintf(stderr,"Total number of bad elements %-d\n\n",nbad);

      fprintf(stderr,"Tagging bad elements in GMV file %s\n\n",gmvfname);
      MESH_ExportToGMV(mesh,gmvfname,0,NULL,NULL);
      exit(-1);
    }
  }

}


  
