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

  int MESH_InitFromFile(Mesh_ptr mesh, const char *filename, MSTK_Comm comm) {
  FILE *fp;
  char inp_rtype[16], temp_str[256], fltype_str[16], rltype_str[16];
  char attname[256], atttype_str[256], attent_str[256];
  int i, j, found, NV=0, NE=0, NF=0, NR=0, nav, nar, gdim, gid;
  int id, dim, vid1, vid2, eid, fid, rid, adjvid, adjrid, adjv_flag;
  int nfv, max_nfv=0, nfe, max_nfe=0, nrv, max_nrv=0, nrf, max_nrf=0;
  int *fedirs, *rfdirs, ival, ncomp, nent, done, status;
  int processed_vertices=0, processed_adjv=0, processed_edges=0;
  int processed_faces=0, processed_regions=0, processed_adjr=0;
  double ver, xyz[3], rval;
  double *rval_arr = NULL;
  MVertex_ptr mv, ev1, ev2, adjv, *fverts, *rverts;
  MEdge_ptr me, *fedges;
  MFace_ptr mf, *rfaces;
  MRegion_ptr mr, adjr;
  MEntity_ptr ent;
  MAttrib_ptr attrib;
  MType attent;
  MAttType atttype;
  RepType file_reptype;

  if (!(fp = fopen(filename,"r"))) {
    sprintf(temp_str,"Cannot open file %s",filename);
    MSTK_Report("MESH_InitFromFile",temp_str,MSTK_ERROR);
    return 0;
  }

  status = fscanf(fp,"%s %lf",temp_str,&ver);
  if (strcmp(temp_str,"MSTK") != 0) {
    MSTK_Report("MESH_InitFromFile","Not a MSTK file",MSTK_ERROR);
    fclose(fp);
    return 0;
  }
  if (status == EOF)
    MSTK_Report("MESH_InitFromFile",
		"Premature end of file before any mesh data is read",MSTK_FATAL);

  if (ver != MSTK_VER) {
    MSTK_Report("MESH_InitFromFile","Version mismatch",MSTK_WARN);
  }

  status = fscanf(fp,"%s %d %d %d %d\n",inp_rtype,
		  &(mesh->nv),&(mesh->ne),&(mesh->nf),&(mesh->nr));
  if (status == EOF)
    MSTK_Report("MESH_InitFromFile",
		"Premature end of file before any mesh data is read",MSTK_FATAL);


  found = 0;
  for (i = 0; i < MSTK_MAXREP; i++) {
    if (strncmp(inp_rtype,MESH_rtype_str[i],2) == 0) {
      file_reptype = MESH_rtype[i];
      found = 1;

      if (mesh->reptype == UNKNOWN_REP)
	MESH_SetRepType(mesh,file_reptype); /* function does additional things
					       that a simple assignement to
					       mesh->reptype does not */
      break;
    }
  }

  if (!found) {
    MSTK_Report("MESH_InitFromFile","Unrecognized representation type",MSTK_ERROR);
    fclose(fp);
    return 0;
  }


  /* For now, the reduced representations only allow region and
     face elements - no edge elements are allowed */

  if (file_reptype >= R1 && file_reptype <= R4) {
    if (mesh->ne)
      MSTK_Report("Mesh_InitFromFile",
		  "Representation does not allow edges to be specified",MSTK_WARN);
  }


  status = fscanf(fp,"%s",temp_str);
  if (status == EOF)
    MSTK_Report("MESH_InitFromFile",
		"Premature end of file while looking for vertex data",MSTK_FATAL);
  else if (status == 0)
    MSTK_Report("MESH_InitFromFile",
		"Error in reading vertex data",MSTK_FATAL);

  if (strncmp(temp_str,"vertices",8) == 0) {

    NV = mesh->nv;
    mesh->mvertex = List_New(NV);

    for (i = 0; i < NV; i++) {
      status = fscanf(fp,"%lf %lf %lf %d %d",&xyz[0],&xyz[1],&xyz[2],&gdim,&gid);
      if (status == EOF)
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file while reading vertices",MSTK_FATAL);
      else if (status == 0)
	MSTK_Report("MESH_InitFromFile",
		    "Error in reading vertex data",MSTK_FATAL);

      mv = MV_New(mesh);
      MV_Set_Coords(mv,xyz);
      
      MV_Set_GEntID(mv,gid);
      MV_Set_GEntDim(mv,gdim);
      
      MV_Set_ID(mv,(i+1));
    }

    processed_vertices = 1;
  }
  else {
    MSTK_Report("MESH_InitFromFile",
		"Vertex information should be listed first",MSTK_ERROR);
    fclose(fp);
    return 0;
  }


  status = fscanf(fp,"%s",temp_str);
  if (status == EOF) {
    if (mesh->ne || mesh->nf || mesh->nr)
      MSTK_Report("MESH_InitFromFile",
		  "Premature end of file after vertex data",MSTK_FATAL);
    else
      return 0;
  }



  /* ADJACENT VERTEX DATA */

  if (strncmp(temp_str,"adjvertices",11) == 0) {
    adjv_flag = 1;
    if (mesh->reptype == R2 || mesh->reptype == R4) {

      /* Mesh representation supports storing of adjacent vertex info */

      for (i = 0; i < NV; i++) {
	mv = List_Entry(mesh->mvertex, i);
	status = fscanf(fp,"%d",&nav);
	if (status == EOF)
	  MSTK_Report("MESH_InitFromFile",
		      "Premature end of file while reading adjacent vertices",
		      MSTK_FATAL);
	else if (status == 0)
	  MSTK_Report("MESH_InitFromFile",
		      "Error in reading adjacent vertex data",MSTK_FATAL);

	for (j = 0; j < nav; j++) {
	  status = fscanf(fp,"%d",&adjvid);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading adjacent vertices"
			,MSTK_FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading adjacent vertex data",MSTK_FATAL);

	  adjv = List_Entry(mesh->mvertex,adjvid-1);
#ifdef DEBUG
	  if (MV_ID(adjv) != adjvid)
	    MSTK_Report("MESH_InitFromFile","Adjacent vertex ID mismatch",
			MSTK_ERROR);
#endif
	  
	  MV_Add_AdjVertex(mv,adjv);
	}
      }
    }
    else {

      /* Mesh representation does not support storing of adjacent
	 vertex information. Just read and discard */

      for (i = 0; i < NV; i++) {
	status = fscanf(fp,"%d",&nav);
	if (status == EOF)
	  MSTK_Report("MESH_InitFromFile",
		      "Premature end of file while reading adjacent vertices",
		      MSTK_FATAL);
	else if (status == 0)
	  MSTK_Report("MESH_InitFromFile",
		      "Error in reading adjacent vertex data",MSTK_FATAL);
	
	for (j = 0; j < nav; j++)
	  status = fscanf(fp,"%d",&adjvid);
      }
    }

    processed_adjv = 1;
  }
  else {
    adjv_flag = 0;
    if (file_reptype == R2 || file_reptype == R4) {
      MSTK_Report("MESH_InitFromFile",
		  "Expected adjacent vertex information",MSTK_ERROR);
    }
  }
    

  if (processed_adjv) {
    status = fscanf(fp,"%s",temp_str);
    if (status == EOF) {
      if (mesh->ne || mesh->nf || mesh->nr)
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file after adjacent vertex data",MSTK_FATAL);
      else {
        fclose(fp);
	return 1;
      }
    }
  }



  /* EDGE DATA */

  if (strncmp(temp_str,"edges",5) == 0) {

    NE = mesh->ne;
    mesh->medge = List_New(NE);

    if (mesh->reptype >= F1 && mesh->reptype <= F4) {

      /* Mesh representation supports edges */

      for (i = 0; i < NE; i++) {
	status = fscanf(fp,"%d %d %d %d",&vid1,&vid2,&gdim,&gid);
	if (status == EOF)
	  MSTK_Report("MESH_InitFromFile",
		      "Premature end of file while reading edges",MSTK_FATAL);
	else if (status == 0)
	  MSTK_Report("MESH_InitFromFile",
		      "Error in reading edge data",MSTK_FATAL);

	ev1 = List_Entry(mesh->mvertex,vid1-1);
	ev2 = List_Entry(mesh->mvertex,vid2-1);
#ifdef DEBUG
	if (MV_ID(ev1) != vid1 || MV_ID(ev2) != vid2)
	  MSTK_Report("MESH_InitFromFile","Mesh vertex ID mismatch",MSTK_ERROR);
#endif

	me = ME_New(mesh);

	ME_Set_Vertex(me,0,ev1);
	ME_Set_Vertex(me,1,ev2);

	ME_Set_GEntID(me,gid);
	ME_Set_GEntDim(me,gdim);

	ME_Set_ID(me, i+1);
      }
    }
    else {

      /* Mesh representation does not support edges. Read and discard */

      for (i = 0; i < NE; i++) {
	status = fscanf(fp,"%d %d %d %d",&vid1,&vid2,&gdim,&gid);
	if (status == EOF)
	  MSTK_Report("MESH_InitFromFile",
		      "Premature end of file while reading edges",MSTK_FATAL);
	else if (status == 0)
	  MSTK_Report("MESH_InitFromFile",
		      "Error in reading edge data",MSTK_FATAL);
      }
    }

    processed_edges = 1;
  }
  else {
    if (file_reptype >= F1 && file_reptype <= F4) {
      MSTK_Report("MESH_InitFromFile","Expected edge information",MSTK_ERROR);
      return 0;
    }
  }
  
  if (processed_edges) {
    status = fscanf(fp,"%s",temp_str);
    if (status == EOF) {
      if (mesh->nf || mesh->nr) {
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file after edge data",MSTK_FATAL);
      }
      else {
        fclose(fp);
	return 1;
      }
    }
  }



  /* FACE DATA */

  /* It is possible for someone to specify faces for a solid mesh specified 
   in the R1 format - do we have to idiot-proof? */

  if (strncmp(temp_str,"face",4) == 0) {
    if (strncmp(temp_str,"faces",5) != 0) 
      MSTK_Report("MESH_InitFromFile","Expected keyword \"faces\"",MSTK_ERROR);
    
    NF = mesh->nf;
    mesh->mface = List_New(NF);

    status = fscanf(fp,"%s",fltype_str);
    if (status == EOF)
      MSTK_Report("MESH_InitFromFile",
		  "Premature end of file while reading faces",MSTK_FATAL);

    if (file_reptype >= R1 && file_reptype <= R4) {

      if (strncmp(fltype_str,"vertex",6) == 0) {      

	/* Mesh representation supports faces for surface meshes */

	fverts = NULL;
	for (i = 0; i < NF; i++) {
	  mf = MF_New(mesh);

	  status = fscanf(fp,"%d",&nfv);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading faces",MSTK_FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading face data",MSTK_FATAL);

	  if (fverts) {
	    if (nfv > max_nfv) {
	      max_nfv = nfv;
	      fverts = realloc(fverts,max_nfv*sizeof(MVertex_ptr));
	    }
	  }
	  else {
	    max_nfv = nfv;
	    fverts = malloc(nfv*sizeof(MVertex_ptr));
	  }
	  
	  for (j = 0; j < nfv; j++) {
	    status = fscanf(fp,"%d",&vid1);
	    if (status == EOF)
	      MSTK_Report("MESH_InitFromFile",
			  "Premature end of file while reading face data",
			  MSTK_FATAL);
	    else if (status == 0)
	      MSTK_Report("MESH_InitFromFile",
			  "Error in reading edge data",MSTK_FATAL);

	    fverts[j] = List_Entry(mesh->mvertex,vid1-1);
#ifdef DEBUG
	    if (MV_ID(fverts[j]) != vid1)
	      MSTK_Report("MESH_InitFromFile","Mesh vertex ID mismatch",MSTK_ERROR);
#endif
	  }
	  MF_Set_Vertices(mf,nfv,fverts);  

	  status = fscanf(fp,"%d %d",&gdim,&gid);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading face data",MSTK_FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading face data",MSTK_FATAL);
	  
	  MF_Set_GEntID(mf,gid);
	  MF_Set_GEntDim(mf,gdim);
	  
	  MF_Set_ID(mf,i+1);
	}
	if (fverts)
	  free(fverts);

	processed_faces = 1;
      }
      else {
	MSTK_Report("MESH_InitFromFile",
		    "Expected face description in terms of vertices",MSTK_ERROR);
        fclose(fp);
	return 0;
      }

    }
    else if (file_reptype >= F1 && file_reptype <= F4) { 

      if (strncmp(fltype_str,"edge",4) == 0) {

	fedges = NULL;	fedirs = NULL;
	for (i = 0; i < NF; i++) {
	  mf = MF_New(mesh);

	  status = fscanf(fp,"%d",&nfe);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading face edge data",
			MSTK_FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading face data",MSTK_FATAL);

	  if (fedges) {
	    if (nfe > max_nfe) {
	      max_nfe = nfe;
	      fedges = realloc(fedges,max_nfe*sizeof(MVertex_ptr));
	      fedirs = realloc(fedirs,max_nfe*sizeof(int));
	    }
	  }
	  else {
	    max_nfe = nfe;
	    fedges = malloc(max_nfe*sizeof(MVertex_ptr));
	    fedirs = malloc(max_nfe*sizeof(MVertex_ptr));
	  }
	  
	  for (j = 0; j < nfe; j++) {
	    status = fscanf(fp,"%d",&eid);
	    if (status == EOF)
	      MSTK_Report("MESH_InitFromFile",
			  "Premature end of file while reading face edge data",
			  MSTK_FATAL);
	    else if (status == 0)
	      MSTK_Report("MESH_InitFromFile",
			  "Error in reading face data",MSTK_FATAL);

	    fedirs[j] = eid > 0 ? 1 : 0;
	    fedges[j] = List_Entry(mesh->medge,abs(eid)-1);
#ifdef DEBUG
	    if (ME_ID(fedges[j]) != abs(eid))
	      MSTK_Report("MESH_InitFromFile","Mesh edge ID mismatch",MSTK_ERROR);
#endif
	  }
	  
	  if (mesh->reptype >= F1 && mesh->reptype <= F4) {

	    MF_Set_Edges(mf,nfe,fedges,fedirs);

	  }
	  else {

	    MSTK_Report("MESH_InitFromFile","Need to convert edge list to vertex list",MSTK_ERROR);
	    fclose(fp);
	    return 0;

	  }
	  
	  status = fscanf(fp,"%d %d",&gdim,&gid);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading faces",MSTK_FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading face data",MSTK_FATAL);
	  
	  MF_Set_GEntDim(mf,gdim);
	  MF_Set_GEntID(mf,gid);
	  
	  MF_Set_ID(mf,i+1);
	}
	if (fedges) {
	  free(fedges);
	  free(fedirs);
	}
	
	processed_faces = 1;
      }
      else {
	MSTK_Report("MESH_InitFromFile",
		    "Expect face description in terms of edges",MSTK_ERROR);
	return 0;
      }
    }
  }
  else {
    if (file_reptype >= F1 && file_reptype <= F4) {
      MSTK_Report("MESH_InitFromFile","Expected face information",MSTK_ERROR);
      fclose(fp);
      return 0;
    }
  }
  
  
  if (processed_faces) {
    status = fscanf(fp,"%s",temp_str);
    if (status == EOF) {
      if (mesh->nr != 0)
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file after face data",MSTK_FATAL);
      else {
        fclose(fp);
	return 1;
      }
    }
    else if (status == 0)
      MSTK_Report("MESH_InitFromFile",
		  "Error in reading region data",MSTK_FATAL);
  }


  /* REGION DATA */

  if (strncmp(temp_str,"region",6) == 0) {
    if (strncmp(temp_str,"regions",7) != 0) 
      MSTK_Report("MESH_InitFromFile","Expected keyword \"regions\"",MSTK_ERROR);

    NR = mesh->nr;
    mesh->mregion = List_New(NR);

    status = fscanf(fp,"%s",rltype_str);
    if (status == EOF)
      MSTK_Report("MESH_InitFromFile",
		  "Premature end of file while reading regions",MSTK_FATAL);
    else if (status == 0)
      MSTK_Report("MESH_InitFromFile",
		  "Error in reading region data",MSTK_FATAL);

    if (file_reptype == R1 || file_reptype == R2) {
      if (strncmp(rltype_str,"vertex",6) == 0) {
	rverts = NULL;
	for (i = 0; i < NR; i++) {
	  mr = MR_New(mesh);
	  
	  status = fscanf(fp,"%d",&nrv);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading region data"
			,MSTK_FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading region data",MSTK_FATAL);

	  if (rverts) {
	    if (nrv > max_nrv) {
	      max_nrv = nrv;
	      rverts = realloc(rverts,max_nrv*sizeof(MVertex_ptr));
	    }
	  }
	  else {
	    max_nrv = nrv;
	    rverts = malloc(nrv*sizeof(MVertex_ptr));
	  }
	  
	  for (j = 0; j < nrv; j++) {
	    status = fscanf(fp,"%d",&vid1);
	    if (status == EOF)
	      MSTK_Report("MESH_InitFromFile",
		      "Premature end of file while reading region data",
			  MSTK_FATAL);
	    else if (status == 0)
	      MSTK_Report("MESH_InitFromFile",
			  "Error in reading region data",MSTK_FATAL);

	    rverts[j] = List_Entry(mesh->mvertex,vid1-1);
#ifdef DEBUG
	    if (MV_ID(rverts[j]) != vid1)
	      MSTK_Report("MESH_InitFromFile","Mesh vertex ID mismatch",MSTK_ERROR);
#endif
	  }
	  MR_Set_Vertices(mr,nrv,rverts,0,NULL);  
	  
	  status = fscanf(fp,"%d %d",&gdim,&gid);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading regions",MSTK_FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading region data",MSTK_FATAL);
	  
	  MR_Set_GEntDim(mr,gdim);

	  MR_Set_GEntID(mr,gid);
	  
	  MR_Set_ID(mr,i+1);
	}
	if (rverts)
	  free(rverts);

	processed_regions = 1;
      }
      else {
	MSTK_Report("MESH_InitFromFile",
		    "Expected region description in terms of faces",MSTK_FATAL);
        fclose(fp);
	return 0;
      }
    }
    else {
      if (strncmp(rltype_str,"face",4) == 0) {
	rfaces = NULL;
	rfdirs = NULL;
	for (i = 0; i < NR; i++) {
	  mr = MR_New(mesh);
	  
	  status = fscanf(fp,"%d",&nrf);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading region data",
			MSTK_FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading region data",MSTK_FATAL);
	  
	  if (rfaces) {
	    if (nrf > max_nrf) {
	      max_nrf = nrf;
	      rfaces = realloc(rfaces,max_nrf*sizeof(MFace_ptr));
	      rfdirs = realloc(rfdirs,max_nrf*sizeof(int));
	    }
	  }
	  else {
	    max_nrf = nrf;
	    rfaces = malloc(nrf*sizeof(MFace_ptr));
	    rfdirs = malloc(nrf*sizeof(int));
	  }
	  
	  for (j = 0; j < nrf; j++) {
	    status = fscanf(fp,"%d",&fid);
	    if (status == EOF)
	      MSTK_Report("MESH_InitFromFile",
			  "Premature end of file while reading region data",
			  MSTK_FATAL);
	    else if (status == 0)
	      MSTK_Report("MESH_InitFromFile",
			  "Error in reading region data",MSTK_FATAL);
	    
	    rfdirs[j] = fid > 0 ? 1 : 0;
	    rfaces[j] = List_Entry(mesh->mface,abs(fid)-1);
#ifdef DEBUG
	    if (MF_ID(rfaces[j]) != abs(fid))
	      MSTK_Report("MESH_InitFromFile","Mesh face ID mismatch",MSTK_ERROR);
#endif
	  }

	  if (mesh->reptype != R1 && mesh->reptype != R2) {

	    MR_Set_Faces(mr,nrf,rfaces,rfdirs);

	  }
	  else {

	    MSTK_Report("MESH_InitFromFile","Need to convert face list to vertex list", MSTK_ERROR);
	    fclose(fp);
	    return 0;

	  }
	  
	  status = fscanf(fp,"%d %d",&gdim,&gid);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			  "Premature end of file while reading region",MSTK_FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading region data",MSTK_FATAL);
	  
	  MR_Set_GEntDim(mr,gdim);
	  
	  MR_Set_GEntID(mr,gid);
	  
	  MR_Set_ID(mr,i+1);
	}
	if (rfaces) {
	  free(rfaces);
	  free(rfdirs);
	}
	
	processed_regions = 1;
      }
      else {
	MSTK_Report("MESH_InitFromFile",
		    "Expected region description in terms of vertices",MSTK_ERROR);
	fclose(fp);
	return 0;
      }
    }
  }

  if (processed_regions) {
    status = fscanf(fp,"%s",temp_str);
    if (status == EOF) {
      if (mesh->reptype == R2)
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file after reading regions",MSTK_FATAL);
      else {
        fclose(fp);
	return 1;
      }
    }
  }


  /* ADJACENT REGION DATA */

  if (strncmp(temp_str,"adjregions",10) == 0) {
    if (mesh->reptype == R2) {
      for (i = 0; i < NR; i++) {
	mr = List_Entry(mesh->mregion,i);

	status = fscanf(fp,"%d",&nar);
	if (status == EOF)
	  MSTK_Report("MESH_InitFromFile",
		      "Premature end of file while reading adjacent regions",
		      MSTK_FATAL);
	else if (status == 0)
	  MSTK_Report("MESH_InitFromFile",
		      "Error in reading adjacent region data",MSTK_FATAL);

	for (j = 0; j < nar; j++) {
	  status = fscanf(fp,"%d",&adjrid);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading ajdacent regions",
			MSTK_FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading adjacent region data",MSTK_FATAL);
	  
	  if (adjrid == 0)
	    continue;
	  adjr = List_Entry(mesh->mregion,adjrid-1);
#ifdef DEBUG
	  if (MR_ID(adjr) != adjrid)
	    MSTK_Report("MESH_InitFromFile",
			"Adjacent region ID mismatch",MSTK_ERROR);
#endif

	  MR_Add_AdjRegion(mr,j,adjr);
	}
      }

      processed_adjr = 1;
    }
    else {
      for (i = 0; i < NR; i++) {
	status = fscanf(fp,"%d",&nar);
	if (status == EOF)
	  MSTK_Report("MESH_InitFromFile",
		      "Premature end of file while reading adjacent regions",
		      MSTK_FATAL);
	else if (status == 0)
	  MSTK_Report("MESH_InitFromFile",
		      "Error in reading adjacent region data",MSTK_FATAL);
	
	for (j = 0; j < nar; j++) {
	  status = fscanf(fp,"%d",&rid);
	  if (status == EOF)
	    MSTK_Report("MESH_InitFromFile",
			"Premature end of file while reading adjacent regions",
			MSTK_FATAL);
	  else if (status == 0)
	    MSTK_Report("MESH_InitFromFile",
			"Error in reading adjacent region data",MSTK_FATAL);
	}
      }

      processed_adjr = 1;
    }
  }

  if (processed_adjr) {
    status = fscanf(fp,"%s",temp_str);
    if (status == EOF) {
      if (mesh->reptype == R2)
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file after reading regions",MSTK_FATAL);
      else {
        fclose(fp);
	return 1;
      }
    }
  }

  /* ATTRIBUTE DATA */

  if (strncmp(temp_str,"attributes",10) == 0) {

    done = 0;
    while (!done) {
      status = fscanf(fp,"%s",attname);
      if (status == EOF) {
	done = 1;
	continue;
      }
      else if (status == 0)
	MSTK_Report("MESH_InitFromFile",
		    "Error in reading attribute data",MSTK_FATAL);

      status = fscanf(fp,"%s",atttype_str);
      if (status == EOF)
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file while reading attributes",MSTK_FATAL);
      else if (status == 0)
	MSTK_Report("MESH_InitFromFile",
		    "Error in reading attribute data",MSTK_FATAL);

      atttype = INT;
      if (strncmp(atttype_str,"INT",3) == 0)
	atttype = INT;
      else if (strncmp(atttype_str,"DOUBLE",6) == 0)
	atttype = DOUBLE;
      else if (strncmp(atttype_str,"POINTER",7) == 0)
	MSTK_Report("MESH_InitFromFile",
		    "Cannot specify POINTER attributes in file",MSTK_FATAL);
      else if (strncmp(atttype_str,"VECTOR",6) == 0)
	atttype = VECTOR;
      else if (strncmp(atttype_str,"TENSOR",6) == 0)
	atttype = TENSOR;
      else {
	sprintf(temp_str,"%-s not a recognized attribute type",atttype_str);
	MSTK_Report("MESH_InitFromFile",temp_str,MSTK_FATAL);
      }

      status = fscanf(fp,"%d",&ncomp);
      if (status == EOF)
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file while reading attributes",MSTK_FATAL);
      else if (status == 0)
	MSTK_Report("MESH_InitFromFile",
		    "Error in reading attribute data",MSTK_FATAL);

      if ((atttype == INT || atttype == DOUBLE) && (ncomp != 1)) 
	MSTK_Report("MESH_InitFromFile","Number of components should be 1 for attributes of type INT or DOUBLE",MSTK_WARN);
      else if ((atttype == VECTOR || atttype == TENSOR) && (ncomp == 0))
	MSTK_Report("MESH_InitFromFile","Number of components should be non-zero for attributes of type VECTOR or TENSOR",MSTK_FATAL);



      status = fscanf(fp,"%s",attent_str);
      if (status == EOF)
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file while reading attributes",MSTK_FATAL);
      else if (status == 0)
	MSTK_Report("MESH_InitFromFile",
		    "Error in reading attribute data",MSTK_FATAL);
      attent = MALLTYPE;
      if (strncmp(attent_str,"MVERTEX",7) == 0) 
	attent = MVERTEX;
      else if (strncmp(attent_str,"MEDGE",5) == 0)
	attent = MEDGE;
      else if (strncmp(attent_str,"MFACE",5) == 0)
	attent = MFACE;
      else if (strncmp(attent_str,"MREGION",7) == 0)
	attent = MREGION;
      else if (strncmp(attent_str,"MALLTYPE",8) == 0)
	attent = MALLTYPE;
      else {
	sprintf(temp_str,"%s not a recognized entity type for attributes",
		attent_str);
	MSTK_Report("MESH_InitFromFile",temp_str,MSTK_FATAL);
      }

      status = fscanf(fp,"%d",&nent);
      if (status == EOF)
	MSTK_Report("MESH_InitFromFile",
		    "Premature end of file while reading attributes",MSTK_FATAL);
      else if (status == 0)
	MSTK_Report("MESH_InitFromFile",
		    "Error in reading attribute data",MSTK_FATAL);

      if (nent < 1) 
	MSTK_Report("MESH_InitFromFile",
		    "Attribute applied on no entities?",MSTK_FATAL);

      if (atttype == INT || atttype == DOUBLE || atttype == POINTER)
	attrib = MAttrib_New(mesh,attname,atttype,attent);
      else
	attrib = MAttrib_New(mesh,attname,atttype,attent,ncomp);
      
      for (i = 0; i < nent; i++) {
	status = fscanf(fp,"%d %d",&dim,&id);
	if (status == EOF)
	  MSTK_Report("MESH_InitFromFile",
		      "Premature end of file while reading attributes",MSTK_FATAL);
	else if (status == 0)
	  MSTK_Report("MESH_InitFromFile",
		      "Error in reading attribute data",MSTK_FATAL);
	
	if (attent != dim && attent != MALLTYPE) {
	  MSTK_Report("MESH_InitFromFile",
		      "Attribute not applicable to this type of entity",MSTK_WARN);
	  if (atttype == INT)
	    status = fscanf(fp,"%d",&ival);
	  else
	    for (j = 0; j < ncomp; j++)
	      status = fscanf(fp,"%lf",&rval);
	}
	else {

	  if (atttype == INT)
	    status = fscanf(fp,"%d",&ival);
	  else if (atttype == DOUBLE)
	    status = fscanf(fp,"%lf",&rval);
	  else if (atttype == VECTOR || atttype == TENSOR) {
	    ival = ncomp;
	    rval_arr = (double *) malloc(ncomp*sizeof(double));
	    for (j = 0; j < ncomp; j++) {
	      status = fscanf(fp,"%lf",&(rval_arr[j]));
	    }
	  }

	  switch (dim) {
	  case MVERTEX:
	    ent = MESH_VertexFromID(mesh,id);
	    break;
	  case MEDGE:
	    ent = MESH_EdgeFromID(mesh,id);
	    break;
	  case MFACE:
	    ent = MESH_FaceFromID(mesh,id);
	    break;
	  case MREGION:
	    ent = MESH_RegionFromID(mesh,id);
	    break;
	  default:
	    ent = NULL;
	    MSTK_Report("MESH_InitFromFile","Invalid entity type",MSTK_FATAL);
	  }
	  
	  MEnt_Set_AttVal(ent,attrib,ival,rval,(void *)rval_arr);
	  
	}

      } /* for (i = 0; i < nent; i++) */

    } /* while (!done) */

  }

  fclose(fp);

  return 1;
}


#ifdef __cplusplus
}
#endif
