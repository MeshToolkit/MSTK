#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif

int MESH_ImportFromFLAGX3D(Mesh_ptr mesh, const char *filename, const int rank, const int numprocs) {

  char funcname[256] = "MESH_ImportFromFLAGX3D";
  char mesg[256], temp_str[1028], keyword[1028];

  FILE *fp;

  int i, j, status, endheader, nfe, nrf, done, value, matid;
  int ndim, nmats, nverts, nfaces, nelems, num_ghost_nodes, num_slaved_nodes;
  int nodes_per_slave, nodes_per_face, faces_per_cell, node_data_fields, cell_data_fields;
  int cellid, nodeid, edgeid, faceid, *medir, *mfdir;

  MVertex_ptr *mv;
  MEdge_ptr *me;
  MFace_ptr *mf, mface;
  MRegion_ptr mregion;

  char (*matnames)[256], (*mateos)[256], (*matopc)[256];

  
  if (!(fp = fopen(filename,"r"))) {
    sprintf(mesg,"Cannot open file %s for reading\n",filename);
    MSTK_Report(funcname,mesg,MSTK_FATAL);
  }

  

  /* Confirm identifying string */

  status = fscanf(fp,"%s %s",keyword,temp_str);

  if (status == EOF)
    MSTK_Report(funcname,"Premature end of file",MSTK_FATAL);
    
  if (strncmp(keyword,"x3dtoflag",9) != 0)
    MSTK_Report(funcname,"Not an X3D file",MSTK_FATAL);


  /* Read header keyword for meta-data */
    
  status = fscanf(fp,"%s",keyword);
  if (status == EOF)
    MSTK_Report(funcname,"Premature end of file",MSTK_FATAL);

  if (strncmp(keyword,"header",6) != 0)
    MSTK_Report(funcname,"Expected keyword 'header'",MSTK_FATAL);

  done = 0; 
  i = 0;
  while (!done) {
 
    /* Read the meta data */

    status = fscanf(fp,"%s",keyword);
    if (status == EOF)
      MSTK_Report(funcname,"Premature end of file",MSTK_FATAL);
  
    if (strncmp(keyword,"end_header",10) == 0)
      done = 1;
    else {

      status = fscanf(fp,"%d",&value);
      if (status == EOF)
        MSTK_Report(funcname,"Premature end of file",MSTK_FATAL);

      if (strcmp(keyword,"process") == 0) {
	if (value != rank+1) {
	  sprintf(mesg,"Reading file for processor %d on processor %d",value-1,rank);
	  MSTK_Report(funcname,mesg,MSTK_FATAL);
	}
      }
      else if (strcmp(keyword,"numdim") == 0) {
	ndim = value;
      }
      else if (strcmp(keyword,"materials") == 0) {
	nmats = value;
        matnames = (char (*)[256]) malloc(nmats*sizeof(char [256]));
        mateos   = (char (*)[256]) malloc(nmats*sizeof(char [256]));
        matopc   = (char (*)[256]) malloc(nmats*sizeof(char [256]));
      }
      else if (strcmp(keyword,"nodes") == 0) {
	nverts = value;
        mv = (MVertex_ptr *) malloc(nverts*sizeof(MVertex_ptr));
      }
      else if (strcmp(keyword,"faces") == 0) {
	nfaces = value;
        if (ndim == 2) {
          me = (MEdge_ptr *) malloc(nfaces*sizeof(MEdge_ptr));
          medir = (int *) malloc(nfaces*sizeof(int));
        }
        else if (ndim == 3) {
          mf = (MFace_ptr *) malloc(nfaces*sizeof(MFace_ptr));
          mfdir = (int *) malloc(nfaces*sizeof(int));
        }
      }
      else if (strcmp(keyword,"elements") == 0) {
	nelems = value;
      }
      else if (strcmp(keyword,"ghost_nodes") == 0) {
	num_ghost_nodes = value;
      }
      else if (strcmp(keyword,"slaved_nodes") == 0) {
	num_slaved_nodes = value;
      }
      else if (strcmp(keyword,"faces_per_cell") == 0) {
	faces_per_cell = value;
      }
      else if (strcmp(keyword,"node_data_fields") == 0) {
	node_data_fields = value;
      }
      else if (strcmp(keyword,"cell_data_fields") == 0) {
	cell_data_fields = value;
      }

      i++;

      if (!done && i == 15) {
	MSTK_Report(funcname,"Cannot find keyword 'end_header'",MSTK_ERROR);
        break;
      }
      
    }

  } /* while (!end_header) */

  
  /* Read the rest of the file */

  done = 0;
  while (!done) {

    status = fscanf(fp,"%s",keyword);
    if (status == EOF)
      MSTK_Report(funcname,"Premature end of file",MSTK_FATAL);

    if (strncmp(keyword,"matnames",8) == 0) {
      for (i = 0; i < nmats; i++) {
	status = fscanf(fp,"%d %s",&matid,temp_str);
        if (status == EOF)
          MSTK_Report(funcname,"Premature end of file",MSTK_FATAL);

	if (matid != -1)
	  strcpy(matnames[matid-1],temp_str);
      }      
    }
    else if (strncmp(keyword,"mateos",6) == 0) {
      for (i = 0; i < nmats; i++) {
	status = fscanf(fp,"%d %s",&matid,temp_str);
        if (status == EOF)
          MSTK_Report(funcname,"Premature end of file",MSTK_FATAL);

	if (matid != -1)
	  strcpy(mateos[matid-1],temp_str);
      }      
    }
    else if (strncmp(keyword,"matopc",6) == 0) {
      for (i = 0; i < nmats; i++) {
	status = fscanf(fp,"%d %s",&matid,temp_str);
        if (status == EOF)
          MSTK_Report(funcname,"Premature end of file",MSTK_FATAL);

	if (matid != -1)
	  strcpy(matopc[matid-1],temp_str);
      }      
    }
    else if (strncmp(keyword,"nodes",5) == 0) {

      for (i = 0; i < nverts; i++) {
        double xyz[3];
	status = fscanf(fp,"%d %lf %lf %lf",&nodeid,&(xyz[0]),&(xyz[1]),&(xyz[2]));
        if (status == EOF)
          MSTK_Report(funcname,"Premature end of file",MSTK_FATAL);

	mv[nodeid-1] = MV_New(mesh);
	MV_Set_Coords(mv[nodeid-1],xyz);
      }

    }
    else if (strncmp(keyword,"end_nodes",9) == 0) {
      /* nothing to do */
    }
    else if (strncmp(keyword,"faces",5) == 0) {

      int owner_pid, fnbr_pid, fnbr_lid, dummy;
      
      if (ndim == 2) {
        int nev, evids[2];
        MVertex_ptr everts[2];

	for (i = 0; i < nfaces; i++) {
	  status = fscanf(fp,"%d %d",&edgeid,&nev); /* nev = 2 */
          if (status == EOF)
            MSTK_Report(funcname,"Premature end of file",MSTK_FATAL);
	  
	  for (j = 0; j < nev; j++) {
	    status = fscanf(fp,"%d",&(evids[j]));	    
	    everts[j] = mv[evids[j]-1];
	  }
	  
	  status = fscanf(fp,"%d %d %d",&owner_pid, &fnbr_pid, &fnbr_lid);
          if (status == EOF)
            MSTK_Report(funcname,"Premature end of file",MSTK_FATAL);

	  status = fscanf(fp,"%d %d %d %d %d",&dummy,&dummy,&dummy,&dummy,&dummy);
          if (status == EOF)
            MSTK_Report(funcname,"Premature end of file",MSTK_FATAL);

          /* check if the duplicate of this edge was already
             created. If it was, then just use it in the opposite
             direction */

          if (owner_pid == rank+1 && fnbr_pid == rank+1 && edgeid > fnbr_lid) {
            me[edgeid-1] = me[fnbr_lid-1];
            medir[edgeid-1] = 0;
          }
          else {
            me[edgeid-1] = ME_New(mesh);
            ME_Set_Vertex(me[edgeid-1],0,everts[0]);
            ME_Set_Vertex(me[edgeid-1],1,everts[1]);
            medir[edgeid-1] = 1;
          }
	}
      }
      else {
        int nfv, fvids[MAXPV2];
        MVertex_ptr fverts[MAXPV2];

	for (i = 0; i < nfaces; i++) {
	  status = fscanf(fp,"%d %d",&faceid,&nfv);
          if (status == EOF)
            MSTK_Report(funcname,"Premature end of file",MSTK_FATAL);
	  
	  for (j = 0; j < nfv; j++) {
	    status = fscanf(fp,"%d",&(fvids[j]));
            fverts[j] = mv[fvids[j]-1];
          }


          status = fscanf(fp,"%d %d %d",&owner_pid, &fnbr_pid, &fnbr_lid);
          if (status == EOF)
            MSTK_Report(funcname,"Premature end of file",MSTK_FATAL);
          
          status = fscanf(fp,"%d %d %d %d %d",&dummy,&dummy,&dummy,&dummy,&dummy);
          if (status == EOF)
            MSTK_Report(funcname,"Premature end of file",MSTK_FATAL);
          
          /* check if the duplicate of this face was already created. If it was, then
             just use it in the opposite direction */
          
          if (owner_pid == rank+1 && fnbr_pid == rank+1 && faceid > fnbr_lid) {
            mf[faceid-1] = mf[fnbr_lid-1];
            mfdir[faceid-1] = 0;
          }
          else {
            mf[faceid-1] = MF_New(mesh);
            MF_Set_Vertices(mf[faceid-1],nfv,fverts);
            mfdir[faceid-1] = 1;
          }
        }
      }
    }
    else if (strncmp(keyword,"end_faces",9) == 0) {
      /* nothing to do */
    }
    else if (strncmp(keyword,"cells",5) == 0) {

      /* faces are duplicated in the x3d specification. So a cell does
         not share a face with another cell but has its own copy - so,
         we create the cells and follow up by merging the duplicate
         faces */

      if (ndim == 2) {

        int feids[MAXPF3], fedirs[MAXPF3];
        MEdge_ptr fedges[MAXPF3];        
        
        for (i = 0; i < nelems; i++) {
          status = fscanf(fp,"%d %d",&cellid,&nfe);
          if (status == EOF)
            MSTK_Report(funcname,"Premature end of file",MSTK_FATAL);

          for (j = 0; j < nfe; j++) {
            fscanf(fp,"%d",&(feids[j]));
            fedges[j] = me[feids[j]-1];
            fedirs[j] = medir[feids[j]-1];              
          }
          
          mface = MF_New(mesh);
          MF_Set_Edges(mface,nfe,fedges,fedirs);          
        }

        

      }
      else if (ndim == 3) {
        
        int rfids[MAXPF3], rfdirs[MAXPF3];
        MFace_ptr rfaces[MAXPF3];        
        
        for (i = 0; i < nelems; i++) {
          status = fscanf(fp,"%d %d",&cellid,&nrf);
          if (status == EOF)
            MSTK_Report(funcname,"Premature end of file",MSTK_FATAL);

          for (j = 0; j < nrf; j++) {
            fscanf(fp,"%d",&(rfids[j]));
            rfaces[j] = mf[rfids[j]-1];
            rfdirs[j] = mfdir[rfids[j]-1];
          }
          
          mregion = MR_New(mesh);
          MR_Set_Faces(mregion,nrf,rfaces,rfdirs);          
        }
      }

    }
    else if (strncmp(keyword,"end_cells",9) == 0) {
      /* nothing to do */
    }
    else if (strncmp(keyword,"slaved_nodes",12) == 0) {

    }
    else if (strncmp(keyword,"ghost_nodes",11) == 0) {

    }
    else if (strncmp(keyword,"cell_data",9) == 0) {

    }
    else if (strncmp(keyword,"end_cell_data",13) == 0) {
    }
    else if (strncmp(keyword,"node_data",9) == 0) {

    }
    else if (strncmp(keyword,"end_node_data",13) == 0) {
    }
    else if (strncmp(keyword,"end_dump",8) == 0) {
      done = 1;
    }

  }


  if (ndim == 2)
    free(medir);
  else if (ndim == 3)
    free(mfdir);

  free(matnames);
  free(mateos);
  free(matopc);
 
  return 1;
}

#ifdef __cplusplus
}
#endif
