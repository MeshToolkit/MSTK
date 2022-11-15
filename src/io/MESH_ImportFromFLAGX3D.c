/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "MSTK.h"
#include "MSTK_private.h"


#ifdef __cplusplus
extern "C" {
#endif

  /* Import a single or distributed mesh from a FLAG X3D format file 

   NOTE: For parallel import, code does not guarantee that of a group of
   coincident vertices on processor boundaries, the vertex tagged as
   master in the FLAG X3D file will remain the master. Rather, the
   vertex on the lowest rank processor will be tagged as the master */

  /* The input arguments parallel_opts indicate if we should build
     parallel adjacencies and some options for controlling that 
     
     parallel_opts[0] = 1/0  ---- Weave/Don't weave the distributed meshes
                                  together to form parallel connections
     parallel_opts[1] = N    ---- Number of ghost layers around mesh
  */                            


  void trimLeadingWhiteSpace(char *string) {
    int b = 0, l = strlen(string);
    while (b < l && string[b] == ' ')
      b++;

    int lnew = l-b;
    memcpy(string, string+b, lnew);
    string[lnew] = '\0';
  }

  int MESH_ImportFromFLAGX3D(Mesh_ptr mesh, const char *filename, int *parallel_opts, MSTK_Comm comm) {

  char funcname[32] = "MESH_ImportFromFLAGX3D";
  char mesg[256], temp_str[1028], keyword[1028], modfilename[256];

  FILE *fp;

  int i, j, status, endheader, nfe, nrf, done, value, matid;
  int ndim=0, nmats=0, nverts=0, nfaces=0, nelems=0, num_ghost_nodes=0, num_slaved_nodes=0;
  int nodes_per_slave=0, nodes_per_face=0, faces_per_cell=0, node_data_fields=0, cell_data_fields=0;
  int cellid, nodeid, edgeid, faceid, *medir=NULL, *mfdir=NULL;

  MVertex_ptr *mv=NULL;
  MEdge_ptr *me=NULL;
  MFace_ptr *mf=NULL, mface;
  MRegion_ptr mregion;

  char (*matnames)[256]=NULL, (*mateos)[256]=NULL, (*matopc)[256]=NULL;

  int rank, numprocs;
  int distributed = 0;

  numprocs = 1;
  rank = 0;

#ifdef MSTK_HAVE_MPI
  MPI_Comm_size(comm,&numprocs);
  MPI_Comm_rank(comm,&rank);
#endif

  char *ext = strstr(filename,".x3d"); /* Search for the x3d extension */
  if (!ext)
    MSTK_Report(funcname,"FLAG X3D files must have .x3d extension",MSTK_FATAL);


  if (!(fp = fopen(filename,"r"))) {
    sprintf(mesg,"Cannot open file %s for reading", filename);
    MSTK_Report(funcname,mesg,MSTK_FATAL);
  }

#ifdef MSTK_HAVE_MPI
  /* Also figure out if we are really doing a parallel read */

  if (numprocs > 1) {
    char rankext[6];
    sprintf(rankext,"%05d",rank+1);
    if (strstr(filename,rankext))  /* does the filename have the rank in it? */
      distributed = 1;
  }

  if (distributed)
    MESH_Set_Prtn(mesh,rank,numprocs); /* necessary for allocation of 
                                          parallel adjacency arrays */
#endif


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
#ifdef MSTK_HAVE_MPI
	if (value != rank+1) {
	  sprintf(mesg,"Reading file for processor %d on processor %d",value-1,rank);
	  MSTK_Report(funcname,mesg,MSTK_FATAL);
	}
#endif
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
            status = fscanf(fp,"%d",&(feids[j]));
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
            status = fscanf(fp,"%d",&(rfids[j]));
            if (status == EOF)
              MSTK_Report(funcname,"Premature end of file",MSTK_FATAL);
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
      int dum, vid, pid, gvid;
      status = fscanf(fp,"%d",&dum);
      for (i = 0; i < num_ghost_nodes; i++) {
	status = fscanf(fp,"%d %d %d %d",&vid, &pid, &dum, &dum);
        if (status == EOF)
          MSTK_Report(funcname,"Premature end of file",MSTK_FATAL);
	pid--;

#ifdef MSTK_HAVE_MPI
        MVertex_ptr mv = MESH_Vertex(mesh,vid-1);
        if (pid != rank) {
          MV_Set_PType(mv, PGHOST);
          MV_Set_MasterParID(mv, pid);
          MESH_Flag_Has_Ghosts_From_Prtn(mesh,pid,MVERTEX);
          MESH_Flag_Has_Ghosts_From_Prtn(mesh,pid,MEDGE);
          if(ndim == 3) {
            MESH_Flag_Has_Ghosts_From_Prtn(mesh,pid,MFACE);
          }
        } else {
          MV_Set_PType(mv, POVERLAP);
          MV_Set_MasterParID(mv, rank);
        }
#endif

      }

#ifdef MSTK_HAVE_MPI
      if (distributed) {
        /* Compute which processors must receive overlap info from this
         * processor from the reverse info */
        MESH_Get_OverlapAdj_From_GhostAdj(mesh,comm);
      }
#endif

    }
    else if (strncmp(keyword,"cell_data",9) == 0) {
      /* In X3D files, cell data is always a scalar */
      int cellvarsdone = 0;
      while (!cellvarsdone) {
        status = fscanf(fp,"%s",temp_str);
        if (status == EOF) {
          MSTK_Report("MESH_ImportFromFLAGX3D",
                      "Premature end of file while reading cell data",MSTK_ERROR);
          return 0;
        }

        if (strncmp(temp_str,"end_cell_data",13) == 0)
          cellvarsdone = 1;
        else if (strncmp(temp_str,"matid",5) == 0) {
          /* matids are special cell data */
          /* there are 10 strings in a row, each taking up 10
           * characters, except in the last row where there are as
           * many as needed to complete the count. So we have to use
           * fgets instead of fscanf */
          int matids_done = 0;
          while (!matids_done) {
            if (ndim == 2) {
              int idx = 0;
              MFace_ptr mf;
              while ((mf = MESH_Next_Face(mesh, &idx))) {
                if (fgets(temp_str, 11, fp) == NULL || feof(fp)) {
                  MSTK_Report("MESH_ImportFromFLAGX3D",
                              "Premature end of file while reading matids",MSTK_ERROR);
                  return 0;
                }
                if (strlen(temp_str) == 1 && temp_str[0] == '\n') {
                  /* ignore newline and read again */
                  if (fgets(temp_str, 11, fp) == NULL || feof(fp)) {
                    MSTK_Report("MESH_ImportFromFLAGX3D",
                                "Premature end of file while reading matids",MSTK_ERROR);
                    return 0;
                  }
                }

                
                trimLeadingWhiteSpace(temp_str);
                
                /* Linear search but hopefully nmats is small */
                for (i = 0; i < nmats; i++)
                  if (strcmp(matnames[i],temp_str) == 0) {
                    MF_Set_GEntDim(mf, 2);
                    MF_Set_GEntID(mf, i+1);
                    break;
                  }
              }
            } else if (ndim == 3) {
              int idx = 0;
              MRegion_ptr mr;
              while ((mr = MESH_Next_Region(mesh, &idx))) {              
                if (fgets(temp_str, 11, fp) == NULL || feof(fp)) {
                  MSTK_Report("MESH_ImportFromFLAGX3D",
                              "Premature end of file while reading matids",MSTK_ERROR);
                  return 0;
                }
                if (strlen(temp_str) == 1 && temp_str[0] == '\n') {
                  /* ignore newline and read again */
                  if (fgets(temp_str, 11, fp) == NULL || feof(fp)) {
                    MSTK_Report("MESH_ImportFromFLAGX3D",
                                "Premature end of file while reading matids",MSTK_ERROR);
                    return 0;
                  }
                }
                                
                trimLeadingWhiteSpace(temp_str);
                
                /* Linear search but hopefully nmats is small */
                for (i = 0; i < nmats; i++)
                  if (strcmp(matnames[i],temp_str) == 0) {
                    MR_Set_GEntDim(mr, 3);
                    MR_Set_GEntID(mr, i+1);
                    break;
                  }
              }
            }
            status = fscanf(fp,"%s",temp_str);
            if (status == EOF) {
              MSTK_Report("MESH_ImportFromFLAGX3D",
                          "Premature end of file while reading matids",MSTK_ERROR);
              return 0;
            }
            
            if (strncmp(temp_str,"end_matid",9) == 0)
              matids_done = 1;
            else {
              char temp_str2[256];
              sprintf(temp_str2, "Expected end_matid but got %s", temp_str);
              MSTK_Report("MESH_ImportFromFLAGX3D",
                          temp_str2, MSTK_ERROR);
            }
          }
        }
        else {
          MAttrib_ptr matt;
          double rval;
          if (ndim == 2) {
            matt = MAttrib_New(mesh, temp_str, DOUBLE, MFACE);
            int idx = 0;
            MFace_ptr mf;
            while ((mf = MESH_Next_Face(mesh, &idx))) {
              status = fscanf(fp,"%lf",&rval);
              if (status == EOF) {
                MSTK_Report("MESH_ImportFromFLAGX3D",
                            "Premature end of file while reading cell data",
                            MSTK_ERROR);
                return 0;
              }
              MEnt_Set_AttVal(mf,matt,0,rval,NULL);
            }
          }
          else if (ndim == 3) {
            matt = MAttrib_New(mesh, temp_str, DOUBLE, MREGION);
            int idx = 0;
            MRegion_ptr mr;
            while ((mr = MESH_Next_Region(mesh, &idx))) {
              status = fscanf(fp,"%lf",&rval);
              if (status == EOF) {
                MSTK_Report("MESH_ImportFromFLAGX3D",
                            "Premature end of file while reading cell data",
                            MSTK_ERROR);
                return 0;
              }
              MEnt_Set_AttVal(mr,matt,0,rval,NULL);
            }
          }
          char temp_str1[256], temp_str2[256];
          strcpy(temp_str1, "end_");
          strcat(temp_str1, temp_str);
          status = fscanf(fp,"%s",temp_str2);
          if (status == EOF) {
            MSTK_Report("MESH_ImportFromFLAGX3D",
                        "Premature end of file while reading cell data",
                        MSTK_ERROR);
            return 0;
          } else if (strcmp(temp_str1, temp_str2) != 0) {
            char temp_str3[256];
            sprintf(temp_str3,"Expected %s but got %s while reading cell data",
                    temp_str1, temp_str2);
            MSTK_Report("MESH_ImportFromFLAGX3D", temp_str3, MSTK_FATAL);
          }
        }
      }
    }
    else if (strncmp(keyword,"end_cell_data",13) == 0) {
    }
    else if (strncmp(keyword,"node_data",9) == 0) {
      /* In X3D files, node data is always a vector */
      int nodevarsdone = 0;
      while (!nodevarsdone) {
        status = fscanf(fp,"%s",temp_str);
        if (status == EOF) {
          MSTK_Report("MESH_ImportFromFLAGX3D",
                      "Premature end of file while reading cell data",MSTK_ERROR);
          return 0;
        }

        if (strncmp(temp_str,"end_node_data",13) == 0)
          nodevarsdone = 1;
        else {
          MAttrib_ptr matt;
          matt = MAttrib_New(mesh, temp_str, VECTOR, MVERTEX, 3);
          double vec[3];
          int idx = 0;
          MVertex_ptr mv;
          while ((mv = MESH_Next_Vertex(mesh, &idx))) {
            status = fscanf(fp,"%lf %lf %lf",&(vec[0]), &(vec[1]), &(vec[2]));
            if (status == EOF) {
              MSTK_Report("MESH_ImportFromFLAGX3D",
                          "Premature end of file while reading cell data",
                          MSTK_ERROR);
              return 0;
            }
            MEnt_Set_AttVal(mv,matt,0,0.0,vec);
          }
          char temp_str1[256], temp_str2[256];
          strcpy(temp_str1, "end_");
          strcat(temp_str1, temp_str);
          status = fscanf(fp,"%s",temp_str2);
          if (status == EOF) {
            MSTK_Report("MESH_ImportFromFLAGX3D",
                        "Premature end of file while reading cell data",
                        MSTK_ERROR);
            return 0;
          } else if (strcmp(temp_str1, temp_str2) != 0) {
            char temp_str3[256];
            sprintf(temp_str3,"Expected %s but got %s while reading node data",
                    temp_str1, temp_str2);
            MSTK_Report("MESH_ImportFromFLAGX3D", temp_str3, MSTK_FATAL);
          }
        }
      }
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


#ifdef MSTK_HAVE_MPI  

  if (distributed && parallel_opts && parallel_opts[0]) {
    /* Must weave distributed meshes to create correct ghost links */
  
    int num_ghost_layers = parallel_opts[1];
    int input_type = 2;
    int topodim = ndim;
  
    int weavestatus = MSTK_Weave_DistributedMeshes(mesh, topodim,
                                          num_ghost_layers, input_type, comm);
  
    if (!weavestatus)
      MSTK_Report(funcname,
                  "Could not weave distributed meshes together correctly",
                  MSTK_FATAL);

    /* Run a parallel check to make sure all connectivity is consistent */

    int parallel_check = MESH_Parallel_Check(mesh,comm);        
    
    if (!parallel_check)
      MSTK_Report(funcname, "Parallel mesh checks failed", MSTK_FATAL);
  }

#endif
 
  return 1;
}

#ifdef __cplusplus
}
#endif
