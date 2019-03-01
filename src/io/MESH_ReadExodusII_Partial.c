/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#include "MSTK.h"
#include "MSTK_private.h"
#include "MSTK_VecFuncs.h"

#include "exodusII.h"
#ifdef EXODUSII_4
#include "exodusII_ext.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif


  /* Read a list of elements from a single Exodus II file on one processor into MSTK */
  /* Nemesis files are the distributed versions of Exodus II files
     with the global IDs of elements and nodes encoded as element maps
     and node maps. Nemesis files also contain some additional info
     that can be read by the Nemesis API but we are choosing to ignore
     that for now. */

  /* Right now we are creating an attribute for material sets, side
     sets and node sets. We are ALSO creating meshsets for each of
     these entity sets. We are also setting mesh geometric entity IDs
     - Should we pick one of the first two? */

  /* From the documentation of netcdf 4.1.2 - "One of the goals of
     netCDF is to support efficient access to small subsets of large
     datasets.  To support this goal, netCDF uses direct access rather than
     sequential access. This can be much more efficient when the order in
     which data is read is different from the order in which it was
     written, or when it must be read in different orders for different
     applications." */

  
#define DEF_MAXFACES 20


  int MESH_ReadExodusII_Partial(Mesh_ptr mesh, const char *filename,
				const int rank, const int nelems,
				int *elem_ids) {

    char mesg[256], funcname[32]="MESH_ReadExodusII_Partial";
    char title[256], sidesetname[256], nodesetname[256];
    char face_type[256];
    char matsetname[256];
    char **elblock_names;
    int comp_ws = sizeof(double), io_ws = 0;
    int exoid=0, status;
    int *elblock_ids, *connect, *node_map, *elem_map, *nnpe;
    int max_el_in_blk=0;
    int *sideset_ids, *ss_elem_list, *ss_side_list, *nodeset_ids, *ns_node_list;
    int num_nodes_in_set, num_sides_in_set, num_df_in_set;
    RepType reptype = MESH_RepType(mesh);

    double *xvals, *yvals, *zvals, xyz[3], rval;
    float version;
    void *pval;

    int exo_nrf[3] = {4,5,6};
    int exo_nrfverts[3][6] =
      {{3,3,3,3,0,0},{4,4,4,3,3,0},{4,4,4,4,4,4}};

    int exo_rfverts[3][6][4] =
      {{{0,1,3,-1},{1,2,3,-1},{2,0,3,-1},{0,2,1,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}},
       {{0,1,4,3},{1,2,5,4},{2,0,3,5},{0,2,1,-1},{3,4,5,-1},{-1,-1,-1,-1}},
       {{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7},{0,3,2,1},{4,5,6,7}}};
    int exo_rfdirs[3][6] =
      {{1,1,1,1,-99,-99},
       {1,1,1,1,1,-99},
       {1,1,1,1,1,1}};

    List_ptr fedges, rfaces;
    MVertex_ptr mv;
    MEdge_ptr me;
    MAttrib_ptr nmapatt=NULL, nodesetatt=NULL, sidesetatt=NULL, dirflagatt=NULL;
    MSet_ptr faceset=NULL, nodeset=NULL, sideset=NULL;
    int distributed=0;
    
  
    ex_init_params exopar;

    FILE *fp;

    exoid = ex_open(filename, EX_READ, &comp_ws, &io_ws, &version);

    if (exoid < 0) {
      sprintf(mesg,"Cannot open/read Exodus II file %s\n",filename);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }

    status = ex_get_init_ext(exoid, &exopar);
    if (status < 0) {
      sprintf(mesg,"Error reading Exodus II file %s\n",filename);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }
  
  
    strcpy(title,exopar.title);
    int ndim = exopar.num_dim;
    int nnodes_total = exopar.num_nodes;
    int nedges_total = exopar.num_edge;
    int nedgeblock_total = exopar.num_edge_blk;
    int nfaces_total = exopar.num_face;
    int nfaceblock_total = exopar.num_face_blk;
    int nelems_total = exopar.num_elem;
    int nelblock_total = exopar.num_elem_blk;
    int nnodesets_total = exopar.num_node_sets;
    int nedgesets_total = exopar.num_edge_sets;
    int nfacesets_total = exopar.num_face_sets;
    int nelemsets_total = exopar.num_elem_sets;
    int nsidesets_total = exopar.num_side_sets;
    int nnodemaps_total = exopar.num_node_maps;
    int nedgemaps_total = exopar.num_edge_maps;
    int nfacemaps_total = exopar.num_face_maps;
    int nelemmaps_total = exopar.num_elem_maps;

    if (ndim == 1)
      MSTK_Report("MESH_ReadExodusII_Partial.c","Cannot handle 1D Exodus meshes",
		  MSTK_FATAL);


    /* sort the elem IDs - so that we can specify ranges of elements */
    qsort(elem_ids, nelems, sizeof(int), compareINT);
    int min_elem_id = elem_ids[0];
    int max_elem_id = elem_ids[nelems-1];


    /* maps from elements in full mesh to partial mesh and vice versa */
    int nalloc = (max_elem_id > nelems_total) ? max_elem_id+1 : nelems_total+1;
    int *global2local_elem_map = (int *) calloc(nalloc, sizeof(int));
    int *local2global_elem_map = (int *) malloc(nelems*sizeof(int));

    /* Get element block IDs and element block info */
    
    elblock_ids = (int *) malloc(nelblock_total*sizeof(int));
    status = ex_get_ids(exoid, EX_ELEM_BLOCK, elblock_ids);
    if (status < 0) {
      sprintf(mesg,
	      "Error reading element block ids in Exodus II file %s\n",
	      filename);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }
  
    elblock_names = (char **) malloc(nelblock_total*sizeof(char *));
    for (int i = 0; i < nelblock_total; i++)
      elblock_names[i] = (char *) malloc(256*sizeof(char));
    status = ex_get_names(exoid, EX_ELEM_BLOCK, elblock_names);
    if (status < 0) {
      sprintf(mesg,"Error reading element block ids in Exodus II file %s\n",filename);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }

    
    int *elblock_min = (int *) malloc(nelblock_total*sizeof(int));
    int *elblock_max = (int *) malloc(nelblock_total*sizeof(int));
    int *elblock_nelems = (int *) malloc(nelblock_total*sizeof(int));
    char (*elblock_type)[256] = (char (*)[256]) malloc(nelblock_total*sizeof(char[256]));
    int *elblock_itype = (int *) malloc(nelblock_total*sizeof(int));
    int *elblock_nelnodes = (int *) malloc(nelblock_total*sizeof(int));
    int *elblock_neledges = (int *) malloc(nelblock_total*sizeof(int));
    int *elblock_nelfaces = (int *) malloc(nelblock_total*sizeof(int));
    int *elblock_natts = (int *) malloc(nelblock_total*sizeof(int));

    /* Parse the element blocks */

    int mesh_type = (ndim == 2) ? 1 : 0;  /* 0: unknown, 1: surface, 2: solid */
    int surf_elems = 0, solid_elems = 0;
    
    int elstart = 1;
    for (int b = 0; b < nelblock_total; b++) {
      status = ex_get_block(exoid, EX_ELEM_BLOCK, elblock_ids[b], 
			    elblock_type[b], &(elblock_nelems[b]),
			    &(elblock_nelnodes[b]), &(elblock_neledges[b]),
			    &(elblock_nelfaces[b]), &(elblock_natts[b]));

      if (status < 0) {
	sprintf(mesg, "Error reading element block %s in Exodus II file %s\n",
		elblock_names[b],filename);
	MSTK_Report(funcname,mesg,MSTK_FATAL);
      }

      elblock_min[b] = elstart;
      elblock_max[b] = elstart + elblock_nelems[b] - 1;
      elstart += elblock_nelems[b];

      elblock_itype[b] = -1;  /* need this only for tets, prisms and hexes */
      if (ndim == 3) {
	if (strncasecmp(elblock_type[b],"NSIDED",6) == 0 ||
	    strncasecmp(elblock_type[b],"TRI",3) == 0 ||
	    strncasecmp(elblock_type[b],"QUAD",4) == 0 ||
	    strncasecmp(elblock_type[b],"SHELL",5) == 0)
	  mesh_type = 1;  /* surface mesh in 3D */
	else {
	  mesh_type = 2;  /* must be a solid element*/

          if (strncasecmp(elblock_type[b],"TET",3) == 0) {
	      elblock_itype[b] = 0;
	      if (elblock_nelnodes[b] > 4)
		MSTK_Report(funcname,"Higher order tets not supported",
			    MSTK_WARN);
	  } else if (strncasecmp(elblock_type[b],"WEDGE",5) == 0) {
	    elblock_itype[b] = 1;
	    if (elblock_nelnodes[b] > 6)
	      MSTK_Report(funcname,"Higher order prisms/wedges not supported",
			  MSTK_WARN);
	  } else if (strncasecmp(elblock_type[b],"HEX",3) == 0) {
	    elblock_itype[b] = 2;
	    if (elblock_nelnodes[b] > 8)
	      MSTK_Report(funcname,"Higher order hexes not supported",
			  MSTK_WARN);
	  } else {
            sprintf(mesg,"Unrecognized or unsupported solid element type: %s",
                    elblock_type[b]);
            MSTK_Report(funcname,mesg,MSTK_FATAL);
            continue;
          }
	}
      }
    }  /* parse element blocks */

    MSet_ptr *matset = (MSet_ptr *) malloc(nelblock_total*sizeof(MSet_ptr));
    MType mstk_elem_type = (mesh_type == 1) ? MFACE : MREGION;
    MAttrib_ptr elblockatt = MAttrib_New(mesh, "elemblock", INT,
					 mstk_elem_type);
    for (int b = 0; b < nelblock_total; b++) {
      sprintf(matsetname,"matset_%-d",elblock_ids[b]);
      matset[b] = MSet_New(mesh, matsetname, mstk_elem_type);
    }
    

    int nelfaces_sum = 0;  /* total number of face entries over all elems */
    int nelfaces_alloc = 0;
    int *nelfaces = NULL;
    int *elfaces = NULL;
    int nfnodes_sum = 0;  /* total number of nodes entries over all faces */
    int nfnodes_alloc = 0;
    int nfaces2read = 0;
    int *faces2read = NULL;
    int *fnodes = NULL;
    int *nfnodes = NULL;
    int nnodes2read = 0;
    int *nodes2read = NULL;

    int *elem_blockidx = (int *) malloc(nelems*sizeof(int));

    int elcounter = 0;

    
    /* Gather info needed for any polyhedral elements */

    int reliable_rfdirs = -1;
#ifdef MSTK_USE_MARKERS
      reliable_rfdirs = MSTK_GetMarker();  /* region knows its face dirs? */
#else
      dirflagatt = MAttrib_New(mesh, "dirflagatt", INT, MREGION);
#endif
      
    
    if (ndim == 3 && nfaceblock_total) {

      nelfaces = (int *) malloc(nelems*sizeof(int));  /* num faces for each element */
      nelfaces_alloc = nelems*5;  /* guess on what nelfaces_sum will be */
      elfaces = (int *) malloc(nelfaces_alloc*sizeof(int));  /* faces of each element */
      
      int ibeg = 0;  /* ibeg, iend are indices into elem_ids array
			(not the global list of element ids) */
      while (ibeg < nelems) {
	
	/* find a contiguous range of elem_ids starting from ibeg */
	int iend = ibeg;
	while (elem_ids[iend+1] == elem_ids[iend]+1) iend++;

	int begel = elem_ids[ibeg];
	int endel = elem_ids[iend];
	for (int b = 0; b < nelblock_total; b++) {

	  if (elblock_min[b] <= begel && begel <= elblock_max[b]) {

	    /* found block containing starting elem 'begel'*/
	    
	    int jbeg = begel - elblock_min[b];  /* jbeg is index in global
						   list of elements */
	    if (elblock_max[b] < endel)
	      endel = elblock_max[b];    /* range spans multiple blocks */
	    int nelems_cur = endel - begel + 1;
	    
	    if (strncasecmp(elblock_type[b],"NFACED",6) == 0) {

	      /* cannot get face counts for a partial set of polyhedra
		 - so get all and filter out the ones not needed */
	      
	      int *nelfaces_cur = (int *) malloc(elblock_nelems[b]*sizeof(int));
	      status = ex_get_entity_count_per_polyhedra(exoid, EX_ELEM_BLOCK,
							 elblock_ids[b],
							 nelfaces_cur);
	      if (status < 0) {
		sprintf(mesg,
			"Error reading elem block info in Exodus II file %s\n",
			filename);
		MSTK_Report(funcname,mesg,MSTK_FATAL);
	      }

	      /* Calculate the offset where the faces of this element
		 will start in the element-face connectivity array */
	      int elfaces_beg = 0;
	      int j;
	      for (j = 0; j < jbeg; j++)
		elfaces_beg += nelfaces_cur[j];
	      
	      int nelfaces_cur_sum = 0;
	      for (j = 0; j < nelems_cur; j++) {
		nelfaces_cur_sum += nelfaces_cur[jbeg+j];
		nelfaces[elcounter] = nelfaces_cur[jbeg+j];
		elem_blockidx[elcounter] = b;
		global2local_elem_map[begel+j] = elcounter+1;
		local2global_elem_map[elcounter] = begel+j;
		elcounter++;
	      }

	      if (nelfaces_alloc < nelfaces_sum + nelfaces_cur_sum) {
		nelfaces_alloc = nelfaces_sum + nelfaces_cur_sum;
		elfaces = (int *) realloc(elfaces, nelfaces_alloc*sizeof(int));
	      }
	      
	      /* ex_get_partial_conn for NFACED elements is different
		 than for standard elements; for NFACED elements you
		 have to give it the offset to the first face of the
		 starting element (begel) and the number of faces to
		 be read */

	      status = ex_get_partial_conn(exoid, EX_ELEM_BLOCK,
					   elblock_ids[b], elfaces_beg+1,
					   nelfaces_cur_sum,
					   NULL, NULL,
					   elfaces + nelfaces_sum);
	      if (status < 0) {
		sprintf(mesg,
			"Error reading elem block info in Exodus II file %s\n",
			filename);
		MSTK_Report(funcname,mesg,MSTK_FATAL);
	      }

	      nelfaces_sum += nelfaces_cur_sum;
	      
	      free(nelfaces_cur);
	    }  /* if elem_type[i] == "nfaced" */
	  }
	}
	ibeg = iend+1;
      }  /* while (ibeg < nelems) */
	  
	  
      /* Deduplicate faces in elfaces list and compile a list of faces
	 that we actually have to read */

      int *faces2read_tmp = (int *) malloc(nelfaces_sum*sizeof(int));
      memcpy(faces2read_tmp, elfaces, nelfaces_sum*sizeof(int));
      qsort(faces2read_tmp, nelfaces_sum, sizeof(int), compareINT);
      faces2read = (int *) malloc(nelfaces_sum*sizeof(int));
      faces2read[0] = faces2read_tmp[0];
      int nfaces2read = 1;
      int curf = 0, nxtf = 1;
      while (nxtf < nfaces2read) {
	if (faces2read_tmp[nxtf] != faces2read_tmp[curf]) {
	  faces2read[nfaces2read++] = faces2read_tmp[nxtf];
	  curf = nxtf;
	}
	nxtf++;
      }
      free(faces2read_tmp);

      
      /* Now read the faces from the face block and collect the nodes
	 they reference */

      int *faceblock_ids = (int *) malloc(nfaceblock_total*sizeof(int));
      status = ex_get_ids(exoid, EX_FACE_BLOCK, faceblock_ids);
      if (status < 0) {
        sprintf(mesg,
		"Error reading face block ids in Exodus II file %s\n",
		filename);
        MSTK_Report(funcname,mesg,MSTK_FATAL);
      }

      int *faceblock_min = (int *) malloc(nfaceblock_total*sizeof(int));
      int *faceblock_max = (int *) malloc(nfaceblock_total*sizeof(int));
      int *faceblock_nfaces = (int *) malloc(nfaceblock_total*sizeof(int));

      int nfnodes_sum_max = 0;
      for (int i = 0; i < nfaceblock_total; i++) {
	int nfnodes_cur_sum;
	status = ex_get_block(exoid, EX_FACE_BLOCK, faceblock_ids[i], face_type,
			      &(faceblock_nfaces[i]), &nfnodes_cur_sum,
			      NULL, NULL, NULL);
	faceblock_min[i] = (i == 0) ? 0 : faceblock_max[i-1]+1;
	faceblock_max[i] = faceblock_min[i] + faceblock_nfaces[i];
	nfnodes_sum_max += nfnodes_cur_sum;
      }

      nfnodes = (int *) malloc(nfaces2read*sizeof(int));
      fnodes = (int *) malloc(nfnodes_sum_max*sizeof(int));

      nfnodes_sum = 0;  /* reset to use as counter */
      ibeg = 0;  /* ibeg, iend are indices into faces2read array */
      while (ibeg < nfaces2read) {
	/* Find range of contiguous face numbers */
	int iend = ibeg;
	while (faces2read[iend+1] == faces2read[iend]+1) iend++;

	int begf = faces2read[ibeg];
	int endf = faces2read[iend];
	
	for (int i = 0; i < nfaceblock_total; i++) {
	  if (faceblock_min[i] <= begf && begf <= faceblock_max[i]) {

	    if (endf > faceblock_max[i])
	      endf = faceblock_max[i];   /* range spans multiple face blocks */

	    int nface_cur = endf - begf + 1;
	    int jbeg = begf - faceblock_min[i];
	    
	    /* Node count for each face so we know how to step through
	       connectivity array */

	    int *nfnodes_cur = (int *) malloc(faceblock_nfaces[i]*sizeof(int));
	    status = ex_get_entity_count_per_polyhedra(exoid, EX_FACE_BLOCK,
						       faceblock_ids[i],
						       nfnodes_cur);
	    
	    if (status < 0) {
	      sprintf(mesg,
		      "Error reading face block info in Exodus II file %s\n",
		      filename);
	      MSTK_Report(funcname,mesg,MSTK_FATAL);
	    }

	    /* Calculate the offset for the first node of the first face that is
	       to be read */
	    int fnodes_beg = 0;
	    int j;
	    for (j = 0; j < jbeg; j++)
	      fnodes_beg += nfnodes_cur[0];

	    int nfnodes_cur_sum = 0;
	    for (int j = 0; j < nface_cur; j++) {
	      nfnodes_cur_sum += nfnodes_cur[jbeg+j];
	      nfnodes[nfnodes_sum+j] = nfnodes_cur[jbeg+j];
	    }
	    
	    /* ex_get_partial_conn requires us to provde the offset of
	       the first node of the starting face and the number of
	       face nodes to retrieve */
	    status = ex_get_partial_conn(exoid, EX_FACE_BLOCK, faceblock_ids[i],
					 fnodes_beg+1, nfnodes_cur_sum,
					 fnodes + nfnodes_sum, NULL, NULL);
	    if (status < 0) {
	      sprintf(mesg,
		      "Error reading face block info in Exodus II file %s\n",
		      filename);
	      MSTK_Report(funcname,mesg,MSTK_FATAL);
	    }

	    nfnodes_sum += nfnodes_cur_sum;
	    free(nfnodes_cur);
	  }
	}
	ibeg = iend+1;
      }  /* while (ibeg < nfaces2read) */

      free(faceblock_ids);
    }  /* ndim == 3 */



    int nelnodes_sum = 0;
    int nelnodes_alloc = nnodes_total;
    int *elnodes = (int *) malloc(nelnodes_alloc*sizeof(int));
    int *nelnodes = (int *) malloc(nelems_total*sizeof(int));

    int ibeg = 0;
    while (ibeg < nelems) {
      
      /* find a contiguous range of elem_ids starting from begel */
      int iend = ibeg;
      while (iend < nelems-1 && elem_ids[iend+1] == elem_ids[iend]+1) iend++;

      int begel = elem_ids[ibeg];
      
      for (int b = 0; b < nelblock_total; b++) {
	
	if (elblock_min[b] <= begel && begel <= elblock_max[b]) {
	  
	  /* found block containing starting elem 'begel'*/
	  
	  int jbeg = begel - elblock_min[b];
	  int endel = elem_ids[iend];
	  if (elblock_max[b] < endel) {
	    endel = elblock_max[b];  /* range spans multiple blocks */
	    iend = ibeg + endel - begel;
	  }
	  int nelems_cur = endel - begel + 1;

	  if (strncasecmp(elblock_type[b], "NSIDED", 6) == 0) {
	    int *nelnodes_cur = (int *) malloc(elblock_nelems[b]*sizeof(int));

	    /* Number of elements per node for every element in block */
	    status = ex_get_entity_count_per_polyhedra(exoid, EX_ELEM_BLOCK,
						       elblock_ids[b],
						       nelnodes_cur);
	    if (status < 0) {
	      sprintf(mesg,
		      "Error reading elem block info in Exodus II file %s\n",
		      filename);
	      MSTK_Report(funcname,mesg,MSTK_FATAL);
	    }

	    /* calculate offset where elem nodes of begel start - need
	       for extracting partial connectivity of NSIDED
	       elements */
	    
	    int elnodes_beg = 0;
	    int j;
	    for (j = 0; j < jbeg; j++)
	      elnodes_beg += nelnodes_cur[j];

	    /* calculate number of nodes that will be read; also, populate maps */
	    int nelnodes_cur_sum = 0;
	    for (j = 0; j < nelems_cur; j++) {
	      nelnodes_cur_sum += nelnodes_cur[jbeg+j];
	      nelnodes[elcounter] = nelnodes_cur[jbeg+j];
	      elem_blockidx[elcounter] = b;
	      global2local_elem_map[begel+j] = elcounter+1;
	      local2global_elem_map[elcounter] = begel+j;
	      elcounter++;
	    }

	    if (nelnodes_sum + nelnodes_cur_sum > nelnodes_alloc) {
	      nelnodes_alloc = 1.5*(nelnodes_alloc + nelnodes_cur_sum);
	      elnodes = (int *) realloc(elnodes, nelnodes_alloc*sizeof(int));
	    }
	    
	    /* ex_get_partial_conn for NSIDED elements is different
	       from the call for standard elements - you have to give
	       it the starting offset (+1) for where the first node of
	       the first element you want is listed and then tell it
	       the total number of nodes you want read */
	    status = ex_get_partial_conn(exoid, EX_ELEM_BLOCK, elblock_ids[b],
					 elnodes_beg+1, nelnodes_cur_sum,
					 elnodes + nelnodes_sum, NULL, NULL);

	    if (status < 0) {
	      sprintf(mesg,
		      "Error reading face block info in Exodus II file %s\n",
		      filename);
	      MSTK_Report(funcname,mesg,MSTK_FATAL);
	    }

	    nelnodes_sum += nelnodes_cur_sum;
	    
	  } else if (strncasecmp(elblock_type[b],"TRI",3) == 0 ||
		     strncasecmp(elblock_type[b],"QUAD",4) == 0 ||
		     strncasecmp(elblock_type[b],"SHELL",5) == 0 ||
		     strncasecmp(elblock_type[b],"TET",3) == 0 ||
		     strncasecmp(elblock_type[b],"WEDGE",5) == 0 ||
		     strncasecmp(elblock_type[b],"HEX",3) == 0) {
	    
	    int nelnodes_cur_sum = nelems_cur*elblock_nelnodes[b];
	    
	    if (nelnodes_sum + nelnodes_cur_sum > nelnodes_alloc) {
	      nelnodes_alloc = 1.5*(nelnodes_alloc + nelnodes_cur_sum);
	      elnodes = (int *) realloc(elnodes, nelnodes_alloc*sizeof(int));
	    }	    
	    
	    status = ex_get_partial_conn(exoid, EX_ELEM_BLOCK, elblock_ids[b],
					 jbeg+1, nelems_cur,
					 elnodes + nelnodes_sum, NULL, NULL);
	    if (status < 0) {
	      sprintf(mesg,
		      "Error reading element block %s in Exodus II file %s\n",
		      elblock_names[b],filename);
	      MSTK_Report(funcname,mesg,MSTK_FATAL);
	    }
	    
	    
	    /* This is a block of standard elements - each element has
	       the same number of nodes */
	    int j;
	    for (j = 0; j < nelems_cur; j++) {
	      nelnodes[elcounter] = elblock_nelnodes[b];
	      elem_blockidx[elcounter] = b;
	      global2local_elem_map[begel+j] = elcounter+1;
	      local2global_elem_map[elcounter] = begel+j;
	      elcounter++;
	    }

	    nelnodes_sum += elblock_nelnodes[b]*nelems_cur;
	  } else if (strncasecmp(elblock_type[b],"NFACED",5) == 0) {
	    /* we already collected data from this block */
	  } else {
	    if (status < 0) {
	      sprintf(mesg,
		      "UNKNOWN element type (%s) info in Exodus II file %s\n",
		      elblock_type[b], filename);
	      MSTK_Report(funcname,mesg,MSTK_FATAL);
	    }
	  }
	}  /* if (elblock_min[b] <= begel ....etc */
      }  /* for b = 0, nelblock_total-1 */

      ibeg = iend+1;
    }  /* while ibeg < nelems */

    
    /* collect face nodes (for faces of polyhedral elements) and
       element nodes (for non-polyhedral elements) and filter out any
       duplicates */ 
    
    int maxnodes2read = nfnodes_sum + nelnodes_sum;
    nodes2read = (int *) malloc(maxnodes2read*sizeof(int));

    int *nodes2read_tmp = (int *) malloc(maxnodes2read*sizeof(int));
    memcpy(nodes2read_tmp, fnodes, nfnodes_sum*sizeof(int));
    memcpy(nodes2read_tmp+nfnodes_sum, elnodes, nelnodes_sum*sizeof(int));

    qsort(nodes2read_tmp, maxnodes2read, sizeof(int), compareINT);

    nodes2read[0] = nodes2read_tmp[0];
    nnodes2read = 1;
    int curn = 0, nxtn = 1;
    while (nxtn < maxnodes2read) {
      if (nodes2read_tmp[nxtn] != nodes2read_tmp[curn]) {
	nodes2read[nnodes2read++] = nodes2read_tmp[nxtn];
	curn = nxtn;
      }
      nxtn++;
    }
    free(nodes2read_tmp);

    
    /* Create the nodes required by this partition */

    int min_node_id = nodes2read[0];
    int max_node_id = nodes2read[nnodes2read-1];
    MVertex_ptr *mverts = (MVertex_ptr *) malloc(nnodes2read*sizeof(MVertex_ptr));

    nalloc = (max_node_id > nnodes_total) ? max_node_id+1 : nnodes_total+1;
    int *global2local_node_map = (int *) calloc(nalloc, sizeof(int));
    int *local2global_node_map = (int *) malloc(nnodes2read*sizeof(int));
    xvals = (double *) malloc(nnodes2read*sizeof(double));
    yvals = (double *) malloc(nnodes2read*sizeof(double));
    zvals = (double *) malloc(nnodes2read*sizeof(double));

    ibeg = 0;
    while (ibeg < nnodes2read) {
      int iend = ibeg;
      while (iend < nnodes2read-1 && nodes2read[iend+1] == nodes2read[iend]+1) iend++;
      int nnodes_cur = iend - ibeg + 1;
      int begnode = nodes2read[ibeg];
      
      status = ex_get_partial_coord(exoid, begnode, nnodes_cur, xvals, yvals, zvals);
      if (status < 0) {
	sprintf(mesg, "Error reading node coordinates in Exodus II file %s\n",
		filename);
	MSTK_Report(funcname,mesg,MSTK_FATAL);
      }

      for (int n = 0; n < nnodes_cur; n++) {
	MVertex_ptr mv = MV_New(mesh);
    
	xyz[0] = xvals[n]; xyz[1] = yvals[n]; xyz[2] = zvals[n];
	MV_Set_Coords(mv, xyz);
#ifdef MSTK_HAVE_MPI
	MV_Set_GlobalID(mv, nodes2read[ibeg+n]);
#endif

	mverts[ibeg+n] = mv;
	global2local_node_map[nodes2read[ibeg+n]] = ibeg+n+1;
	local2global_node_map[ibeg+n] = nodes2read[ibeg+n];
      }

      ibeg = iend+1;
    }


    /* now subsitute the node numbers in the fnodes and elnodes array with the
       local node numbers */

    /* we will use a hash table later so that lookup time can be close
       to O(1) instead of O(logN) */
  
    for (int n = 0; n < nfnodes_sum; n++) {
      int *p = bsearch(&(fnodes[n]), nodes2read, nnodes2read, sizeof(int),
		       compareINT);
      int loc = p - nodes2read;
      fnodes[n] = loc;
    }

    for (int n = 0; n < nelnodes_sum; n++) {
      int *p = bsearch(&(elnodes[n]), nodes2read, nnodes2read, sizeof(int),
		       compareINT);
      int loc = p - nodes2read;
      elnodes[n] = loc;
    }


    /* Now make the elements */

    MVertex_ptr fverts[MAXPV2];
    int nfaces = 0;
    if (ndim == 2) {

      int offset = 0;
      for (int e = 0; e < nelems; e++) {
	int ib = elem_blockidx[e];
	int b = elblock_ids[ib];
	MFace_ptr mf = MF_New(mesh);

	for (int n = 0; n < nelnodes[e]; n++)
	  fverts[n] = mverts[elnodes[offset+n]];
	offset += nelnodes[e];

	MF_Set_Vertices(mf, nelnodes[e], fverts);
#ifdef MSTK_HAVE_MPI
	MF_Set_GlobalID(mf, elem_ids[e]);
#endif
	MF_Set_GEntID(mf, b);
	MF_Set_GEntID(mf, 2);

	MSet_Add(matset[ib], mf);
	MEnt_Set_AttVal(mf, elblockatt, b, 0.0, NULL);  /* necessary? */
      }

    } else {  /* ndim == 3 */

      if (mesh_type == 1) {  /* surface mesh */

	int noffset = 0;
	for (int e = 0; e < nelems; e++) {
	  int ib = elem_blockidx[e];
	  int b = elblock_ids[ib];
	
	  MFace_ptr mf = MF_New(mesh);
	
	  for (int n = 0; n < nelnodes[e]; n++)
	    fverts[n] = mverts[elnodes[noffset+n]];
	  noffset += nelnodes[e];
	
	  MF_Set_Vertices(mf, nelnodes[e], fverts);
#ifdef MSTK_HAVE_MPI
	  MF_Set_GlobalID(mf, elem_ids[e]);
#endif
	  MF_Set_GEntID(mf, b);
	  MF_Set_GEntDim(mf, 3);
	
	  MSet_Add(matset[ib], mf);
	  MEnt_Set_AttVal(mf, elblockatt, b, 0.0, NULL);
	}
      } else {  /* solid mesh */

	MFace_ptr *mfaces = NULL;
	if (nfaces2read) {

	  /* create faces of any polyhedral elements */

	  mfaces = (MFace_ptr *) malloc(nfaces2read*sizeof(MFace_ptr));

	  for (int f = 0; f < nfaces2read; f++) {
	    mfaces[f] = MF_New(mesh);
	    for (int n = 0; n < nfnodes[f]; n++)
	      fverts[n] = mverts[fnodes[n]];
	    MF_Set_Vertices(mfaces[f], nfnodes[f], fverts);
	    MSet_Add(faceset, mfaces[f]);
	  }

	  /* replace face numbers in the elfaces array with local face numbers */
       
	  for (int f = 0; f < nelfaces_sum; f++) {
	    int *p = bsearch(&(elfaces[f]), faces2read, nfaces2read,
			     sizeof(int), compareINT);
	    int loc = p - faces2read;
	    elfaces[f] = loc;
	  }
	}  /* if nfaces2read */

      
	/* now make the elements */

	MVertex_ptr rverts[MAXPV3];
	int foffset = 0, noffset = 0;
	for (int e = 0; e < nelems; e++) {
	  int ib = elem_blockidx[e];
	  int b = elblock_ids[ib];

	  int eltype = elblock_itype[ib];
	  
	  MRegion_ptr mr = MR_New(mesh);

	  MFace_ptr rfaces_arr[MAXPF3];
	  int rfdirs[MAXPF3];
	  int nrf = 0;
	  if (eltype == -1) {  /* polyhedron */

	    nrf = nelfaces[e];

	    for (int f = 0; f < nrf; f++) {
	      rfaces_arr[f] = mfaces[elfaces[foffset+f]];
	      List_ptr fregs = MF_Regions(rfaces_arr[f]);
	      if (fregs) {
		if (List_Num_Entries(fregs) == 1) {
		  MRegion_ptr adjreg = List_Entry(fregs, 0);
		  int oppfdir = MR_FaceDir(adjreg, rfaces_arr[f]);
		  rfdirs[f] = !oppfdir;
		} else
		  MSTK_Report(funcname, "Face already connected to two regions",
			      MSTK_FATAL);
		List_Delete(fregs);
	      } else
		rfdirs[f] = 1;
	    }
	    foffset += nrf;

	  } else {  /* standard element */
	  
	    for (int n = 0; n < nelnodes[e]; n++)
	      rverts[n] = mverts[elnodes[noffset+n]];
	    noffset += nelnodes[e];

	    nrf = exo_nrf[eltype];
	    if (reptype == F1 || reptype == F4) {                        
	      for (int f = 0; f < nrf; f++) {
		MFace_ptr face = NULL;

		int nfv = exo_nrfverts[eltype][f];              
		for (int k1 = 0; k1 < nfv; k1++)
		  fverts[k1] = rverts[exo_rfverts[eltype][f][k1]];

		/* Use a hash to find faces? */
		List_ptr vfaces = MV_Faces(fverts[0]);
		int nvfaces = vfaces ? List_Num_Entries(vfaces) : 0;
		int ii;
		for (ii = 0; ii < nvfaces; ii++) {
		  MFace_ptr vface = List_Entry(vfaces,ii);
	      
		  int jj, has_all_verts = 1;
		  for (jj = 0; jj < nfv; jj++)
		    if (!MF_UsesEntity(vface,fverts[jj],MVERTEX)) {
		      has_all_verts = 0;
		      break;
		    }
	      
		  if (has_all_verts) {
		    face = vface;
		    break;
		  }
		}
		if (nvfaces) List_Delete(vfaces);
            
		if (face) {
		  List_ptr fregs;
	      
		  rfaces_arr[f] = face;
		  fregs = MF_Regions(face);
		  if (!fregs) {
		    /* Comes here when the other region attached to this
		       face will be a polyhedral element that will be
		       created later. Because polyhedral elements need to
		       be described in terms of their faces, the shared
		       face explicitly described may have a different
		       vertex ordering than the implicit ordering required
		       by this standard element. In this case, let the
		       standard element ordering rule */

		    MF_Set_Vertices(face, nfv, fverts);
		    rfdirs[f] = exo_rfdirs[eltype][f]; 

		  } else {

		    if (List_Num_Entries(fregs) == 1) {
		      int dir;
		      MRegion_ptr freg;
		  
		      freg = List_Entry(fregs,0);
		      dir = MR_FaceDir(freg,face);

		      /* since this is a standard type element, it has to
		       * follow its face direction according to the
		       * template and not derive its direction from the
		       * other region. If the other region is a standard
		       * element too, it should be consistent in a valid
		       * mesh. If its a polyhedral element, its face dir
		       * will be determined from this standard element
		       * later in the code but we will fix it here to be
		       * consistent with the standard element */

		      rfdirs[f] = exo_rfdirs[eltype][f];
		      if (dir != !rfdirs[f]) {
			MF_Set_Vertices(face, nfv, fverts);
			MR_Set_FaceDir(freg, face, !dir);
		      }

		    }
		    else if (List_Num_Entries(fregs) == 2)
		      MSTK_Report(funcname,"Face already connected two faces",
				  MSTK_FATAL);
		    List_Delete(fregs);
		  }  /* if (fregs) {} else {} */
		} else {
		  face = MF_New(mesh);
	      
		  MF_Set_Vertices(face, exo_nrfverts[eltype][f], fverts);
		  rfaces_arr[f] = face;
		  rfdirs[f] = exo_rfdirs[eltype][f];
		}
	      }
	      
	      MR_Set_Faces(mr,nrf,rfaces_arr,rfdirs);
#ifdef MSTK_USE_MARKERS
	      MEnt_Mark(mr, reliable_rfdirs);
#else
	      MEnt_Set_AttVal(mr, dirflagatt, 1, 0.0, NULL);
#endif
	    }
	    else if (reptype == R1 || reptype == R2)
	      MR_Set_Vertices(mr, nelnodes[e], rverts, 0, NULL);
	  }

#ifdef MSTK_HAVE_MPI
	  MR_Set_GlobalID(mr, elem_ids[e]);
#endif
	  MR_Set_GEntID(mr, b);
	  MR_Set_GEntDim(mr, 3);
	
	  MSet_Add(matset[ib], mr);
	  MEnt_Set_AttVal(mr, elblockatt, b, 0.0, NULL);
	}
      }
    

      /* Now we have to go through some painful steps to reliably
       * get all face directions of the regions. If we were sure
       * that no element would be degenerate (0 volume or having 0
       * area faces), we could compare the normal of the face with
       * the vector from the center of the region to the center of
       * the face for each cell. If, however, we expect degenerate
       * faces we have to eliminate or at least minimize such
       * comparisons, which is what is done here. NOTE THAT DURING
       * THIS PROCESS THE MESH MAY BECOME TEMPORARILY TOPOLOGICALLY
       * INCONSISTENT*/
    
      /* Find all regions whose region face directions have to be finalized */
    
      int rfdirs_good;
      int nbad = 0;
      int idx = 0;
      MRegion_ptr mr;
      while ((mr = MESH_Next_Region(mesh, &idx))) {
#ifdef MSTK_USE_MARKERS
	rfdirs_good = MEnt_IsMarked(mr, reliable_rfdirs);
#else
	MEnt_Get_AttVal(mr, dirflagatt, &rfdirs_good, &rval, &pval);
#endif
	if (!rfdirs_good)
	  nbad++;
      }
    
      if (nbad) {

#ifdef MSTK_USE_MARKERS
	int inreglist_mark = MSTK_GetMarker();
#else
	MAttrib_ptr reglistatt = MAttrib_New(mesh, "reglist", INT, MREGION);
#endif

	int iter = 0;
	int done = 0;
	while (!done && iter < 1000) { 
        
	  /* Number of iterations determined by number of disjoint
	   * blocks containing elements with unknown directions */

	  int nfixed = 0;
	  int found_start = 0;
	  List_ptr reglist = List_New(0);

	  // Skip this part if 'nbad' == number of regions in the mesh,
	  // i.e., the direction of faces of no region is known

	  if (nbad != MESH_Num_Regions(mesh)) {

	    /* Find a bad region which is next to a reliable region */
	    idx = 0;
	    while ((mr = MESH_Next_Region(mesh, &idx))) {
#ifdef MSTK_USE_MARKERS
	      rfdirs_good = MEnt_IsMarked(mr, reliable_rfdirs);
#else              
	      double rval;
	      void *pval;
	      MEnt_Get_AttVal(mr, dirflagatt, &rfdirs_good, &rval, &pval);
#endif
	      if (rfdirs_good) continue;  /* region-face dirs are reliable */
            
	      List_ptr rfaces = MR_Faces(mr);
	      int idx2 = 0;
	      int found = 0;
	      MFace_ptr rf;
	      while ((rf = List_Next_Entry(rfaces, &idx2))) {
		List_ptr fregs = MF_Regions(rf);
		if (List_Num_Entries(fregs) == 2) {
		  MRegion_ptr oppr = List_Entry(fregs, 0);
		  if (mr == oppr)
		    oppr = List_Entry(fregs, 1);
#ifdef MSTK_USE_MARKERS
		  rfdirs_good = MEnt_IsMarked(oppr, reliable_rfdirs);
#else
		  MEnt_Get_AttVal(oppr, dirflagatt, &rfdirs_good, &rval, &pval);
#endif
		  if (rfdirs_good) {
		    List_Add(reglist, mr);
		    found = 1;
#ifdef MSTK_USE_MARKERS
		    MEnt_Mark(mr, inreglist_mark);
#else
		    MEnt_Set_AttVal(mr, reglistatt, 1, 0.0, NULL);
#endif
		  }
		}
		List_Delete(fregs);
		if (found) break;
	      }
	      List_Delete(rfaces);
	    }
	  }  // if nbad != MESH_Num_Regions(mesh)

	  if (!List_Num_Entries(reglist)) {
	    /* We did not find a single region with unreliable face
	     * directions adjacent to a region with reliable face
	     * directions - This may happen in the first iteration
	     * through a disjointed block of polyhedral elements. We
	     * have to do geometric checks to find such a pair */
          
	    MRegion_ptr mrgood = NULL;
	    idx = 0;
	    while (!mrgood && ((mr = MESH_Next_Region(mesh, &idx)))) {
            
#ifdef MSTK_USE_MARKERS
	      rfdirs_good = MEnt_IsMarked(mr, reliable_rfdirs);
#else              
	      double rval;
	      void *pval;
	      MEnt_Get_AttVal(mr, dirflagatt, &rfdirs_good, &rval, &pval);
#endif
	      if (rfdirs_good) continue;  /* region-face dirs are reliable */
            
	      /* Found a region whose face directions could not be
		 determined topologically - so determine it
		 geometrically */

	      double (*fxyz)[3] =
		(double (*)[3]) malloc(MAXPV2*sizeof(double [3]));
            
	      /* Geometric center of region */
	      /* Average normal and geometric center of faces */
            
	      List_ptr rfaces = MR_Faces(mr);
	      int nrf = List_Num_Entries(rfaces);
            
	      List_ptr rvlist = List_New(10);
	      int nrv = 0;
            
	      double rcen[3];
	      rcen[0] = rcen[1] = rcen[2] = 0.0;
	      double fcen[MAXPF3][3], fnormal[MAXPF3][3];
            
	      for (int k = 0; k < nrf; k++) {
		MFace_ptr rf = List_Entry(rfaces, k);
              
		List_ptr fvlist = MF_Vertices(rf, 1, 0);
		int nfv = List_Num_Entries(fvlist);
              
		fcen[k][0] = fcen[k][1] = fcen[k][2] = 0.0;
		int k1;
		for (k1 = 0; k1 < nfv; k1++) {
		  int k2;
		  MVertex_ptr fv = List_Entry(fvlist, k1);
		  MV_Coords(fv, fxyz[k1]);
                
		  if (!List_Contains(rvlist, fv)) {
		    List_Add(rvlist, fv);
		    for (k2 = 0; k2 < 3; k2++) rcen[k2] += fxyz[k1][k2];
		    nrv++;
		  }
                
		  for (k2 = 0; k2 < 3; k2++) fcen[k][k2] += fxyz[k1][k2];
		}
		for (k1 = 0; k1 < 3; k1++) fcen[k][k1] /= nfv;
		List_Delete(fvlist);
		List_Delete(rvlist);
              
		fnormal[k][0] = fnormal[k][1] = fnormal[k][2] = 0.0;
		for (k1 = 0; k1 < nfv; k1++) {
		  double vec1[3], vec2[3], normal[3];
		  MSTK_VDiff3(fxyz[(k1+1)%nfv], fxyz[k1], vec1);
		  MSTK_VDiff3(fxyz[(k1+nfv-1)%nfv], fxyz[k1], vec2);
		  MSTK_VCross3(vec1, vec2, normal);
		  double len2 = MSTK_VLenSqr3(normal);
		  if (len2 > 0.0) MSTK_VNormalize3(normal);
		  MSTK_VSum3(fnormal[k], normal, fnormal[k]);
		}
	      }
	      for (int k = 0; k < 3; k++) rcen[k] /= nrv;
            
	      int nfgood = 0;
	      for (int k = 0; k < nrf; k++) {
		/* @todo: do checks to make sure face is not too warped */
              
		if (MSTK_VLenSqr3(fnormal[k]) < 1.0e-20) continue;  /* zero area face */
              
		double outvec[3];              
		MSTK_VDiff3(fcen[k], rcen, outvec);
		if (MSTK_VLenSqr3(outvec) < 1.0e-20) continue;  /* zero vol element */
              
		MSTK_VNormalize3(outvec);
		MSTK_VNormalize3(fnormal[k]);
		double dp = MSTK_VDot3(outvec, fnormal[k]);
		int dir = (dp > 0) ? 1 : 0;
		if (dir != MR_FaceDir_i(mr, k))
		  MR_Rev_FaceDir_i(mr, k);
		nfgood++;
	      }            
            
	      if (nfgood == nrf) {
		/* Found a region where we could reliably orient all the faces. */
		/* Take note of this region and mark an adjacent region as the 
		   start region for fixing process */

		mrgood = mr;
#ifdef MSTK_USE_MARKERS
		MEnt_Mark(mrgood, reliable_rfdirs);
#else
		MEnt_Set_AttVal(mrgood, dirflagatt, 1, 0.0, NULL);
#endif
		nfixed++;

		for (int k = 0; k < nrf; k++) {
		  MFace_ptr rf = List_Entry(rfaces, k);
		  MRegion_ptr adjr = NULL;
		  List_ptr fregs = MF_Regions(rf);
		  if (List_Num_Entries(fregs) == 2) {
		    adjr = List_Entry(fregs, 0);
		    if (adjr == mr)
		      adjr = List_Entry(fregs, 1);
		  }
		  List_Delete(fregs);
		  if (adjr) {
		    List_Add(reglist, adjr);
#ifdef MSTK_USE_MARKERS
		    MEnt_Mark(adjr, inreglist_mark);
#else
		    MEnt_Set_AttVal(adjr, reglistatt, 1, 0.0, NULL);
#endif
		    break;
		  }
		}
	      }
	      List_Delete(rfaces);
	      free(fxyz);
	    }  /* while (!mrgood && (mr = MESH_Next_Entry(mesh, &idx))) */

	    if (!mrgood)
	      MSTK_Report(funcname, "Could not find a single region whose faces could be reliably oriented", MSTK_FATAL);
	  }  /* if (!List_Num_Entries(reglist)) */


	  /* Starting from the reliable direction region, walk through
	   * all other connected mesh regions in the mesh or in this
	   * connected block of mesh regions */

	  idx = 0;
	  while ((mr = List_Next_Entry(reglist, &idx))) {
	    List_ptr rfaces = MR_Faces(mr);
	    int nrf = List_Num_Entries(rfaces);

	    /* if we know the direction in which an adjacent reliable
	     * region uses the common face, we can determine the
	     * direction of use of the face in this region */             

	    MFace_ptr fstart = NULL;
	    for (int k = 0; k < nrf; k++) {
	      MFace_ptr rf = List_Entry(rfaces, k);
	      List_ptr fregs = MF_Regions(rf);
	      if (List_Num_Entries(fregs) == 2) {
		MRegion_ptr oppr = List_Entry(fregs, 0);
		if (oppr == mr)
		  oppr = List_Entry(fregs, 1);
#ifdef MSTK_USE_MARKERS
		rfdirs_good = MEnt_IsMarked(oppr, reliable_rfdirs);
#else
		MEnt_Get_AttVal(oppr, dirflagatt, &rfdirs_good, &rval, &pval);
#endif
		if (rfdirs_good) {
		  int oppfdir = MR_FaceDir(oppr, rf);
		  int fdir = !oppfdir;
		  if (fdir != MR_FaceDir_i(mr, k))
		    MR_Set_FaceDir_i(mr, k, fdir);
		  fstart = rf;                
		}
	      }
	      List_Delete(fregs);
	      if (fstart) break;
	    }
          

	    /* Using the known face direction in this region, find the
	     * direction of use of all the other faces in the region */

	    List_ptr rfaces2 = List_New(nrf);
	    List_Add(rfaces2, fstart);

	    int idx2 = 0;
	    MFace_ptr curface;
	    while ((curface = List_Next_Entry(rfaces2, &idx2))) {
	      int curdir = MR_FaceDir(mr, curface);
                
	      List_ptr fedges = MF_Edges(curface, 1, 0);
	      int nfe = List_Num_Entries(fedges);
                
	      /* For each of the edges, find an adjacent face of the
	       * region sharing that edge and set its direction
	       * based on directional consistency rules for a
	       * 2-manifold object (if two faces of a region sharing
	       * an edge use the edge in opposite senses, then the
	       * two faces are used by the region in the same
	       * sense) */
            
	      int i2;
	      for (i2 = 0; i2 < nfe; ++i2) {
		MEdge_ptr cmnedge = List_Entry(fedges, i2);
		int fedir = MF_EdgeDir_i(curface, i2);
              
		/* Find another face in this region that uses this edge
		 * and whose direction of use by the region is not
		 * known */
              
		int j2;
		for (j2 = 0; j2 < nrf; ++j2) {
		  MFace_ptr adjface = List_Entry(rfaces, j2);
		  if (adjface == curface) continue;
		  if (MF_UsesEntity(adjface, cmnedge, MEDGE)) {
		    /* Found adjacent face in region */
		    if (!List_Contains(rfaces2, adjface)) { /* small list - no need of marking */
		      int adjfedir = MF_EdgeDir(adjface, cmnedge);
		      int adjdir = (fedir != adjfedir) ? curdir : !curdir;
		      if (adjdir != MR_FaceDir_i(mr, j2))
			MR_Rev_FaceDir_i(mr, j2);
		      List_Add(rfaces2, adjface);
		    }
		    break;
		  }
		}  /* for j2 = 0, nrf-1 */
	      }  /* for i2 = 0, nfe-1 */
	      List_Delete(fedges);
	    }  /* while (curface = List_Next_Entry(rfaces2, &idx2)) */
	    List_Delete(rfaces2);

	    /* This region's face directions are fixed */
#ifdef MSTK_USE_MARKERS
	    MEnt_Mark(mr, reliable_rfdirs);
#else
	    MEnt_Set_AttVal(mr, dirflagatt, 1, 0.0, NULL);
#endif
	    nfixed++;

	    /* Add its neighbors into the list if their face directions
	     * are not reliably known */

	    for (int k = 0; k < nrf; ++k) {
	      MFace_ptr rf = List_Entry(rfaces, k);
	      List_ptr fregs = MF_Regions(rf);
	      if (List_Num_Entries(fregs) == 2) {
		MRegion_ptr oppr = List_Entry(fregs, 0);
		if (oppr == mr)
		  oppr = List_Entry(fregs, 1);

		int inrlist;
#ifdef MSTK_USE_MARKERS
		inrlist = MEnt_IsMarked(oppr, inreglist_mark);
		rfdirs_good = MEnt_IsMarked(oppr, reliable_rfdirs);
#else
		MEnt_Get_AttVal(oppr, reglistatt, &inrlist, &rval, &pval);
		MEnt_Get_AttVal(oppr, dirflagatt, &rfdirs_good, &rval, &pval);
#endif
		if (!inrlist && !rfdirs_good) {
		  List_Add(reglist, oppr);
#ifdef MSTK_USE_MARKERS
		  MEnt_Mark(oppr, inreglist_mark);
#else
		  MEnt_Set_AttVal(oppr, reglistatt, 1, 0.0, NULL);
#endif
		}
	      }
	      List_Delete(fregs);
	    }
	    List_Delete(rfaces);

	  }  /* while (mr = List_Next_Entry(reglist, &idx)) */   
        
#ifdef MSTK_USE_MARKERS
	  List_Unmark(reglist, inreglist_mark);
#else
	  idx = 0;
	  while ((mr = List_Next_Entry(reglist, &idx)))
	    MEnt_Rem_AttVal(mr, reglistatt);
#endif
	  List_Delete(reglist);

	  nbad -= nfixed;
	  if (!nbad) done = 1;
	  iter++;
	}  /* while (!done && iter < 100) */

#ifdef MSTK_USE_MARKERS
	MSTK_FreeMarker(inreglist_mark);
#else
	MAttrib_Delete(reglistatt);
#endif
      }  /* if (nbad) */


#ifdef MSTK_USE_MARKERS
      idx = 0;
      while ((mr = MESH_Next_Region(mesh, &idx)))
        MEnt_Unmark(mr, reliable_rfdirs);
      MSTK_FreeMarker(reliable_rfdirs);
#else
      MAttrib_Delete(dirflagatt);
#endif

      if (nbad)
	MSTK_Report(funcname, "Could not fix face dirs of some regions",
		    MSTK_FATAL);

      /* Fix classifications of interior mesh faces only - if this mesh
	 is part of a distributed mesh we will get boundary info wrong */

      idx = 0;
      MFace_ptr mf;
      while ((mf = MESH_Next_Face(mesh,&idx))) {
	List_ptr fregs = MF_Regions(mf);
	if (fregs) {
	  int nregs = List_Num_Entries(fregs);
	  if (nregs == 2) {
	    MRegion_ptr reg0, reg1;
	    int greg0, greg1;
	    reg0 = List_Entry(fregs,0);
	    reg1 = List_Entry(fregs,1);
	    greg0 = MR_GEntID(reg0);
	    greg1 = MR_GEntID(reg1);

	    MF_Set_GEntDim(mf,3);
	    if (greg0 == greg1)
	      MF_Set_GEntID(mf,greg0);
	  }
	  List_Delete(fregs);
	}
      }

    }  /* ndim == 3 */
  


    
    /* read node number map - store it as an attribute to spit out later
       if necessary */
  
    node_map = (int *) malloc(nnodes_total*sizeof(int));

#ifdef EXODUS_6_DEPRECATED
    status = ex_get_node_num_map(exoid, node_map);
#else
    status = ex_get_id_map(exoid, EX_NODE_MAP, node_map);
#endif
    if (status < 0) {
      sprintf(mesg,"Error reading node map in Exodus II file %s\n",filename);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }
  
    nmapatt = MAttrib_New(mesh, "node_map", INT, MVERTEX);
    
    for (int i = 0; i < nnodes2read; i++) {
      mv = MESH_Vertex(mesh, i);
      int global_node_id = local2global_node_map[i];
      MEnt_Set_AttVal(mv, nmapatt, node_map[global_node_id-1], 0.0, NULL);
      
#ifdef MSTK_HAVE_MPI
        MV_Set_GlobalID(mv, node_map[global_node_id-1]);
        MV_Set_MasterParID(mv, rank); /* This might get modified later */
#endif
    }

    free(node_map);
  

    /* Read node sets */

    if (nnodesets_total) {
      nodeset_ids = (int *) malloc(nnodesets_total*sizeof(int));

#ifdef EXODUS_6_DEPRECATED
      status = ex_get_node_set_ids(exoid,nodeset_ids);
#else
      status = ex_get_ids(exoid, EX_NODE_SET, nodeset_ids);
#endif
      if (status < 0) {
        sprintf(mesg,
                "Error reading nodeset IDs in Exodus II file %s\n",
                filename);
        MSTK_Report(funcname,mesg,MSTK_FATAL);
      }
    
    
      for (int i = 0; i < nnodesets_total; i++) {

        sprintf(nodesetname,"nodeset_%-d",nodeset_ids[i]);

        nodesetatt = MAttrib_New(mesh,nodesetname,INT,MVERTEX);

        nodeset = MSet_New(mesh,nodesetname,MVERTEX);

#ifdef EXODUS_6_DEPRECATED    
        status = ex_get_node_set_param(exoid,nodeset_ids[i],&num_nodes_in_set,
                                       &num_df_in_set);
#else
        status = ex_get_set_param(exoid, EX_NODE_SET, nodeset_ids[i],
                                  &num_nodes_in_set, &num_df_in_set);
#endif
      
        ns_node_list = (int *) malloc(num_nodes_in_set*sizeof(int));
      
#ifdef EXODUS_6_DEPRECATED
        status = ex_get_node_set(exoid,nodeset_ids[i],ns_node_list);
#else
        status = ex_get_set(exoid, EX_NODE_SET, nodeset_ids[i],
                            ns_node_list, NULL);
#endif
        if (status < 0) {
          sprintf(mesg,
                  "Error reading nodes in nodeset %-d in Exodus II file %s\n",nodeset_ids[i],filename);
          MSTK_Report(funcname,mesg,MSTK_FATAL);
        }
      
      
        for (int j = 0; j < num_nodes_in_set; j++) {
	  int global_node_id = ns_node_list[j];
	  int local_node_id = global2local_node_map[global_node_id];
	  if (!local_node_id) continue;  /* this node was not read in */
	  
          mv = MESH_VertexFromID(mesh,local_node_id);
	
          /* Set attribute value for this node */
	
          MEnt_Set_AttVal(mv,nodesetatt,nodeset_ids[i],0.0,NULL);

          /* Add node to a nodeset */

          MSet_Add(nodeset,mv);
        }
      
        free(ns_node_list);
      }
    
      free(nodeset_ids);
      
    }

    

    /* Read side sets */
    
    if (nsidesets_total) {
      
      sideset_ids = (int *) malloc(nsidesets_total*sizeof(int));
     
#ifdef EXODUS_6_DEPRECATED 
      status = ex_get_side_set_ids(exoid,sideset_ids);
#else
      status = ex_get_ids(exoid, EX_SIDE_SET, sideset_ids);
#endif
      if (status < 0) {
	sprintf(mesg,
		"Error reading sideset IDs in Exodus II file %s\n",
		filename);
	MSTK_Report(funcname,mesg,MSTK_FATAL);
      }
      
      
      for (int i = 0; i < nsidesets_total; i++) {

	sprintf(sidesetname,"sideset_%-d",sideset_ids[i]);
	
	sidesetatt = MAttrib_New(mesh,sidesetname,INT,mstk_elem_type-1);
	sideset = MSet_New(mesh,sidesetname,mstk_elem_type-1);
      
#ifdef EXODUS_6_DEPRECATED
	status = ex_get_side_set_param(exoid,sideset_ids[i],&num_sides_in_set,
				       &num_df_in_set);
#else
	status = ex_get_set_param(exoid, EX_SIDE_SET, sideset_ids[i],
				  &num_sides_in_set, &num_df_in_set);
#endif
	
	ss_elem_list = (int *) malloc(num_sides_in_set*sizeof(int));
	ss_side_list = (int *) malloc(num_sides_in_set*sizeof(int));

#ifdef EXODUS_6_DEPRECATED	
	status = ex_get_side_set(exoid,sideset_ids[i],ss_elem_list,ss_side_list);
#else
	status = ex_get_set(exoid, EX_SIDE_SET, sideset_ids[i], ss_elem_list,
			    ss_side_list);
#endif
	if (status < 0) {
	  sprintf(mesg,
		  "Error reading elements in sideset %-d in Exodus II file %s\n",sideset_ids[i],filename);
	  MSTK_Report(funcname,mesg,MSTK_FATAL);
	}
	
	  
	for (int j = 0; j < num_sides_in_set; j++) {
	  int global_elem_id = ss_elem_list[j];
	  int local_elem_id = global2local_elem_map[global_elem_id];
	  if (!local_elem_id) continue;  /* this elem was not read in */

	  if (mstk_elem_type == MFACE) {
	    MFace_ptr mf = MESH_FaceFromID(mesh,local_elem_id);
	      
	    List_ptr fedges = MF_Edges(mf,1,0);
	    MEdge_ptr me = List_Entry(fedges,ss_side_list[j]-1);
	    List_Delete(fedges);
	      
	    /* Set attribute value for this edge */	      
	    MEnt_Set_AttVal(me,sidesetatt,sideset_ids[i],0.0,NULL);
	      
	    /* Add the edge to a set */	      
	    MSet_Add(sideset,me);
	      
	  } else {
	    MRegion_ptr mr = MESH_RegionFromID(mesh,local_elem_id);
	      
	    List_ptr rfaces = MR_Faces(mr);
	    MFace_ptr mf = List_Entry(rfaces,ss_side_list[j]-1);
	    List_Delete(rfaces);
	      
	    /* Set attribute value for this edge */
	    MEnt_Set_AttVal(mf,sidesetatt,sideset_ids[i],0.0,NULL);
	      
	    /* Add the face to a set */
	    MSet_Add(sideset,mf);
	  }
	}
	
	free(ss_elem_list);
	free(ss_side_list);
      }

      free(sideset_ids);
    }


    /* read element sets */


    if (nelemsets_total) {

      int *elemset_ids = (int *) malloc(nelemsets_total*sizeof(int));
      status = ex_get_ids(exoid,EX_ELEM_SET,elemset_ids);
      if (status < 0) {
	sprintf(mesg,
		"Error reading element set IDs in Exodus II file %s\n",
		filename);
	MSTK_Report(funcname,mesg,MSTK_FATAL);
      }
      
      
      for (int i = 0; i < nelemsets_total; i++) {

	int num_elems_in_set;
	char elemsetname[256];
	sprintf(elemsetname,"elemset_%-d",elemset_ids[i]);
	
	MSet_ptr elemset = MSet_New(mesh,elemsetname,mstk_elem_type);
      
	status = ex_get_set_param(exoid, EX_ELEM_SET, elemset_ids[i],
				  &num_elems_in_set, &num_df_in_set);
	
	int *es_elem_list = (int *) malloc(num_elems_in_set*sizeof(int));
	
	status = ex_get_set(exoid, EX_ELEM_SET, elemset_ids[i], es_elem_list,
			    NULL);
	if (status < 0) {
	  sprintf(mesg,
		  "Error reading elements in elemset %-d in Exodus II file %s\n",elemset_ids[i],filename);
	  MSTK_Report(funcname,mesg,MSTK_FATAL);
	}
	
		  
	/* Add the faces to the set */

	for (int j = 0; j < num_elems_in_set; j++) {
	  int global_elem_id = es_elem_list[j];
	  int local_elem_id = global2local_elem_map[global_elem_id];
	  if (!local_elem_id) continue;  /* this element was not read in */
	  
	  if (mstk_elem_type == MFACE) {
	    MFace_ptr mf = MESH_FaceFromID(mesh,local_elem_id);
	    MSet_Add(elemset,mf);
	  } else {
	    MRegion_ptr mr = MESH_RegionFromID(mesh,local_elem_id);
	    MSet_Add(elemset,mr);
	  }
	}	
	free(es_elem_list);
      }
      free(elemset_ids);
    }
    
    
    /* read element number map - interpret it as the GLOBAL ID of the element 
       for multi-processor runs */
    
    elem_map = (int *) malloc(nelems_total*sizeof(int));
    
#ifdef EXODUS_6_DEPRECATED
    status = ex_get_elem_num_map(exoid, elem_map);
#else
    status = ex_get_id_map(exoid, EX_ELEM_MAP, elem_map);
#endif
    if (status < 0) {
      sprintf(mesg,"Error reading element map in Exodus II file %s\n",filename);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }
    
    nmapatt = MAttrib_New(mesh, "elem_map", INT, mstk_elem_type);
    
    if (mstk_elem_type == MFACE) {
      for (int i = 0; i < nelems; i++) {
	MFace_ptr mf = MESH_Face(mesh, i);
	int global_elem_id = local2global_elem_map[i];
	MEnt_Set_AttVal(mf, nmapatt, elem_map[global_elem_id-1], 0.0, NULL);
	
#ifdef MSTK_HAVE_MPI
	MF_Set_GlobalID(mf,elem_map[global_elem_id-1]);
	MF_Set_MasterParID(mf,rank);
#endif
      }
    } else {
      for (int i = 0; i < nelems; i++) {
	MRegion_ptr mr = MESH_Region(mesh, i);
	int global_elem_id = local2global_elem_map[i];
	MEnt_Set_AttVal(mr, nmapatt, elem_map[global_elem_id-1], 0.0, NULL);
	
#ifdef MSTK_HAVE_MPI
	MR_Set_GlobalID(mr,elem_map[global_elem_id-1]);
	MR_Set_MasterParID(mr,rank);
#endif
      }
    }

    free(elem_map);


    /* Read in fields on elements and nodes */
    int nelemvars=0, nnodevars=0;
    int time_step = 1;
    char veckey[16] = "_veccomp";
    int keylen = strlen(veckey);


    /* How many variables are there on elements */

    status = ex_get_variable_param(exoid, EX_ELEM_BLOCK, &nelemvars);
    if (status < 0) {
      sprintf(mesg, "Error reading element variables in Exodus II file %s\n",
	      filename);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }

    if (nelemvars) {

      /* What are their names */

      char **elvarnames = (char **) malloc(nelemvars*sizeof(char *));
      for (int i = 0; i < nelemvars; i++)
        elvarnames[i] = (char *) malloc(256*sizeof(char));

      status = ex_get_variable_names(exoid, EX_ELEM_BLOCK, nelemvars,
				     elvarnames);
      if (status < 0) {
        sprintf(mesg,
		"Error reading element variable names in Exodus II file %s\n",
		filename);
        MSTK_Report(funcname,mesg,MSTK_FATAL);
      }

      /* Work through all the variables */
  
      int varindex;
      for (varindex = 0; varindex < nelemvars; varindex++) {
        char varname[256];
        int namelen;

        strcpy(varname,elvarnames[varindex]);
        namelen = strlen(varname);

        /* Check if this variable name contains the string "_veccomp".
	   This is a special keyword MSTK uses to understand that this
	   variable is part of a set of vector components and that they
	   should be imported together as a vector or tensor
	   attribute. The components will have "_veccomp0", "_veccomp1",
	   "_veccomp2" etc. */

        char *keystring = strstr(varname,veckey);
        if (keystring) {

	  /* Check if this is the first component of a vector field by
	     checking if the last character is "0"; if so, discover the
	     rest of the components and aggregate them into a vector or
	     tensor attribute.  If the last digit is something other than
	     0, we can ignore it since it would have been taken care of
	     when processing component '0'. THE ASSUMPTION HERE IS THAT
	     THE REMAINING COMPONENTS ARE IN THE IMMEDIATELY SUCCEEDING
	     INDICES. */

	  if (keystring[keylen] == '0') {
	    int ncomp = 1;
	    int varindex2;

	    /* If all but the last character are the same, this is a
	       component of the same vector */
	    for (varindex2 = varindex+1; varindex2 < nelemvars; varindex2++)
	      if (strncmp(varname,elvarnames[varindex2],namelen-1) == 0)
		ncomp++;
	      else
		break; /* chain is broken */

	    char attname[256];
	    strcpy(attname,varname);
	    char *tempstr = strstr(attname,veckey);
	    tempstr[0] = '\0';

	    MAttType atttype = (ncomp == ndim) ? VECTOR : TENSOR;
	    MAttrib_ptr mattrib = MAttrib_New(mesh,attname,atttype,
					      mstk_elem_type,ncomp);

	    double **elem_var_vals = (double **) malloc(ncomp*sizeof(double *));
	    for (int j = 0; j < ncomp; j++)
	      elem_var_vals[j] = malloc(nelems*sizeof(double));

	    int ibeg = 0;
	    while (ibeg < nelems-1) {
	      int iend = ibeg;
	      while (elem_ids[iend+1] == elem_ids[iend]+1) iend++;
	      int nelems_cur = iend - ibeg + 1;

	      int begel = elem_ids[ibeg];
	      
	      for (int b = 0; b < nelblock_total; b++) {

		if (begel >= elblock_min[b] && begel <= elblock_max[b]) {

		  /* found block containing starting elem 'begel'*/
	  
		  int jbeg = begel - elblock_min[b];
		  int endel = elem_ids[iend];
		  if (elblock_max[b] < endel)
		    endel = elblock_max[b];  /* range spans multiple blocks */
		  int nelems_cur = endel - begel + 1;

		  for (int k = 0; k < ncomp; k++) {
		    varindex2 = varindex + k;
		    /* ex_get_partial_var is using indices starting from 1 */
		    status = ex_get_partial_var(exoid, time_step, EX_ELEM_BLOCK,
						varindex2+1, elblock_ids[b],
						jbeg+1, nelems_cur,
						elem_var_vals[k] + begel);
		    if (status < 0) {
		      sprintf(mesg,
			      "Error reading element variables in Exodus II file %s\n",
			      filename);
		      MSTK_Report(funcname,mesg,MSTK_FATAL);
		    }
		  }

		  for (int j = 0; j < nelems_cur; j++) {
		    double *pval = malloc(ncomp*sizeof(double)); // freed by MESH_Delete
		    for (int k = 0; k < ncomp; k++)
		      pval[k] = elem_var_vals[k][jbeg+j];

		    int localid = global2local_elem_map[begel+j];
		    MEntity_ptr ment = (mstk_elem_type == MREGION) ?
		      MESH_Region(mesh,localid) :
		      MESH_Face(mesh,localid);
		
		    MEnt_Set_AttVal(ment,mattrib,0,0.0,pval);
		  }
		}
	      }
	    }
	    
	    for (int k = 0; k < ncomp; k++)
	      free(elem_var_vals[k]);
	    free(elem_var_vals);
	  } /* if its component 0 of a vector */
	  
	} /* if its a vector */
	else { /* scalar variable */
	  
	  MAttrib_ptr mattrib =
	    MAttrib_New(mesh,varname,DOUBLE,mstk_elem_type);
	  
	  double *elem_var_vals = (double *) malloc(max_el_in_blk*sizeof(double));
	    
	  int ibeg = 0;
	  while (ibeg < nelems-1) {
	    int iend = ibeg;
	    while (elem_ids[iend+1] == elem_ids[iend]+1) iend++;
	    int nelems_cur = iend - ibeg + 1;
	    
	    int begel = elem_ids[ibeg];
	    
	    for (int b = 0; b < nelblock_total; b++) {
	      if (begel >= elblock_min[b] && begel <= elblock_max[b]) {
		
		/* found block containing starting elem 'begel'*/
		
		int jbeg = begel - elblock_min[b];
		int endel = elem_ids[iend];
		if (elblock_max[b] < endel)
		  endel = elblock_max[b];  /* range spans multiple blocks */
		int nelems_cur = endel - begel + 1;
		
		/* ex_get_partial_var is using indices starting from 1 */
		status = ex_get_partial_var(exoid, time_step, EX_ELEM_BLOCK,
					    varindex+1, elblock_ids[b], jbeg+1,
					    nelems_cur, elem_var_vals+begel);
		if (status < 0) {
		  sprintf(mesg,
			  "Error reading element variables in Exodus II file %s\n",
			  filename);
		  MSTK_Report(funcname,mesg,MSTK_FATAL);
		}
		
		for (int j = 0; j < nelems_cur; j++) {
		  int localid = global2local_elem_map[begel+j];
		  MEntity_ptr ment = (mstk_elem_type == MREGION) ?
		    MESH_Region(mesh,localid) : MESH_Face(mesh,localid);
		  
		  MEnt_Set_AttVal(ment,mattrib,0,elem_var_vals[jbeg+j],NULL);
		}
	      }
	    }
	  }
	    
	  free(elem_var_vals);	  
        }  /* If its a scalar variable */
	
      } /* for each element variable */

      for (int i = 0; i < nelemvars; i++)
        free(elvarnames[i]);
      free(elvarnames);
    } /* if (nelemvars) */

    
    /* Read in variables associated with nodes */
    
    /* How many variables are there on nodes */

    status = ex_get_variable_param(exoid, EX_NODAL, &nnodevars);
    if (status < 0) {
      sprintf(mesg, "Error reading node variables in Exodus II file %s\n",
	      filename);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }

    if (nnodevars) {
      /* What are their names */

      char **nodevarnames = (char **) malloc(nnodevars*sizeof(char *));
      for (int i = 0; i < nnodevars; i++)
        nodevarnames[i] = (char *) malloc(256*sizeof(char));

      status = ex_get_variable_names(exoid, EX_NODAL, nnodevars, nodevarnames);
      if (status < 0) {
        sprintf(mesg,
		"Error reading node variable names in Exodus II file %s\n",
		filename);
        MSTK_Report(funcname,mesg,MSTK_FATAL);
      }

      /* Work through all the variables */

      int varindex;
      for (varindex = 0; varindex < nnodevars; varindex++) {
        char varname[256];
        int namelen;

        strcpy(varname,nodevarnames[varindex]);
        namelen = strlen(varname);

        /* Check if this variable name contains the string "_veccomp" at
	   the end.  This is a special keyword MSTK uses to understand
	   that this variable is part of a set of vector components and
	   that they should be imported together as vector or tensor
	   attribute. The components will have "_veccomp0", "_veccomp1",
	   "_veccomp2" etc. */

        char *keystring = strstr(varname,veckey);
        if (keystring) {

	  /* Check if this is the first component of a vector field by
	     checking if the last character is "0"; if so, discover the
	     rest and aggregate them into a vector or tensor attribute.
	     If the last digit is something other than 0, we can ignore it
	     since it would have been taken care of when processing
	     component '0'. THE ASSUMPTION HERE IS THAT THE REMAINING
	     COMPONENTS ARE IN THE IMMEDIATELY SUCCEEDING INDICES. */

	  if (keystring[keylen] == '0') {
	    int ncomp = 1;
	    int varindex2;

	    /* If all but the last character are the same, this is a
	       component of the same vector */
	    for (varindex2 = varindex+1; varindex2 < nnodevars; varindex2++)
	      if (strncmp(varname,nodevarnames[varindex2],namelen-1) == 0)
		ncomp++;
	      else
		break; /* chain is broken */

	    char attname[256];
	    strcpy(attname,varname);
	    char *tempstr = strstr(attname,veckey);
	    tempstr[0] = '\0';

	    MAttType atttype = (ncomp == ndim) ? VECTOR : TENSOR;
	    MAttrib_ptr mattrib = MAttrib_New(mesh,attname,atttype,MVERTEX,ncomp);

	    int nodeid = 0;

	    double **node_var_vals = (double **) malloc(ncomp*sizeof(double *));
	    for (int k = 0; k < ncomp; k++)
	      node_var_vals[k] = (double *) malloc(nnodes2read*sizeof(double));

	    int ibeg = 0;
	    while (ibeg < nnodes2read) {
	      int iend = ibeg;
	      while (nodes2read[iend+1] == nodes2read[iend]+1) iend++;
	      int nnodes_cur = iend - ibeg + 1;

	      for (int k = 0; k < ncomp; k++) {
		varindex2 = varindex + k;
		
		/* ex_get_partial_var is using indices starting from 1 */
		status = ex_get_partial_var(exoid, time_step, EX_NODAL,
					    varindex2+1, 1, ibeg+1, nnodes_cur,
					    node_var_vals[k] + ibeg);
		if (status < 0) {
		  sprintf(mesg,
			  "Error reading node variables in Exodus II file %s\n",
			  filename);
		  MSTK_Report(funcname,mesg,MSTK_FATAL);
		}
	      }
	    }

	    for (int j = 0; j < nnodes2read; j++) {
	      double *pval = malloc(ncomp*sizeof(double)); // freed by MESH_DELETE
	      for (int k = 0; k < ncomp; k++)
		pval[k] = node_var_vals[k][j];
	      MEntity_ptr ment = MESH_Vertex(mesh,j);
	      MEnt_Set_AttVal(ment,mattrib,0,0.0,pval);
	    }

	    for (int k = 0; k < ncomp; k++)
	      free(node_var_vals[k]);
	    free(node_var_vals);

	  } /* if its component 0 of a vector */

        } /* if its a vector */
        else { /* scalar variable */

	  MAttrib_ptr mattrib = MAttrib_New(mesh,varname,DOUBLE,MVERTEX);

	  double *node_var_vals = (double *) malloc(nnodes2read*sizeof(double));

	  int ibeg = 0;
	  while (ibeg < nnodes2read) {
	    int iend = ibeg;
	    while (nodes2read[iend+1] == nodes2read[iend]+1) iend++;
	    int nnodes_cur = iend - ibeg + 1;
	    
	    /* ex_get_partial_conn is using indices starting from 1 */
	    status = ex_get_partial_var(exoid, time_step, EX_NODAL,
					varindex+1, 1, ibeg+1, nnodes_cur,
					node_var_vals+ibeg);
            if (status < 0) {
              sprintf(mesg, "Error reading node variables in Exodus II file %s\n",
                      filename);
              MSTK_Report(funcname,mesg,MSTK_FATAL);
            }

            for (int j = 0; j < nnodes2read; j++) {
              MEntity_ptr ment = MESH_Vertex(mesh,j);
              MEnt_Set_AttVal(ment,mattrib,0,node_var_vals[j],NULL);
            }

            free(node_var_vals);
	  }
	} /* If its a scalar variable */
      } /* for each node variable */
  
      for (int i = 0; i < nnodevars; i++)
	free(nodevarnames[i]);
      free(nodevarnames);
    } /* if (nnodevars) */

    ex_close(exoid);

    free(elblock_ids);
    for (int i = 0; i < nelblock_total; i++) free(elblock_names[i]);
    free(elblock_names);
    free(elblock_nelems);
    
    return 1;
  }


#ifdef __cplusplus
}
#endif
