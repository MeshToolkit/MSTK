#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MSTK.h"
#include "MSTK_private.h"

#include "exodusII.h"
#ifdef EXODUSII_4
#include "exodusII_ext.h"
#endif

#include "uthash.h"


#ifdef __cplusplus
extern "C" {
#endif


  /* Read a single Exodus II file on one processor and return the
     element connectivity graph */

  /* create a hash function from n positive integers. Since these
     integers are the IDs of face vertices, we know that if the same
     set of integers is passed into the function repeatedly they will
     be in the same order or reverse order as the first time,
     i.e. given numbers a,b,c,d the first time, we will get a,b,c,d or
     b,c,d,a or c,d,a,b or d,a,b,c or a,d,c,b or b,a,d,c or c,b,a,d or
     d,c,b,a but never a,c,d,b. So to hash we start at the minimum
     number and then go in the direction of the smaller next
     number. For example, if we are given 22, 10, 24, 220, we order
     them as 10, 22, 220, 24 - this lets us get away with having just
     the first two numbers in ascending order instead of having to
     sort all the numbers. Hash table size should be of medium size
     according to the size of the mesh to avoid having a very large
     but sparsely populated or very small but very densely populated
     structure. It should also preferably be a prime number to avoid
     too many conflicts. 

     Can also consider hash functions shown here:

     https://cstheory.stackexchange.com/questions/3390/is-there-a-hash-function-for-a-collection-i-e-multi-set-of-integers-that-has

     and here:

     https://en.wikipedia.org/wiki/Universal_hashing#Hashing_vectors

*/

  unsigned int hash_ints(int n, unsigned int *nums, unsigned int hash_table_size) {
    /* Doesn't work well for 2 integers */
    int i, k=0;
    unsigned int p1 = 3;

    unsigned int h = 1;
    unsigned int min = 1e+9;
    for (i = 0; i < n; i++)
      if (nums[i] < min) {
        min = nums[i];
        k = i;
      }
    if (nums[(k+1)%n] < nums[(k-1+n)%n]) {
      for (i = 0; i < n; i++)
        h = h*p1 + nums[(k+i)%n];
    }
    else {
      for (i = 0; i < n; i++)
        h = h*p1 + nums[(k-i+n)%n];
    }
    
    /* h = (int) ((double)h/(n-1) + 0.5); */

    return h % hash_table_size;
  }

  /* https://stackoverflow.com/questions/1536393/good-hash-function-for-permutations */
  
  unsigned int hash_ints2(int n, unsigned int *nums, unsigned int hash_table_size) {
    int i, k=0;
    unsigned int R = 1234567891;

    unsigned int h = 1;
    unsigned int min = 1e+9;
    for (i = 0; i < n; i++)
      if (nums[i] < min) {
        min = nums[i];
        k = i;
      }
    if (nums[(k+1)%n] < nums[(k-1+n)%n]) {
      for (i = 0; i < n; i++)
        h *= R + 2*nums[(k+i)%n];
    }
    else {
      for (i = 0; i < n; i++)
        h *= R + 2*nums[(k-i+n)%n];
    }
    h /= 2;


    return h % hash_table_size;
  }


#define MAX_COLLISIONS 3
  
  typedef struct {
    int key;     /* hash key from sorted list of vertices*/ 
    int id;      /* actual id in list of faces */
    int elem[2]; /* ids of elements connected to face - can be -1 */
    int nv;      /* number of vertices */
    int *verts;  /* SORTED list of vertices IDs (we don't need to care
		  about the actual face topology or geometry here so
		  we can sort them) */
  } face_t;

  typedef struct {
    int key;           /* hash key that all faces in the list share */
    int nf;            /* number of faces in list - typically small 2 or so */
    face_t *faces[MAX_COLLISIONS];     /* list of faces */

    UT_hash_handle hh;
  } facelist_t;
  

  /* find a face defined by these vertices - assumption is that verts
     is sorted in increasing order */
  
  face_t *find_face(facelist_t *face_hashtable, int hash_table_size,
		    int nv, int *verts) {
    int hash_key = hash_ints2(nv, verts, hash_table_size);

    facelist_t *flist;
    HASH_FIND_INT(face_hashtable, &hash_key, flist);

    if (!flist)
      return NULL;
    else {
      if (flist->nf == 1)
	return flist->faces[0];
      else {
	int i;
	face_t *face = NULL;
	for (i = 0; i < flist->nf; i++) {
	  face = flist->faces[i];
	  if (face->nv != nv) continue;
	  int allmatch = 1;
	  int j;
	  for (j = 0; j < nv; j++)
	    if (face->verts[j] != verts[j]) {
	      allmatch = 0;
	      break;
	    }
	  if (allmatch) break;
	}
	return face;  /* could be null if a match was not found */
      }
    }

    return NULL;
  }

  /* Need address of pointer to face_hashtable since the pointer could
     be modified */
  
  int add_face(facelist_t **face_hashtable, int hash_table_size, face_t *face) {

    if (!find_face(*face_hashtable, hash_table_size, face->nv, face->verts)) {
      face->key = hash_ints2(face->nv, face->verts, hash_table_size);

      facelist_t *flist = (facelist_t *) malloc(sizeof(facelist_t));
      flist->key = face->key;
      flist->nf = 1;
      flist->faces[0] = face;
      
      HASH_ADD_INT(*face_hashtable, key, flist);
      return 1;
    }
    else {
      MSTK_Report("add_face","Face structure with vertices already in table",
		  MSTK_FATAL);
      return 0;
    }
  }


  
  int ExodusII_GetElementGraph(const char *filename, int *nelems, int **adjbeg,
			       int **adjelems) {

    char mesg[256], funcname[32]="ExodusII_GetElementGraph";
    char title[256];
    char elem_type[256], face_type[256];
    char **elem_blknames;
    int i, j, k, k1;
    int comp_ws = sizeof(double), io_ws = 0;
    int exoid=0, status;
    int ndim, nnodes, nelblock;
    int nedges, nedge_blk, nfaces, nface_blk;
    int *elem_blk_ids, *connect, *nnpe;
    int nelnodes, neledges, nelfaces;
    float version;

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

    ex_init_params exopar;

    FILE *fp;

    exoid = ex_open(filename, EX_READ, &comp_ws, &io_ws, &version);

    if (exoid < 0) {
      sprintf(mesg,"Cannot open/read Exodus II file %s\n",filename);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }

    status = ex_get_init_ext(exoid, &exopar);
    if (status < 0) {
      sprintf(mesg,"Error while reading Exodus II file %s\n",filename);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }
  
  
    strcpy(title,exopar.title);
    ndim = exopar.num_dim;
    nnodes = exopar.num_nodes;
    nfaces = exopar.num_face;
    nface_blk = exopar.num_face_blk;
    *nelems = exopar.num_elem;
    nelblock = exopar.num_elem_blk;

    if (ndim == 1)
      MSTK_Report(funcname,"Cannot read 1D meshes",MSTK_FATAL);


    /* Get element block IDs and names */
  
  
    elem_blk_ids = (int *) malloc(nelblock*sizeof(int));
    status = ex_get_ids(exoid, EX_ELEM_BLOCK, elem_blk_ids);
    if (status < 0) {
      sprintf(mesg,"Error while reading element block ids in Exodus II file %s\n",filename);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }
  
    elem_blknames = (char **) malloc(nelblock*sizeof(char *));
    for (i = 0; i < nelblock; i++)
      elem_blknames[i] = (char *) malloc(256*sizeof(char));
    status = ex_get_names(exoid, EX_ELEM_BLOCK, elem_blknames);
    if (status < 0) {
      sprintf(mesg,"Error while reading element block ids in Exodus II file %s\n",filename);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }

    
    /* Get an upper bound on the number of faces in the mesh by adding
       the explicitly specified number of faces and the number of
       faces of each element (without deduplication) */

    int max_faces = nfaces;
    int surf_elems = 0, solid_elems = 0;
    for (i = 0; i < nelblock; i++) {
      int nel_in_blk, nelnodes, neledges, nelfaces, natts;
	
      status = ex_get_block(exoid, EX_ELEM_BLOCK, elem_blk_ids[i], 
			    elem_type, &nel_in_blk, &nelnodes, 
			    &neledges, &nelfaces, &natts);
      if (status < 0) {
	sprintf(mesg,
		"Error while reading element block %s in Exodus II file %s\n",
		elem_blknames[i],filename);
	MSTK_Report(funcname,mesg,MSTK_FATAL);
      }
      
      if (strcmp(elem_type,"NFACED") == 0 || 
	  strcmp(elem_type,"nfaced") == 0) {  /* Polyhedral block */
	
	nnpe = (int *) malloc(nel_in_blk*sizeof(int));
	
	status = ex_get_entity_count_per_polyhedra(exoid, EX_ELEM_BLOCK,
						   elem_blk_ids[i], nnpe);
	if (status < 0) {
	  sprintf(mesg,
		  "Error while reading elem block info in Exodus II file %s\n",
		  filename);
	  MSTK_Report(funcname,mesg,MSTK_FATAL);
	}

	for (j = 0; j < nel_in_blk; j++)
	  max_faces += nnpe[j];
	free(nnpe);

	solid_elems = 1;
	
      } else if (strncasecmp(elem_type,"TET",3) == 0 ||
		 strncasecmp(elem_type,"WEDGE",5) == 0 ||
		 strncasecmp(elem_type,"HEX",3) == 0) {

	max_faces += 6*nel_in_blk;  /* overestimate to num faces of hex */

	solid_elems = 1;
      } 
      else if (strncasecmp(elem_type,"NSIDED",6) == 0) {

	nnpe = (int *) malloc(nel_in_blk*sizeof(int));

	status = ex_get_entity_count_per_polyhedra(exoid, EX_ELEM_BLOCK,
						   elem_blk_ids[i], nnpe);
	if (status < 0) {
	  sprintf(mesg,
		  "Error while reading element block info in Exodus II file %s\n",
		  filename);
	  MSTK_Report(funcname,mesg,MSTK_FATAL);
	}

	for (j = 0; j < nel_in_blk; j++)
	  max_faces += nnpe[j];
	free(nnpe);

	surf_elems = 1;
      }
      else if (strncasecmp(elem_type,"TRI",3) == 0 ||
	       strncasecmp(elem_type,"QUAD",4) == 0 ||
	       strncasecmp(elem_type,"SHELL3",6) == 0 ||
	       strncasecmp(elem_type,"SHELL4",6) == 0) {
	
	max_faces += 4*nel_in_blk;  /* overestimate to num faces of quad */

	surf_elems = 1;
      }
    }

    if (surf_elems && solid_elems)
      MSTK_Report("ExodusII_GetElementGraph",
		  "Cannot build element graph for non-manifold meshes",
		  MSTK_FATAL);


    
    /* Create face pointers array for easy look up by ID */

    face_t *faces = (face_t *) malloc(max_faces*sizeof(face_t));

    /* Create a hash table of face pointer lists for easy look up and
       deduplication using a hash key */
    
    int hash_table_size = 100*max_faces + 1;
    facelist_t *face_hashtable = NULL;    

    
    
    /* Read face blocks if they are there - these will be used by
       polyhedral elements - There should only be one per EXODUS II
       specification but this generalizes it */

    int max_face_id = 0;
    
    if (nface_blk) {
      int nfblock, ntotnodes, *face_blk_ids;    

      /* Face block IDs */
      
      face_blk_ids = (int *) malloc(nface_blk*sizeof(int));

      status = ex_get_ids(exoid, EX_FACE_BLOCK, face_blk_ids);
      if (status < 0) {
        sprintf(mesg,
		"Error while reading face block ids in Exodus II file %s\n",
		filename);
        MSTK_Report(funcname,mesg,MSTK_FATAL);
      }
      for (i = 0; i < nface_blk; i++) {

	/* Some info about this face block like the number of faces in
	   the block and the total number of node entries used to
	   describe the faces. Useful for sizing */
	
        status = ex_get_block(exoid, EX_FACE_BLOCK, face_blk_ids[i], face_type,
                              &nfblock, &ntotnodes, NULL, NULL, NULL);

        if (status < 0) {
          sprintf(mesg,
                  "Error while reading face block info in Exodus II file %s\n",
                  filename);
          MSTK_Report(funcname,mesg,MSTK_FATAL);
        }

	/* Face descriptions in terms of their nodes */
	
        connect = (int *) malloc(ntotnodes*sizeof(int));

        status = ex_get_conn(exoid, EX_FACE_BLOCK, face_blk_ids[i], connect,
                             NULL, NULL);

        if (status < 0) {
          sprintf(mesg,
                  "Error while reading face block info in Exodus II file %s\n",
                  filename);
          MSTK_Report(funcname,mesg,MSTK_FATAL);
        }

	/* Node count for each face so we know how to step through
	   connectivity array */
        nnpe = (int *) malloc(nfblock*sizeof(int));

        status = ex_get_entity_count_per_polyhedra(exoid, EX_FACE_BLOCK,
                                                   face_blk_ids[i], nnpe);
        if (status < 0) {
          sprintf(mesg,
                  "Error while reading face block info in Exodus II file %s\n",
                  filename);
          MSTK_Report(funcname,mesg,MSTK_FATAL);
        }

        int offset = 0;
        for (j = 0; j < nfblock; j++) {
	  if (max_face_id == max_faces)
	    MSTK_Report("ExodusII_GetElementGraph",
			"Face estimate fell short", MSTK_FATAL);

	  // add face to hash
	  face_t *newface = &(faces[max_face_id]);
	  newface->id = max_face_id++;
	  newface->nv = nnpe[j];
	  newface->verts = (int *) malloc(nnpe[j]*sizeof(int));
	  memcpy(newface->verts, connect+offset, nnpe[j]*sizeof(int));
	  qsort(newface->verts, newface->nv, sizeof(int), compareINT);
	  newface->elem[0] = newface->elem[1] = -1;

	  add_face(&face_hashtable, hash_table_size, newface);
	  
          offset += nnpe[j];
        }

        free(connect);
        free(nnpe);
      }
      free(face_blk_ids);
    }



    /* Now actually process the element blocks */
      
    int *nelemfaces = (int *) calloc(*nelems, sizeof(int));
    int **elem_faces = (int **) calloc(*nelems, sizeof(int *));
    int max_elem_id = 0;
    
    for (i = 0; i < nelblock; i++) {
      int nel_in_blk, nelnodes, neledges, nelfaces, natts;
	
      status = ex_get_block(exoid, EX_ELEM_BLOCK, elem_blk_ids[i], 
			    elem_type, &nel_in_blk, &nelnodes, 
			    &neledges, &nelfaces, &natts);
      if (status < 0) {
	sprintf(mesg,
		"Error while reading element block %s in Exodus II file %s\n",
		elem_blknames[i],filename);
	MSTK_Report(funcname,mesg,MSTK_FATAL);
      }
      
      if (strcmp(elem_type,"NFACED") == 0 || 
	  strcmp(elem_type,"nfaced") == 0) {  /* Polyhedral block */
	
	connect = (int *) malloc(nelfaces*sizeof(int));

	status = ex_get_conn(exoid, EX_ELEM_BLOCK, elem_blk_ids[i], NULL,
			     NULL, connect);
	if (status < 0) {
	  sprintf(mesg,
		  "Error while reading element block info in Exodus II file %s\n",
		  filename);
	  MSTK_Report(funcname,mesg,MSTK_FATAL);
	}

      
	nnpe = (int *) malloc(nel_in_blk*sizeof(int));
	
	status = ex_get_entity_count_per_polyhedra(exoid, EX_ELEM_BLOCK,
						   elem_blk_ids[i], nnpe);
	if (status < 0) {
	  sprintf(mesg,
		  "Error while reading elem block info in Exodus II file %s\n",
		  filename);
	  MSTK_Report(funcname,mesg,MSTK_FATAL);
	}

	for (j = 0; j < nel_in_blk; j++) {
	  nelemfaces[j] = nnpe[j];
	  elem_faces[max_elem_id+j] = (int *) malloc(nnpe[j]*sizeof(int));
	}

	int offset = 0;
	for (j = 0; j < nel_in_blk; j++) {

	  for (k = 0; k < nnpe[j]; k++) {
	    int fid = connect[offset+k]-1;
	      
	    elem_faces[max_elem_id][k] = fid;
	    if (faces[fid].elem[0] == -1)
	      faces[fid].elem[0] = max_elem_id;
	    else
	      faces[fid].elem[1] = max_elem_id;
	  }
	  max_elem_id++;
	  offset += nnpe[j];

	}  /* for (j = 0, nel_in_blk) */

	free(connect);
	free(nnpe);

      } else if (strncasecmp(elem_type,"TET",3) == 0 ||
		 strncasecmp(elem_type,"WEDGE",5) == 0 ||
		 strncasecmp(elem_type,"HEX",3) == 0) {
	int nrf, eltype;
	
	/* Get the connectivity of all elements in this block */
	
	connect = (int *) calloc(nelnodes*nel_in_blk,sizeof(int));
	
	status = ex_get_conn(exoid, EX_ELEM_BLOCK, elem_blk_ids[i], connect,
			     NULL, NULL);
	if (status < 0) {
	  sprintf(mesg,
		  "Error while reading element block %s in Exodus II file %s\n",
		  elem_blknames[i],filename);
	  MSTK_Report(funcname,mesg,MSTK_FATAL);
	}
	
	if (strncasecmp(elem_type,"TET",3) == 0) {
	  eltype = 0;
	  nrf = 4;
	  if (nelnodes > 4) {
	    MSTK_Report(funcname,"Higher order tets not supported",MSTK_WARN);
	    continue;
	  }
	}
	else if (strncasecmp(elem_type,"WEDGE",5) == 0) {
	  eltype = 1;
	  nrf = 5;
	  if (nelnodes > 6) {
	    MSTK_Report(funcname,"Higher order prisms/wedges not supported",MSTK_WARN);
	    continue;
	  }
	}
	else if (strncasecmp(elem_type,"HEX",3) == 0) {
	  eltype = 2;
	  nrf = 6;
	  if (nelnodes > 8) {
	    MSTK_Report(funcname,"Higher order hexes not supported",MSTK_WARN);
	    continue;
	  }
	}
	else {
	  sprintf(mesg,"Unrecognized or unsupported solid element type: %s",
		  elem_type);
	  MSTK_Report(funcname,mesg,MSTK_WARN);
	  continue;
	}
	
	for (j = 0; j < nel_in_blk; j++) {

	  elem_faces[max_elem_id] = (int *) malloc(nrf*sizeof(int));
	  nelemfaces[max_elem_id] = nrf;
	  
	  for (k = 0; k < nrf; k++) {
	    int fverts[4];

	    int nfv = exo_nrfverts[eltype][k];              
	    for (k1 = 0; k1 < nfv; k1++) {
	      int vloc = exo_rfverts[eltype][k][k1];
	      fverts[k1] = connect[nelnodes*j+vloc];
	    }
	    qsort(fverts, nfv, sizeof(int), compareINT);

	    face_t *face = find_face(face_hashtable, hash_table_size, nfv,
				     fverts);
	    if (!face) {  // no such "face"

	      if (max_face_id == max_faces)
		MSTK_Report("ExodusII_GetElementGraph",
			    "Face estimate fell short", MSTK_FATAL);

	      // add face to hash
	      face_t *newface = &(faces[max_face_id]);
	      newface->id = max_face_id++;
	      newface->nv = nfv;
	      newface->verts = (int *) malloc(nfv*sizeof(int));
	      memcpy(newface->verts, fverts, nfv*sizeof(int));

	      newface->elem[0] = max_elem_id;
	      newface->elem[1] = -1;
	      
	      add_face(&face_hashtable, hash_table_size, newface);
		
	      // update element-to-"face" list
	      elem_faces[max_elem_id][k] = newface->id;

	    } else {

	      if (face->elem[0] >= 0 && face->elem[1] >= 0)
		MSTK_Report(funcname,
			    "Hash table clash. Increase hash table size or Use a different hashing function",
			    MSTK_FATAL);
	      else {

		// update element-to-"face" list
		elem_faces[max_elem_id][k] = face->id;

		// update "face"-to-element info
		face->elem[1] = max_elem_id;
		
	      }
	    }
	  }
	  max_elem_id++;
	}
	
	free(connect);
	
      } 
      else if (strncasecmp(elem_type,"NSIDED",6) == 0) {
	
	/* In this case the nelnodes parameter is actually total number of 
	   nodes referenced by all the polygons (counting duplicates) */

	connect = (int *) malloc(nelnodes*sizeof(int));

	status = ex_get_conn(exoid, EX_ELEM_BLOCK, elem_blk_ids[i], connect,
			     NULL, NULL);
	if (status < 0) {
	  sprintf(mesg,
		  "Error while reading face block info in Exodus II file %s\n",
		  filename);
	  MSTK_Report(funcname,mesg,MSTK_FATAL);
	}

      
	nnpe = (int *) malloc(nel_in_blk*sizeof(int));

	status = ex_get_entity_count_per_polyhedra(exoid, EX_ELEM_BLOCK,
						   elem_blk_ids[i], nnpe);
	if (status < 0) {
	  sprintf(mesg,
		  "Error while reading element block info in Exodus II file %s\n",
		  filename);
	  MSTK_Report(funcname,mesg,MSTK_FATAL);
	}

	for (j = 0; j < nel_in_blk; j++) {
	  nelemfaces[j] = nnpe[j];
	  elem_faces[max_elem_id+j] = (int *) malloc(nnpe[j]*sizeof(int));
	}

	int offset = 0;
	for (j = 0; j < nel_in_blk; j++) {

	  for (k = 0; k < nnpe[j]; k++) {
	    int nodeids[2];
	    nodeids[0] = connect[offset + k];
	    nodeids[1] = connect[offset + (k+1)%nnpe[j]];
	    if (nodeids[1] < nodeids[0]) {
	      nodeids[1] = connect[offset + k];
	      nodeids[0] = connect[offset + (k+1)%nnpe[j]];
	    }
	    
	    face_t *face = find_face(face_hashtable, hash_table_size,
				     2, nodeids);
	    if (!face) {  // no such "face"

	      // create face and add to hash
	      if (max_face_id == max_faces)
		MSTK_Report("ExodusII_GetElementGraph",
			    "Face estimate fell short", MSTK_FATAL);
	  
	      face_t *newface = &(faces[max_face_id]);
	      newface->id = max_face_id++;
	      newface->nv = 2;
	      newface->verts = (int *) malloc(2*sizeof(int));
	      newface->verts[0] = nodeids[0];
	      newface->verts[1] = nodeids[2];

	      newface->elem[0] = max_elem_id;
	      newface->elem[1] = -1;
		
	      add_face(&face_hashtable, hash_table_size, newface);

	      // update element-to-"face" list
	      elem_faces[max_elem_id][k] = newface->id;

	    } else {

	      if (face->elem[0] >= 0 && face->elem[1] >= 0)
		MSTK_Report(funcname,
			    "Hash table clash. Increase hash table size or use a different hashing function",
			    MSTK_FATAL);
	      else {

		// update element-to-"face" list
		elem_faces[max_elem_id][k] = face->id;

		// update "face"-to-element info
		face->elem[1] = max_elem_id;
		
	      }
	    }
	  }
	  offset += nnpe[j];
	  max_elem_id++;
	}
	
	free(connect);
	free(nnpe);

      }
      else if (strncasecmp(elem_type,"TRI",3) == 0 ||
	       strncasecmp(elem_type,"QUAD",4) == 0 ||
	       strncasecmp(elem_type,"SHELL3",6) == 0 ||
	       strncasecmp(elem_type,"SHELL4",6) == 0) {
	
	/* Get the connectivity of all elements in this block */
	
	connect = (int *) calloc(nelnodes*nel_in_blk,sizeof(int));

	status = ex_get_conn(exoid, EX_ELEM_BLOCK, elem_blk_ids[i], connect,
			     NULL, NULL);
	if (status < 0) {
	  sprintf(mesg,
		  "Error while reading element block %s in Exodus II file %s\n",
		  elem_blknames[i],filename);
	  MSTK_Report(funcname,mesg,MSTK_FATAL);
	}
	
	int offset = 0;
	for (j = 0; j < nel_in_blk; j++) {

	  elem_faces[max_elem_id] = (int *) malloc(nelnodes*sizeof(int));
	  
	  // Add each "face" (edge for surface meshes) to hash table
	  // (if not already there). Update face-element connectivity
	  // info

	  for (k = 0; k < nelnodes; k++) {
	    int nodeids[2];
	    nodeids[0] = connect[offset + k];
	    nodeids[1] = connect[offset + (k+1)%nelnodes];
	    if (nodeids[1] < nodeids[0]) {
	      nodeids[1] = connect[offset + k];
	      nodeids[0] = connect[offset + (k+1)%nelnodes];
	    }

	    face_t *face =  find_face(face_hashtable, hash_table_size,
				      2, nodeids);

	    if (!face) {  // no such "face"

	      if (max_face_id == max_faces)
		MSTK_Report("ExodusII_GetElementGraph",
			    "Face estimate fell short", MSTK_FATAL);
	      
	      // add face to hash
	      face_t *newface = &(faces[max_face_id]);
	      newface->id = max_face_id++;
	      newface->nv = 2;
	      newface->verts = (int *) malloc(2*sizeof(int));
	      newface->verts[0] = nodeids[0];
	      newface->verts[1] = nodeids[1];
	      
	      // update "face"-to-element info
	      newface->elem[0] = max_elem_id;
	      newface->elem[1] = -1;
	      
	      add_face(&face_hashtable, hash_table_size, newface);

	      // update element-to-"face" list
	      elem_faces[max_elem_id][k] = newface->id;

	    } else {

	      if (face->elem[0] >= 0 && face->elem[1] >= 0)
		MSTK_Report(funcname,
			    "Hash table clash. Increase hash table size or Use a different hashing function",
			    MSTK_FATAL);
	      else {

		// update element-to-"face" list
		elem_faces[max_elem_id][k] = face->id;

		// update "face"-to-element info
		face->elem[1] = max_elem_id;
		
	      }
	    }
	  }
	  nelemfaces[max_elem_id] = nelnodes;
	  offset += nelnodes;
	  max_elem_id++;
	}
	
	free(connect);
      } /* if (strcmp(elem_type,"NFACED") ... else */
    } /* for (i = 0; i < nelblock; i++) */
    ex_close(exoid);

    free(elem_blk_ids);
    for (i = 0; i < nelblock; i++) free(elem_blknames[i]);
    free(elem_blknames);

    *nelems = max_elem_id;

    
    /* Now build the graph */
    *adjbeg = (int *) malloc((*nelems+1)*sizeof(int));
    int nalloc = 4*(*nelems);  /* start with assuming 4 nbrs per entity */
    *adjelems = (int *) malloc(nalloc*sizeof(int));

    int pos = 0;
    for (i = 0; i < *nelems; i++) {
      (*adjbeg)[i] = pos;
      
      int numf = nelemfaces[i];
      for (j = 0; j < numf; j++) {
	int fid = elem_faces[i][j];
	face_t *face = &(faces[fid]);
	int oppelem = -1;
	if (face->elem[0] == i)
	  oppelem = face->elem[1];
	else
	  oppelem = face->elem[0];
	if (oppelem != -1) {
	  if (pos == nalloc) {
	    nalloc *= 2;
	    *adjelems = (int *) realloc(*adjelems, nalloc*sizeof(int));
	  }
	  (*adjelems)[pos] = oppelem;
	  pos++;
	}
      }
    }
    (*adjbeg)[*nelems] = pos;
    
    for (i = 0; i < *nelems; i++)
      free(elem_faces[i]);
    free(elem_faces);
    free(nelemfaces);

    /* Free the hash table */
    facelist_t *current_flist, *tmp_flist;
    HASH_ITER(hh, face_hashtable, current_flist, tmp_flist) {
      HASH_DEL(face_hashtable, current_flist);
      free(current_flist);
    }

    for (i = 0; i < max_face_id; i++)
      free(faces[i].verts);
    free(faces);

    return 1;

  }


#ifdef __cplusplus
}
#endif
