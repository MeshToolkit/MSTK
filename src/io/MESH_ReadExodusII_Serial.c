#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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


  /* Read a single Exodus II (or Nemesis) file on one processor into MSTK */
  /* Nemesis files are the distributed versions of Exodus II files
     with the global IDs of elements and nodes encoded as element maps
     and node maps. Nemesis files also contain some additional info
     that can be read by the Nemesis API but we are choosing to ignore
     that for now. */

  /* Right now we are creating an attribute for material sets, side
     sets and node sets. We are ALSO creating meshsets for each of
     these entity sets. We are also setting mesh geometric entity IDs
     - Should we pick one of the first two? */

#define DEF_MAXFACES 20

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
     too many conflicts. */

  unsigned int hash_ints(int n, unsigned int *nums, unsigned int hash_table_size) {
    int i, k=0;
    unsigned int p1 = 3;

    unsigned int h = 0;
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


  int MESH_ReadExodusII_Serial(Mesh_ptr mesh, const char *filename, const int rank) {

    char mesg[256], funcname[32]="MESH_ImportFromExodusII";
    char title[256], sidesetname[256], nodesetname[256];
    char elem_type[256], face_type[256];
    char matsetname[256];
    char **elem_blknames;
    int i, j, k, k1;
    int comp_ws = sizeof(double), io_ws = 0;
    int exoid=0, status;
    int ndim, nnodes, nelems, nelblock, nnodesets, nsidesets;
    int nedges, nedge_blk, nfaces, nface_blk, nelemsets;
    int nedgesets, nfacesets, nnodemaps, nedgemaps, nfacemaps, nelemmaps;
    int *elem_blk_ids, *connect, *node_map, *elem_map, *nnpe;
    int nelnodes, neledges, nelfaces;
    int *nel_in_blk, max_el_in_blk=0, natts;
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
    MVertex_ptr mv, *fverts, *rverts;
    MEdge_ptr me;
    MFace_ptr mf;
    MRegion_ptr mr;
    MAttrib_ptr nmapatt=NULL, elblockatt=NULL, nodesetatt=NULL, sidesetatt=NULL;
    MSet_ptr faceset=NULL, nodeset=NULL, sideset=NULL, matset=NULL;
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
      sprintf(mesg,"Error while reading Exodus II file %s\n",filename);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }
  
  
    strcpy(title,exopar.title);
    ndim = exopar.num_dim;
    nnodes = exopar.num_nodes;
    nedges = exopar.num_edge;
    nedge_blk = exopar.num_edge_blk;
    nfaces = exopar.num_face;
    nface_blk = exopar.num_face_blk;
    nelems = exopar.num_elem;
    nelblock = exopar.num_elem_blk;
    nnodesets = exopar.num_node_sets;
    nedgesets = exopar.num_edge_sets;
    nfacesets = exopar.num_face_sets;
    nelemsets = exopar.num_elem_sets;
    nsidesets = exopar.num_side_sets;
    nnodemaps = exopar.num_node_maps;
    nedgemaps = exopar.num_edge_maps;
    nfacemaps = exopar.num_face_maps;
    nelemmaps = exopar.num_elem_maps;
  
  
  
  
  
    /* read the vertex information */
  
    xvals = (double *) malloc(nnodes*sizeof(double));
    yvals = (double *) malloc(nnodes*sizeof(double));  
    if (ndim == 2)
      zvals = (double *) calloc(nnodes,sizeof(double));
    else
      zvals = (double *) malloc(nnodes*sizeof(double));
  
    status = ex_get_coord(exoid, xvals, yvals, zvals);
    if (status < 0) {
      sprintf(mesg,"Error while reading vertex info in Exodus II file %s\n",filename);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }


    for (i = 0; i < nnodes; i++) {
    
      mv = MV_New(mesh);
    
      xyz[0] = xvals[i]; xyz[1] = yvals[i]; xyz[2] = zvals[i];
      MV_Set_Coords(mv,xyz);
    }

    free(xvals);
    free(yvals);
    free(zvals);
  

    /* read node number map - store it as an attribute to spit out later
       if necessary */
  
    node_map = (int *) malloc(nnodes*sizeof(int));
  
    status = ex_get_node_num_map(exoid, node_map);
    if (status < 0) {
      sprintf(mesg,"Error while reading node map in Exodus II file %s\n",filename);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }
  
    if (status == 0) {
    
      nmapatt = MAttrib_New(mesh, "node_map", INT, MVERTEX);
    
      for (i = 0; i < nnodes; i++) {
        mv = MESH_Vertex(mesh, i);
        MEnt_Set_AttVal(mv, nmapatt, node_map[i], 0.0, NULL);

#ifdef MSTK_HAVE_MPI
        MV_Set_GlobalID(mv, node_map[i]);
        MV_Set_MasterParID(mv, rank); /* This might get modified later */
#endif
      }
    
    }

    free(node_map);
  

    /* Read node sets */

    if (nnodesets) {
      nodeset_ids = (int *) malloc(nnodesets*sizeof(int));

      status = ex_get_node_set_ids(exoid,nodeset_ids);
      if (status < 0) {
        sprintf(mesg,
                "Error while reading nodeset IDs in Exodus II file %s\n",
                filename);
        MSTK_Report(funcname,mesg,MSTK_FATAL);
      }
    
    
      for (i = 0; i < nnodesets; i++) {

        sprintf(nodesetname,"nodeset_%-d",nodeset_ids[i]);

        nodesetatt = MAttrib_New(mesh,nodesetname,INT,MVERTEX);

        nodeset = MSet_New(mesh,nodesetname,MVERTEX);
    
        status = ex_get_node_set_param(exoid,nodeset_ids[i],&num_nodes_in_set,
                                       &num_df_in_set);
      
        ns_node_list = (int *) malloc(num_nodes_in_set*sizeof(int));
      
        status = ex_get_node_set(exoid,nodeset_ids[i],ns_node_list);
        if (status < 0) {
          sprintf(mesg,
                  "Error while reading nodes in nodeset %-d in Exodus II file %s\n",nodeset_ids[i],filename);
          MSTK_Report(funcname,mesg,MSTK_FATAL);
        }
      
      
        for (j = 0; j < num_nodes_in_set; j++) {
          mv = MESH_VertexFromID(mesh,ns_node_list[j]);
	
          /* Set attribute value for this node */
	
          MEnt_Set_AttVal(mv,nodesetatt,nodeset_ids[i],0.0,NULL);

          /* Add node to a nodeset */

          MSet_Add(nodeset,mv);
        }
      
        free(ns_node_list);
      }
    
      free(nodeset_ids);
      
    }



    /* Read face blocks if they are there - these will be used by
       polyhedral elements - There should only be one per EXODUS II
       specification and the following coding */

    if (nface_blk) {
      int nfblock, ntotnodes, *face_blk_ids;    

      face_blk_ids = (int *) malloc(nface_blk*sizeof(int));

      status = ex_get_ids(exoid, EX_FACE_BLOCK, face_blk_ids);
      if (status < 0) {
        sprintf(mesg,"Error while reading element block ids in Exodus II file %s\n",filename);
        MSTK_Report(funcname,mesg,MSTK_FATAL);
      }
      for (i = 0; i < nface_blk; i++) {

        status = ex_get_block(exoid, EX_FACE_BLOCK, face_blk_ids[i], face_type,
                              &nfblock, &ntotnodes, NULL, NULL, NULL);

        if (status < 0) {
          sprintf(mesg,
                  "Error while reading face block info in Exodus II file %s\n",
                  filename);
          MSTK_Report(funcname,mesg,MSTK_FATAL);
        }

        connect = (int *) malloc(ntotnodes*sizeof(int));

        status = ex_get_conn(exoid, EX_FACE_BLOCK, face_blk_ids[i], connect,
                             NULL, NULL);

        if (status < 0) {
          sprintf(mesg,
                  "Error while reading face block info in Exodus II file %s\n",
                  filename);
          MSTK_Report(funcname,mesg,MSTK_FATAL);
        }

      
        nnpe = (int *) malloc(nfblock*sizeof(int));

        status = ex_get_entity_count_per_polyhedra(exoid, EX_FACE_BLOCK,
                                                   face_blk_ids[i], nnpe);
        if (status < 0) {
          sprintf(mesg,
                  "Error while reading face block info in Exodus II file %s\n",
                  filename);
          MSTK_Report(funcname,mesg,MSTK_FATAL);
        }

        int max = 0;
        for (j = 0; j < nfblock; j++) 
          if (nnpe[j] > max) max = nnpe[j];


        fverts = (MVertex_ptr *) malloc(max*sizeof(MVertex_ptr));

        faceset = MSet_New(mesh,"faceset",MFACE);

        int offset = 0;
        for (j = 0; j < nfblock; j++) {
          mf = MF_New(mesh);

          for (k = 0; k < nnpe[j]; k++) 
            fverts[k] = MESH_VertexFromID(mesh,connect[offset+k]);

          MF_Set_Vertices(mf,nnpe[j],fverts);
          offset += nnpe[j];

          MSet_Add(faceset,mf);
        }

        free(fverts);
        free(connect);
        free(nnpe);
      }

    }



  
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

    nel_in_blk = (int *) malloc(nelblock*sizeof(int));
  
  

    int surf_elems=0, solid_elems=0;
    int mesh_type=0;  /* 1 - Surface, 2 - Solid, 3 - mixed */

    if (ndim == 1) {
    
      MSTK_Report(funcname,"Not implemented",MSTK_FATAL);
    
    }
    else if (ndim == 2) {

      mesh_type = 1;
    
      elblockatt = MAttrib_New(mesh, "elemblock", INT, MFACE);
    
    
      /* Read each element block */
    
      for (i = 0; i < nelblock; i++) {
      
        status = ex_get_block(exoid, EX_ELEM_BLOCK, elem_blk_ids[i], 
                              elem_type, &(nel_in_blk[i]), &nelnodes, 
                              &neledges, &nelfaces, &natts);
        if (status < 0) {
          sprintf(mesg,
                  "Error while reading element block %s in Exodus II file %s\n",
                  elem_blknames[i],filename);
          MSTK_Report(funcname,mesg,MSTK_FATAL);
        }

        if (max_el_in_blk < nel_in_blk[i])
          max_el_in_blk = nel_in_blk[i];

        sprintf(matsetname,"matset_%-d",elem_blk_ids[i]);
        matset = MSet_New(mesh,matsetname,MFACE);
      
        if (strcmp(elem_type,"NSIDED") == 0 || 
            strcmp(elem_type,"nsided") == 0) {

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

      
          nnpe = (int *) malloc(nel_in_blk[i]*sizeof(int));

          status = ex_get_entity_count_per_polyhedra(exoid, EX_ELEM_BLOCK,
                                                     elem_blk_ids[i], nnpe);
          if (status < 0) {
            sprintf(mesg,
                    "Error while reading face block info in Exodus II file %s\n",
                    filename);
            MSTK_Report(funcname,mesg,MSTK_FATAL);
          }

          int max = 0;
          for (j = 0; j < nel_in_blk[i]; j++) 
            if (nnpe[j] > max) max = nnpe[j];

          fverts = (MVertex_ptr *) malloc(max*sizeof(MVertex_ptr));

          int offset = 0;
          for (j = 0; j < nel_in_blk[i]; j++) {
            mf = MF_New(mesh);

            for (k = 0; k < nnpe[j]; k++) 
              fverts[k] = MESH_VertexFromID(mesh,connect[offset+k]);
	  
            MF_Set_Vertices(mf,nnpe[j],fverts);

            MF_Set_GEntID(mf, elem_blk_ids[i]);
            MF_Set_GEntDim(mf, 2);

            MSet_Add(matset,mf);
            MEnt_Set_AttVal(mf,elblockatt,elem_blk_ids[i],0.0,NULL);

            offset += nnpe[j];
          }
	
          free(fverts);
          free(connect);
          free(nnpe);

        }
        else {
	
          /* Get the connectivity of all elements in this block */

          connect = (int *) calloc(nelnodes*nel_in_blk[i],sizeof(int));
	
          status = ex_get_elem_conn(exoid, elem_blk_ids[i], connect);
          if (status < 0) {
            sprintf(mesg,"Error while reading element block %s in Exodus II file %s\n",elem_blknames[i],filename);
            MSTK_Report(funcname,mesg,MSTK_FATAL);
          }
	

          /* Create the MSTK faces */
	
          fverts = (MVertex_ptr *) calloc(nelnodes,sizeof(MVertex_ptr));
	
          for (j = 0; j < nel_in_blk[i]; j++) {
	  
            mf = MF_New(mesh);
	  
            for (k = 0; k < nelnodes; k++)
              fverts[k] = MESH_VertexFromID(mesh,connect[nelnodes*j+k]);
	  
            MF_Set_Vertices(mf,nelnodes,fverts);

            MF_Set_GEntID(mf, elem_blk_ids[i]);
            MF_Set_GEntDim(mf, 2);

            MSet_Add(matset,mf);
            MEnt_Set_AttVal(mf,elblockatt,elem_blk_ids[i],0.0,NULL);

          }
	
          free(fverts);
          free(connect);
	
        } /* if (strcmp(elem_type,"NSIDED") ... else */
      
      } /* for (i = 0; i < nelblock; i++) */
    

      /* Read side sets */

      if (nsidesets) {

        sideset_ids = (int *) malloc(nsidesets*sizeof(int));
      
        status = ex_get_side_set_ids(exoid,sideset_ids);
        if (status < 0) {
          sprintf(mesg,
                  "Error while reading sideset IDs in Exodus II file %s\n",
                  filename);
          MSTK_Report(funcname,mesg,MSTK_FATAL);
        }
      
      
        for (i = 0; i < nsidesets; i++) {

          sprintf(sidesetname,"sideset_%-d",sideset_ids[i]);
	
          sidesetatt = MAttrib_New(mesh,sidesetname,INT,MEDGE);
          sideset = MSet_New(mesh,sidesetname,MEDGE);
      
          status = ex_get_side_set_param(exoid,sideset_ids[i],&num_sides_in_set,
                                         &num_df_in_set);
	
          ss_elem_list = (int *) malloc(num_sides_in_set*sizeof(int));
          ss_side_list = (int *) malloc(num_sides_in_set*sizeof(int));
	
          status = ex_get_side_set(exoid,sideset_ids[i],ss_elem_list,ss_side_list);
          if (status < 0) {
            sprintf(mesg,
                    "Error while reading elements in sideset %-d in Exodus II file %s\n",sideset_ids[i],filename);
            MSTK_Report(funcname,mesg,MSTK_FATAL);
          }
	
	
          for (j = 0; j < num_sides_in_set; j++) {
            mf = MESH_FaceFromID(mesh,ss_elem_list[j]);
	  
            fedges = MF_Edges(mf,1,0);
            me = List_Entry(fedges,ss_side_list[j]-1);
            List_Delete(fedges);
	  
            /* Set attribute value for this edge */
	  
            MEnt_Set_AttVal(me,sidesetatt,sideset_ids[i],0.0,NULL);

            /* Add the edge to a set */

            MSet_Add(sideset,me);
	  
            /* Interpret sideset attribute as classification info for the edge */
	  
            /* TURN OFF SINCE SIDESETS THAT ARE INTERNAL TO DOMAINS TRIGGER
               SPURIOUS CLASSIFICATION - 05/04/2015 */
            /* ME_Set_GEntDim(me,1); */
            /* ME_Set_GEntID(me,sideset_ids[i]); */
          }
	
          free(ss_elem_list);
          free(ss_side_list);
	
        }

        free(sideset_ids);
      }


      /* read element sets */


      if (nelemsets) {

        int *elemset_ids = (int *) malloc(nelemsets*sizeof(int));
        status = ex_get_ids(exoid,EX_ELEM_SET,elemset_ids);
        if (status < 0) {
          sprintf(mesg,
                  "Error while reading element set IDs in Exodus II file %s\n",
                  filename);
          MSTK_Report(funcname,mesg,MSTK_FATAL);
        }
      
      
        for (i = 0; i < nelemsets; i++) {

          int num_elems_in_set;
          char elemsetname[256];
          sprintf(elemsetname,"elemset_%-d",elemset_ids[i]);
	
          MSet_ptr elemset = MSet_New(mesh,elemsetname,MFACE);
      
          status = ex_get_set_param(exoid,EX_ELEM_SET,elemset_ids[i],
                                    &num_elems_in_set,
                                    &num_df_in_set);
	
          int *es_elem_list = (int *) malloc(num_elems_in_set*sizeof(int));
	
          status = ex_get_set(exoid,EX_ELEM_SET,elemset_ids[i],es_elem_list,NULL);
          if (status < 0) {
            sprintf(mesg,
                    "Error while reading elements in elemset %-d in Exodus II file %s\n",elemset_ids[i],filename);
            MSTK_Report(funcname,mesg,MSTK_FATAL);
          }
	
		  
          /* Add the faces to the set */

          for (j = 0; j < num_elems_in_set; j++) {
            mf = MESH_FaceFromID(mesh,es_elem_list[j]);
            MSet_Add(elemset,mf);	  
          }	
          free(es_elem_list);
        }
        free(elemset_ids);
      }
    
    
      /* read element number map - interpret it as the GLOBAL ID of the element 
         for multi-processor runs */
    
      elem_map = (int *) malloc(nelems*sizeof(int));
    
      status = ex_get_elem_num_map(exoid, elem_map);
      if (status < 0) {
        sprintf(mesg,"Error while reading element map in Exodus II file %s\n",filename);
        MSTK_Report(funcname,mesg,MSTK_FATAL);
      }
    
      if (status == 0) {
      
        nmapatt = MAttrib_New(mesh, "elem_map", INT, MFACE);
      
        for (i = 0; i < nelems; i++) {
          mf = MESH_Face(mesh, i);
          MEnt_Set_AttVal(mf, nmapatt, elem_map[i], 0.0, NULL);

#ifdef MSTK_HAVE_MPI
          MF_Set_GlobalID(mf,elem_map[i]);
          MF_Set_MasterParID(mf,rank);
#endif
        }
      
      }

      free(elem_map);
      surf_elems = 1;
    
    }
    else if (ndim == 3) {


      /* Check if this is a strictly solid mesh, strictly surface mesh
         or mixed */

      int nface_est = 0;
      int max_nfv = 4; /* max number of vertices per face - for faces that
                          have more than 4 vertices, the faces are created
                          explicitly and we don't need to build a hash function
                          to check for their existence */
      for (i = 0; i < nelblock; i++) {
      
        status = ex_get_block(exoid, EX_ELEM_BLOCK, elem_blk_ids[i], 
                              elem_type, &(nel_in_blk[i]), &nelnodes, 
                              &neledges, &nelfaces, &natts);
        if (status < 0) {
          sprintf(mesg,
                  "Error while reading element block %s in Exodus II file %s\n",
                  elem_blknames[i],filename);
          MSTK_Report(funcname,mesg,MSTK_FATAL);
        }

        if (max_el_in_blk < nel_in_blk[i])
          max_el_in_blk = nel_in_blk[i];

        if (strncasecmp(elem_type,"NFACED",6) == 0 ||
            strncasecmp(elem_type,"TET",3) == 0 ||
            strncasecmp(elem_type,"WEDGE",5) == 0 ||
            strncasecmp(elem_type,"HEX",3) == 0) {

          solid_elems = 1;

          if (strncasecmp(elem_type,"NFACED",6) == 0)
            nface_est += nel_in_blk[i]*nelfaces;
          else
            nface_est += nel_in_blk[i]*6;

        }
        else if (strncasecmp(elem_type,"NSIDED",6) == 0 ||
                 strncasecmp(elem_type,"TRI",3) == 0 ||
                 strncasecmp(elem_type,"QUAD",4) == 0 ||
                 strncasecmp(elem_type,"SHELL",5) == 0)
          surf_elems = 1;      
      }

      if (surf_elems && !solid_elems)
        mesh_type = 1;
      else if (!surf_elems && solid_elems)
        mesh_type = 2;
      else if (surf_elems && solid_elems)
        mesh_type = 3;
      else {
        MSTK_Report(funcname,"Mesh has neither surface nor solid elements?",MSTK_FATAL);
      }
    
      if (mesh_type == 1) 
        elblockatt = MAttrib_New(mesh, "elemblock", INT, MFACE);
      else if (mesh_type == 2)
        elblockatt = MAttrib_New(mesh, "elemblock", INT, MREGION);
      else if (mesh_type == 3)
        elblockatt = MAttrib_New(mesh, "elemblock", INT, MALLTYPE);

#ifdef MSTK_USE_MARKERS
      int reliable_rfdirs = MSTK_GetMarker();  /* region knows its face dirs? */
#else
      MAttrib_ptr dirflagatt = MAttrib_New(mesh, "dirflag", INT, MREGION);
#endif

      for (i = 0; i < nelblock; i++) {
      
        status = ex_get_block(exoid, EX_ELEM_BLOCK, elem_blk_ids[i], 
                              elem_type, &(nel_in_blk[i]), &nelnodes, 
                              &neledges, &nelfaces, &natts);
        if (status < 0) {
          sprintf(mesg,
                  "Error while reading element block %s in Exodus II file %s\n",
                  elem_blknames[i],filename);
          MSTK_Report(funcname,mesg,MSTK_FATAL);
        }
      
        sprintf(matsetname,"matset_%-d",elem_blk_ids[i]);
        if (mesh_type == 1)
          matset = MSet_New(mesh,matsetname,MFACE);
        else if (mesh_type == 2)
          matset = MSet_New(mesh,matsetname,MREGION);
        else if (mesh_type == 3)
          matset = MSet_New(mesh,matsetname,MALLTYPE);

      
        if (strcmp(elem_type,"NFACED") == 0 || 
            strcmp(elem_type,"nfaced") == 0) {
	
          /* Exodus II format for polyhedra does not include face
           * directions so we have to try and divine them by various
           * means - hence this section is painfully long */

          /* In this case the nelnodes parameter is actually total number of 
             nodes referenced by all the polyhedra (counting duplicates) */

          connect = (int *) malloc(nelfaces*sizeof(int));

          status = ex_get_conn(exoid, EX_ELEM_BLOCK, elem_blk_ids[i], NULL, NULL,
                               connect);

          if (status < 0) {
            sprintf(mesg,
                    "Error while reading elem block info in Exodus II file %s\n",
                    filename);
            MSTK_Report(funcname,mesg,MSTK_FATAL);
          }

      
          nnpe = (int *) malloc(nel_in_blk[i]*sizeof(int));
	
          status = ex_get_entity_count_per_polyhedra(exoid, EX_ELEM_BLOCK,
                                                     elem_blk_ids[i], nnpe);
          if (status < 0) {
            sprintf(mesg,
                    "Error while reading elem block info in Exodus II file %s\n",
                    filename);
            MSTK_Report(funcname,mesg,MSTK_FATAL);
          }

          int max = 0;
          for (j = 0; j < nel_in_blk[i]; j++) 
            if (nnpe[j] > max) max = nnpe[j];

          MFace_ptr *rfarr = (MFace_ptr *) malloc(max*sizeof(MFace_ptr));
          int *rfdirs = (int *) malloc(max*sizeof(int *));

          int offset = 0;
          for (j = 0; j < nel_in_blk[i]; j++) {

            for (k = 0; k < nnpe[j]; k++) {
              /* Face of region */
              rfarr[k] = MSet_Entry(faceset,connect[offset+k]-1);
            
              /* Exodus II doesn't tell us the direction in which face
                 is used by region. */
              rfdirs[k] = 1;
              List_ptr fregs = MF_Regions(rfarr[k]);
              if (fregs) {
                if (List_Num_Entries(fregs) == 1) {
                  MRegion_ptr adjreg = List_Entry(fregs, 0);
                  int oppfdir = MR_FaceDir(adjreg, rfarr[k]);
                  rfdirs[k] = !oppfdir;
                }
                List_Delete(fregs);
              }
            }
          
            /* create the region (possibly with erroneous face
             * directions) and mark it as being in this block */
          
            mr = MR_New(mesh);            
            MR_Set_Faces(mr, nnpe[j], rfarr, rfdirs);            
            MR_Set_GEntID(mr, elem_blk_ids[i]);	  
            MR_Set_GEntDim(mr, 3);
          
            MEnt_Set_AttVal(mr,elblockatt,elem_blk_ids[i],0.0,NULL);
            MSet_Add(matset,mr);
          
            offset += nnpe[j];
          }  /* for (j = 0, nel_in_blk[i]) */

          free(rfarr); free(rfdirs);
          free(connect);
          free(nnpe);

        }
        else if (strncasecmp(elem_type,"TET",3) == 0 ||
                 strncasecmp(elem_type,"WEDGE",5) == 0 ||
                 strncasecmp(elem_type,"HEX",3) == 0) {
          int nrf, eltype;
          MFace_ptr face;
	
          /* Get the connectivity of all elements in this block */
	
          connect = (int *) calloc(nelnodes*nel_in_blk[i],sizeof(int));
	
          status = ex_get_elem_conn(exoid, elem_blk_ids[i], connect);
          if (status < 0) {
            sprintf(mesg,"Error while reading element block %s in Exodus II file %s\n",elem_blknames[i],filename);
            MSTK_Report(funcname,mesg,MSTK_FATAL);
          }
	
	
          /* Create the MSTK regions */
	
          rverts = (MVertex_ptr *) calloc(nelnodes,sizeof(MVertex_ptr));
	
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
	
          fverts = (MVertex_ptr *) malloc(4*sizeof(MVertex_ptr));

          MFace_ptr *rfarr = (MFace_ptr *) malloc(nrf*sizeof(MFace_ptr));
          int *rfdirs = (int *) malloc(nrf*sizeof(int *));

          for (j = 0; j < nel_in_blk[i]; j++) {
	  
            mr = MR_New(mesh);
	  
            /* Exodus II and MSTK node conventions are the same but the
               face conventions are not - so for type R1 we will use
               MR_Set_Vertices but for F1 instead of using
               MR_Set_Vertices we will create the faces individually and
               then do MR_Set_Faces */
	  
            for (k = 0; k < nelnodes; k++)
              rverts[k] = MESH_VertexFromID(mesh,connect[nelnodes*j+k]);
            
            if (reptype == F1 || reptype == F4) {
                        
              for (k = 0; k < exo_nrf[eltype]; k++) {
                face = NULL;

                int nfv = exo_nrfverts[eltype][k];              
                for (k1 = 0; k1 < nfv; k1++)
                  fverts[k1] = rverts[exo_rfverts[eltype][k][k1]];

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

                  rfarr[k] = face;
                  fregs = MF_Regions(face);
                  if (!fregs || !List_Num_Entries(fregs)) {
                    /* Comes here when the other region attached to this
                       face will be a polyhedral element that will be
                       created later. Because polyhedral elements need
                       to be described in terms of their faces, the
                       shared face explicitly described may have a
                       different vertex ordering than the implicit
                       ordering required by this standard element. In
                       this case, let the standard element ordering
                       rule */
                    MF_Set_Vertices(face, nfv, fverts);
                    rfdirs[k] = exo_rfdirs[eltype][k]; 
                  }
                  else {
                    if (List_Num_Entries(fregs) == 1) {
                      int dir;
                      MRegion_ptr freg;

                      freg = List_Entry(fregs,0);
                      dir = MR_FaceDir(freg,face);
                      /* since this is a standard type element, it has
                       * to follow its face direction according to the
                       * template and not derive its direction from the
                       * other region. If the other region is a standard
                       * element too, it should be consistent in a valid
                       * mesh. If its a polyhedral element, its face dir
                       * will be determined from this standard element
                       * later in the code but we will fix it here to be
                       * consistent with the standard element */
                      rfdirs[k] = exo_rfdirs[eltype][k];
                      if (dir != !rfdirs[k]) {
                        MF_Set_Vertices(face, nfv, fverts);
                        MR_Set_FaceDir(freg, face, !dir);
                      }
                    }
                    else if (List_Num_Entries(fregs) == 2) {
                      MSTK_Report(funcname,"Face already connected two faces",MSTK_FATAL);
                    }
                    List_Delete(fregs);
                  }
                }
                else {
                  face = MF_New(mesh);

                  MF_Set_Vertices(face,exo_nrfverts[eltype][k],fverts);
                  rfarr[k] = face;
                  rfdirs[k] = exo_rfdirs[eltype][k];
                }
              }
	  
              MR_Set_Faces(mr,nrf,rfarr,rfdirs);
#ifdef MSTK_USE_MARKERS
              MEnt_Mark(mr, reliable_rfdirs);
#else
              MEnt_Set_AttVal(mr, dirflagatt, 1, 0.0, NULL);
#endif
            }
            else if (reptype == R1 || reptype == R2) {

              MR_Set_Vertices(mr,nelnodes,rverts,0,NULL);

            }

            MR_Set_GEntID(mr, elem_blk_ids[i]);
            MR_Set_GEntDim(mr, 3);
          
            MEnt_Set_AttVal(mr,elblockatt,elem_blk_ids[i],0.0,NULL);
            MSet_Add(matset,mr);
          
          }
	
          free(fverts);
          free(rverts);
          free(rfarr);
          free(rfdirs);
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

      
          nnpe = (int *) malloc(nel_in_blk[i]*sizeof(int));

          status = ex_get_entity_count_per_polyhedra(exoid, EX_ELEM_BLOCK,
                                                     elem_blk_ids[i], nnpe);
          if (status < 0) {
            sprintf(mesg,
                    "Error while reading face block info in Exodus II file %s\n",
                    filename);
            MSTK_Report(funcname,mesg,MSTK_FATAL);
          }

          int max = 0;
          for (j = 0; j < nel_in_blk[i]; j++) 
            if (nnpe[j] > max) max = nnpe[j];

          fverts = (MVertex_ptr *) malloc(max*sizeof(MVertex_ptr));

          int offset = 0;
          for (j = 0; j < nel_in_blk[i]; j++) {
            mf = MF_New(mesh);

            for (k = 0; k < nnpe[j]; k++) 
              fverts[k] = MESH_VertexFromID(mesh,connect[offset+k]);
	  
            MF_Set_Vertices(mf,nnpe[j],fverts);

            MF_Set_GEntID(mf, elem_blk_ids[i]);
            MF_Set_GEntID(mf, 2);

            if (mesh_type == 1 || mesh_type == 3) {
              MEnt_Set_AttVal(mf,elblockatt,elem_blk_ids[i],0.0,NULL);
              MSet_Add(matset,mf);
            }

            offset += nnpe[j];
          }
	
          free(fverts);
          free(connect);
          free(nnpe);

        }
        else if (strncasecmp(elem_type,"TRI",3) == 0 ||
                 strncasecmp(elem_type,"QUAD",4) == 0 ||
                 strncasecmp(elem_type,"SHELL3",6) == 0 ||
                 strncasecmp(elem_type,"SHELL4",6) == 0) {
	
          /* Get the connectivity of all elements in this block */
	
          connect = (int *) calloc(nelnodes*nel_in_blk[i],sizeof(int));
	
          status = ex_get_elem_conn(exoid, elem_blk_ids[i], connect);
          if (status < 0) {
            sprintf(mesg,"Error while reading element block %s in Exodus II file %s\n",elem_blknames[i],filename);
            MSTK_Report(funcname,mesg,MSTK_FATAL);
          }
	
	
          /* Create the MSTK faces */
	
          fverts = (MVertex_ptr *) calloc(nelnodes,sizeof(MVertex_ptr));
	
          for (j = 0; j < nel_in_blk[i]; j++) {
	  
            mf = MF_New(mesh);
	  
            for (k = 0; k < nelnodes; k++)
              fverts[k] = MESH_VertexFromID(mesh,connect[nelnodes*j+k]);
	  
            MF_Set_Vertices(mf,nelnodes,fverts);

            MF_Set_GEntID(mf, elem_blk_ids[i]);
            MF_Set_GEntID(mf, 2);

            if (mesh_type == 1 || mesh_type == 3) {
              MEnt_Set_AttVal(mf,elblockatt,elem_blk_ids[i],0.0,NULL);
              MSet_Add(matset,mf);
            }

          }
	
          free(fverts);
          free(connect);
	
        } /* if (strcmp(elem_type,"NFACED") ... else */
      
      
      } /* for (i = 0; i < nelblock; i++) */

        
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
             * directions. We have to do geometric checks to find such a
             * pair */
          
            MRegion_ptr mrgood = NULL;
            idx = 0;
            while (!mrgood && ((mr = MESH_Next_Region(mesh, &idx)))) {
            
              double (*fxyz)[3] =
                  (double (*)[3]) malloc(MAXPV2*sizeof(double [3]));
            
              /* Geometric center of region */
              /* Average normal and geometric center of faces */
            
              List_ptr rfaces = MR_Faces(mr);
              int nrf = List_Num_Entries(rfaces);
            
              List_ptr rverts = List_New(10);
              int nrv = 0;
            
              double rcen[3];
              rcen[0] = rcen[1] = rcen[2] = 0.0;
              double fcen[MAXPF3][3], fnormal[MAXPF3][3];
            
              for (k = 0; k < nrf; k++) {
                MFace_ptr rf = List_Entry(rfaces, k);
              
                List_ptr fverts = MF_Vertices(rf, 1, 0);
                int nfv = List_Num_Entries(fverts);
              
                fcen[k][0] = fcen[k][1] = fcen[k][2] = 0.0;
                int k1;
                for (k1 = 0; k1 < nfv; k1++) {
                  int k2;
                  MVertex_ptr fv = List_Entry(fverts, k1);
                  MV_Coords(fv, fxyz[k1]);
                
                  if (!List_Contains(rverts, fv)) {
                    List_Add(rverts, fv);
                    for (k2 = 0; k2 < 3; k2++) rcen[k2] += fxyz[k1][k2];
                    nrv++;
                  }
                
                  for (k2 = 0; k2 < 3; k2++) fcen[k][k2] += fxyz[k1][k2];
                }
                for (k1 = 0; k1 < 3; k1++) fcen[k][k1] /= nfv;
                List_Delete(fverts);
              
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
              for (k = 0; k < 3; k++) rcen[k] /= nrv;
            
              int nfgood = 0;
              for (k = 0; k < nrf; k++) {
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

                for (k = 0; k < nrf; k++) {
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
            }  /* while (!mrgood && (mr = MESH_Next_Entry(mesh, &idx))) */

            if (!mrgood)
              MSTK_Report(funcname, "Could not find a single region whose faces could be reliably oriented", MSTK_FATAL);
          }  /* if (!List_Num_Entries(reglist)) */


          /* Starting from the reliable direction region, walk through
           * all other connected regions in the mesh */

          idx = 0;
          while ((mr = List_Next_Entry(reglist, &idx))) {
            List_ptr rfaces = MR_Faces(mr);
            int nrf = List_Num_Entries(rfaces);

            /* if we know the direction in which an adjacent reliable
             * region uses the common face, we can determine the
             * direction of use of the face in this region */             

            MFace_ptr fstart = NULL;
            for (k = 0; k < nrf; k++) {
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

            for (k = 0; k < nrf; ++k) {
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
    
    
      /* Read side sets */

      if (nsidesets) {

        if (mesh_type == 3)
          MSTK_Report(funcname,"Cannot handle sidesets in meshes with surface and solid elements",MSTK_FATAL);

        sideset_ids = (int *) malloc(nsidesets*sizeof(int));

        status = ex_get_side_set_ids(exoid,sideset_ids);
        if (status < 0) {
          sprintf(mesg,
                  "Error while reading sideset IDs %s in Exodus II file %s\n",
                  elem_blknames[i],filename);
          MSTK_Report(funcname,mesg,MSTK_FATAL);
        }
      
      
        for (i = 0; i < nsidesets; i++) {

          sprintf(sidesetname,"sideset_%-d",sideset_ids[i]);
	
          status = ex_get_side_set_param(exoid,sideset_ids[i],&num_sides_in_set,
                                         &num_df_in_set);
	
          ss_elem_list = (int *) malloc(num_sides_in_set*sizeof(int));
          ss_side_list = (int *) malloc(num_sides_in_set*sizeof(int));

          status = ex_get_side_set(exoid,sideset_ids[i],ss_elem_list,ss_side_list);
          if (status < 0) {
            sprintf(mesg,
                    "Error while reading elements in sideset %-d in Exodus II file %s\n",sideset_ids[i],filename);
            MSTK_Report(funcname,mesg,MSTK_FATAL);
          }
	
	

          if (mesh_type == 1) {

            sidesetatt = MAttrib_New(mesh,sidesetname,INT,MEDGE);
            sideset = MSet_New(mesh,sidesetname,MEDGE);          
      
            for (j = 0; j < num_sides_in_set; j++) {
              int eltype;
              List_ptr fedges;
            
              mf = MESH_FaceFromID(mesh,ss_elem_list[j]);
              if (!mf)
                MSTK_Report(funcname,"Cannot find element in sideset",MSTK_FATAL);
            
              fedges = MF_Edges(mf,1,0);
            
              int lfnum = ss_side_list[j]-1;
              
              /* Since we made the faces of the region according the
                 convention that Exodus II uses, the local face number
                 of the Exodus II element matches that of the MSTK element */
            
              me = List_Entry(fedges,lfnum);

              if (!me)
                MSTK_Report(funcname,"Could not find sideset edge",MSTK_FATAL);
              else {
                List_ptr efaces;
                efaces = ME_Faces(me);
                if (List_Num_Entries(efaces) != 1) {
                  if (MF_GEntID(List_Entry(efaces,0)) == 
                      MF_GEntID(List_Entry(efaces,1)))
                    MSTK_Report(funcname,
                                "Interior edge identified as being on the boundary",
                                MSTK_ERROR);
                }
                if (efaces) List_Delete(efaces);
              }
              List_Delete(fedges);
            
              /* Set attribute value for this face */
            
              MEnt_Set_AttVal(me,sidesetatt,sideset_ids[i],0.0,NULL);
            
              /* Add the face to a set */
            
              MSet_Add(sideset,me);
            
              /* Interpret sideset attribute as classification info for the edge */
              /* TURN OFF BECAUSE SIDESETS THAT ARE INTERNAL TO THE
                 DOMAIN TRIGGER SPURIOUS CLASSIFICATION - 05/04/2015 */
            
              /* ME_Set_GEntDim(me,1); */
              /* ME_Set_GEntID(me,sideset_ids[i]); */
            }
          }
          else if (mesh_type == 2) {

            sidesetatt = MAttrib_New(mesh,sidesetname,INT,MFACE);
            sideset = MSet_New(mesh,sidesetname,MFACE);
      
            for (j = 0; j < num_sides_in_set; j++) {
              int eltype;
            
              mr = MESH_RegionFromID(mesh,ss_elem_list[j]);
              if (!mr)
                MSTK_Report(funcname,"Could not find element in sideset",MSTK_FATAL);

              rfaces = MR_Faces(mr);       /* No particular order */
            
              int lfnum = ss_side_list[j]-1;
            
              /* Since we made the faces of the region according the
                 convention that Exodus II uses, the local face number
                 of the Exodus II element matches that of the MSTK element */
            
              mf = List_Entry(rfaces,lfnum);

              if (!mf)
                MSTK_Report(funcname,"Could not find sideset face",MSTK_FATAL);
              else {
                List_ptr fregs;
                fregs = MF_Regions(mf);
                if (List_Num_Entries(fregs) != 1) {
#ifdef DEBUG
                  if (MF_GEntID(List_Entry(fregs,0)) == 
                      MF_GEntID(List_Entry(fregs,1))) {
                    char mesg[256];
                    sprintf(mesg,"Sideset %-d contains non-boundary faces",sideset_ids[i]);
                    MSTK_Report(funcname,mesg,MSTK_WARN);
                  }
#endif
                }
                if (fregs) List_Delete(fregs);
              }
              List_Delete(rfaces);
            
              /* Set attribute value for this face */
            
              MEnt_Set_AttVal(mf,sidesetatt,sideset_ids[i],0.0,NULL);
            
              /* Add the face to a set */
            
              MSet_Add(sideset,mf);
            
              /* Interpret sideset attribute as classification info for the edge */
              /* TURN OFF BECAUSE SIDESETS THAT ARE INTERNAL TO THE
                 DOMAIN TRIGGER SPURIOUS CLASSIFICATION - 05/04/2015 */

              /* MF_Set_GEntDim(mf,2); */
              /* MF_Set_GEntID(mf,sideset_ids[i]); */
            }
          }
	
          free(ss_elem_list);
          free(ss_side_list);
	
        }

        free(sideset_ids);
      }

      /* Read element sets */

      if (nelemsets) {

        int *elemset_ids = (int *) malloc(nelemsets*sizeof(int));
        status = ex_get_ids(exoid,EX_ELEM_SET,elemset_ids);
        if (status < 0) {
          sprintf(mesg,
                  "Error while reading element set IDs in Exodus II file %s\n",
                  filename);
          MSTK_Report(funcname,mesg,MSTK_FATAL);
        }
      
      
        for (i = 0; i < nelemsets; i++) {

          int num_elems_in_set;
          char elemsetname[256];
          sprintf(elemsetname,"elemset_%-d",elemset_ids[i]);
	
          MSet_ptr elemset = MSet_New(mesh,elemsetname,MREGION);
      
          status = ex_get_set_param(exoid,EX_ELEM_SET,elemset_ids[i],
                                    &num_elems_in_set,
                                    &num_df_in_set);
	
          int *es_elem_list = (int *) malloc(num_elems_in_set*sizeof(int));
	
          status = ex_get_set(exoid,EX_ELEM_SET,elemset_ids[i],es_elem_list,NULL);
          if (status < 0) {
            sprintf(mesg,
                    "Error while reading elements in elemset %-d in Exodus II file %s\n",elemset_ids[i],filename);
            MSTK_Report(funcname,mesg,MSTK_FATAL);
          }
	
		  
          /* Add the elements to the set */

          if (mesh_type == 1) {
            for (j = 0; j < num_elems_in_set; j++) {
              mf = MESH_FaceFromID(mesh,es_elem_list[j]);
              MSet_Add(elemset,mf);	  
            }	
          }
          else if (mesh_type == 2) {
            for (j = 0; j < num_elems_in_set; j++) {
              mr = MESH_RegionFromID(mesh,es_elem_list[j]);
              MSet_Add(elemset,mr);	  
            }	
          }
          else {
            MSTK_Report(funcname,"Cannot read element sets in mixed surface/solid meshes",MSTK_WARN);
          }

          free(es_elem_list);
        }
        free(elemset_ids);
      }
    
    
      /* read element number map - store it as an attribute to spit out
         later if necessary */
    
      elem_map = (int *) malloc(nelems*sizeof(int));
    
      status = ex_get_elem_num_map(exoid, elem_map);
      if (status < 0) {
        sprintf(mesg,"Error while reading element map in Exodus II file %s\n",filename);
        MSTK_Report(funcname,mesg,MSTK_FATAL);
      }
    
      if (status == 0) {
      
        if (mesh_type == 1) {
          nmapatt = MAttrib_New(mesh, "elem_map", INT, MFACE);
	
          for (i = 0; i < nelems; i++) {
            mf = MESH_Face(mesh, i);
            MEnt_Set_AttVal(mf, nmapatt, elem_map[i], 0.0, NULL);
#ifdef MSTK_HAVE_MPI
            MF_Set_GlobalID(mf,elem_map[i]);
            MF_Set_MasterParID(mf,rank);
#endif
          }
        }
        else if (mesh_type == 2) {
          nmapatt = MAttrib_New(mesh, "elem_map", INT, MREGION);
	
          for (i = 0; i < nelems; i++) {
            mr = MESH_Region(mesh, i);
            MEnt_Set_AttVal(mr, nmapatt, elem_map[i], 0.0, NULL);
#ifdef MSTK_HAVE_MPI
            MR_Set_GlobalID(mr,elem_map[i]);
            MR_Set_MasterParID(mr,rank);
#endif
          }
        }
        else {
          MSTK_Report(funcname,
                      "Element maps not supported for mixed surface-solid meshes",
                      MSTK_WARN);
        }
      
      }

      free(elem_map);
    } /* ndim = 3 */



    /* Read in fields on elements and nodes */
    int nelemvars=0, nnodevars=0;
    int time_step = 1;
    char veckey[16] = "_veccomp";
    int keylen = strlen(veckey);


    /* How many variables are there on elements */

    status = ex_get_var_param(exoid, "e", &nelemvars);
    if (status < 0) {
      sprintf(mesg, 
              "Error while reading element variables in Exodus II file %s\n",
              filename);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }

    if (nelemvars) {

      /* What are their names */

      char **elvarnames = (char **) malloc(nelemvars*sizeof(char *));
      for (i = 0; i < nelemvars; i++)
        elvarnames[i] = (char *) malloc(256*sizeof(char));

      status = ex_get_var_names(exoid, "e", nelemvars, elvarnames);
      if (status < 0) {
        sprintf(mesg, 
                "Error while reading element variable names in Exodus II file %s\n",
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

            MAttrib_ptr mattrib;
            if (surf_elems) {
              MAttType atttype = (ncomp == 2) ? VECTOR : TENSOR;
              mattrib = MAttrib_New(mesh,attname,atttype,MFACE,ncomp);
            }
            else if (solid_elems) {
              MAttType atttype = (ncomp == 3) ? VECTOR : TENSOR;
              mattrib = MAttrib_New(mesh,attname,atttype,MREGION,ncomp);
            }
            else {
              MSTK_Report("MESH_ReadExodusII_Serial",
                          "Attribute export on wire meshes and node meshes not implemented",
                          MSTK_ERROR);
              continue;
            }

            int elid = 0;
            double **elem_var_vals = (double **) malloc(ncomp*sizeof(double *));
            for (j = 0; j < ncomp; j++)
              elem_var_vals[j] = malloc(max_el_in_blk*sizeof(double));

            for (i = 0; i < nelblock; i++) {

              for (j = 0; j < ncomp; j++) {
                varindex2 = varindex + j;
                status = ex_get_elem_var(exoid, time_step, varindex2+1,
                                         elem_blk_ids[i], nel_in_blk[i], 
                                         elem_var_vals[j]);
                if (status < 0) {
                  sprintf(mesg, 
                          "Error while reading element variables in Exodus II file %s\n",
                          filename);
                  MSTK_Report(funcname,mesg,MSTK_FATAL);
                }
              }

              for (k = 0; k < nel_in_blk[i]; k++) {
                double *pval = malloc(ncomp*sizeof(double)); // freed by MESH_Delete
                for (j = 0; j < ncomp; j++)
                  pval[j] = elem_var_vals[j][k];
                MEntity_ptr ment = solid_elems ? MESH_Region(mesh,elid) :
                    MESH_Face(mesh,elid);
                elid++;

                MEnt_Set_AttVal(ment,mattrib,0,0.0,pval);
              }

            } // for each element block

            for (j = 0; j < ncomp; j++)
              free(elem_var_vals[j]);
            free(elem_var_vals);

          } /* if its component 0 of a vector */

        } /* if its a vector */
        else { /* scalar variable */

          MAttrib_ptr mattrib;
          if (surf_elems)
            mattrib = MAttrib_New(mesh,varname,DOUBLE,MFACE);
          else if (solid_elems)
            mattrib = MAttrib_New(mesh,varname,DOUBLE,MREGION);
          else {
            MSTK_Report("MESH_ReadExodusII_Serial",
                        "Attribute export on wire meshes and node meshes not implemented",
                        MSTK_ERROR);
            continue;
          }

          int elid = 0;
          double *elem_var_vals = (double *) malloc(max_el_in_blk*sizeof(double));

          for (i = 0; i < nelblock; i++) {
            status = ex_get_elem_var(exoid, time_step, varindex+1, elem_blk_ids[i], 
                                     nel_in_blk[i], elem_var_vals);
            if (status < 0) {
              sprintf(mesg, 
                      "Error while reading element variables in Exodus II file %s\n",
                      filename);
              MSTK_Report(funcname,mesg,MSTK_FATAL);
            }

            for (k = 0; k < nel_in_blk[i]; k++) {
              MEntity_ptr ment = solid_elems ? MESH_Region(mesh,elid) :
                  MESH_Face(mesh,elid);

              MEnt_Set_AttVal(ment,mattrib,0,elem_var_vals[k],NULL);
              elid++;
            }

          } /* for each element block */

          free(elem_var_vals);

        } /* If its a scalar variable */
      } /* for each element variable */

      free(elem_blk_ids);
      for (i = 0; i < nelblock; i++) free(elem_blknames[i]);
      free(elem_blknames);

      for (i = 0; i < nelemvars; i++)
        free(elvarnames[i]);
      free(elvarnames);
    } /* if (nelemvars) */


    /* Read in variables associated with nodes */

    /* How many variables are there on nodes */

    status = ex_get_var_param(exoid, "n", &nnodevars);
    if (status < 0) {
      sprintf(mesg, "Error while reading node variables in Exodus II file %s\n",
              filename);
      MSTK_Report(funcname,mesg,MSTK_FATAL);
    }

    if (nnodevars) {
      /* What are their names */

      char **nodevarnames = (char **) malloc(nnodevars*sizeof(char *));
      for (i = 0; i < nnodevars; i++)
        nodevarnames[i] = (char *) malloc(256*sizeof(char));

      status = ex_get_var_names(exoid, "n", nnodevars, nodevarnames);
      if (status < 0) {
        sprintf(mesg, "Error while reading node variable names in Exodus II file %s\n",filename);
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

            MAttrib_ptr mattrib;        
            if (surf_elems) {
              MAttType atttype = (ncomp == 2) ? VECTOR : TENSOR;
              mattrib = MAttrib_New(mesh,attname,atttype,MVERTEX,ncomp);
            }
            else if (solid_elems) {
              MAttType atttype = (ncomp == 3) ? VECTOR : TENSOR;
              mattrib = MAttrib_New(mesh,attname,atttype,MVERTEX,ncomp);
            }
            else {
              MSTK_Report("MESH_ReadExodusII_Serial",
                          "Attribute export on wire meshes and node meshes not implemented",
                          MSTK_ERROR);
              continue;
            }

            int nodeid = 0;

            double **node_var_vals = (double **) malloc(ncomp*sizeof(double *));
            for (j = 0; j < ncomp; j++)
              node_var_vals[j] = malloc(nnodes*sizeof(double));

            for (j = 0; j < ncomp; j++) {
              varindex2 = varindex + j;
              status = ex_get_nodal_var(exoid, time_step, varindex2+1,
                                        nnodes, node_var_vals[j]);
              if (status < 0) {
                sprintf(mesg, 
                        "Error while reading node variables in Exodus II file %s\n",
                        filename);
                MSTK_Report(funcname,mesg,MSTK_FATAL);
              }
            }

            for (k = 0; k < nnodes; k++) {
              double *pval = malloc(ncomp*sizeof(double)); // freed by MESH_DELETE
              for (j = 0; j < ncomp; j++)
                pval[j] = node_var_vals[j][k];
              MEntity_ptr ment = MESH_Vertex(mesh,k);
              MEnt_Set_AttVal(ment,mattrib,0,0.0,pval);
            }

            for (j = 0; j < ncomp; j++)
              free(node_var_vals[j]);
            free(node_var_vals);

          } /* if its component 0 of a vector */

        } /* if its a vector */
        else { /* scalar variable */

          MAttrib_ptr mattrib = MAttrib_New(mesh,varname,DOUBLE,MVERTEX);

          double *node_var_vals = (double *) malloc(nnodes*sizeof(double));

          status = ex_get_nodal_var(exoid, time_step, varindex+1, nnodes, 
                                    node_var_vals);
          if (status < 0) {
            sprintf(mesg, 
                    "Error while reading node variables in Exodus II file %s\n",
                    filename);
            MSTK_Report(funcname,mesg,MSTK_FATAL);
          }

          for (k = 0; k < nnodes; k++) {
            MEntity_ptr ment = MESH_Vertex(mesh,k);
            MEnt_Set_AttVal(ment,mattrib,0,node_var_vals[k],NULL);
          }

          free(node_var_vals);

        } /* If its a scalar variable */
      } /* for each element variable */
  
      for (i = 0; i < nnodevars; i++)
        free(nodevarnames[i]);
      free(nodevarnames);
    } /* if (nnodevars) */

    ex_close(exoid);

    return 1;

  }


#ifdef __cplusplus
}
#endif
