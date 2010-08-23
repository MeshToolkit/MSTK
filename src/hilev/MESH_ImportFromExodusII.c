#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MSTK.h"

#include "exodusII.h"
#include "exodusII_ext.h"


#ifdef __cplusplus
extern "C" {
#endif


int MESH_ImportFromExodusII(Mesh_ptr mesh, const char *filename) {

  char mesg[256], funcname[256]="MESH_ImportFromExodusII";
  char title[256], elem_type[256];
  char (*elem_blknames)[256];
  int i, j, k, k1;
  int comp_ws = sizeof(double), io_ws = 0;
  int exoid, status;
  int ndim, nnodes, nelems, nelblock, nnodesets, nsidesets;
  int nedges, nedge_blk, nfaces, nface_blk, nelemsets;
  int nedgesets, nfacesets, nnodemaps, nedgemaps, nfacemaps, nelemmaps;
  int *elem_blk_ids, *connect, *node_map, *elem_map;
  int nelem_i, nelnodes, neledges, nelfaces, natts;
  double *xvals, *yvals, *zvals, xyz[3];
  float version;
  
  MVertex_ptr mv, *fverts, *rverts;
  MFace_ptr mf;
  MRegion_ptr mr;
  MAttrib_ptr nmapatt, nelmapatt, elblockatt;
  
  
  ex_init_params exopar;
  
  
  exoid = ex_open(filename, EX_READ, &comp_ws, &io_ws, &version);
  
  if (exoid < 0) {
    sprintf(mesg,"Cannot open/read Exodus II file %s\n",filename);
    MSTK_Report(funcname,mesg,FATAL);
  }
  
  
  status = ex_get_init_ext(exoid, &exopar);
  if (status < 0) {
    sprintf(mesg,"Error while reading Exodus II file %s\n",filename);
    MSTK_Report(funcname,mesg,FATAL);
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
  nnodemaps = exopar.num_node_maps;
  nedgemaps = exopar.num_edge_maps;
  nfacemaps = exopar.num_face_maps;
  nelemmaps = exopar.num_elem_maps;
  
  
  
  
  
  /* read the vertex information */
  
  xvals = (double *) malloc(nnodes*sizeof(double));
  yvals = (double *) malloc(nnodes*sizeof(double));
  zvals = (double *) malloc(nnodes*sizeof(double));
  
  status = ex_get_coord(exoid, xvals, yvals, zvals);
  if (status < 0) {
    sprintf(mesg,"Error while reading vertex info in Exodus II file %s\n",filename);
    MSTK_Report(funcname,mesg,FATAL);
  }

  for (i = 0; i < nnodes; i++) {
    
    mv = MV_New(mesh);
    
    xyz[0] = xvals[i]; xyz[1] = yvals[i]; xyz[2] = zvals[i];
    MV_Set_Coords(mv,xyz);
  }
  
  

  /* read node number map - store it as an attribute to spit out later
     if necessary */
  
  node_map = (int *) malloc(nnodes*sizeof(int));
  
  status = ex_get_node_num_map(exoid, node_map);
  if (status < 0) {
    sprintf(mesg,"Error while reading node map in Exodus II file %s\n",filename);
    MSTK_Report(funcname,mesg,FATAL);
  }
  
  if (status == 0) {
    
    nmapatt = MAttrib_New(mesh, "node_map", INT, MVERTEX);
    
    for (i = 0; i < nnodes; i++) {
      mv = MESH_Vertex(mesh, i);
      MEnt_Set_AttVal(mv, nmapatt, node_map[i], 0.0, NULL);
    }
    
  }
  
  
  
  /* Get element block IDs and names */
  
  
  elem_blk_ids = (int *) malloc(nelblock*sizeof(int));
  status = ex_get_elem_blk_ids(exoid, elem_blk_ids);
  if (status < 0) {
    sprintf(mesg,"Error while reading element block ids in Exodus II file %s\n",filename);
    MSTK_Report(funcname,mesg,FATAL);
  }
  
  elem_blknames = (char (*)[256]) malloc(nelblock*sizeof(char [256]));
  status = ex_get_names(exoid, EX_ELEM_BLOCK, elem_blknames);
  if (status < 0) {
    sprintf(mesg,"Error while reading element block ids in Exodus II file %s\n",filename);
    MSTK_Report(funcname,mesg,FATAL);
  }
  
  
  
  
  
  if (ndim == 1) {
    
    MSTK_Report("MESH_ImportFromExodusII","Not implemented",FATAL);
    
  }
  else if (ndim == 2) {
    
    elblockatt = MAttrib_New(mesh, "elblock_id", INT, MFACE);
    
    
    /* Read each element block */
    
    for (i = 0; i < nelblock; i++) {
      
      status = ex_get_block(exoid, EX_ELEM_BLOCK, elem_blk_ids[i], 
			    elem_type, &nelem_i, &nelnodes, 
			    &neledges, &nelfaces, &natts);
      if (status < 0) {
	sprintf(mesg,"Error while reading element block %s in Exodus II file %s\n",elem_blknames[i],filename);
	MSTK_Report(funcname,mesg,FATAL);
      }
      
      
      if (strcmp(elem_type,"NSIDED") == 0 || 
	  strcmp(elem_type,"nsided") == 0) {
	
	MSTK_Report("MESH_ImportFromExodusII","Reading of polyhedral elements not yet implemented",WARN);
	
      }
      else {
	
	/* Get the connectivity of all elements in this block */
	
	connect = (int *) calloc(nelnodes*nelem_i,sizeof(int));
	
	status = ex_get_elem_conn(exoid, elem_blk_ids[i], connect);
	if (status < 0) {
	  sprintf(mesg,"Error while reading element block %s in Exodus II file %s\n",elem_blknames[i],filename);
	  MSTK_Report(funcname,mesg,FATAL);
	}
	
	
	/* Create the MSTK faces */
	
	fverts = (MVertex_ptr *) calloc(nelnodes,sizeof(MVertex_ptr));
	
	for (j = 0; j < nelem_i; j++) {
	  
	  mf = MF_New(mesh);
	  
	  for (k = 0; k < nelnodes; k++)
	    fverts[k] = MESH_VertexFromID(mesh,connect[nelnodes*j+k]);
	  
	  MF_Set_Vertices(mf,nelnodes,fverts);
	}
	
	free(fverts);
	
      } /* if (strcmp(elem_type,"NSIDED") ... else */
      
    } /* for (i = 0; i < nelblock; i++) */
    
    
    
    /* read element number map - store it as an attribute to spit out
       later if necessary */
    
    elem_map = (int *) malloc(nelems*sizeof(int));
    
    status = ex_get_node_num_map(exoid, elem_map);
    if (status < 0) {
      sprintf(mesg,"Error while reading element map in Exodus II file %s\n",filename);
      MSTK_Report(funcname,mesg,FATAL);
    }
    
    if (status == 0) {
      
      nmapatt = MAttrib_New(mesh, "elem_map", INT, MFACE);
      
      for (i = 0; i < nelems; i++) {
	mf = MESH_Face(mesh, i);
	MEnt_Set_AttVal(mf, nmapatt, elem_map[i], 0.0, NULL);
      }
      
    }
    
  }
  else if (ndim == 3) {
    
    elblockatt = MAttrib_New(mesh, "elblock_id", INT, MREGION);
    
    
    for (i = 0; i < nelblock; i++) {
      
      status = ex_get_block(exoid, EX_ELEM_BLOCK, elem_blk_ids[i], 
			    elem_type, &nelem_i, &nelnodes, 
			    &neledges, &nelfaces, &natts);
      if (status < 0) {
	sprintf(mesg,
		"Error while reading element block %s in Exodus II file %s\n",
		elem_blknames[i],filename);
	MSTK_Report(funcname,mesg,FATAL);
      }
      
      
      if (strcmp(elem_type,"NFACED") == 0 || 
	  strcmp(elem_type,"nfaced") == 0) {
	
	MSTK_Report("MESH_ImportFromExodusII","Reading of polyhedral elements not yet implemented",WARN);
	
      }
      else {
	
	/* Get the connectivity of all elements in this block */
	
	connect = (int *) calloc(nelnodes*nelem_i,sizeof(int));
	
	status = ex_get_elem_conn(exoid, elem_blk_ids[i], connect);
	if (status < 0) {
	  sprintf(mesg,"Error while reading element block %s in Exodus II file %s\n",elem_blknames[i],filename);
	  MSTK_Report(funcname,mesg,FATAL);
	}
	
	
	/* Create the MSTK regions */
	
	rverts = (MVertex_ptr *) calloc(nelnodes,sizeof(MVertex_ptr));
	
	if (strcmp(elem_type,"TETRA") == 0) {
	  if (nelnodes > 4) {
	    MSTK_Report(funcname,"Higher order tets not supported",WARN);
	    continue;
	  }
	}
	else if (strcmp(elem_type,"WEDGE") == 0) {
	  if (nelnodes > 6) {
	    MSTK_Report(funcname,"Higher order prisms/wedges not supported",WARN);
	    continue;
	  }
	}
	else if (strcmp(elem_type,"HEX") == 0) {
	  if (nelnodes > 8) {
	    MSTK_Report(funcname,"Higher order hexes not supported",WARN);
	    continue;
	  }
	}
	else {
	  sprintf(mesg,"Unrecognized or unsupported solid element type: %s",elem_type);
	  MSTK_Report(funcname,mesg,WARN);
	  continue;
	}
	
	
	for (j = 0; j < nelem_i; j++) {
	  
	  mr = MR_New(mesh);
	  
	  /* Exodus II and MSTK node conventions are the same */
	  
	  for (k = 0; k < nelnodes; k++)
	    rverts[k] = MESH_VertexFromID(mesh,connect[nelnodes*j+k]);
	  
	  MR_Set_Vertices(mr,nelnodes,rverts,0,NULL);
	}
	
	free(rverts);
	
      } /* if (strcmp(elem_type,"NSIDED") ... else */
      
    } /* for (i = 0; i < nelblock; i++) */
    
    
    

    
    
    /* read element number map - store it as an attribute to spit out
       later if necessary */
    
    elem_map = (int *) malloc(nelems*sizeof(int));
    
    status = ex_get_node_num_map(exoid, elem_map);
    if (status < 0) {
      sprintf(mesg,"Error while reading element map in Exodus II file %s\n",filename);
      MSTK_Report(funcname,mesg,FATAL);
    }
    
    if (status == 0) {
      
      nmapatt = MAttrib_New(mesh, "elem_map", INT, MREGION);
      
      for (i = 0; i < nelems; i++) {
	mr = MESH_Region(mesh, i);
	MEnt_Set_AttVal(mr, nmapatt, elem_map[i], 0.0, NULL);
      }
      
    }
  }

  ex_close(exoid);

  return 1;

}


#ifdef __cplusplus
}
#endif
