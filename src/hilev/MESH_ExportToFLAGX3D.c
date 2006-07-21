#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif



  /* Function to export MSTK mesh to FLAG X3D format (LA-UR-04-9033) */

  /* if natt = 0, all attributes are written out if natt = -1, no
     attributes are written out if natt > 0, only attributes specified
     in attnames are written out Except for cell material IDs,
     attributes are ignored by FLAG; still we will write them out for
     completeness */

  /* opts is an array of flags that controls how mesh is exported to FLAG 
     Currently, it is a dummy argument
  */

  /* WE ARE NOT HANDLING DISTRIBUTED MESHES */



int MESH_ExportToFLAGX3D(Mesh_ptr mesh, const char *filename, const int natt, 
			 const char **attnames, int *opts) {
  int			gentid, *gentities;
  List_ptr	        fverts, rfaces, fregs, efaces, fedges;
  MVertex_ptr           vertex;
  MEdge_ptr             edge;
  MFace_ptr	        face;
  MRegion_ptr           region;
  MAttrib_ptr           attrib, *nodatts=NULL, *cellatts=NULL, oppatt;
  MAttType              atttype;
  char                  attname[256], matname[256], tmpstr[256];
  int                   jv, je, jr, jf;
  int                   nv, ne, nf, nr, nrf, nfv, nef, nfe, nfr;
  int			i, found, k, nmeshatt, ival;
  int                   nalloc, ngent;
  int                   attentdim, j, idx;
  int                   ndup, max_nrf, max_nfe;
  int                   oppfid, oppeid, nf2, ne2;
  int                   nnodatt, ncellatt;
  double		vxyz[3], rval;
  void                 *pval;
  FILE		        *fp;


  if (!(fp = fopen(filename,"w"))) {
    fprintf(stderr,"MESH_ExportToFLAGX3D: Couldn't open output file %s\n",
	    filename);
    exit(2);
  }
  
  
  nv = MESH_Num_Vertices(mesh);
  ne = MESH_Num_Edges(mesh);
  nf = MESH_Num_Faces(mesh);
  nr = MESH_Num_Regions(mesh);
  if (!nr && !nf) {
    fprintf(stderr,
	    "Format does not support meshes with only edges or nodes\n");
    return 0;
  }
  
  /***********************************************************************/
  /* Opening line                                                        */
  /***********************************************************************/
  
  
  fprintf(fp,"x3dtoflag ascii\n");
  
  
  /***********************************************************************/
  /* Header information */
  /***********************************************************************/

  fprintf(fp,"header\n");

  /* Processor number - For now we will not bother with distributed I/O */

  fprintf(fp,"   %-22s %10d\n","process",1);

  /* Problem dimension - 2 if we only have faces, 3 if we have regions */
  
  if (nr)
    fprintf(fp,"   %-22s %10d\n","numdim",3);
  else 
    fprintf(fp,"   %-22s %10d\n","numdim",2);
  
  
  /* Number of materials */
  
  ngent = 0;
  nalloc = 10;
  gentities = (int *) MSTK_malloc(nalloc*sizeof(int));
  
  if (nr) {
    idx = 0;
    while ((region = MESH_Next_Region(mesh,&idx))) {
      gentid = MR_GEntID(region);
      
      found = 0;
      for (i = 0; i < ngent; i++)
	if (gentities[i] == gentid) {
	  found = 1;
	  break;
	}
      
      if (!found) {
	if (ngent+1 >= nalloc) {
	  nalloc *= 2;
	  gentities = (int *) MSTK_realloc(gentities,nalloc*sizeof(int));
	}
	gentities[ngent] = gentid;
	ngent++;
      }
    }
  }
  else {
    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {
      gentid = MF_GEntID(face);
      
      found = 0;
      for (i = 0; i < ngent; i++)
	if (gentities[i] == gentid) {
	  found = 1;
	  break;
	}
      
      if (!found) {
	if (ngent+1 >= nalloc) {
	  nalloc *= 2;
	  gentities = (int *) MSTK_realloc(gentities,nalloc*sizeof(int));
	}
	gentities[ngent] = gentid;
	ngent++;
      }
    }
  }
  
  if (!ngent) {
    /* No materMake all cells of material type 1 */
    gentities[ngent] = 1;
    ngent++;
  }
  
  fprintf(fp,"   %-22s %10d\n","materials",ngent);
  
  /* Number of nodes */
  
  fprintf(fp,"   %-22s %10d\n","nodes",nv);
  
  /* Number of "faces"; "faces" is the unfortunate term used in the
     FLAG X3D file to mean boundaries of elements/cells/zones in FLAG
     X3D, i.e., in MSTK parlance, 'faces' are edges in a surface mesh,
     faces in a volume mesh */
  
  if (nr) {
    
    oppatt = MAttrib_New(mesh,"oppfaceID",INT,MFACE); 
    
    ndup = 0;
    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {      
      
      /* Each face in FLAG X3D must be referenced by only one
	 cell. A face description in FLAG X3D is supposed to be such
	 that normal of 'face' points out of its referring
	 cell. Therefore, if a face is connected to two regions,
	 then we need to write out the face and its duplicate,
	 pointing in the opposite direction. If a face is connected
	 to only one region, we need to write out the face in its
	 natural direction if its normal points out of the owning
	 region; we need to write out the face in the reverse
	 direction if its normal points into the owning region */
      
      fregs = MF_Regions(face);
      if (fregs && List_Num_Entries(fregs) == 2) {
	
	ndup++;
	MEnt_Set_AttVal(face,oppatt,nf+ndup,0,NULL);
      }
      if (fregs)
	List_Delete(fregs);
    }
    
    nf2 = nf + ndup;
    
    fprintf(fp,"   %-22s %10d\n","faces",nf2);
    fprintf(fp,"   %-22s %10d\n","elements",nr);
  }
  else {
    oppatt = MAttrib_New(mesh,"oppedgeID",INT,MEDGE); 
    
    ndup = 0;
    idx = 0;
    while ((edge = MESH_Next_Edge(mesh,&idx))) {      
      
      /* Each edge in FLAG X3D must be referenced by only one cell. An
	 edge description in FLAG X3D is supposed to be such that the
	 cell uses the edge in the +ve direction. Therefore, if an
	 edge is connected to two cells, then we need to write out the
	 edge and its duplicate, pointing in the opposite
	 direction. If an edge is connected to only one cell, we need
	 to write out the edge in its natural direction if the cell
	 uses the edge in the +ve sense; we need to write out the edge
	 in the reverse direction if the cell uses it in the -ve
	 sense */
      
      efaces = ME_Faces(edge);
      if (efaces && List_Num_Entries(efaces) == 2) {	  
	ndup++;
	MEnt_Set_AttVal(edge,oppatt,ne+ndup,0,NULL);
      }
      
#ifdef DEBUG
      if (efaces && List_Num_Entries(efaces) > 2) {
	fprintf(stderr,"Non-manifold model - more than two faces connected to edge in surface mesh.\n");
      }
#endif
      
      if (efaces)
	List_Delete(efaces);
    }
    
    ne2 = ne + ndup;
    
    fprintf(fp,"   %-22s %10d\n","faces",ne2);
    fprintf(fp,"   %-22s %10d\n","elements",nf);
  }
  
  /* We are not yet supporting distributed meshes - so number of ghost
     nodes is 0 */
  
  fprintf(fp,"   %-22s %10d\n","ghost_nodes",0);
  
  /* No slave (constrained) node info - there are no mechanisms in
     place to do this in MSTK and therefore, this info cannot be
     written out */
  
  fprintf(fp,"   %-22s %10d\n","slave_nodes",0);
  fprintf(fp,"   %-22s %10d\n","nodes_per_slave",2);  /* required default */
  
  /* Nodes per "face", required for array dimensioning in FLAG. If
     this is a surface mesh, this will be 2; if it is a volume mesh,
     it is the maximum number of nodes that any face in the mesh
     has */
  
  if (nr) {
    idx = 0; max_nfe = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {
      nfe = MF_Num_Edges(face);
      if (nfe > max_nfe)
	max_nfe = nfe;
    }
    
    fprintf(fp,"   %-22s %10d\n","nodes_per_face",max_nfe); /* real faces */
  }
  else {
    fprintf(fp,"   %-22s %10d\n","nodes_per_face",2); /* edges */
  }
  
  /* "faces" per cell, required for array dimensioning in FLAG. If
     this is a surface mesh, this the maximum number of edges any face
     has; if it is a volume mesh, this is the maximum number of faces
     any region has */
  
  if (nr) {
    idx = 0; max_nrf = 0;
    while ((region = MESH_Next_Region(mesh,&idx))) {
      nrf = MR_Num_Faces(region);
      if (nrf > max_nrf)
	max_nrf = nrf;
    }
    
    fprintf(fp,"   %-22s %10d\n","faces_per_cell",max_nrf);
  }
  else {
    idx = 0; max_nfe = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {
      nfe = MF_Num_Edges(face);
      if (nfe > max_nfe)
	max_nfe = nfe;
    }
    
    fprintf(fp,"   %-22s %10d\n","faces_per_cell",max_nfe);
  }
  
  /* Number of vertex and element centered attributes */
  
  nmeshatt = MESH_Num_Attribs(mesh);
  
  /* if (natt < 0), no attributes are to be written out */
  
  nnodatt = ncellatt = 0;
  if (natt >= 0 && nmeshatt) {
    
    /* Collect the attributes */
    
    cellatts = (MAttrib_ptr *) MSTK_malloc(nmeshatt*sizeof(int));
    nodatts = (MAttrib_ptr *) MSTK_malloc(nmeshatt*sizeof(int));
    
    for (i = 0; i < nmeshatt; i++) {
      attrib = MESH_Attrib(mesh,i);

      /* No need to write out the temporary array we created in this routine */
      if (attrib == oppatt)
        continue;
      
      /* If the attribute is not a INT or a DOUBLE we cannot write it out */
      atttype = MAttrib_Get_Type(attrib);      
      if (atttype != INT && atttype != DOUBLE)
	continue;

      /* If natt == 0, all qualifying attributes are to be written out */
      /* If natt > 0, only specified ones are to be exported */
      
      if (natt > 0) {
	MAttrib_Get_Name(attrib,attname);
	
	found = 0; j = 0;
	while (!found && j < natt) {
	  if (strcmp(attname,attnames[j]) == 0)
	    found = 1;
	  else
	    j++;
	}
	
	if (!found)
	  continue;
      }
      
      
      attentdim = MAttrib_Get_EntDim(attrib);
      if ((attentdim == MVERTEX) || (attentdim == MALLTYPE))
	nodatts[nnodatt++] = attrib;
      
      
      if (nr) {
	if ((attentdim == MREGION) || (attentdim == MALLTYPE)) {
	  cellatts[ncellatt++] = attrib;
	}
      }
      else if (nf) {
	if ((attentdim == MFACE) || (attentdim == MALLTYPE)) {
	  cellatts[ncellatt++] = attrib;
	}
      }
    }
  }
  
  fprintf(fp,"   %-22s %10d\n","node_data_fields",nnodatt);
  
  /* must write out matid and partelm and then add additional cell data */
  fprintf(fp,"   %-22s %10d\n","cell_data_fields",ncellatt+2); 
  
  
  /* End of header information */
  fprintf(fp,"end_header\n");
  
  
  /***********************************************************************/
  /* Material info */
  /***********************************************************************/
  
  /* Material names */
  fprintf(fp,"matnames\n");
  for (i = 0; i < ngent; i++) {
    sprintf(matname,"mat%-d",i+1);
    fprintf(fp,"   %10d   %s\n",i+1,matname);
  }
  fprintf(fp,"end_matnames\n");
  
  /* Material EOS - will not specify */
  fprintf(fp,"mateos\n");
  for (i = 0; i < ngent; i++)
    fprintf(fp,"   %10d   -1\n",i+1);
  fprintf(fp,"end_mateos\n");
  
  /* Material opacity - will not specify */
  fprintf(fp,"matopc\n");
  for (i = 0; i < ngent; i++)
    fprintf(fp,"   %10d   -1\n",i+1);
  fprintf(fp,"end_matopc\n");
  
  
  
  /***********************************************************************/
  /* Node coordinates                                                    */
  /* FLAG X3D manual says the format must be (I10, 3(1X, 1PE22.14E3))    */
  /* I believe PE means that instead of 12.3 is printed as 1.23E+01      */
  /* instead of 0.123E+02. I don't know how to force that in C. I also   */
  /* don't know how to force C to print E+xx or E-yy for the exponent    */
  /***********************************************************************/
  
  fprintf(fp,"nodes\n");
  for (jv = 0; jv < nv; jv++) {
    vertex = MESH_Vertex(mesh,jv);
    MV_Coords(vertex,vxyz);
    fprintf(fp,"% 10d % 22.14E % 22.14E % 22.14E\n",MV_ID(vertex),
	    vxyz[0],vxyz[1],vxyz[2]);
    MV_Set_ID(vertex,(jv+1));
  }
  fprintf(fp,"end_nodes\n");

  /***********************************************************************/
  /* Faces Data Block                                                    */
  /* Faces is the unfortunate term that the FLAG X3D file uses to refer  */
  /* to the boundaries of elements/cells, i.e., "faces" are mesh faces   */
  /* in a volume mesh, and are mesh edges in a surface mesh              */
  /* Contrary to what the FLAG X3D documentation says, the format that   */
  /* works seems to be 13I10                                             */
  /***********************************************************************/

  fprintf(fp,"faces\n");
  if (nr) {

    /***********************************************************/
    /* Solid mesh                                              */
    /***********************************************************/

    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {

      fregs = MF_Regions(face);
      nfr = fregs ? List_Num_Entries(fregs) : 0;

      fprintf(fp,"% 10d",MF_ID(face));
      k = 1;

      if (nfr == 2) {
	/* Two regions connected to face - write face out in its
	   natural orientation here, write out the duplicate face in
	   the opposite orientation later */

	fverts = MF_Vertices(face,1,0);
	nfv = List_Num_Entries(fverts);
	
	fprintf(fp,"% 10d",nfv);
	k++;
	
	for (jv = 0; jv < nfv; jv++) {
	  vertex = List_Entry(fverts,jv);
	  fprintf(fp,"% 10d",MV_ID(vertex));
	  if ((++k)%13 == 0)
	    fprintf(fp,"\n");
	}

	List_Delete(fverts);
	List_Delete(fregs);


	/* ID of the processor owning this face */
	/* We are not yet handling distributed meshes */

	fprintf(fp,"% 10d",1);
	if ((++k)%13 == 0)
	  fprintf(fp,"\n");
	
	MEnt_Get_AttVal(face,oppatt,&oppfid,&rval,&pval);
#ifdef DEBUG
	if (!oppfid)
	  fprintf(stderr,"Internal face has no duplicate face info?\n");
#endif

	/* ID of the processor owning the duplicate face */
	/* We are not yet handling distributed meshes */
      
	fprintf(fp,"% 10d",1);
	if ((++k)%13 == 0)
	  fprintf(fp,"\n");
	
	/* ID of the duplicate face */
	
	fprintf(fp,"% 10d",oppfid);
	if ((++k)%13 == 0)
	  fprintf(fp,"\n");

	/* Five dummy arguments required by FLAG X3D format specification */
	
	for (i = 0; i < 5; i++) {
	  fprintf(fp,"% 10d",1);
	  if (i <= 3 && (++k)%13 == 0)
	    fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
      }
      else if (nfr == 1) {
	/* Only one region connected to face - If the face normal
	   points out of region, then write the face out in its
	   natural orientation; if not, write it out in the opposite
	   orientation */

	if (MR_FaceDir(List_Entry(fregs,0),face) == 1)
	  fverts = MF_Vertices(face,1,0);
	else
	  fverts = MF_Vertices(face,0,0);

	nfv = List_Num_Entries(fverts);
	
	fprintf(fp,"% 10d",nfv);
	k++;
	
	for (jv = 0; jv < nfv; jv++) {
	  vertex = List_Entry(fverts,jv);
	  fprintf(fp,"% 10d",MV_ID(vertex));
	  if ((++k)%13 == 0)
	    fprintf(fp,"\n");
	}	

	List_Delete(fverts);
	List_Delete(fregs);

	/* ID of the processor owning this face */
	/* We are not yet handling distributed meshes */
	
	fprintf(fp,"% 10d",1);
	if ((++k)%13 == 0)
	  fprintf(fp,"\n");
	
	/* Boundary face - no duplicate face info */
	
	fprintf(fp,"% 10d",0); /* processor owning duplicate face */
	if ((++k)%13 == 0)
	  fprintf(fp,"\n");
	fprintf(fp,"% 10d",0); /* duplicate face ID */
	if ((++k)%13 == 0)
	  fprintf(fp,"\n");

	/* Five dummy arguments required by FLAG X3D format specification */
	
	for (i = 0; i < 5; i++) {
	  fprintf(fp,"% 10d",1);
	  if (i <= 3 && (++k)%13 == 0)
	    fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
      }
#ifdef DEBUG
      else {
	if (nfr == 0)
	  fprintf(stderr,
		  "Meshes of non-manifold objects not supported in FLAG X3D\n");
	else 
	  fprintf(stderr,
		  "Invalid mesh - too many regions connected to face\n");
      }
#endif
    }


    jf = 0;
    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {
      MEnt_Get_AttVal(face,oppatt,&oppfid,&rval,&pval);
      if (!oppfid) continue;

      /* Write out the duplicate of this internal face */

      fprintf(fp,"% 10d",oppfid);
      k = 1;

      /* Write out the vertices in the opposite direction to the
	 definition of the face */

      fverts = MF_Vertices(face,0,0);
      nfv = List_Num_Entries(fverts);

      fprintf(fp,"% 10d",nfv);
      k++;

      for (jv = 0; jv < nfv; jv++) {
	vertex = List_Entry(fverts,jv);
	fprintf(fp,"% 10d",MV_ID(vertex));
	if ((++k)%13 == 0)
	  fprintf(fp,"\n");
      }
      
      /* ID of the processor owning this face */
      /* We are not yet handling distributed meshes */
      
      fprintf(fp,"% 10d",1);
      if ((++k)%13 == 0)
	fprintf(fp,"\n");
      
      /* ID of the processor owning the duplicate face */
      /* We are not yet handling distributed meshes */
      
      fprintf(fp,"% 10d",1);
      if ((++k)%13 == 0)
	fprintf(fp,"\n");
      
      /* ID of the duplicate (in this case, original) face */
      
      fprintf(fp,"% 10d",MF_ID(face));
      if ((++k)%13 == 0)
	fprintf(fp,"\n");

      /* Five dummy arguments required by FLAG X3D format specification */
      
      for (i = 0; i < 5; i++) {
	fprintf(fp,"% 10d",1);
	if (i <= 3 && (++k)%13 == 0)
	  fprintf(fp,"\n");
      }
      fprintf(fp,"\n");
    }
  }
  else {

    /***********************************************************/
    /* Surface mesh (maybe planar)                             */
    /* The info line for each "face" in a 2D mesh (i.e. edge)  */
    /* will only be 12 items long; so we needn't worry about   */
    /* special care to satisfy the 13I10 format                */
    /***********************************************************/

    idx = 0;
    while ((edge = MESH_Next_Edge(mesh,&idx))) {

      efaces = ME_Faces(edge);
      nef = efaces ? List_Num_Entries(efaces) : 0;

      fprintf(fp,"% 10d",ME_ID(edge));

      if (nef == 2) {
	/* Two faces connected to edge - write edge out in its
	   natural orientation here, write out the duplicate edge in
	   the opposite orientation later */

	fprintf(fp,"% 10d",2);
	
	for (jv = 0; jv < 2; jv++) {
	  vertex = ME_Vertex(edge,jv);
	  fprintf(fp,"% 10d",MV_ID(vertex));
	}

	List_Delete(efaces);


	/* ID of the processor owning this edge */
	/* We are not yet handling distributed meshes */
	
	fprintf(fp,"% 10d",1);
	
	MEnt_Get_AttVal(edge,oppatt,&oppeid,&rval,&pval);
#ifdef DEBUG
	if (!oppeid)
	  fprintf(stderr,"Internal edge has no duplicate edge info?\n");
#endif

	/* ID of the processor owning the duplicate edge */
	/* We are not yet handling distributed meshes */
      
	fprintf(fp,"% 10d",1);
	
	/* ID of the duplicate edge */
	
	fprintf(fp,"% 10d",oppeid);

	/* Five dummy arguments required by FLAG X3D format specification */
	
	for (i = 0; i < 5; i++) fprintf(fp,"% 10d",1);
	fprintf(fp,"\n");
      }
      else if (nef == 1) {
	/* Only one face connected to edge - If the edge is used in
	   the +ve sense by the face, then write the edge out in its
	   natural orientation; if not, write it out in the opposite
	   orientation */

	fprintf(fp,"% 10d",2);

	if (MF_EdgeDir(List_Entry(efaces,0),edge) == 1) {
	  for (jv = 0; jv < 2; jv++) {
	    vertex = ME_Vertex(edge,jv);
	    fprintf(fp,"% 10d",MV_ID(vertex));
	  }
	}
	else {
	  for (jv = 0; jv < 2; jv++) {
	    vertex = ME_Vertex(edge,!jv);
	    fprintf(fp,"% 10d",MV_ID(vertex));
	  }
	}

	List_Delete(efaces);

	/* ID of the processor owning this edge */
	/* We are not yet handling distributed meshes */
	
	fprintf(fp,"% 10d",1);
	
	/* Boundary face - no duplicate edge info */
	
	fprintf(fp,"% 10d",0); /* id of processor owning duplicate edge */
	fprintf(fp,"% 10d",0); /* duplicate edge id */

	/* Five dummy arguments required by FLAG X3D format specification */
	
	for (i = 0; i < 5; i++) fprintf(fp,"% 10d",1);
	fprintf(fp,"\n");
      }
#ifdef DEBUG
      else {
	if (nfr == 0)
	  fprintf(stderr,
		  "Meshes of non-manifold objects not supported in FLAG X3D\n");
	else 
	  fprintf(stderr,
		  "Invalid mesh - too many regions connected to face\n");
      }
#endif
    }


    idx = 0;
    while ((edge = MESH_Next_Edge(mesh,&idx))) {
      MEnt_Get_AttVal(edge,oppatt,&oppeid,&rval,&pval);
      if (!oppeid) continue;

      fprintf(fp,"% 10d",oppeid);
      
      fprintf(fp,"% 10d",2);

      for (jv = 0; jv < 2; jv++) {
	vertex = ME_Vertex(edge,!jv);
	fprintf(fp,"% 10d",MV_ID(vertex));
      }
      
      /* ID of the processor owning this edge */
      /* We are not yet handling distributed meshes */
      
      fprintf(fp,"% 10d",1);
      
      /* ID of the processor owning the duplicate edge */
      /* We are not yet handling distributed meshes */
      
      fprintf(fp,"% 10d",1);
      
      /* ID of the duplicate (in this case, original) face */
      
      fprintf(fp,"% 10d",ME_ID(edge));

      /* Five dummy arguments required by FLAG X3D format specification */
      
      for (i = 0; i < 5; i++) fprintf(fp,"% 10d",1);
      fprintf(fp,"\n");
    }

  }
  fprintf(fp,"end_faces\n");
  

  /***********************************************************************/
  /* Cell Data Block                                                     */
  /***********************************************************************/

  fprintf(fp,"cells\n");
  if (nr) {

    idx = 0;    
    while ((region = MESH_Next_Region(mesh,&idx))) {
      fprintf(fp,"% 10d",MR_ID(region));
      k = 1;

      rfaces = MR_Faces(region);
      nrf = List_Num_Entries(rfaces);
      
      fprintf(fp,"% 10d",nrf);
      k++;

      for (jf = 0; jf < nrf; jf++) {
	face = List_Entry(rfaces,jf);

	MEnt_Get_AttVal(face,oppatt,&oppfid,&rval,&pval);
	if (oppfid) {
	  /* Interior face */
	  /* Write out the face id, if the region uses the face in the
	     +ve sense (the face normal points out of the region); write
	     out the duplicate face id, if not */
	  
	  if (MR_FaceDir_i(region,jf))
	    fprintf(fp,"% 10d",MF_ID(face));
	  else
	    fprintf(fp,"% 10d",oppfid);
	}
	else {
	  /* Boundary face */
	  /* Write out the face id - any necessary reordered writing
	     of nodes so that the face points out of the region was
	     done earlier */

	  fprintf(fp,"% 10d",MF_ID(face));
	}

	k++;
	if (k%14 == 0 && jf != nrf-1)
	  fprintf(fp,"\n");
	  
      }

      List_Delete(rfaces);
      fprintf(fp,"\n");
    }

    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx)))
      MEnt_Rem_AttVal(face,oppatt);
    MAttrib_Delete(oppatt);
  }
  else {

    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {
      fprintf(fp,"% 10d",MF_ID(face));
      k = 1;

      fedges = MF_Edges(face,1,0);
      nfe = List_Num_Entries(fedges);

      fprintf(fp,"% 10d",nfe);
      k++;

      for (je = 0; je < nfe; je++) {
	edge = List_Entry(fedges,je);

	MEnt_Get_AttVal(edge,oppatt,&oppeid,&rval,&pval);
	if (oppeid) {
	  /* Interior edge */
	  /* Write out the edge id, if the face uses the edge in the
	     +ve sense; write out the duplicate edge id, if not */

	  if (MF_EdgeDir_i(face,je))
	    fprintf(fp,"% 10d",ME_ID(edge));
	  else
	    fprintf(fp,"% 10d",oppeid);
	}
	else {
	  /* Boundary edge - any necessary reordered writing of nodes
	     so that the face uses the edge in the +ve sense was done
	     earlier */

	  fprintf(fp,"% 10d",ME_ID(edge));
	}

	k++;
	if (k%14 == 0 && je != nfe-1)
	  fprintf(fp,"\n");
      }

      List_Delete(fedges);

      fprintf(fp,"\n");
    }

    idx = 0;
    while ((edge = MESH_Next_Edge(mesh,&idx)))
      MEnt_Rem_AttVal(edge,oppatt);
    MAttrib_Delete(oppatt);

  }

  fprintf(fp,"end_cells\n");


  /***********************************************************************/
  /* Slaved Nodes Data Block - not written out                           */
  /***********************************************************************/

  fprintf(fp,"slaved_nodes       0\n");
  fprintf(fp,"end_slaved_nodes\n");

  /***********************************************************************/
  /* Ghost Nodes Data Block - not written out                            */
  /***********************************************************************/

  fprintf(fp,"ghost_nodes        0\n");
  fprintf(fp,"end_ghost_nodes\n");


  /***********************************************************************/
  /* Cell centered data                                                  */
  /***********************************************************************/

  fprintf(fp,"cell_data\n");

  /***********************************************************************/
  /* Material ID data block                                              */
  /***********************************************************************/

  fprintf(fp,"matid\n");
  if (nr) {
    k = 0;
    idx = 0;
    while ((region = MESH_Next_Region(mesh,&idx))) {
      gentid = MR_GEntID(region);
      if (!gentid) gentid = 1;
      fprintf(fp,"% 10d",gentid);
      k++;
      if (k%10 == 0 || k == nr)
	fprintf(fp,"\n");
    }
  }
  else {
    k = 0;
    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {
      gentid = MF_GEntID(face);
      if (!gentid) gentid = 1;
      fprintf(fp,"% 10d",gentid);
      k++;
      if (k%10 == 0 || k == nf)
	fprintf(fp,"\n");
    }
  }
  fprintf(fp,"end_matid\n");
 
  /***********************************************************************/
  /* Partition/Processor Data for cells - Not used but must be present   */
  /* Since we are not handling distributed meshes, processor ID for all  */
  /* cells is 1                                                          */
  /***********************************************************************/

  fprintf(fp,"partelm\n");
  if (nr) {
    for (jr = 0; jr < nr; jr++) {
      fprintf(fp,"% 10d",1);
      if ((jr+1)%10 == 0 || jr == nr-1)
	fprintf(fp,"\n");
    }
  }
  else {
    for (jf = 0; jf < nf; jf++) {
      fprintf(fp,"% 10d",1);
      if ((jf+1)%10 == 0 || jf == nf-1)
	fprintf(fp,"\n");
    }
  }
  fprintf(fp,"end_partelm\n");


  /* Write out any other cell based attributes, if requested */

  for (i = 0; i < ncellatt; i++) {
    attrib = cellatts[i];
    atttype = MAttrib_Get_Type(attrib);
    
    MAttrib_Get_Name(attrib,attname);
    fprintf(fp,"%s ",attname);
    
    if (nr) {
      idx = 0;
      while ((region = MESH_Next_Region(mesh,&idx))) {
	
	MEnt_Get_AttVal(region,attrib,&ival,&rval,&pval);
	
	if (atttype == INT)
	  fprintf(fp,"% 20.12E\n",(double)ival);
	else if (atttype == DOUBLE)
	  fprintf(fp,"% 20.12E\n", rval);
      }
    }
    else {
      idx = 0;
      while ((face = MESH_Next_Face(mesh,&idx))) {
	
	MEnt_Get_AttVal(face,attrib,&ival,&rval,&pval);
	
	if (atttype == INT)
	  fprintf(fp,"% 20.12E\n",(double)ival);
	else if (atttype == DOUBLE)
	  fprintf(fp,"% 20.12E\n", rval);
      }
    }

    strcpy(tmpstr,"end_");
    strcat(tmpstr,attname);
    fprintf(fp,"%s\n",tmpstr);
  }

  if (cellatts) MSTK_free(cellatts);
  
  fprintf(fp,"end_cell_data\n");
  
  
  /***********************************************************************/
  /* Node centered data                                                  */
  /***********************************************************************/

  fprintf(fp,"node_data\n");

  for (i = 0; i < ncellatt; i++) {
    attrib = nodatts[i];
    atttype = MAttrib_Get_Type(attrib);

    MAttrib_Get_Name(attrib,attname);
    fprintf(fp,"%s ",attname);
    
    idx = 0;
    while ((vertex = MESH_Next_Vertex(mesh,&idx))) {
      
      MEnt_Get_AttVal(vertex,attrib,&ival,&rval,&pval);
      
      if (atttype == INT)
	fprintf(fp,"% 20.12E\n",(double)ival);
      else if (atttype == DOUBLE)
	fprintf(fp,"% 20.12E\n", rval);
    }

    strcpy(tmpstr,"end_");
    strcat(tmpstr,attname);
    fprintf(fp,"%s\n",tmpstr);
  }

  if (nodatts) MSTK_free(nodatts);
  
  fprintf(fp,"end_node_data\n");

  /* Finish Export */

  fprintf(fp,"end_dump\n");

  fclose(fp);



  /* Clean up */

  MSTK_free(gentities);

  return 1;
}



#ifdef __cplusplus
  }
#endif


