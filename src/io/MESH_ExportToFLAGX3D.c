#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"
#include "MSTK_private.h"

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


int MESH_ExportToFLAGX3D(Mesh_ptr mesh, const char *filename, const int natt, 
			 const char **attnames, const int *opts, MSTK_Comm comm) {
  int			gentid, *gentities;
  List_ptr	        fverts, rfaces, fregs, efaces, fedges;
  MVertex_ptr           vertex;
  MEdge_ptr             edge;
  MFace_ptr	        face;
  MRegion_ptr           region;
  MAttrib_ptr           attrib, *nodatts=NULL, *cellatts=NULL;
  MAttrib_ptr           vidatt, eidatt, fidatt, ridatt, oppatt, opppidatt;
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
  int                   vid, eid, fid, rid;
  double		vxyz[3], rval;
  void                 *pval;
  FILE		        *fp;

#ifdef MSTK_HAVE_MPI
  int pid, numprocs;
  MPI_Comm_size(comm,&numprocs);
  MPI_Comm_rank(comm,&pid);
  pid += 1;  /* FLAG X3D counts processor IDs from 1 */
#endif

  vidatt = MAttrib_New(mesh,"vidatt",INT,MVERTEX);
  eidatt = MAttrib_New(mesh,"eidatt",INT,MEDGE);
  fidatt = MAttrib_New(mesh,"fidatt",INT,MFACE);
  ridatt = MAttrib_New(mesh,"ridatt",INT,MREGION);
  
  idx = 0; i = 0;
  while ((vertex = MESH_Next_Vertex(mesh,&idx)))
    MEnt_Set_AttVal(vertex,vidatt,++i,0.0,NULL);
  
  idx = 0; i = 0;
  while ((edge = MESH_Next_Edge(mesh,&idx)))
    MEnt_Set_AttVal(edge,eidatt,++i,0.0,NULL);
  
  idx = 0; i = 0;
  while ((face = MESH_Next_Face(mesh,&idx)))
    MEnt_Set_AttVal(face,fidatt,++i,0.0,NULL);
  
  idx = 0; i = 0;
  while ((region = MESH_Next_Region(mesh,&idx)))
    MEnt_Set_AttVal(region,ridatt,++i,0.0,NULL);
  
  
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
  
  fprintf(fp,"   %-22s %10d\n","process",pid);
  
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
    /* No materials. Make all cells of material type 1 */
    gentities[ngent] = 1;
    ngent++;
  }
  
  fprintf(fp,"   %-22s %10d\n","materials",ngent);
  
  /* Number of nodes */
  
  fprintf(fp,"   %-22s %10d\n","nodes",nv);
  
  /* Number of "faces"; "faces" is the term used in the FLAG X3D file
     to mean boundaries of elements/cells/zones in FLAG X3D, i.e., in
     MSTK parlance, 'faces' are edges in a surface mesh, faces in a
     volume mesh */
  
  /* Each 'face' in FLAG X3D must be referenced by only one cell. A
     'face' description in FLAG X3D is supposed to be such that normal
     of 'face' points out of its referring cell. Therefore, if a
     'face' is connected to two cells, then we need to write out the
     face and its duplicate, pointing in the opposite direction. If
     'face' is on a processor boundary, this duplicate ID has to be
     the local ID of the corresponding 'face' from the other
     processor. If a 'face' is connected to only one region, we need
     to write out the 'face' in its natural direction if its normal
     points out of the owning region; we need to write out the 'face'
     in the reverse direction if its normal points into the owning
     region */
  
  
  if (nr) {
    
    oppatt = MAttrib_New(mesh,"oppfaceID",INT,MFACE); 
    opppidatt = MAttrib_New(mesh,"oppfacePID",INT,MFACE); 
    
    /* First count the number of faces that are not on a partition boundary */
    
    int nfbdry = 0;
    int nfint = 0;
    int nfparbdry = 0;
    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {      
      
      fregs = MF_Regions(face);
      if (!fregs) 
	MSTK_Report("MESH_ExportToFLAGX3D","FLAGX3D cannot handle this type of non-manifold mesh",MSTK_FATAL);
      
      if (List_Num_Entries(fregs) == 2) { /* face not on domain boundary */
	
	MRegion_ptr freg0 = List_Entry(fregs,0);
	MRegion_ptr freg1 = List_Entry(fregs,1);
        
	/* count interior faces that are completely in the partition
	 * (nfint) and on partition boundaries (nfparbdry) but avoid
	 * counting interior to the ghost layer */
	
	if (MR_PType(freg0) != PGHOST && MR_PType(freg1) != PGHOST) 
	  nfint++;
	else
	  nfparbdry++;
      }
      else {
	/* count boundary faces but avoid counting faces that are on
	 * the outer boundary of the ghost layer */
	
	MRegion_ptr freg0 = List_Entry(fregs,0);
	if (MR_PType(freg0) != PGHOST) 
	  nfbdry++;
      }
      
      if (fregs)
	List_Delete(fregs);
    }
    
    
    int ndup_local = 0; /* Number of duplicate faces from this partition */
    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {      
      
      fregs = MF_Regions(face);
      if (fregs && List_Num_Entries(fregs) == 2) {
	
	int nghostreg = 0;
	MRegion_ptr freg0 = List_Entry(fregs,0);
	if (MR_PType(freg0) == PGHOST) nghostreg++;
	MRegion_ptr freg1 = List_Entry(fregs,1);
	if (MR_PType(freg1) == PGHOST) nghostreg++;
	
	if (nghostreg == 0) {
	  
	  /* face is interior to this partition - we will pretend
	     there are duplicate faces at the end of the regular face list */
          
	  ndup_local++;
	  int dupfaceid = nfint+nfbdry+ndup_local;
	  MEnt_Set_AttVal(face,oppatt,dupfaceid,0,NULL);
          
	}
	else if (nghostreg == 1) {
	  
	  /* face is on partition boundary - the duplicate face ID has
	     to be the local ID of the face from the adjacent
	     partition. Set the oppatt value to the local face ID for
	     now.  After all faces are processed, do communicate the
	     local ID to the opposite processor */
          
	  MEnt_Set_AttVal(face,oppatt,fid,0,NULL);
	  MEnt_Set_AttVal(face,opppidatt,pid,0,NULL);
	}
      }
      if (fregs)
	List_Delete(fregs);
    }
    
    int nftot = nfbdry + 2*nfint + nfparbdry;
    fprintf(fp,"   %-22s %10d\n","faces",nftot);
    
    int nrtot = 0;
    idx = 0;
    while ((region = MESH_Next_Region(mesh,&idx))) 
      if (MR_PType(region) != PGHOST)
	nrtot++;
    
    fprintf(fp,"   %-22s %10d\n","elements",nrtot);
    
  }
  else {
    oppatt = MAttrib_New(mesh,"oppedgeID",INT,MEDGE); 
    opppidatt = MAttrib_New(mesh,"oppedgePID",INT,MEDGE); 
    
    /* First count the number of faces that are not on a partition boundary */
    
    int nebdry = 0;
    int neint = 0;
    int neparbdry = 0;
    idx = 0;
    while ((edge = MESH_Next_Edge(mesh,&idx))) {      
      
      efaces = ME_Faces(edge);
      if (!efaces || List_Num_Entries(efaces) > 2) 
	MSTK_Report("MESH_ExportToFLAGX3D","FLAGX3D cannot handle this type of non-manifold mesh",MSTK_FATAL);
      
      if (List_Num_Entries(efaces) == 2) { /* edge not on domain boundary */
	
	MFace_ptr eface0 = List_Entry(efaces,0);
	MRegion_ptr eface1 = List_Entry(efaces,1);
        
	/* Count interior edges that are completely in the partition
	 * (neint) and on partition boundaries (neparbdry) but avoid
	 * counting edges interior to the ghost layer */
	
	if (MF_PType(eface0) != PGHOST && MF_PType(eface1) != PGHOST) 
	  neint++;
	else
	  neparbdry++;
      }
      else {
	/* Count boundary edges but avoid counting edges that are on
	 * the outer boundary of the ghost layer */
	
	MFace_ptr eface0 = List_Entry(efaces,0);
	if (MF_PType(eface0) != PGHOST) 
	  nebdry++;
      }
      
      if (efaces)
	List_Delete(efaces);
    }
    
    int ndup_local = 0; /* Number of duplicate faces from this partition */
    idx = 0;
    while ((edge = MESH_Next_Edge(mesh,&idx))) {      
      
      efaces = ME_Faces(edge);
      if (efaces && List_Num_Entries(efaces) == 2) {	  
	
	int nghostface = 0;
	MFace_ptr eface0 = List_Entry(efaces,0);
	if (MR_PType(eface0) == PGHOST) nghostface++;
	MRegion_ptr eface1 = List_Entry(efaces,1);
	if (MR_PType(eface1) == PGHOST) nghostface++;
	
	if (nghostface == 0) {
	  
	  /* edge is interior to this partition - we will pretend
	     there is a duplicate edge at the end of the regular edge list */
	  
	  ndup_local++;
	  int dupedgeid = neint+nebdry+ndup_local;
	  MEnt_Set_AttVal(edge,oppatt,dupedgeid,0,NULL);
	  
	}
	else if (nghostface == 1) {
	  
	  /* edge is on partition boundary - the duplicate edge ID has
	     to be the local ID of the edge from the adjacent
	     partition. Attach the edge ID as an attribute to the edge
	     and exchange to transmit the info */
	  
	  MEnt_Set_AttVal(edge,oppatt,eid,0,NULL);
	  MEnt_Set_AttVal(edge,opppidatt,pid,0,NULL);
	}
      }
      
      if (efaces)
	List_Delete(efaces);
    }
    
    int netot = nebdry + 2*neint + neparbdry;
    fprintf(fp,"   %-22s %10d\n","faces",netot);
    
    int nftot = 0;
    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx)))
      if (MF_PType(face) != PGHOST)
	nftot++;
    
    fprintf(fp,"   %-22s %10d\n","elements",nftot);
  }
  
  /* Exchange the local edge/face IDs between master and slave edges/faces */
  
#ifdef MSTK_HAVE_MPI
  MESH_XchngEdgeFaceAttrib(mesh,oppatt,comm);
  MESH_XchngEdgeFaceAttrib(mesh,opppidatt,comm);
#endif

  
  /* Make a list of vertices on partition boundaries - in doing so we
   * have to account for the fact that some vertices may be entirely
   * in the ghost layer and we have to exclude those */

  MAttrib_ptr vmasteratt = MAttrib_New(mesh,"vmasteratt",INT,MVERTEX);

  List_ptr prtn_bndry_verts = List_New(10);
  if (nr) {
    idx = 0;
    while ((vertex = MESH_Next_Vertex(mesh,&idx))) {
      if (MV_PType(vertex) == PINTERIOR) continue;

      List_ptr vregions = MV_Regions(vertex);
      int foundghost = 0, foundowned = 0;
      int idx2 = 0;
      MRegion_ptr vr;
      while ((vr = List_Next_Entry(vregions,&idx))) {
	if (MR_PType(vr) != PGHOST)
	  foundowned = 1;
	else 
	  foundghost = 1;
	if (foundowned && foundghost)
	  break;
      }
      List_Delete(vregions);
      
      if (foundowned && foundghost) { /* Here is a vertex on a partition boundary */
	List_Add(prtn_bndry_verts,vertex);
	if (MV_PType(vertex) == POVERLAP)  /* vertex is master on prtn bndry */
	  MEnt_Set_AttVal(vertex,vmasteratt,MV_ID(vertex),0.0,NULL);
      }
    }
  }
  else {
    idx = 0;
    while ((vertex = MESH_Next_Vertex(mesh,&idx))) {
      if (MV_PType(vertex) == PINTERIOR) continue;

      List_ptr vfaces = MV_Faces(vertex);
      int foundghost = 0, foundowned = 0;
      int idx2 = 0;
      MFace_ptr vf;
      while ((vf = List_Next_Entry(vfaces,&idx))) {
	if (MF_PType(vf) != PGHOST)
	  foundowned = 1;
	else 
	  foundghost = 1;
	if (foundowned && foundghost)
	  break;
      }
      List_Delete(vfaces);
      
      if (foundowned && foundghost) { /* Here is a vertex on a partition boundary */	
	List_Add(prtn_bndry_verts,vertex);
	if (MV_PType(vertex) == POVERLAP) /* vertex is master on prtn bndry */
	  MEnt_Set_AttVal(vertex,vmasteratt,MV_ID(vertex),0.0,NULL);
      }
    }
  }

  /* Transmit master vertex IDs to slave IDs */

#ifdef MSTK_HAVE_MPI
  MESH_Update1Attribute(mesh,vmasteratt,comm);
#endif
  
  
  fprintf(fp,"   %-22s %10d\n","ghost_nodes",List_Num_Entries(prtn_bndry_verts));
  
  /* No slave (constrained) node info - there are no mechanisms in
     place to do this in MSTK and therefore, this info cannot be
     written out */
  
  fprintf(fp,"   %-22s %10d\n","slaved_nodes",0);
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

  /* Apparently writing out non-standard node data is causing problems - so skip */
  /*
    fprintf(fp,"   %-22s %10d\n","node_data_fields",nnodatt);
  */
  fprintf(fp,"   %-22s %10d\n","node_data_fields",0);
  
  /* must write out matid and partelm and then add additional cell data */
  /* Apparently writing out non-standard node data is causing problems - so write out only matid and partelm data */
  /*
    fprintf(fp,"   %-22s %10d\n","cell_data_fields",ncellatt+2);
  */
  fprintf(fp,"   %-22s %10d\n","cell_data_fields",2);
  
  
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
    MEnt_Get_AttVal(vertex,vidatt,&vid,&rval,&pval);
    fprintf(fp,"% 10d % 22.14E % 22.14E % 22.14E\n",vid,
	    vxyz[0],vxyz[1],vxyz[2]);
  }
  fprintf(fp,"end_nodes\n");

  /***********************************************************************/
  /* Faces Data Block                                                    */
  /* Faces is the term that the FLAG X3D file uses to refer              */
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

      MEnt_Get_AttVal(face,fidatt,&fid,&rval,&pval);
      fprintf(fp,"% 10d",fid);
      k = 1;

      if (nfr == 2) {
	/* Two regions connected to face - write face out in its
	   natural orientation here, write out the duplicate face in
	   the opposite orientation later */

	MRegion_ptr freg0 = List_Entry(fregs,0);
	MRegion_ptr freg1 = List_Entry(fregs,1);
	List_Delete(fregs);

	if (MR_PType(freg0) == PGHOST && MR_PType(freg1) == PGHOST) {
	  /* This face is not in the partition interior or on the
	   * partition boundary - rather it is a face in the ghost
	   * layer - IGNORE */

	  continue;
	}


	fverts = MF_Vertices(face,1,0);
	nfv = List_Num_Entries(fverts);
	
	fprintf(fp,"% 10d",nfv);
	k++;
	
	for (jv = 0; jv < nfv; jv++) {
	  vertex = List_Entry(fverts,jv);
	  MEnt_Get_AttVal(vertex,vidatt,&vid,&rval,&pval);
	  fprintf(fp,"% 10d",vid);
	  if ((++k)%13 == 0)
	    fprintf(fp,"\n");
	}

	List_Delete(fverts);

	/* ID of the processor owning this face */

	fprintf(fp,"% 10d",pid+1);   
	if ((++k)%13 == 0)
	  fprintf(fp,"\n");
	
	/* ID of the processor owning the duplicate face */

	int oppfpid;
	MEnt_Get_AttVal(face,opppidatt,&oppfpid,&rval,&pval);
	if (oppfpid)
	  fprintf(fp,"% 10d",oppfpid); 
	else
	  fprintf(fp,"% 10d",pid+1);
	if ((++k)%13 == 0)
	  fprintf(fp,"\n");
	
	/* ID of the duplicate face */
	
	MEnt_Get_AttVal(face,oppatt,&oppfid,&rval,&pval);
#ifdef DEBUG
	if (!oppfid)
	  fprintf(stderr,"Non-boundary face has no duplicate face info?\n");
#endif

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

	MRegion_ptr freg0 = List_Entry(fregs,0);
	List_Delete(fregs);

	if (MR_PType(freg0) == PGHOST) {
	  /* face is on the outer boundary of the ghost layer - ignore */
	  continue;
	}

	if (MR_FaceDir(freg0,face) == 1)
	  fverts = MF_Vertices(face,1,0);
	else
	  fverts = MF_Vertices(face,0,0);

	nfv = List_Num_Entries(fverts);
	
	fprintf(fp,"% 10d",nfv);
	k++;
	
	for (jv = 0; jv < nfv; jv++) {
	  vertex = List_Entry(fverts,jv);
	  MEnt_Get_AttVal(vertex,vidatt,&vid,&rval,&pval);
	  fprintf(fp,"% 10d",vid);
	  if ((++k)%13 == 0)
	    fprintf(fp,"\n");
	}	

	List_Delete(fverts);

	/* ID of the processor owning this face */
	
	fprintf(fp,"% 10d",pid+1); 
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

      int oppfpid;
      MEnt_Get_AttVal(face,opppidatt,&oppfpid,&rval,&pval);
      if (oppfpid != pid) continue;  /* duplicate face from another processor */

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
	MEnt_Get_AttVal(vertex,vidatt,&vid,&rval,&pval);
	fprintf(fp,"% 10d",vid);
	if ((++k)%13 == 0)
	  fprintf(fp,"\n");
      }
      
      /* ID of the processor owning this face */
      fprintf(fp,"% 10d",oppfpid);
      if ((++k)%13 == 0)
	fprintf(fp,"\n");
      
      /* ID of the processor owning the duplicate (in this case, original) face */
      
      fprintf(fp,"% 10d",pid+1);
      if ((++k)%13 == 0)
	fprintf(fp,"\n");
      
      /* ID of the duplicate (in this case, original) face */
      
      MEnt_Get_AttVal(face,fidatt,&fid,&rval,&pval); 
      fprintf(fp,"% 10d",fid);
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

      MEnt_Get_AttVal(edge,eidatt,&eid,&rval,&pval);
      fprintf(fp,"% 10d",eid);
      k = 1;

      if (nef == 2) {
	/* Two faces connected to edge - write edge out in its
	   natural orientation here, write out the duplicate edge in
	   the opposite orientation later */

	MFace_ptr ef0 = List_Entry(efaces,0);
	MFace_ptr ef1 = List_Entry(efaces,1);
	List_Delete(efaces);

	if (MF_PType(ef0) == PGHOST && MF_PType(ef1) == PGHOST) {
	  /* This edge is not in the partition interior or on the
	   * partition boundary - rather it is a edge in the ghost
	   * layer - IGNORE */

	  continue;
	}

	fprintf(fp,"% 10d",2);
	
	for (jv = 0; jv < 2; jv++) {
	  vertex = ME_Vertex(edge,jv);
	  MEnt_Get_AttVal(vertex,vidatt,&vid,&rval,&pval);
	  fprintf(fp,"% 10d",vid);
	}


	/* ID of the processor owning this edge */
	
	fprintf(fp,"% 10d",pid+1);
	
	/* ID of the processor owning the duplicate edge */

	int oppepid;
	MEnt_Get_AttVal(edge,opppidatt,&oppepid,&rval,&pval);
	if (oppepid) 
	  fprintf(fp,"% 10d",oppepid);
	else
	  fprintf(fp,"% 10d",pid+1);

	/* ID of the duplicate edge */
	
	MEnt_Get_AttVal(edge,oppatt,&oppeid,&rval,&pval);
#ifdef DEBUG
	if (!oppeid)
	  fprintf(stderr,"Non-boundary edge has no duplicate edge info?\n");
#endif

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

	MFace_ptr ef0 = List_Entry(efaces,0);
	List_Delete(efaces);

	if (MF_PType(ef0) == PGHOST) {
	  /* This edge is on the outer boundary of the ghost layer */
	  continue;
	}

	fprintf(fp,"% 10d",2);

	if (MF_EdgeDir(ef0,edge) == 1) {
	  for (jv = 0; jv < 2; jv++) {
	    vertex = ME_Vertex(edge,jv);
	    MEnt_Get_AttVal(vertex,vidatt,&vid,&rval,&pval);
	    fprintf(fp,"% 10d",vid);
	  }
	}
	else {
	  for (jv = 0; jv < 2; jv++) {
	    vertex = ME_Vertex(edge,!jv);
	    MEnt_Get_AttVal(vertex,vidatt,&vid,&rval,&pval);
	    fprintf(fp,"% 10d",vid);
	  }
	}

	/* ID of the processor owning this edge */
	
	fprintf(fp,"% 10d",pid+1);
	
	/* Boundary face - no duplicate edge info */
	
	fprintf(fp,"% 10d",0); /* id of processor owning duplicate edge */
	fprintf(fp,"% 10d",0); /* duplicate edge id */

	/* Five dummy arguments required by FLAG X3D format specification */
	
	for (i = 0; i < 5; i++) fprintf(fp,"% 10d",1);
	fprintf(fp,"\n");
      }
#ifdef DEBUG
      else {
	if (nef == 0)
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

      int oppepid;
      MEnt_Get_AttVal(edge,opppidatt,&oppepid,&rval,&pval);
      if (oppepid != pid) continue; /* duplicate edge from another processor */

      fprintf(fp,"% 10d",oppeid);
      
      fprintf(fp,"% 10d",2);

      for (jv = 0; jv < 2; jv++) {
	vertex = ME_Vertex(edge,!jv);
	MEnt_Get_AttVal(vertex,vidatt,&vid,&rval,&pval);
	fprintf(fp,"% 10d",vid);
      }
      
      /* ID of the processor owning this edge */
      fprintf(fp,"% 10d",oppepid);
      
      /* ID of the processor owning the duplicate edge */
      
      fprintf(fp,"% 10d",pid+1);
      
      /* ID of the duplicate (in this case, original) face */
      
      MEnt_Get_AttVal(edge,eidatt,&eid,&rval,&pval);
      fprintf(fp,"% 10d",eid);

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
      if (MR_PType(region) == PGHOST) continue;

      MEnt_Get_AttVal(region,ridatt,&rid,&rval,&pval);
      fprintf(fp,"% 10d",rid);
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
	  
	  if (MR_FaceDir_i(region,jf)) {
	    MEnt_Get_AttVal(face,fidatt,&fid,&rval,&pval);
	    fprintf(fp,"% 10d",fid);
	  }
	  else
	    fprintf(fp,"% 10d",oppfid);
	}
	else {
	  /* Boundary face */
	  /* Write out the face id - any necessary reordered writing
	     of nodes so that the face points out of the region was
	     done earlier */
	  
	  MEnt_Get_AttVal(face,fidatt,&fid,&rval,&pval);
	  fprintf(fp,"% 10d",fid);
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
      if (MF_PType(face) == PGHOST) continue;

      MEnt_Get_AttVal(face,fidatt,&fid,&rval,&pval);
      fprintf(fp,"% 10d",fid);
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

	  if (MF_EdgeDir_i(face,je)) {
	    MEnt_Get_AttVal(edge,eidatt,&eid,&rval,&pval);
	    fprintf(fp,"% 10d",eid);
	  }
	  else
	    fprintf(fp,"% 10d",oppeid);
	}
	else {
	  /* Boundary edge - any necessary reordered writing of nodes
	     so that the face uses the edge in the +ve sense was done
	     earlier */

	  MEnt_Get_AttVal(edge,eidatt,&eid,&rval,&pval);
	  fprintf(fp,"% 10d",eid);
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

  fprintf(fp,"ghost_nodes        %10d\n",List_Num_Entries(prtn_bndry_verts));
  idx = 0;
  while ((vertex = List_Next_Entry(prtn_bndry_verts,&idx))) {
    int masterpid = MV_MasterParID(vertex);
    if (masterpid == pid) /* This is the master node */
      fprintf(fp,"% 10d % 10d % 10d % 10d\n",MV_ID(vertex),pid+1,MV_ID(vertex),1);
    else {
      int mastervid;
      MEnt_Get_AttVal(vertex,vmasteratt,&mastervid,&rval,&pval);
      fprintf(fp,"% 10d % 10d % 10d % 10d\n",MV_ID(vertex),masterpid,mastervid,1);
    }
  }
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
  k = 0;
  if (nr) {
    idx = 0; 
    while ((region = MESH_Next_Region(mesh,&idx))) 
      if (MR_PType(region) != PGHOST) {
	fprintf(fp,"% 10d",pid);
	if ((k+1)%10 == 0 || k == nr-1)
	  fprintf(fp,"\n");
        k++;
      }
  }
  else {
    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx)))
      if (MF_PType(face) != PGHOST) {
	fprintf(fp,"% 10d",pid);
	if ((k+1)%10 == 0 || k == nf-1)
	  fprintf(fp,"\n");
        k++;
      }

  }
  fprintf(fp,"end_partelm\n");


  /* Write out any other cell based attributes, if requested */

  /* Apparently these are causing problems so we will comment it out */

  /*
    for (i = 0; i < ncellatt; i++) {
    attrib = cellatts[i];
    atttype = MAttrib_Get_Type(attrib);
    
    MAttrib_Get_Name(attrib,attname);
    fprintf(fp,"%s ",attname);
    
    if (nr) {
    idx = 0;
    while ((region = MESH_Next_Region(mesh,&idx))) {
    if (MR_PType(region) == PGHOST) continue;

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
    if (MF_PType(face) == PGHOST) continue;

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
  */

  if (cellatts) free(cellatts);
  
  fprintf(fp,"end_cell_data\n");
  
  
  /***********************************************************************/
  /* Node centered data                                                  */
  /***********************************************************************/

  fprintf(fp,"node_data\n");

  /*
    for (i = 0; i < ncellatt; i++) {
      attrib = nodatts[i];
      atttype = MAttrib_Get_Type(attrib);

      MAttrib_Get_Name(attrib,attname);
      fprintf(fp,"%s ",attname);
    
      idx = 0;
      while ((vertex = MESH_Next_Vertex(mesh,&idx))) {
        \/\* have to check that this is not a pure ghost node surrounded by only ghost elements \*\/
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
    */
  
  fprintf(fp,"end_node_data\n");

  /* Finish Export */
  
  fprintf(fp,"end_dump\n");

  fclose(fp);

  
  
  /* Clean up */
  
  free(gentities);
  
  idx = 0;
  while ((vertex = MESH_Next_Vertex(mesh,&idx)))
    MEnt_Rem_AttVal(vertex,vidatt);
  
  idx = 0;
  while ((edge = MESH_Next_Edge(mesh,&idx)))
    MEnt_Rem_AttVal(edge,eidatt);
  
  idx = 0;
  while ((face = MESH_Next_Face(mesh,&idx)))
    MEnt_Rem_AttVal(face,fidatt);
  
  idx = 0;
  while ((region = MESH_Next_Region(mesh,&idx)))
    MEnt_Rem_AttVal(region,ridatt);
  
  
  MAttrib_Delete(vidatt);
  MAttrib_Delete(eidatt);
  MAttrib_Delete(fidatt);
  MAttrib_Delete(ridatt);
  
  return 1;
}


#ifdef __cplusplus
  }
#endif


