#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"
#include "MSTK_private.h"


#ifdef __cplusplus
extern "C" {
#endif



  /* Function to export a partitioned surface mesh to the Parallel FLAG
     X3D format (LA-UR-04-9033) */

  /* if natt = 0, all attributes are written out if natt = -1, no
     attributes are written out if natt > 0, only attributes specified
     in attnames are written out Except for cell material IDs,
     attributes are ignored by FLAG; still we will write them out for
     completeness */

  /* opts is an array of flags that controls how mesh is exported to FLAG 
     Currently, it is a dummy argument
  */

  /* It is assumed that the mesh has been partitioned and the
     processor ID and Local Number have been assigned for each mesh
     element of the highest order (faces in this case) */

int MESH_Surf_ExportToFLAGX3D_Par(Mesh_ptr mesh, const char *filename,  
			     const int nparts, const int natt, 
			     const char **attnames, int *opts, int *fprocids) {
  int			gentid, *gentities;
  List_ptr	        efaces, fedges, vfaces;
  List_ptr              localcells, localedges, localnodes;
  MVertex_ptr           mv;
  MEdge_ptr             me;
  MFace_ptr	        mf;
  MAttrib_ptr           attrib, *nodatts=NULL, *cellatts=NULL, oppatt;
  MAttType              atttype;
  char                  attname[256], matname[256], tmpstr[256];
  char                  partfilename[256], basename[256], mapfilename[256];
  int                   jv, je, len, suff, fid;
  int                   nef, nfe;
  int			i, ip, found, k, nmeshatt, ival, npb;
  int                   nalloc, ngent;
  int                   attentdim, j, idx, idx2;
  int                   ndup, max_nfe;
  int                   oppeid, ne2;
  int                   nnodatt, ncellatt, pbdry, *lnumarr;
  int                   *mvprocs, *meprocs, *mfprocs, pnum, minprocid;
  int                   neloc, nvloc, nfloc, nfp, nep, nvp, lnum;
  double		vxyz[3], rval;
  void                 *pval;
  FILE		        *fp;

  if (nparts < 2) {
    MESH_ExportToFLAGX3D(mesh,filename,natt,attnames,opts);
    return 1;
  }


  /* Extract the basename from the filename */

  strcpy(basename,filename);
  len = strlen(basename);
  suff = 0;
  if (len > 4) {
    k = len-4;
    while (!suff && k > 1) {
      if (strncmp(&(basename[k]),".x3d",4) == 0) {
        suff = 1;
      }
      else
        k--;
    }
  }
  if (suff)
    basename[len-4] = '\0';

  fprintf(stderr,"File Basename: %s\n",basename);




  /* Catalog the materials */
  
  
  ngent = 0;
  nalloc = 10;
  gentities = (int *) MSTK_malloc(nalloc*sizeof(int));
  
  idx = 0;
  while ((mf = MESH_Next_Face(mesh,&idx))) {
    gentid = MF_GEntID(mf);
    
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
  
  if (!ngent) {
    /* No materials. Make all cells of material type 1 */
    gentities[ngent] = 1;
    ngent++;
  }




  /* Initialize some attributes for writing out partitioned meshes */

  MESH_Init_ParAtts(mesh);


  /* Attach a local number and processor ID for each cell */

  lnumarr = (int *) calloc(nparts,sizeof(int));

  idx = 0;
  while ((mf = MESH_Next_Face(mesh,&idx))) {
    fid = MF_ID(mf);
    
    pnum = fprocids[fid-1];
    (lnumarr[pnum])++;
    
    MEnt_Set_ProcIDs(mf,1,&pnum); /* Resides on only one processor */
    MEnt_Set_LocalID(mf,pnum,lnumarr[pnum]);
  }



  pbdry = MSTK_GetMarker();
  oppatt = MAttrib_New(mesh,"oppedgeID",INT,MEDGE); 
      

  /***********************************************************************/
  /* Attach local number and processor data for each node                */
  /***********************************************************************/

  mfprocs = (int *) malloc(1*sizeof(int));
  mvprocs = (int *) malloc(nparts*sizeof(int));

  for (i = 0; i < nparts; i++) lnumarr[i] = 0;

  idx = 0;
  while ((mv = MESH_Next_Vertex(mesh,&idx))) {
    
    vfaces = MV_Faces(mv);
    
    nvp = 0;

    idx2 = 0;
    while ((mf = List_Next_Entry(vfaces,&idx2))) {

      MEnt_ProcIDs(mf,&nfp,mfprocs);

      if (nfp > 1) {
	MSTK_Report("MESH_Surf_ExportToFlagX3D","Cell of surface mesh is on more than two processors",MSTK_FATAL);
      }
      else if (nfp == 0) {
	MSTK_Report("MESH_Surf_ExportToFlagX3D","Cell of surface mesh is not on any processors",MSTK_FATAL);
      }

      pnum = mfprocs[0];             

      /* If there are cells from two different processors connected to
	 the node, then this a node on a partition boundary */

      found = 0;
      for (i = 0; i < nvp; i++) {
	if (mvprocs[i] == pnum) {
	  found = 1;
	  break;
	}
      }

      if (!found) {
	mvprocs[nvp] = pnum;
	nvp++;
      }
    }
    List_Delete(vfaces);
    

    /* If more than one processor is connected to node, Mark it as
	 being connected to a boundary node */

    if (nvp > 1)
      MEnt_Mark(mv,pbdry);


    /* Attach processor and local number data to node */

    MEnt_Set_ProcIDs(mv,nvp,mvprocs);

    for (i = 0; i < nvp; i++) {
      pnum = mvprocs[i];
      lnumarr[pnum]++;
      MEnt_Set_LocalID(mv,mvprocs[i],lnumarr[pnum]);
    }	
  }





  /***********************************************************************/
  /* Assign local number and processor data to all edges                 */
  /***********************************************************************/
  

  for (i = 0; i < nparts; i++) lnumarr[i] = 0;
  
  meprocs = (int *) malloc(2*sizeof(int));
  
  idx = 0;
  while ((me = MESH_Next_Edge(mesh,&idx))) {

    efaces = ME_Faces(me);
    nef = List_Num_Entries(efaces);

    if (nef > 2) {
      fprintf(stderr,"FLAG X3D format does not support non-manifold edges\n");
      return 0;
    }

    /* Attach processor and local number data to edge */

    for (i = 0; i < nef; i++) {

      mf = List_Entry(efaces,i);
      
      MEnt_ProcIDs(mf,&nfp,mfprocs);
      pnum = mfprocs[0];

      if (i == 0) {
	meprocs[0] = pnum;
	nep = 1;
      }
      else {
	/* if the two faces connected to the edge are on different
	   processors, mark edge as being on a processor boundary */

	if (pnum != meprocs[0]) {
	  
	  meprocs[1] = pnum;
	  nep = 2;

	  MEnt_Mark(me,pbdry);
	}	  
      }
    }
    
    List_Delete(efaces);


    /* Set processor and local number data for edge */

    MEnt_Set_ProcIDs(me,nep,meprocs);
    
    for (i = 0; i < nep; i++) {
      pnum = meprocs[i];
      lnumarr[pnum]++;
      MEnt_Set_LocalID(me,pnum,lnumarr[pnum]);
    }

  }





  /***********************************************************************/
  /* Now write out each partition                                        */
  /***********************************************************************/


  for (ip = 0; ip < nparts; ip++) {

    /* Open the file */

    sprintf(partfilename,"%s.x3d.%05d",basename,ip+1);
    if (!(fp = fopen(partfilename,"w"))) {
      fprintf(stderr,"MESH_ExportToFLAGX3D: Couldn't open output file %s\n",
	      partfilename);
      exit(2);
    }
    


    /* Compile a list of cells that are local to this processor */
  
    localcells = List_New(10);

    max_nfe = 0;
    idx = 0;
    while ((mf = MESH_Next_Face(mesh,&idx))) {

      MEnt_ProcIDs(mf,&nfp,mfprocs);

      if (nfp != 1)
	MSTK_Report("MESH_Surf_ExportToFlagX3D","Cell of surface mesh must be on one processor only",MSTK_FATAL);
       
      if (mfprocs[0] == ip) {
	List_Add(localcells,mf);

	
	nfe = MF_Num_Edges(mf);
	if (nfe > max_nfe)
	  max_nfe = nfe;
      }
    }

    nfloc = List_Num_Entries(localcells);




    /* Compile a list of edges (or "faces" according to FLAG
       terminology) local to this processor */
    /* While it would be more efficient to compile this list by going
       through the local cells and adding their edges to the list,
       this generates non-sequential numbering in the X3D files. This
       causes problems for FLAG which ignores the explicit number
       given to an edge or a node and uses the implicit number of the
       edge or node based on where it occurs in the list*/

    localedges = List_New(10);
    idx = 0;
    while ((me = MESH_Next_Edge(mesh,&idx))) {

      efaces = ME_Faces(me);

      idx2 = 0;
      while ((mf = List_Next_Entry(efaces,&idx2))) {

	MEnt_ProcIDs(mf,&nfp,mfprocs);

	if (mfprocs[0] == ip) { /* At least one face is on processor ip */
	  List_Add(localedges,me);
	  break;
	}

      }

      List_Delete(efaces);

    }

    neloc = List_Num_Entries(localedges);


    /* Compile a list of nodes local to this processor */

    npb = 0;           /* Number of vertices on the processor boundary */

    localnodes = List_New(10);

    idx = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx))) {

      vfaces = MV_Faces(mv);

      idx2 = 0;
      while ((mf = List_Next_Entry(vfaces,&idx2))) {

	MEnt_ProcIDs(mf,&nfp,mfprocs);

	if (mfprocs[0] == ip) { /* At least one face is on processor ip */

	  List_Add(localnodes,mv);

	  if (MEnt_IsMarked(mv,pbdry))
	    npb++;

	  break;
	}
      }

      List_Delete(vfaces);

    }

    nvloc = List_Num_Entries(localnodes);


    /***********************************************************************/
    /* Now write the file out                                              */
    /***********************************************************************/
    /* Opening line                                                        */
    /***********************************************************************/
    
    
    fprintf(fp,"x3dtoflag ascii\n");
  
  
    /***********************************************************************/
    /* Header information */
    /***********************************************************************/

    fprintf(fp,"header\n");

    /* Processor number - For now we will not bother with distributed I/O */
    
    fprintf(fp,"   %-22s %10d\n","process",ip+1);

    /* Problem dimension - 2 since we only have faces */
  
    fprintf(fp,"   %-22s %10d\n","numdim",2);
  
  
    /* Number of materials */
    
    fprintf(fp,"   %-22s %10d\n","materials",ngent);
  
    /* Number of nodes */
  
    fprintf(fp,"   %-22s %10d\n","nodes",nvloc);
  
    /* Number of "faces"; "faces" is the unfortunate term used in the
       FLAG X3D file to mean boundaries of elements/cells/zones in FLAG
       X3D, i.e., in MSTK parlance, 'faces' are edges in a surface mesh */
    
    ndup = 0;
    idx = 0;
    while ((me = List_Next_Entry(localedges,&idx))) {      
	
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
	
      efaces = ME_Faces(me);
      nef = List_Num_Entries(efaces);

      if (nef == 2 && !MEnt_IsMarked(me,pbdry)) { 
	ndup++;
	MEnt_Set_AttVal(me,oppatt,neloc+ndup,0,NULL);
      }
      
      if (efaces)
	List_Delete(efaces);
    }
    
    ne2 = neloc + ndup;
    
    fprintf(fp,"   %-22s %10d\n","faces",ne2);
    fprintf(fp,"   %-22s %10d\n","elements",nfloc);
  
    /* number of ghost nodes, i.e., number of nodes on partition boundaries */
  
    fprintf(fp,"   %-22s %10d\n","ghost_nodes",npb);
  
    /* No slave (constrained) node info - there are no mechanisms in
       place to do this in MSTK and therefore, this info cannot be
       written out */
    
    fprintf(fp,"   %-22s %10d\n","slaved_nodes",0);
    fprintf(fp,"   %-22s %10d\n","nodes_per_slave",2);  /* required default */
    
    /* Nodes per "face", required for array dimensioning in FLAG. If
       this is a surface mesh, this will be 2; if it is a volume mesh,
       it is the maximum number of nodes that any face in the mesh
       has */
    
    fprintf(fp,"   %-22s %10d\n","nodes_per_face",2); /* edges */
    
    /* "faces" per cell, required for array dimensioning in
     FLAG. Since this is a surface mesh, this the maximum number of
     edges any face has */
  
    fprintf(fp,"   %-22s %10d\n","faces_per_cell",max_nfe);
  
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
	
	
	if ((attentdim == MFACE) || (attentdim == MALLTYPE)) {
	  cellatts[ncellatt++] = attrib;
	}
      }
    }

    /* Apparently writing out non-standard node data is causing problems - so skip */
    /*
    fprintf(fp,"   %-22s %10d\n","node_data_fields",nnodatt);
    */
    fprintf(fp,"   %-22s %10d\n","node_data_fields",0);
    
    /* must write out matid and partelm and then add additional cell data */
    /* Apparently writing out non-standard cell data is causing problems - so write out only matid and partelm data */
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

    idx = 0;
    while ((mv = List_Next_Entry(localnodes,&idx))) {
      MV_Coords(mv,vxyz);
      fprintf(fp,"% 10d % 22.14E % 22.14E % 22.14E\n",MEnt_LocalID(mv,ip),
	      vxyz[0],vxyz[1],vxyz[2]);
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

    /***********************************************************/
    /* Surface mesh (maybe planar)                             */
    /* The info line for each "face" in a 2D mesh (i.e. edge)  */
    /* will only be 12 items long; so we needn't worry about   */
    /* special care to satisfy the 13I10 format                */
    /***********************************************************/
    
    idx = 0;
    while ((me = List_Next_Entry(localedges,&idx))) {
      
      efaces = ME_Faces(me);
      nef = List_Num_Entries(efaces);

      lnum = MEnt_LocalID(me,ip);
      fprintf(fp,"% 10d",lnum);
      
      if (nef == 2) {

	if (!MEnt_IsMarked(me,pbdry)) {

	  /* duplicate edge is on processor */

	  /* Both faces connected to edge are on processore - Write
	     edge out in its natural orientation here, write out the
	     duplicate edge in the opposite orientation later */
	  
	  fprintf(fp,"% 10d",2);
	  
	  for (jv = 0; jv < 2; jv++) {
	    mv = ME_Vertex(me,jv);
	    fprintf(fp,"% 10d",MEnt_LocalID(mv,ip));
	  }
	  
	  MEnt_Get_AttVal(me,oppatt,&oppeid,&rval,&pval);
#ifdef DEBUG
	  if (!oppeid)
	    MSTK_Report("MESH_Surf_ExportToFLAGX3D",
			"Internal edge has no duplicate edge info?",MSTK_ERROR);
#endif

	  /* ID of the processor owning this edge */
	
	  fprintf(fp,"% 10d",ip+1);

	  /* ID of the processor owning the duplicate edge */
      
	  fprintf(fp,"% 10d",ip+1);
	
	  /* ID of the duplicate edge */
	
	  fprintf(fp,"% 10d",oppeid);
	}
	else {
	  
	  /* duplicate edge is off processor */

	  /* Only one face connected to edge is on the processor, so
	     it is like a boundary edge - If the edge is used in the
	     +ve sense by the face, then write the edge out in its
	     natural orientation; if not, write it out in the opposite
	     orientation */

	  fprintf(fp,"% 10d",2);

	  /* Find which of the connected faces is on processor */

	  mf = List_Entry(efaces,0);
	  MEnt_ProcIDs(mf,&nfp,mfprocs);

	  if (mfprocs[0] != ip)
	    mf = List_Entry(efaces,1);
	  
	  if (MF_EdgeDir(mf,me) == 1) {
	    for (jv = 0; jv < 2; jv++) {
	      mv = ME_Vertex(me,jv);
	      fprintf(fp,"% 10d",MEnt_LocalID(mv,ip));
	    }
	  }
	  else {
	    for (jv = 0; jv < 2; jv++) {
	      mv = ME_Vertex(me,!jv);
	      fprintf(fp,"% 10d",MEnt_LocalID(mv,ip));
	    }
	  }
	  

	  /* ID of the processor owning this edge */
	
	  fprintf(fp,"% 10d",ip+1);


	  MEnt_ProcIDs(me,&nep,meprocs);
#ifdef DEBUG
	  if (nep > 2)
	    MSTK_Report("MESH_Surf_ExportToFlagX3D",
			"Non-manifold edges not supported in X3D format",
			MSTK_ERROR);
#endif	  
	
	  if (meprocs[0] != ip) {
	    
	    /* ID of the processor owning the duplicate edge on the
	       neighboring processor */
      
	    fprintf(fp,"% 10d",meprocs[0]+1);
	
	    /* Local ID of the duplicate edge on the neighboring processor */
	
	    lnum = MEnt_LocalID(me,meprocs[0]);

	    fprintf(fp,"% 10d",lnum);

	  }
	  else {

	    /* ID of the processor owning the duplicate edge on the
	       neighboring processor */
      
	    fprintf(fp,"% 10d",meprocs[1]+1);
	
	    /* Local ID of the duplicate edge on the neighboring processor */
	
	    lnum = MEnt_LocalID(me,meprocs[1]);

	    fprintf(fp,"% 10d",lnum);

	  }
	}

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

	if (MF_EdgeDir(List_Entry(efaces,0),me) == 1) {
	  for (jv = 0; jv < 2; jv++) {
	    mv = ME_Vertex(me,jv);
	    fprintf(fp,"% 10d",MEnt_LocalID(mv,ip));
	  }
	}
	else {
	  for (jv = 0; jv < 2; jv++) {
	    mv = ME_Vertex(me,!jv);
	    fprintf(fp,"% 10d",MEnt_LocalID(mv,ip));
	  }
	}


	/* ID of the processor owning this edge */
	
	fprintf(fp,"% 10d",ip+1);
	
	/* Boundary face - no duplicate edge info */
	
	fprintf(fp,"% 10d",0); /* id of processor owning duplicate edge */
	fprintf(fp,"% 10d",0); /* duplicate edge id */

	/* Five dummy arguments required by FLAG X3D format specification */
	
	for (i = 0; i < 5; i++) fprintf(fp,"% 10d",1);
	fprintf(fp,"\n");
      }

      List_Delete(efaces);

    }


    idx = 0;
    while ((me = List_Next_Entry(localedges,&idx))) {

      MEnt_Get_AttVal(me,oppatt,&oppeid,&rval,&pval);

      if (!oppeid) continue;

      fprintf(fp,"% 10d",oppeid);
      
      fprintf(fp,"% 10d",2);

      for (jv = 0; jv < 2; jv++) {
	mv = ME_Vertex(me,!jv);

	fprintf(fp,"% 10d",MEnt_LocalID(mv,ip));
      }
      
      /* ID of the processor owning this edge */
      
      fprintf(fp,"% 10d",ip+1);
      
      /* ID of the processor owning the duplicate edge */
      
      fprintf(fp,"% 10d",ip+1);
      
      /* ID of the duplicate (in this case, original) face */

      lnum = MEnt_LocalID(me,ip);
      fprintf(fp,"% 10d",lnum);

      /* Five dummy arguments required by FLAG X3D format specification */
      
      for (i = 0; i < 5; i++) fprintf(fp,"% 10d",1);
      fprintf(fp,"\n");
    }

    fprintf(fp,"end_faces\n");
  

    /***********************************************************************/
    /* Cell Data Block                                                     */
    /***********************************************************************/
    
    fprintf(fp,"cells\n");

    idx = 0;
    while ((mf = List_Next_Entry(localcells,&idx))) {
      fprintf(fp,"% 10d",MEnt_LocalID(mf,ip));
      k = 1;

      fedges = MF_Edges(mf,1,0);
      nfe = List_Num_Entries(fedges);

      fprintf(fp,"% 10d",nfe);
      k++;

      for (je = 0; je < nfe; je++) {
	me = List_Entry(fedges,je);

	MEnt_Get_AttVal(me,oppatt,&oppeid,&rval,&pval);
	if (oppeid) {
	  /* Interior edge */
	  /* Write out the edge id, if the face uses the edge in the
	     +ve sense; write out the duplicate edge id, if not */

	  if (MF_EdgeDir_i(mf,je))
	    fprintf(fp,"% 10d",MEnt_LocalID(me,ip));
	  else
	    fprintf(fp,"% 10d",oppeid);
	}
	else {
	  /* Boundary edge - any necessary reordered writing of nodes
	     so that the face uses the edge in the +ve sense was done
	     earlier */

	  fprintf(fp,"% 10d",MEnt_LocalID(me,ip));
	}

	k++;
	if (k%14 == 0 && je != nfe-1)
	  fprintf(fp,"\n");
      }

      List_Delete(fedges);

      fprintf(fp,"\n");
    }

    fprintf(fp,"end_cells\n");


    /***********************************************************************/
    /* Slaved Nodes Data Block - not written out                           */
    /***********************************************************************/
    
    fprintf(fp,"slaved_nodes       0\n");
    fprintf(fp,"end_slaved_nodes\n");
    
    /***********************************************************************/
    /* Ghost Nodes Data Block - Assumes that for a node on a processor     */
    /* boundary, the image of the node on the processor with the smallest  */
    /* ID is the master and the others are slaves                          */
    /***********************************************************************/
    
    fprintf(fp,"ghost_nodes %10d\n",npb);

    idx = 0;
    while ((mv = List_Next_Entry(localnodes,&idx))) {

      if (!MEnt_IsMarked(mv,pbdry)) continue;

      MEnt_ProcIDs(mv,&nvp,mvprocs);

      minprocid = nparts;
      for (i = 0; i < nvp; i++)
	if (mvprocs[i] < minprocid)
	  minprocid = mvprocs[i];

      lnum = MEnt_LocalID(mv,minprocid);

      fprintf(fp,"% 10d",MEnt_LocalID(mv,ip));
      fprintf(fp,"% 10d",minprocid+1);
      fprintf(fp,"% 10d",lnum);
      fprintf(fp,"% 10d\n",1);
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

    k = 0;
    idx = 0;
    while ((mf = List_Next_Entry(localcells,&idx))) {

      gentid = MF_GEntID(mf);
      if (!gentid) gentid = 1;

      fprintf(fp,"% 10d",gentid);

      k++;
      if (k%10 == 0 || k == nfloc)
	fprintf(fp,"\n");
    }
    fprintf(fp,"end_matid\n");
 
    /***********************************************************************/
    /* Partition/Processor Data for cells - Not used but must be present   */
    /* Since we are not handling distributed meshes, processor ID for      */
    /* all cells is 1                                                      */
    /***********************************************************************/

    fprintf(fp,"partelm\n");

    idx = 0; i = 0;
    while ((mf = List_Next_Entry(localcells,&idx))) {

      MEnt_ProcIDs(mf,&nfp,mfprocs);      

      fprintf(fp,"% 10d",mfprocs[0]+1);

      if ((i+1)%10 == 0 || i == nfloc-1)
	fprintf(fp,"\n");

      i++;
    }
    fprintf(fp,"end_partelm\n");


    /***********************************************************************/
    /* Write out any other cell based attributes, if requested             */
    /***********************************************************************/
    
    /* Causing problems */
    /*
    for (i = 0; i < ncellatt; i++) {
      attrib = cellatts[i];
      atttype = MAttrib_Get_Type(attrib);
      
      MAttrib_Get_Name(attrib,attname);
      fprintf(fp,"%s ",attname);
    
      idx = 0;
      while ((mf = List_Next_Entry(localcells,&idx))) {
	
	MEnt_Get_AttVal(mf,attrib,&ival,&rval,&pval);
	
	if (atttype == INT)
	  fprintf(fp,"% 20.12E\n",(double)ival);
	else if (atttype == DOUBLE)
	  fprintf(fp,"% 20.12E\n", rval);
      }

      strcpy(tmpstr,"end_");
      strcat(tmpstr,attname);
      fprintf(fp,"%s\n",tmpstr);
    }
    */
    if (cellatts) MSTK_free(cellatts);
  
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
      while ((mv = List_Next_Entry(localnodes,&idx))) {
	
	MEnt_Get_AttVal(mv,attrib,&ival,&rval,&pval);
	
	if (atttype == INT)
	  fprintf(fp,"% 20.12E\n",(double)ival);
	else if (atttype == DOUBLE)
	  fprintf(fp,"% 20.12E\n", rval);
      }
      
      strcpy(tmpstr,"end_");
      strcat(tmpstr,attname);
      fprintf(fp,"%s\n",tmpstr);
    }
    */
    if (nodatts) MSTK_free(nodatts);
    
    fprintf(fp,"end_node_data\n");


    /* Finish Export */

    fprintf(fp,"end_dump\n");
    
    fclose(fp);


   

    /* Write out the mapping from local cell numbers to global cell
       numbers */

    sprintf(mapfilename,"%s.map.%05d",basename,ip+1);
    if (!(fp = fopen(mapfilename,"w"))) {
      fprintf(stderr,"MESH_ExportToFLAGX3D: Couldn't open output file %s\n",
	      mapfilename);
      exit(2);
    }
    

    fprintf(fp,"%d\n",List_Num_Entries(localcells));
    
    idx = 0;
    while ((mf = List_Next_Entry(localcells,&idx)))
      fprintf(fp,"%d %d\n",MEnt_LocalID(mf,ip),MEnt_ID(mf));


    fclose(fp);



    /* Clean up for this processor output */

    idx = 0;
    while ((me = List_Next_Entry(localedges,&idx)))
      MEnt_Rem_AttVal(me,oppatt);

    List_Delete(localnodes);
    List_Delete(localedges);
    List_Delete(localcells);

  } /* End loop on partitions */



  /* Clean up */

  idx = 0;
  while ((mv = MESH_Next_Vertex(mesh,&idx)))
    MEnt_Unmark(mv,pbdry);
  idx = 0;
  while ((me = MESH_Next_Edge(mesh,&idx)))
    MEnt_Unmark(me,pbdry);
  MSTK_FreeMarker(pbdry);
  MAttrib_Delete(oppatt);

  MSTK_free(gentities);


  return 1;
}

#ifdef __cplusplus
  }
#endif

