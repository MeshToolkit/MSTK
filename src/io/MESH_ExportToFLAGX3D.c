/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

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
  MAttrib_ptr           vidatt_tmp, eidatt_tmp, fidatt_tmp, ridatt_tmp, oppatt, opppidatt;
  MAttType              atttype;
  char                  attname[256], matname[256], tmpstr[256];
  char                  modfilename[256];
  int                   jv, je, jr, jf;
  int                   nv, ne, nf, nr, nrf, nfv, nef, nfe, nfr, nvout;
  int			i, found, k, nmeshatt, ival=0, ndim;
  int                   nalloc, ngent;
  int                   attentdim, j, idx;
  int                   ndup, max_nrf, max_nfe;
  int                   oppfid, oppeid, nf2, ne2;
  int                   nnodatt, ncellatt;
  int                   vid, eid, fid, rid;
  double		vxyz[3], rval=0.0;
  void                 *pval=NULL;
  FILE		        *fp;

  int                   pid = 0, numprocs = 1;
#ifdef MSTK_HAVE_MPI
  MPI_Comm_size(comm,&numprocs);
  MPI_Comm_rank(comm,&pid);
  pid += 1;  /* FLAG X3D counts processor IDs from 1 */
#endif

  strcpy(modfilename,filename);
  if (numprocs > 1)
    sprintf(modfilename,"%s.%05d",filename,pid);

  if (!(fp = fopen(modfilename,"w"))) {
    fprintf(stderr,"MESH_ExportToFLAGX3D: Couldn't open output file %s\n",
	    modfilename);
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
  
  vidatt_tmp = MAttrib_New(mesh,"vidatt_tmp",INT,MVERTEX);
  eidatt_tmp = MAttrib_New(mesh,"eidatt_tmp",INT,MEDGE);
  fidatt_tmp = MAttrib_New(mesh,"fidatt_tmp",INT,MFACE);
  ridatt_tmp = MAttrib_New(mesh,"ridatt_tmp",INT,MREGION);
  
  idx = 0; i = 0;
  while ((region = MESH_Next_Region(mesh,&idx)))
    if (MR_PType(region) != PGHOST)
      MEnt_Set_AttVal(region,ridatt_tmp,++i,0.0,NULL);
  
  idx = 0; i = 0;
  while ((face = MESH_Next_Face(mesh,&idx)))
    if (MF_PType(face) != PGHOST || MF_OnParBoundary(face))
      MEnt_Set_AttVal(face,fidatt_tmp,++i,0.0,NULL);
  
  idx = 0; i = 0;
  while ((edge = MESH_Next_Edge(mesh,&idx)))
    if (ME_PType(edge) != PGHOST || ME_OnParBoundary(edge))
      MEnt_Set_AttVal(edge,eidatt_tmp,++i,0.0,NULL);
  
  idx = 0; i = 0;
  while ((vertex = MESH_Next_Vertex(mesh,&idx)))
    if (MV_PType(vertex) != PGHOST || MV_OnParBoundary(vertex))
      MEnt_Set_AttVal(vertex,vidatt_tmp,++i,0.0,NULL);
  
  nvout = i;
  
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
  
  ndim = nr ? 3 : 2;
  fprintf(fp,"   %-22s %10d\n","numdim",ndim);
  
  /* Number of materials */
  
  ngent = 0;
  nalloc = 10;
  gentities = (int *) malloc(nalloc*sizeof(int));
  
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
	  gentities = (int *) realloc(gentities,nalloc*sizeof(int));
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
	  gentities = (int *) realloc(gentities,nalloc*sizeof(int));
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
#ifdef MSTK_HAVE_MPI
  /* make sure each processor sees all materials even if it does not
   * contain elements of that material. Also, in the concatenated
   * list, we may end up with maxngent*numprocs entries if every rank
   * has a unique set of materials */

  int maxngent;
  MPI_Allreduce(&ngent, &maxngent, 1, MPI_INT, MPI_MAX, comm);

  if (nalloc < maxngent) {
    nalloc = maxngent;
    gentities = (int *) realloc(gentities,nalloc*sizeof(int));
  }
  for (i = ngent; i < maxngent; i++) gentities[i] = -1;

  int *allgentities = (int *) calloc(maxngent*numprocs,sizeof(int));
  MPI_Allgather(gentities, maxngent, MPI_INT, allgentities, maxngent, MPI_INT,
                comm);

  
  nalloc = maxngent*numprocs;
  gentities = (int *) realloc(gentities,nalloc*sizeof(int));
  for (i = 0; i < nalloc; i++) gentities[i] = 0;
  ngent = 0;
  for (i = 0; i < maxngent*numprocs; i++) {
    gentid = allgentities[i];
    if (gentid == -1) continue;

    found = 0;
    for (j = 0; j < ngent; j++)
      if (gentities[j] == gentid) {
        found = 1;
        break;
      }      
    if (!found) {
      gentities[ngent] = gentid;
      ngent++;
    }
  }
  qsort(gentities, ngent, sizeof(int), compareINT);
#endif
  
  
  fprintf(fp,"   %-22s %10d\n","materials",ngent);
  
  /* Number of nodes */
  
  fprintf(fp,"   %-22s %10d\n","nodes",nvout);
  
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
    
    int nftot = 0;
    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx)))
      if (MF_PType(face) != PGHOST || MF_OnParBoundary(face))
        nftot++;
        
    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {      
      int fid;
      MEnt_Get_AttVal(face,fidatt_tmp,&fid,&rval,&pval);

      if (MF_OnParBoundary(face)) {
        /* face is on partition boundary - the duplicate face ID has
           to be the local ID of the face from the adjacent
           partition. Set the oppatt value to the local face ID for
           now.  After all faces are processed, communicate the
           local ID to the opposite processor and update the face */
        
        MEnt_Set_AttVal(face,oppatt,fid,0,NULL);
        MEnt_Set_AttVal(face,opppidatt,pid,0,NULL);
      }
      else if (MF_PType(face) != PGHOST) {
        fregs = MF_Regions(face);
        if (fregs && List_Num_Entries(fregs) == 2) {
          /* face is interior to this partition - we will pretend
             there are duplicate faces at the end of the regular face list */
          
          int dupfaceid = ++nftot;
          MEnt_Set_AttVal(face,oppatt,dupfaceid,0,NULL);          
          MEnt_Set_AttVal(face,opppidatt,pid,0,NULL);
        }
        else
          MEnt_Set_AttVal(face,fidatt_tmp,fid,0,NULL);
        if (fregs)
          List_Delete(fregs);
      }
    }
    
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
    
    int netot = 0;
    idx = 0;
    while ((edge = MESH_Next_Edge(mesh,&idx)))
      if (ME_PType(edge) != PGHOST || ME_OnParBoundary(edge))
        netot++;
    
    idx = 0;
    while ((edge = MESH_Next_Edge(mesh,&idx))) {      
      int eid=0;
      MEnt_Get_AttVal(edge,eidatt_tmp,&eid,&rval,&pval);
      
      if (ME_OnParBoundary(edge)) {
        /* edge is on partition boundary - the duplicate edge ID has
           to be the local ID of the edge from the adjacent
           partition. Set the oppatt value to the local edge ID for
           now.  After all edges are processed, communicate the
           local ID to the opposite processor and update the edge */
      
        MEnt_Set_AttVal(edge,oppatt,eid,0,NULL);
        MEnt_Set_AttVal(edge,opppidatt,pid,0,NULL);
      }
      else if (ME_PType(edge) != PGHOST) {
        efaces = ME_Faces(edge);
        if (efaces && List_Num_Entries(efaces) == 2) {
	  /* edge is interior to this partition - we will pretend
	     there are duplicate edges at the end of the regular face list */
          
	  int dupedgeid = ++netot;
	  MEnt_Set_AttVal(edge,oppatt,dupedgeid,0,NULL);          
          MEnt_Set_AttVal(edge,opppidatt,pid,0,NULL);
	}
        else 
          MEnt_Set_AttVal(edge,eidatt_tmp,eid,0,NULL);
        if (efaces)
          List_Delete(efaces);
      }
    }
    
    fprintf(fp,"   %-22s %10d\n","faces",netot);
    
    int nftot = 0;
    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx)))
      if (MF_PType(face) != PGHOST)
	nftot++;
    
    fprintf(fp,"   %-22s %10d\n","elements",nftot);
  }
  
  /* Exchange the local edge/face IDs between master and slave
     edges/faces. After this is done, oppatt and opppidatt on
     partition boundary faces will have the correct information from
     the opposite processors */
  
#ifdef MSTK_HAVE_MPI
  MESH_XchngEdgeFaceAttrib(mesh,oppatt,comm);
  MESH_XchngEdgeFaceAttrib(mesh,opppidatt,comm);
#endif
  
  /* Make a list of vertices on partition boundaries and determine
   * their master node IDs - in doing so we have to account for the
   * fact that some vertices may be entirely in the ghost layer and we
   * have to exclude those */

  MAttrib_ptr vmasteratt = MAttrib_New(mesh,"vmasteratt",INT,MVERTEX);

  List_ptr prtn_bndry_verts = List_New(10);
  idx = 0;
  while ((vertex = MESH_Next_Vertex(mesh,&idx))) {
    if (MV_PType(vertex) == PINTERIOR) continue;
    
    if (MV_OnParBoundary(vertex)) { /* vertex on a partition boundary */
      List_Add(prtn_bndry_verts,vertex);
      if (MV_PType(vertex) == POVERLAP) { /* vertex is master on prtn bndry */
        MEnt_Get_AttVal(vertex,vidatt_tmp,&vid,&rval,&pval);
        MEnt_Set_AttVal(vertex,vmasteratt,vid,0.0,NULL);
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
    
    cellatts = (MAttrib_ptr *) malloc(nmeshatt*sizeof(MAttrib_ptr));
    nodatts = (MAttrib_ptr *) malloc(nmeshatt*sizeof(MAttrib_ptr));
    
    for (i = 0; i < nmeshatt; i++) {
      attrib = MESH_Attrib(mesh,i);

      /* No need to write out the temporary array we created in this routine */
      if (attrib == oppatt || attrib == vidatt_tmp || attrib == eidatt_tmp || attrib == fidatt_tmp || attrib == ridatt_tmp)
        continue;
      
      atttype = MAttrib_Get_Type(attrib);      
      if (atttype != INT && atttype != DOUBLE && atttype != VECTOR && atttype != TENSOR)
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
      if ((attentdim == MVERTEX) || (attentdim == MALLTYPE)) {
        if (atttype == VECTOR) {
          int ncomps = MAttrib_Get_NumComps(attrib);
          if (ncomps == ndim)
            nodatts[nnodatt++] = attrib;
        }
      }
      
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

  /* Writing out non-standard node data is causing problems in FLAG */
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
    /* From the X3D manual - "The matid/end_matid block was meant to
     * specify material ids for each element as given in the matnames
     * entry of the header data block. However, due to a
     * long-established (and hence essentially uncorrectable) error
     * this block actually contains the material names. This error is
     * why the names must be integers even though the data header
     * block identifies them as character strings" */

    sprintf(matname,"%-d",i+1);
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
    if (MV_PType(vertex) != PGHOST || MV_OnParBoundary(vertex)) {
      MV_Coords(vertex,vxyz);
      MEnt_Get_AttVal(vertex,vidatt_tmp,&vid,&rval,&pval);
      fprintf(fp,"% 10d % 22.14E % 22.14E % 22.14E\n",vid,
              vxyz[0],vxyz[1],vxyz[2]);
    }
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
      if (MF_PType(face) == PGHOST && !MF_OnParBoundary(face)) continue;
      
      fregs = MF_Regions(face);
      nfr = fregs ? List_Num_Entries(fregs) : 0;

#ifdef DEBUG
      if (nfr == 0) {
        MSTK_Report("MESH_ExportToFLAGX3D",
                    "Meshes of non-manifold objects not supported in FLAG X3D\n",
                    MSTK_ERROR);
        return -1;
      }
#endif

      MEnt_Get_AttVal(face,fidatt_tmp,&fid,&rval,&pval); 
      fprintf(fp,"% 10d",fid);
      k = 1;

      if (nfr == 1 || MF_OnParBoundary(face)) {
        /* get face verts in outward pointing dir */

        MRegion_ptr freg0 = NULL;
        int idx2 = 0;
        MRegion_ptr freg;
        while ((freg = List_Next_Entry(fregs, &idx2)))
          if (MR_PType(freg) != PGHOST)
            freg0 = freg;
        int rfdir = MR_FaceDir(freg0,face);

        fverts = MF_Vertices(face,rfdir,0);

      } else {
        /* get face vertices in the natural direction */
        fverts = MF_Vertices(face,1,0);
      }
      List_Delete(fregs);

      nfv = List_Num_Entries(fverts);
	
      fprintf(fp,"% 10d",nfv);
      k++;
      
      for (jv = 0; jv < nfv; jv++) {
        vertex = List_Entry(fverts,jv);
        MEnt_Get_AttVal(vertex,vidatt_tmp,&vid,&rval,&pval);
        fprintf(fp,"% 10d",vid);
        if ((++k)%13 == 0)
          fprintf(fp,"\n");
      }
      List_Delete(fverts);


      /* ID of the processor owning this face */

      fprintf(fp,"% 10d",pid);   
      if ((++k)%13 == 0)
        fprintf(fp,"\n");
	

      int oppfid = 0, oppfpid = 0;  /* default for outer boundary faces */
      if (nfr == 2) {
        if (MF_OnParBoundary(face))
          /* processor owning opposite face */
          MEnt_Get_AttVal(face,opppidatt,&oppfpid,&rval,&pval);
        else
          oppfpid = pid;

        /* ID of the duplicate face */
        MEnt_Get_AttVal(face,oppatt,&oppfid,&rval,&pval);
#ifdef DEBUG
        if (!oppfid)
          fprintf(stderr,"Non-boundary face has no duplicate face info?\n");
#endif
      } else {
       if (MF_OnParBoundary(face)) {
          MEnt_Get_AttVal(face,opppidatt,&oppfpid,&rval,&pval);
          MEnt_Get_AttVal(face,oppatt,&oppfid,&rval,&pval);
       }
      }

      fprintf(fp,"% 10d",oppfpid); /* processor owning duplicate face */
      if ((++k)%13 == 0)
        fprintf(fp,"\n");
      
      fprintf(fp,"% 10d",oppfid); /* duplicate face ID (local ID on that proc) */         
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

    
    /* Now write out duplicate faces */
    
    jf = 0;
    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {
      if (MF_PType(face) == PGHOST && !MF_OnParBoundary(face)) continue;

      MEnt_Get_AttVal(face,oppatt,&oppfid,&rval,&pval);
      if (!oppfid) continue;  /* boundary face */

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
        MEnt_Get_AttVal(vertex,vidatt_tmp,&vid,&rval,&pval);
        fprintf(fp,"% 10d",vid);
        if ((++k)%13 == 0)
          fprintf(fp,"\n");
      }
      
      /* ID of the processor owning this face */
      fprintf(fp,"% 10d",oppfpid);
      if ((++k)%13 == 0)
        fprintf(fp,"\n");
      
      /* ID of processor owning the duplicate (in this case, original) face */
      
      fprintf(fp,"% 10d",pid);
      if ((++k)%13 == 0)
        fprintf(fp,"\n");
      
      /* ID of the duplicate (in this case, original) face */
      
      MEnt_Get_AttVal(face,fidatt_tmp,&fid,&rval,&pval); 
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

      List_Delete(fverts);
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
      if (ME_PType(edge) == PGHOST && !ME_OnParBoundary(edge)) continue;
      
      efaces = ME_Faces(edge);
      nef = efaces ? List_Num_Entries(efaces) : 0;
      
#ifdef DEBUG
      if (nef == 0 || nef > 2) {
        MSTK_Report("MESH_ExportToFLAGX3D",
                    "Meshes of non-manifold objects not supported in FLAG X3D",
                    MSTK_ERROR);
        return -1;
      }
#endif


      MEnt_Get_AttVal(edge,eidatt_tmp,&eid,&rval,&pval);
      fprintf(fp,"% 10d",eid);
      k = 1;

      fprintf(fp,"% 10d",2);  /* Number of vertices - always 2 for edge */

      if (nef == 1 || ME_OnParBoundary(edge)) {
        MEdge_ptr eface0 = NULL;
        int idx2 = 0;
        MFace_ptr eface;
        while ((eface = List_Next_Entry(efaces, &idx2)))
          if (MF_PType(eface) != PGHOST)
            eface0 = eface;
        
        int fedir = MF_EdgeDir(eface0,edge);
        for (jv = 0; jv < 2; jv++) {
          vertex = ME_Vertex(edge,!(jv^fedir));
          MEnt_Get_AttVal(vertex,vidatt_tmp,&vid,&rval,&pval);
          fprintf(fp,"% 10d",vid);
        }
      } else {
        for (jv = 0; jv < 2; jv++) {
          vertex = ME_Vertex(edge,jv);
          MEnt_Get_AttVal(vertex,vidatt_tmp,&vid,&rval,&pval);
          fprintf(fp,"% 10d",vid);
        }
      }
      List_Delete(efaces);

      
      /* ID of the processor owning this edge */
	
      fprintf(fp,"% 10d",pid);
	
      
      int oppeid = 0, oppepid = 0;  /* default for outer boundary edges */
      if (nef == 2) {
        if (ME_OnParBoundary(edge))
          /* processor owning the opposite edge */
          MEnt_Get_AttVal(edge,opppidatt,&oppepid,&rval,&pval);
        else
          oppepid = pid;

        /* ID of the duplicate edge */
        MEnt_Get_AttVal(edge,oppatt,&oppeid,&rval,&pval);
#ifdef DEBUG
        if (!oppeid)
          fprintf(stderr,"Non-boundary edge has no duplicate edge info?\n");
#endif
      } else {
        if (ME_OnParBoundary(edge)) {
          MEnt_Get_AttVal(edge,opppidatt,&oppepid,&rval,&pval);
          MEnt_Get_AttVal(edge,oppatt,&oppeid,&rval,&pval);
        }
      }

      fprintf(fp,"% 10d",oppepid);
      fprintf(fp,"% 10d",oppeid);

      /* Five dummy arguments required by FLAG X3D format specification */
	
      for (i = 0; i < 5; i++) fprintf(fp,"% 10d",1);
      fprintf(fp,"\n");
    }

    /* Now write out duplicate edges */
    
    idx = 0;
    while ((edge = MESH_Next_Edge(mesh,&idx))) {
      if (ME_PType(edge) == PGHOST && !ME_OnParBoundary(edge)) continue;
      
      MEnt_Get_AttVal(edge,oppatt,&oppeid,&rval,&pval);
      if (!oppeid) continue;  /* boundary edge */

      int oppepid;
      MEnt_Get_AttVal(edge,opppidatt,&oppepid,&rval,&pval);
      if (oppepid != pid) continue; /* duplicate edge from another processor */

      fprintf(fp,"% 10d",oppeid);
      
      fprintf(fp,"% 10d",2);  /* Number of vertices - always 2 for edges */

      for (jv = 0; jv < 2; jv++) {
        vertex = ME_Vertex(edge,!jv);
        MEnt_Get_AttVal(vertex,vidatt_tmp,&vid,&rval,&pval);
        fprintf(fp,"% 10d",vid);
      }
      
      /* ID of the processor owning this edge */
      fprintf(fp,"% 10d",oppepid);
      
      /* ID of the processor owning the duplicate edge */
      
      fprintf(fp,"% 10d",pid);
      
      /* ID of the duplicate (in this case, original) face */
      
      MEnt_Get_AttVal(edge,eidatt_tmp,&eid,&rval,&pval);
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

      MEnt_Get_AttVal(region,ridatt_tmp,&rid,&rval,&pval);
      fprintf(fp,"% 10d",rid);
      k = 1;

      rfaces = MR_Faces(region);
      nrf = List_Num_Entries(rfaces);
      
      fprintf(fp,"% 10d",nrf);
      k++;

      for (jf = 0; jf < nrf; jf++) {
        face = List_Entry(rfaces,jf);

        fregs = MF_Regions(face);
        if (fregs && List_Num_Entries(fregs) == 2 && !MF_OnParBoundary(face)) {
          /* Interior face */
          /* Write out the face id, if the region uses the face in the
             +ve sense (the face normal points out of the region); write
             out the duplicate face id, if not */

          int rfdir = (MR_FaceDir(region,face) == 1) ? 1 : 0;
          
          if (rfdir == 1) {
            MEnt_Get_AttVal(face,fidatt_tmp,&fid,&rval,&pval);
            fprintf(fp,"% 10d",fid);
          }
          else {
            MEnt_Get_AttVal(face,oppatt,&oppfid,&rval,&pval);
            fprintf(fp,"% 10d",oppfid);
          }
        }
        else {
          /* Boundary face or face on partition boundary */
	  
          MEnt_Get_AttVal(face,fidatt_tmp,&fid,&rval,&pval);
          fprintf(fp,"% 10d",fid);
        }
        if (fregs) List_Delete(fregs);

        k++;
        if (k%14 == 0 && jf != nrf-1)
          fprintf(fp,"\n");
	  
      }

      List_Delete(rfaces);
      fprintf(fp,"\n");
    }

    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {
      MEnt_Rem_AttVal(face,oppatt);
      MEnt_Rem_AttVal(face,opppidatt);
    }
    MAttrib_Delete(oppatt);
    MAttrib_Delete(opppidatt);
  }
  else {

    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {
      if (MF_PType(face) == PGHOST) continue;

      MEnt_Get_AttVal(face,fidatt_tmp,&fid,&rval,&pval);
      fprintf(fp,"% 10d",fid);
      k = 1;

      fedges = MF_Edges(face,1,0);
      nfe = List_Num_Entries(fedges);

      fprintf(fp,"% 10d",nfe);
      k++;

      for (je = 0; je < nfe; je++) {
        edge = List_Entry(fedges,je);

        efaces = ME_Faces(edge);
        if (efaces && List_Num_Entries(efaces) == 2 && !ME_OnParBoundary(edge)) {
          /* Interior edge */
          /* Write out the edge id, if the face uses the edge in the
             +ve sense; write out the duplicate edge id, if not */
          
          if (MF_EdgeDir_i(face,je)) {
            MEnt_Get_AttVal(edge,eidatt_tmp,&eid,&rval,&pval);
            fprintf(fp,"% 10d",eid);
          }
          else {
            MEnt_Get_AttVal(edge,oppatt,&oppeid,&rval,&pval);
            fprintf(fp,"% 10d",oppeid);
          }
        }
        else {
          /* Boundary edge or edge on partition boundary */

          MEnt_Get_AttVal(edge,eidatt_tmp,&eid,&rval,&pval);
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
    while ((edge = MESH_Next_Edge(mesh,&idx))) {
      MEnt_Rem_AttVal(edge,oppatt);
      MEnt_Rem_AttVal(edge,opppidatt);
    }
    MAttrib_Delete(oppatt);
    MAttrib_Delete(opppidatt);

  }

  fprintf(fp,"end_cells\n");


  /***********************************************************************/
  /* Slaved Nodes Data Block - not written out                           */
  /***********************************************************************/

  fprintf(fp,"%-12s% 10d\n","slaved_nodes",0);
  fprintf(fp,"end_slaved_nodes\n");

  /***********************************************************************/
  /* Ghost Nodes Data Block - not written out                            */
  /***********************************************************************/

  fprintf(fp,"ghost_nodes %10d\n",List_Num_Entries(prtn_bndry_verts));
  idx = 0;
  while ((vertex = List_Next_Entry(prtn_bndry_verts,&idx))) {
    MEnt_Get_AttVal(vertex,vidatt_tmp,&vid,&rval,&pval);
    int masterpid = MV_MasterParID(vertex)+1;
    if (masterpid == pid) {/* This is the master node */
      fprintf(fp,"% 10d% 10d% 10d% 10d\n",vid,pid,vid,1);
    }
    else {
      int mastervid;
      MEnt_Get_AttVal(vertex,vmasteratt,&mastervid,&rval,&pval);
      fprintf(fp,"% 10d% 10d% 10d% 10d\n",vid,masterpid,mastervid,1);
    }
  }
  fprintf(fp,"end_ghost_nodes\n");


  /***********************************************************************/
  /* Cell centered data                                                  */
  /***********************************************************************/

  fprintf(fp,"cell_data\n");

  /***********************************************************************/
  /* Material ID data block                                              */
  /*                                                                     */
  /* From the X3D manual - "The matid/end_matid block was meant to
    specify material ids for each element as given in the matnames
    entry of the header data block. However, due to a long-established
    (and hence essentially uncorrectable) error this block actually
    contains the material names. This error is why the names must be
    integers even though the data header block identifies them as
    character strings"                                                   */
  /***********************************************************************/

  fprintf(fp,"matid\n");
  if (nr) {
    k = 0;
    idx = 0;
    while ((region = MESH_Next_Region(mesh,&idx))) {
      if (MR_PType(region) == PGHOST) continue;
      gentid = MR_GEntID(region);
      if (!gentid) gentid = 1;
      fprintf(fp,"% 10d",gentid);
      k++;
      if (k%10 == 0)
        fprintf(fp,"\n");
    }
  }
  else {
    k = 0;
    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {
      if (MF_PType(face) == PGHOST) continue;
      gentid = MF_GEntID(face);
      if (!gentid) gentid = 1;
      fprintf(fp,"% 10d",gentid);
      k++;
      if (k%10 == 0)
        fprintf(fp,"\n");
    }
  }
  if (k%10)
    fprintf(fp,"\n");
  fprintf(fp,"end_matid\n");
 
  /***********************************************************************/
  /* Partition/Processor Data for cells - Not used but must be present   */
  /***********************************************************************/

  fprintf(fp,"partelm\n");
  k = 0;
  if (nr) {
    idx = 0; 
    while ((region = MESH_Next_Region(mesh,&idx))) 
      if (MR_PType(region) != PGHOST) {
        fprintf(fp,"% 10d",pid);
        k++;
        if (k%10 == 0)
          fprintf(fp,"\n");
      }
  }
  else {
    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx)))
      if (MF_PType(face) != PGHOST) {
        fprintf(fp,"% 10d",pid);
        k++;
        if (k%10 == 0)
          fprintf(fp,"\n");
      }

  }
  if (k%10)
    fprintf(fp,"\n");
  fprintf(fp,"end_partelm\n");


  /* Write out any other cell based attributes, if requested */

  for (i = 0; i < ncellatt; i++) {
    attrib = cellatts[i];
    atttype = MAttrib_Get_Type(attrib);
    if (atttype == VECTOR || atttype == TENSOR) {
      MSTK_Report("MESH_ExportToFLAGX3D", "Export of cell vectors and tensors to X3D file not implemented", MSTK_WARN);
      continue;
    }

    MAttrib_Get_Name(attrib,attname);
    if (strcmp("_tmp",attname-4) == 0)
      continue;

    fprintf(fp,"%s\n",attname);
    
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

  if (cellatts) free(cellatts);
  
  fprintf(fp,"end_cell_data\n");
  
  
  /***********************************************************************/
  /* Node centered data                                                  */
  /***********************************************************************/

  fprintf(fp,"node_data\n");

  for (i = 0; i < nnodatt; i++) {
    attrib = nodatts[i];
    atttype = MAttrib_Get_Type(attrib);
    if (atttype != VECTOR) continue;  /* X3D only supports vector attributes for nodes */
    int ncomps = MAttrib_Get_NumComps(attrib);

    MAttrib_Get_Name(attrib,attname);
    if (strcmp("_tmp",attname-4) == 0)
      continue;

    fprintf(fp,"%s\n",attname);
    
    idx = 0;
    while ((vertex = MESH_Next_Vertex(mesh,&idx))) {
      if (MV_PType(vertex) == PGHOST && !MV_OnParBoundary(vertex)) continue;

      double *vval;
      MEnt_Get_AttVal(vertex,attrib,&ival,&rval,&vval);
      
      for (k = 0; k < ncomps; k++)
        fprintf(fp,"% 20.12E", vval[k]);
      for (k = ncomps; k < 3; k++)
        fprintf(fp,"% 20.12E", 0.0);  /* expects 3 components for nodal vecs */
      fprintf(fp,"\n");
    }
    
    strcpy(tmpstr,"end_");
    strcat(tmpstr,attname);
    fprintf(fp,"%s\n",tmpstr);
  }
  
  if (nodatts) free(nodatts);
  
  fprintf(fp,"end_node_data\n");

  /* Finish Export */
  
  fprintf(fp,"end_dump\n");

  fclose(fp);


  /* Write out elements of each material (derived from each element
   * block) as a .reg file */

  int g;
  for (g = 0; g < ngent; ++g) {
    char regfilename[256];
    sprintf(regfilename, "mat.%d.Reg",g+1);
    if (numprocs > 1) {
      char ext[256];
      sprintf(ext, ".%05d",pid);
      strcat(regfilename, ext);
    }
    fp = fopen(regfilename, "w");
    if (nr) {
      idx = 0;
      while ((region = MESH_Next_Region(mesh, &idx))) {
        if (MR_PType(region) != PGHOST && MR_GEntID(region) == gentities[g])
          fprintf(fp, "% 10d\n", MR_ID(region));
      }
      fprintf(fp,"\n");

      idx = 0;
      while ((region = MESH_Next_Region(mesh, &idx))) {
        if (MR_PType(region) != PGHOST && MR_GEntID(region) == gentities[g])
          fprintf(fp, "% 20.18e\n", 1.0);
      }
      
    } else {
      idx = 0;
      while ((face = MESH_Next_Face(mesh, &idx))) {
        if (MF_PType(face) != PGHOST && MF_GEntID(face) == gentities[g])
          fprintf(fp, "% 10d\n", MF_ID(face));
      }
      fprintf(fp,"\n");
      

      idx = 0;
      while ((face = MESH_Next_Face(mesh, &idx))) {
        if (MF_PType(face) != PGHOST && MF_GEntID(face) == gentities[g])
          fprintf(fp, "% 20.18e\n", 1.0);
      }
      
    }
    fclose(fp);
  }
  
  /* Write out each element set as a .reg file; Some or all of these
     will be a duplicate of the mat.N.Reg files */

  MSet_ptr mset;
  idx = 0;
  while ((mset = MESH_Next_MSet(mesh, &idx))) {
    if ((nr > 0 && MSet_EntDim(mset) != MREGION) ||
        (nr == 0 && nf > 0 && MSet_EntDim(mset) != MFACE)) continue;

    char setname[256];
    MSet_Name(mset, setname);

    char regfilename[256];
    strcpy(regfilename, setname);
    strcat(setname, ".Reg");
    if (numprocs > 1) {
      char ext[256];
      sprintf(ext, ".%05d",pid);
      strcat(regfilename, ext);
    }

    fp = fopen(regfilename, "w");
    
    int idx2 = 0;
    MEntity_ptr ment;
    while ((ment = MSet_Next_Entry(mset, &idx2)))
      if (MEnt_PType(ment) != PGHOST)
        fprintf(fp, "% 10d\n", MEnt_ID(ment));
    fprintf(fp,"\n");
    idx2 = 0;
    while ((ment = MSet_Next_Entry(mset, &idx2)))
      if (MEnt_PType(ment) != PGHOST)
        fprintf(fp, "% 20.18e\n",1.0);

    fclose(fp);
  }


  /* Write out each nodeset as a .bdy file */

  idx = 0;
  while ((mset = MESH_Next_MSet(mesh, &idx))) {
    if (MSet_EntDim(mset) != MVERTEX) continue;

    char setname[256];
    MSet_Name(mset, setname);

    char bdyfilename[256];
    strcpy(bdyfilename, setname);
    strcat(bdyfilename, ".bdy");
    if (numprocs > 1) {
      char ext[256];
      sprintf(ext, ".%05d",pid);
      strcat(bdyfilename, ext);
    }

    fp = fopen(bdyfilename, "w");
    
    int idx2 = 0;
    MEntity_ptr ment;
    while ((ment = MSet_Next_Entry(mset, &idx2)))
      if (MEnt_PType(ment) != PGHOST || MEnt_OnParBoundary(ment))
        fprintf(fp, "% 10d\n", MEnt_ID(ment));

    fclose(fp);
  }


  
  
  /* Clean up */
  
  free(gentities);
  
  idx = 0;
  while ((vertex = MESH_Next_Vertex(mesh,&idx)))
    MEnt_Rem_AttVal(vertex,vidatt_tmp);
  
  idx = 0;
  while ((edge = MESH_Next_Edge(mesh,&idx)))
    MEnt_Rem_AttVal(edge,eidatt_tmp);
  
  idx = 0;
  while ((face = MESH_Next_Face(mesh,&idx)))
    MEnt_Rem_AttVal(face,fidatt_tmp);
  
  idx = 0;
  while ((region = MESH_Next_Region(mesh,&idx)))
    MEnt_Rem_AttVal(region,ridatt_tmp);

  MAttrib_Delete(vidatt_tmp);
  MAttrib_Delete(eidatt_tmp);
  MAttrib_Delete(fidatt_tmp);
  MAttrib_Delete(ridatt_tmp);
  
  return 1;
}


#ifdef __cplusplus
  }
#endif


