#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif



  /* Function to export MSTK mesh to GMV format */

  /* if natt = 0, all attributes are written out 
     if natt = -1, no attributes are written out 
     if natt > 0, only attributes specified in attnames are written out */

  /* opts is an array of flags that controls how mesh is exported to GMV 
     
     if opts[0] = 0, nodes are written out with the 'nodev' keyword
     and as triplets of (x,y,z) coordinates. If it is 1, nodes are
     written out with the 'nodes' keyword and all the x coordinates
     are written out, followed by all the y coordinates and then the z
     coordinates 

     if opts[1] = 0, cells are written out with the 'cells' keyword
     and are listed normally as tets, hexes, tris etc and point to
     node ids if opts[1] = 1, cells are written out as 'vface3d' or
     'vface2d' which point to 'vfaces' which in turn point to node ids
  */

  /* If you want default options, call the functions as

     MESH_ExportToGMV(mymesh,mygmvfile,0,NULL,NULL);
  */


int MESH_ExportToGMV(Mesh_ptr mesh, const char *filename, const int natt, 
		     const char **attnames, const int *opts, MSTK_Comm comm) {

  int			gentid, *gentities;
  MFType		ftype;
  MRType                rtype;
  List_ptr		rverts, fverts, rfaces, fregs, efaces;
  List_ptr              vlist=NULL, elist=NULL, flist=NULL, rlist=NULL;
  MVertex_ptr           vertex;
  MEdge_ptr             edge;
  MFace_ptr	        face;
  MRegion_ptr           region, region0, region1, oppreg;
  MAttrib_ptr           attrib, *outattribs, oppatt;
  MAttType              atttype;
  char                  attname[256], matname[256], date_str[256];
  int                   jv, jf, gmodel, nrf, nrv, nfv, dir;
  int			i, found, k, nmeshatt, noutatt, ival, idgaps;
  int                   nalloc, ngent, fnum, fid, vid, icr;
  int                   nv, ne, nf, nr;
  int                   minid, maxid, tempid, prev_tempid;
  int                   attentdim, j, ncells, polygons=0, idx;
  int                   ncells1=0, ncells2=0, ncells3=0;
  int                   ndup, NFACES, rid0, rid1, rid;
  int                   oppfid, opprid;
  int                   sorted;
  double		vxyz[3], rval, *xcoord, *ycoord, *zcoord;
  void                 *pval;
  FILE		        *fp;
  time_t                ctime;

  int			rtmpl[5][8] =   {{0,2,1,3,-1,-1,-1,-1},
					 {0,3,2,1,4,-1,-1,-1},
					 {3,4,5,0,1,2,-1,-1},
					 {-1,-1,-1,-1,-1,-1,-1,-1},
					 {4,5,6,7,0,1,2,3}};
  MAttrib_ptr    vidatt=0,eidatt=0,fidatt=0,ridatt=0;
  char funcname[256] = "MESH_ExportToGMV", mesg[256];

  gmodel = 1;

  char modfilename[256];

  int rank=0, numprocs=1;
#ifdef MSTK_HAVE_MPI
  if (comm != 0 && comm != MPI_COMM_NULL) {
    MPI_Comm_size(comm,&numprocs);
    MPI_Comm_rank(comm,&rank);
  }
  
  if (numprocs > 1) 
    sprintf(modfilename,"%s.%05d",filename,rank);
  else
    strcpy(modfilename,filename);
#else
  strcpy(modfilename,filename);
#endif

  if (!(fp = fopen(modfilename,"w"))) {
    sprintf(mesg,"Cannot open file %s for writing",modfilename);
    MSTK_Report(funcname,mesg,MSTK_FATAL);
  }


  fprintf(fp,"gmvinput ascii\n");
  nv = MESH_Num_Vertices(mesh);

  fprintf(fp,"codename MSTK_V_1.4\n");
  ctime = time(&ctime);
  strftime(date_str,sizeof(date_str),"%m/%d/%Y",localtime(&ctime));
  fprintf(fp,"simdate %s\n",date_str);


  /* Check if there are gaps in the IDs of entities - if there are
     then we cannot use the IDs directly to reference entities */

  idgaps = 0;

  if (nv) {
    maxid = 0;
    minid = 1e+9;
    idx = 0;
    while ((vertex = MESH_Next_Vertex(mesh,&idx))) {
      int id = MV_ID(vertex);
      minid = (id < minid) ? id : minid;
      maxid = (id > maxid) ? id : maxid;
    }
    if (minid != 1 || maxid != nv)
      idgaps = 1;
  }
  ne = MESH_Num_Edges(mesh);
  if (ne && !idgaps) {
    maxid = 0;
    minid = 1e+9;
    idx = 0;
    while ((edge = MESH_Next_Edge(mesh,&idx))) {
      int id = ME_ID(edge);
      minid = (id < minid) ? id : minid;
      maxid = (id > maxid) ? id : maxid;
    }
    if (minid != 1 || maxid != ne)
      idgaps = 1;
   }
  nf = MESH_Num_Faces(mesh);
  if (ne && !idgaps) {
    maxid = 0;
    minid = 1e+9;
    idx = 0;
    while ((face = MESH_Next_Face(mesh,&idx))) {
      int id = MF_ID(face);
      minid = (id < minid) ? id : minid;
      maxid = (id > maxid) ? id : maxid;
    }
    if (minid != 1 || maxid != nf)
      idgaps = 1;
   }
  nr = MESH_Num_Regions(mesh);
  if (nr && !idgaps) {
    maxid = 0;
    minid = 1e+9;
    idx = 0;
    while ((region = MESH_Next_Region(mesh,&idx))) {
      int id = MR_ID(region);
      minid = (id < minid) ? id : minid;
      maxid = (id > maxid) ? id : maxid;
    }
    if (minid != 1 || maxid != nr)
      idgaps = 1;
   }
  

  /* Now build up the vertex list sorted by ID and write them out */

  vidatt = MESH_AttribByName(mesh,"vidatt");
  if (!vidatt)
    vidatt = MAttrib_New(mesh,"vidatt",INT,MVERTEX);

  vlist = List_New(nv);
  idx = 0; i = 1; 
  prev_tempid = 0; sorted = 1;
  while ((vertex = MESH_Next_Vertex(mesh,&idx))) {
    List_Add(vlist,vertex);

    tempid = idgaps ? i : MV_ID(vertex);
    MEnt_Set_AttVal(vertex,vidatt,tempid,0.0,NULL);

    if (tempid != prev_tempid+1) sorted = 0;
    prev_tempid = tempid;

    i++;
  }

  if (!sorted) /* Sort the vertices according to the IDs */
    List_Sort(vlist,nv,sizeof(MVertex_ptr),compareID);


  if (!opts || opts[0] == 0) {
    fprintf(fp,"nodev %d\n",nv);
    idx = 0;
    while ((vertex = List_Next_Entry(vlist,&idx))) {
      MV_Coords(vertex,vxyz);
      fprintf(fp,"% 20.12f % 20.12lf % 20.12lf\n",vxyz[0],vxyz[1],vxyz[2]);
    }
  }
  else {
    fprintf(fp,"nodes %d\n",nv);
    xcoord = (double *) malloc(nv*sizeof(double));
    ycoord = (double *) malloc(nv*sizeof(double));
    zcoord = (double *) malloc(nv*sizeof(double));
    idx = 0; jv = 0;
    while ((vertex = List_Next_Entry(vlist,&idx))) {
      MV_Coords(vertex,vxyz);
      xcoord[jv] = vxyz[0];
      ycoord[jv] = vxyz[1];
      zcoord[jv] = vxyz[2];
      jv++;
    }

    for (jv = 0; jv < nv; jv++) {
      fprintf(fp,"% 20.12lf",xcoord[jv]);
      if (jv%5 == 0 || jv == nv-1)
	fprintf(fp,"\n");
    }
    for (jv = 0; jv < nv; jv++) {
      fprintf(fp,"% 20.12lf",ycoord[jv]);
      if (jv%5 == 0 || jv == nv-1)
	fprintf(fp,"\n");
    }
    for (jv = 0; jv < nv; jv++) {
      fprintf(fp,"% 20.12lf",zcoord[jv]);
      if (jv%5 == 0 || jv == nv-1)
	fprintf(fp,"\n");
    }
    free(xcoord);
    free(ycoord);
    free(zcoord);
  }

#ifdef MSTK_USE_MARKERS
  int cellmk = MSTK_GetMarker();
#else
  MAttrib_ptr cellmkatt = MAttrib_New(mesh, "cellmk", INT, MALLTYPE);
#endif

  /* number of regions */

  ncells = 0;

  if (nr) {    
    ncells3 = nr;
    ncells += nr;

    ridatt = MESH_AttribByName(mesh,"ridatt");
    if (!ridatt)
      ridatt = MAttrib_New(mesh,"ridatt",INT,MREGION);

    rlist = List_New(nr);

    idx = 0; i = 1;
    prev_tempid = 0;
    sorted = 1;
    while ((region = MESH_Next_Region(mesh,&idx))) {
      List_Add(rlist,region);
      
      tempid = idgaps? i : MR_ID(region);
      MEnt_Set_AttVal(region,ridatt,tempid,0.0,NULL);

      if (tempid != prev_tempid+1) sorted = 0;
      prev_tempid = tempid;

      i++;
    }

    if (!sorted) /* sort by region ID */ 
      List_Sort(rlist,ncells,sizeof(MRegion_ptr),compareID);
  }

  /* May have problems if we want to write polygons as 'faces' */
  /* number of faces not connected to a region */  

  if (nf) {

    fidatt = MESH_AttribByName(mesh,"fidatt");
    if (!fidatt)
      fidatt = MAttrib_New(mesh,"fidatt",INT,MFACE);

    flist = List_New(nf);

    ncells2 = 0;
    idx = 0; i = 1;
    prev_tempid = 0;
    sorted = 1;
    while ((face = MESH_Next_Face(mesh,&idx))) {
      if ((fregs = MF_Regions(face)))
	List_Delete(fregs);
      else {
	ncells2++;
#ifdef MSTK_USE_MARKERS
	MEnt_Mark(face,cellmk);
#else
        MEnt_Set_AttVal(face, cellmkatt, 1, 0.0, NULL);
#endif
      }

      List_Add(flist,face);

      tempid = idgaps? i : MF_ID(face);
      MEnt_Set_AttVal(face,fidatt,tempid,0.0,NULL);

      if (tempid != prev_tempid+1) sorted = 0;
      prev_tempid = tempid;

      i++;
    }
    ncells += ncells2;
    
    if (!sorted)
      List_Sort(flist,nf,sizeof(MFace_ptr),compareID);
  }


  /* number of edges not connected to a face */
  if (ne) {

    elist = List_New(ne);

    eidatt = MESH_AttribByName(mesh,"eidatt");
    if (!eidatt)
      eidatt = MAttrib_New(mesh,"eidatt",INT,MEDGE);

    ncells1 = 0;
    idx = 0; i = 1;
    prev_tempid = 0;
    sorted = 1;
    while ((edge = MESH_Next_Edge(mesh,&idx))) {
      if ((efaces = ME_Faces(edge)))
	List_Delete(efaces);
      else {
	ncells1++;
#ifdef MSTK_USE_MARKERS
	MEnt_Mark(edge,cellmk);
#else
        MEnt_Set_AttVal(edge, cellmkatt, 1, 0.0, NULL);
#endif
      }

      List_Add(elist,edge);

      tempid = idgaps ? i : ME_ID(edge);
      MEnt_Set_AttVal(edge,eidatt,tempid,0.0,NULL);

      if (tempid != prev_tempid) sorted = 0;
      prev_tempid = tempid;

      i++;
    }
    ncells += ncells1;

    if (!sorted)
      List_Sort(elist,ne,sizeof(MEdge_ptr),compareID);
  }



  /* write region connectivity for entire mesh */

  ngent = 0;
  nalloc = 10;
  gentities = (int *) malloc(nalloc*sizeof(int));
  
  fprintf(fp,"cells %d\n",ncells);
    
  if (ncells3) {
    if (!opts || opts[1] == 0) { /* write out cells as cells*/
      
      idx = 0;
      while ((region = List_Next_Entry(rlist,&idx))) {

        rverts = MR_Vertices(region);
        nrv = List_Num_Entries(rverts);

        rfaces = MR_Faces(region);
        nrf = List_Num_Entries(rfaces);

        switch (nrv) {
        case 4:
          rtype = TET;
          fprintf(fp,"tet 4 ");
          break;
        case 5:
          if (nrf == 5) {
            rtype = PYRAMID;
            fprintf(fp,"pyramid 5 ");
          }
          else 
            rtype = POLYHED;
          break;
        case 6:
          if (nrf == 5) { /* need additional check to fully verify */
            rtype = PRISM;
            fprintf(fp,"prism 6 ");
          }
          else
            rtype = POLYHED;
          break;
        case 8:
          if (nrf == 6) {
            /* Really verify that its a hex - all faces must be a quad */
            /* It is possible to have a non-hex with 6 faces, 8 vertices
               and 12 edges */

            int allquad = 1;
            for (jf = 0; jf < nrf; jf++) {
              face = List_Entry(rfaces,jf);
              if (MF_Num_Vertices(face) != 4) {
                allquad = 0;
                break;
              }
            }

            if (allquad) {
              rtype = HEX;
              fprintf(fp,"hex 8 ");
            }
            else
              rtype = POLYHED;
          }
          else
            rtype = POLYHED;
          break;
        default:
          rtype = POLYHED;
          break;
        }
        if (rtype != POLYHED) {
          for (jv = 0; jv < nrv; jv++) {
            vertex = List_Entry(rverts,rtmpl[nrv-4][jv]);
            MEnt_Get_AttVal(vertex,vidatt,&vid,&rval,&pval);
            fprintf(fp," % 8d",vid);
          }
          fprintf(fp,"\n");
        }
        else {
          fprintf(fp,"general %d\n",nrf);
          for (jf = 0; jf < nrf; jf++) {
            face = List_Entry(rfaces,jf);
            fprintf(fp,"%d ",MF_Num_Edges(face)); /* assuming linear elements */
          }
          fprintf(fp,"\n");
          for (jf = 0; jf < nrf; jf++) {
            face = List_Entry(rfaces,jf);
            dir = MR_FaceDir_i(region,jf);
	  
            fverts = MF_Vertices(face,dir,0);
            nfv = List_Num_Entries(fverts);
	  
            for (jv = 0; jv < nfv; jv++) {
              vertex = List_Entry(fverts,jv);
              MEnt_Get_AttVal(vertex,vidatt,&vid,&rval,&pval);
              fprintf(fp," %8d",vid);
            }
            fprintf(fp,"\n");
            List_Delete(fverts);
          }
        }
      
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
      
        List_Delete(rfaces);
        List_Delete(rverts);
      }
    }
    else { /* Write out cells as vface3D and faces as vfaces */
      NFACES = nf;

      oppatt = MAttrib_New(mesh,"oppfaceID",INT,MFACE); 

      ndup = 0;
      idx = 0;
      while ((region = List_Next_Entry(rlist,&idx))) {      
        MEnt_Get_AttVal(region,ridatt,&rid,&rval,&pval);

        rfaces = MR_Faces(region);
        nrf = List_Num_Entries(rfaces);
      
        fprintf(fp,"vface3d % d ",nrf);
      
        for (jf = 0; jf < nrf; jf++) {
          face = List_Entry(rfaces,jf);

          /* if face has a region on the other side, we will be writing
             out a duplicate face for the opposite region. We use the
             convention that if current region ID from which we are
             looking at the face is smaller, then the face ID will be
             written out as is, whereas for the opposite region a new
             face with ID = face_id+num_mesh_faces will be written
             out. Vice versa if this region ID is higher */

          fregs = MF_Regions(face);
          if (fregs && List_Num_Entries(fregs) == 2) {
	  
            oppreg = List_Entry(fregs,0);
            if (oppreg == region)
              oppreg = List_Entry(fregs,1);
	  
            MEnt_Get_AttVal(oppreg,ridatt,&opprid,&rval,&pval);
            if (rid < opprid) {
              MEnt_Get_AttVal(face,fidatt,&fid,&rval,&pval);
              fprintf(fp," % d ",fid);
            }
            else {
              ndup++;
              fprintf(fp," % d ",NFACES+ndup);
              MEnt_Set_AttVal(face,oppatt,NFACES+ndup,0,NULL);
            }
          }
          else {
            MEnt_Get_AttVal(face,fidatt,&fid,&rval,&pval);
            fprintf(fp,"% d ",fid);
          }
        }
        fprintf(fp,"\n");

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
    

      /* Now to write out vfaces */

      fprintf(fp,"vfaces % d\n",NFACES+ndup);

      /* Write out all the faces of the mesh first */
      ndup = 0;
      idx = 0;
      while ((face = List_Next_Entry(flist,&idx))) {      

        fregs = MF_Regions(face);
        if (!fregs)	
          continue;

        if (List_Num_Entries(fregs) == 2) {
          ndup++;

          region0 = List_Entry(fregs,0);
          region1 = List_Entry(fregs,1);
	  
          /* This is the face written out with its regular ID. The
             connected region with the lower ID will use this face such
             that its normal points out of the region */

          MEnt_Get_AttVal(region0,ridatt,&rid0,&rval,&pval);
          MEnt_Get_AttVal(region1,ridatt,&rid1,&rval,&pval);
          if (rid0 < rid1) {
            region = region0;
            rid = rid0;
          }
          else {
            region = region1;
            rid = rid1;
          }

          dir = MR_FaceDir(region,face);	  
          fverts = MF_Vertices(face,dir,0);

          nfv = List_Num_Entries(fverts);
	
          MEnt_Get_AttVal(face,fidatt,&fid,&rval,&pval);
          MEnt_Get_AttVal(face,oppatt,&oppfid,&rval,&pval);
	
          fprintf(fp,"% d  1  % d 1 % d ",nfv,oppfid,rid);
          for (i = 0; i < nfv; i++) {
            MEnt_Get_AttVal(List_Entry(fverts,i),vidatt,&vid,&rval,&pval);
            fprintf(fp," % d ",vid);
          }
          fprintf(fp,"\n");
	
          List_Delete(fverts);
        }
        else {
          region = List_Entry(fregs,0);
          MEnt_Get_AttVal(region,ridatt,&rid,&rval,&pval);
          dir = MR_FaceDir(region,face);	  
          fverts = MF_Vertices(face,dir,0);

          nfv = List_Num_Entries(fverts);

          MEnt_Get_AttVal(face,fidatt,&fid,&rval,&pval);
          oppfid = 0;
	
          fprintf(fp,"% d  1  % d 1 % d ",nfv,oppfid,rid);
          for (i = 0; i < nfv; i++) {
            MEnt_Get_AttVal(List_Entry(fverts,i),vidatt,&vid,&rval,&pval);
            fprintf(fp," % d ",vid);
          }
          fprintf(fp,"\n");
	
          List_Delete(fverts);
        }
        List_Delete(fregs);
      }

      /* Write out all the faces that are shared by two regions again.
         Since we didn't store the order which the faces were assigned
         the "duplicate face IDs" (and probably we don't want to store
         this information) we need to loop over regions (cells) again,
         to duplicate the constructive process so that the faces are
         written out in the right order.
      */
    
      idx = 0; 
      ndup = 0;

      while ((region = List_Next_Entry(rlist,&idx))) {      

        MEnt_Get_AttVal(region,ridatt,&rid,&rval,&pval);

        rfaces = MR_Faces(region);
        nrf = List_Num_Entries(rfaces);
      
        for (jf = 0; jf < nrf; jf++) {

          face = List_Entry(rfaces,jf);

          fregs = MF_Regions(face);
 
          if (fregs && List_Num_Entries(fregs) == 2) {

            oppreg = List_Entry(fregs,0);
            if (oppreg == region){
              oppreg = List_Entry(fregs,1);
            }

            MEnt_Get_AttVal(oppreg,ridatt,&opprid,&rval,&pval);
            if ( rid > opprid) {
              MEnt_Get_AttVal(face,fidatt,&fid,&rval,&pval);
              dir = MR_FaceDir(region,face);	  
              fverts = MF_Vertices(face,dir,0);
              nfv = List_Num_Entries(fverts);
              fprintf(fp,"% d  1  % d 1 % d ",nfv,fid,rid);

              for (i = 0; i < nfv; i++) {
                MEnt_Get_AttVal(List_Entry(fverts,i),vidatt,&vid,&rval,&pval);
                fprintf(fp," % d ",vid);
              }
              fprintf(fp,"\n");
              List_Delete(fverts);      
            }  
          }
          List_Delete(fregs);
        }

      }

      idx = 0;
      while ((face = List_Next_Entry(flist,&idx)))
        MEnt_Rem_AttVal(face,oppatt);
      MAttrib_Delete(oppatt);
    }
  }

  
  if (ncells2) {
    if (opts && opts[1] == 1) {
      MSTK_Report(funcname,"vface2d not yet implemented",MSTK_WARN);
      return 0;
    }

    idx = 0;
    while ((face = List_Next_Entry(flist,&idx))) {
      int fmarked;
#ifdef MSTK_USE_MARKERS
      fmarked = MEnt_IsMarked(face,cellmk);
#else
      MEnt_Get_AttVal(face, cellmkatt, &fmarked, &rval, &pval);
#endif
      if (!fmarked)
        continue;
      
      if (MF_Num_Vertices(face) > 4) {
	polygons = 1;
	break;
      }
    }
    
    /* Turn off writing polygons as faces for now */
    polygons = 0;      
    if (polygons) 
      fprintf(fp,"faces %d 0\n",ncells2);
    
    idx = 0;
    while ((face = List_Next_Entry(flist,&idx))) {
      int fmarked;
#ifdef MSTK_USE_MARKERS
      fmarked = MEnt_IsMarked(face,cellmk);
#else
      MEnt_Get_AttVal(face, cellmkatt, &fmarked, &rval, &pval);
#endif
      if (!fmarked)
	continue;
      
      fnum++;
      MEnt_Get_AttVal(face,fidatt,&fid,&rval,&pval);
      
      fverts = MF_Vertices(face,1,0);
      nfv = List_Num_Entries(fverts);
      if (polygons) {
	fprintf(fp,"%-d ",nfv);
	
	for (jv = 0; jv < nfv; jv++) {
	  vertex = List_Entry(fverts,jv);
	  MEnt_Get_AttVal(vertex,vidatt,&vid,&rval,&pval);
	  fprintf(fp," % 8d",vid);
	}
	
	fprintf(fp," 0  0 ");
	
	fprintf(fp,"\n");
	List_Delete(fverts);
      }
      else {
	switch (nfv) {
	case 3:
	  ftype = TRI;
	  fprintf(fp,"tri 3 ");
	  break;
	case 4:
	  ftype = QUAD;
	    fprintf(fp,"quad 4 ");
	    break;
	  default:
	    ftype = POLYGON;
	    fprintf(fp,"general 1\n");
	    fprintf(fp,"%-d ",nfv);
	    break;
	  }
	  
	  for (jv = 0; jv < nfv; jv++) {
	    vertex = List_Entry(fverts,jv);
	    MEnt_Get_AttVal(vertex,vidatt,&vid,&rval,&pval);
	    fprintf(fp," % 8d",vid);
	  }
	  
	  fprintf(fp,"\n");
	List_Delete(fverts);
      }
	
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

  
  /* Write out stand-alone edges (not connected to any face */

  if (ncells1) {
    idx = 0;
    while ((edge = List_Next_Entry(elist,&idx))) {
      int emarked;
#ifdef MSTK_USE_MARKERS
      emarked = MEnt_IsMarked(edge,cellmk);
#else
      MEnt_Get_AttVal(edge, cellmkatt, &emarked, &rval, &pval);
#endif
      if (!emarked)
	continue;
      
      fprintf(fp,"line 2 ");

      for (jv = 0; jv < 2; jv++) {
	vertex = ME_Vertex(edge,jv);
	MEnt_Get_AttVal(vertex,vidatt,&vid,&rval,&pval);
	fprintf(fp," % 8d",vid);
      }
	  
      fprintf(fp,"\n");
	
      gentid = ME_GEntID(edge);
	  
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


  if (gmodel) {
    fprintf(fp,"material %d 0\n",ngent);
    for (i = 0; i < ngent; i++) {
      sprintf(matname,"mat%-d\n",gentities[i]);
      fprintf(fp,"%s",matname);
    }
    
    k = 0;
    if (ncells3) {
      idx = 0;
      while ((region = List_Next_Entry(rlist,&idx))) {
        gentid = MR_GEntID(region);
        
        for (i = 0, found = 0; i < ngent; i++) {
          if (gentities[i] == gentid) {
            found = 1;
            break;
          }
        }
        fprintf(fp,"%d ",(i+1));
        if ((++k)%10 == 0)
          fprintf(fp,"\n");
      }
    }

    if (ncells2) {
      idx = 0;
      while ((face = List_Next_Entry(flist,&idx))) {
        int fmarked;
#ifdef MSTK_USE_MARKERS
	fmarked = MEnt_IsMarked(face,cellmk);
#else
        MEnt_Get_AttVal(face, cellmkatt, &fmarked, &rval, &pval);
#endif
        if (!fmarked)
	  continue;

	gentid = MF_GEntID(face);
	for (i = 0, found = 0; i < ngent; i++) {
	  if (gentities[i] == gentid) {
	    found = 1;
	    break;
	  }
	}
	fprintf(fp,"%d ",(i+1));
	if ((++k)%10 == 0)
	  fprintf(fp,"\n");
      }
    }

    if (ncells1) {
      idx = 0;
      while ((edge = List_Next_Entry(elist,&idx))) {
        int emarked;
#ifdef MSTK_USE_MARKERS
	emarked = MEnt_IsMarked(edge,cellmk);
#else
        MEnt_Get_AttVal(edge, cellmkatt, &emarked, &rval, &pval);
#endif
        if (!emarked)
	  continue;

	gentid = ME_GEntID(edge);
	for (i = 0, found = 0; i < ngent; i++) {
	  if (gentities[i] == gentid) {
	    found = 1;
	    break;
	  }
	}
	fprintf(fp,"%d ",(i+1));
	if ((++k)%10 == 0)
	  fprintf(fp,"\n");
      }
    }
    fprintf(fp,"\n");


    /* Other variables related to geometric classification */
    
    fprintf(fp,"variable \n");
    fprintf(fp,"itetclr  0\n");
    
    k = 0;
    if (ncells3) {
      idx = 0;
      while ((region = List_Next_Entry(rlist,&idx))) {
        gentid = MR_GEntID(region);
        fprintf(fp,"%d ",gentid);      
        if ((++k)%10 == 0)
          fprintf(fp,"\n");
      }
    }

    if (ncells2) {
      while ((face = List_Next_Entry(flist,&idx))) {
        int fmarked;
#ifdef MSTK_USE_MARKERS
	fmarked = MEnt_IsMarked(face,cellmk);
#else
        MEnt_Get_AttVal(face, cellmkatt, &fmarked, &rval, &pval);
#endif
        if (!fmarked)
	  continue;
	else {
	  gentid = MF_GEntID(face);
	  fprintf(fp,"%d ",gentid);
	  if ((++k)%10 == 0)
	    fprintf(fp,"\n");
	}
      }
    }

    if (ncells1) {
      idx = 0;
      while ((edge = List_Next_Entry(elist,&idx))) {
        int emarked;
#ifdef MSTK_USE_MARKERS
	emarked = MEnt_IsMarked(edge,cellmk);
#else
        MEnt_Get_AttVal(edge, cellmkatt, &emarked, &rval, &pval);
#endif
        if (!emarked)
	  continue;
	else {
	  gentid = ME_GEntID(edge);
	  fprintf(fp,"%d ",gentid);
	  if ((++k)%10 == 0)
	    fprintf(fp,"\n");
	}
      }
    }

    fprintf(fp,"\n");

    
    fprintf(fp,"icr1   1 \n"); 
    k = 0;
    idx = 0;
    while ((vertex = List_Next_Entry(vlist,&idx))) {
      icr = 3.0-MV_GEntDim(vertex);
      fprintf(fp,"%f ",(float)icr);
      if ((++k)%10 == 0)
	fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
  }



  nmeshatt = MESH_Num_Attribs(mesh);
      
  /* Must write out other attributes that exist on mesh entities */ 
  /* If natt is 0, all attributes are to be written out */
  
  if (nmeshatt) {

    /* Collect the attributes */
    
    outattribs = (MAttrib_ptr *) malloc(nmeshatt*sizeof(MAttrib_ptr));
    
    for (i = 0, noutatt = 0; i < nmeshatt; i++) {
      attrib = MESH_Attrib(mesh,i);
      
      attentdim = MAttrib_Get_EntDim(attrib);
      if ((attentdim == MVERTEX)  || (attentdim == MREGION) || 
          (ncells2 && attentdim == MFACE) || 
          (ncells1 && attentdim == MEDGE) || (attentdim == MALLTYPE)) {
        
        atttype = MAttrib_Get_Type(attrib);
        if (atttype == INT || atttype == DOUBLE) {
          
          MAttrib_Get_Name(attrib,attname);
          
          if (natt > 0) {
            found = 0;
            for (j = 0; j < natt; j++) {
              if (strcmp(attname,attnames[j]) == 0) {
                found = 1;
              }
	    }
	    if (found) {
              outattribs[noutatt] = attrib;
              noutatt++;
            }
            else
              continue;
	  }
          else {
            if (strncmp(attname+1,"idatt",5) == 0) {
              /* this is a special MSTK attribute created for writing out
                 contiguous entity IDs - ignore it */
              continue;
            }
            else {
              outattribs[noutatt] = attrib;
              noutatt++;
            }
          }
        }
      }
    }

    /* write them out */

    if (noutatt && !gmodel)
      fprintf(fp,"variable \n");

    for (i = 0; i < noutatt; i++) {
      attrib = outattribs[i];
      
      MAttrib_Get_Name(attrib,attname);
      atttype = MAttrib_Get_Type(attrib);
      attentdim = MAttrib_Get_EntDim(attrib);

      k = 0;
      if (attentdim == MVERTEX || attentdim == MALLTYPE) {
	fprintf(fp,"%s ",attname);
	fprintf(fp," 1 \n");
	
	idx = 0;
	while ((vertex = List_Next_Entry(vlist,&idx))) {
	  
	  MEnt_Get_AttVal(vertex,attrib,&ival,&rval,&pval);
	  
	  if (atttype == INT) {
	    fprintf(fp,"%d ",ival);
	    if ((k+1)%10 == 0 && k != nv) fprintf(fp,"\n");
	    k++;
	  }
	  else if (atttype == DOUBLE) {
	    fprintf(fp,"%14.7lf ", rval);
	    if ((k+1)%5 == 0 && k != nv) fprintf(fp,"\n");
	    k++;
	  }
	}
	if (k%5 != 0) fprintf(fp,"\n");
      }
      
      if (attentdim != MVERTEX && attentdim != MUNKNOWNTYPE) {
	fprintf(fp,"%s ",attname);
	fprintf(fp," 0 \n");
	
        if (ncells3) {
          idx = 0;
          while ((region = List_Next_Entry(rlist,&idx))) {
            
            MEnt_Get_AttVal(region,attrib,&ival,&rval,&pval);
            
            if (atttype == INT) {
              fprintf(fp,"%d ",ival);
              if ((k+1)%10 == 0) fprintf(fp,"\n");
              k++;
            }
            else if (atttype == DOUBLE) {
              fprintf(fp,"%lf ", rval);
              if ((k+1)%5 == 0) fprintf(fp,"\n");
              k++;
            }
          }
          if (k%5 != 0) fprintf(fp,"\n");
        }

	if (ncells2) {
	  idx = 0;
	  while ((face = List_Next_Entry(flist,&idx))) {
            int fmarked;
#ifdef MSTK_USE_MARKERS
            fmarked = MEnt_IsMarked(face,cellmk);
#else
            MEnt_Get_AttVal(face, cellmkatt, &fmarked, &rval, &pval);
#endif
            if (!fmarked)
	      continue;

	    MEnt_Get_AttVal(face,attrib,&ival,&rval,&pval);
	    
	    if (atttype == INT) {
	      fprintf(fp,"%d ",ival);
	      if ((k+1)%10 == 0) fprintf(fp,"\n");
	      k++;
	    }
	    else if (atttype == DOUBLE) {
	      fprintf(fp,"%lf ", rval);
	      if ((k+1)%5 == 0) fprintf(fp,"\n");
	      k++;
	    }
	  }
	  if (k%5 != 0) fprintf(fp,"\n");
	}
	
	if (ncells1) {
	  idx = 0;
	  while ((edge = List_Next_Entry(elist,&idx))) {
            int emarked;
#ifdef MSTK_USE_MARKERS
            emarked = MEnt_IsMarked(edge,cellmk);
#else
            MEnt_Get_AttVal(edge, cellmkatt, &emarked, &rval, &pval);
#endif
	    if (!emarked)
	      continue;

	    MEnt_Get_AttVal(edge,attrib,&ival,&rval,&pval);
	    
	    if (atttype == INT) {
	      fprintf(fp,"%d ",ival);
	      if ((k+1)%10 == 0) fprintf(fp,"\n");
	      k++;
	    }
	    else if (atttype == DOUBLE) {
	      fprintf(fp,"%lf ", rval);
	      if ((k+1)%5 == 0) fprintf(fp,"\n");
	      k++;
	    }
	  }
	  if (k%5 != 0) fprintf(fp,"\n");
	}
      }
    }

    if (noutatt)
      fprintf(fp,"\n");
	  
    free(outattribs);
  }


#ifdef MSTK_HAVE_MPI

  /* Write out global ID and Master Processor Rank for each node and
     for each element */

  fprintf(fp,"Node_GlobalID 1\n");
  idx = 0; k = 0;
  while ((vertex = List_Next_Entry(vlist,&idx))) {
    fprintf(fp,"%d ",MEnt_GlobalID((MEntity_ptr)vertex));
    if ((k+1)%5 == 0) fprintf(fp,"\n");
    k++;
  }
  if (k%5 != 0) fprintf(fp,"\n");

  fprintf(fp,"Node_MasterPID 1\n");
  idx = 0; k = 0;
  while ((vertex = List_Next_Entry(vlist,&idx))) {
    fprintf(fp,"%d ",MEnt_MasterParID((MEntity_ptr)vertex));
    if ((k+1)%10 == 0) fprintf(fp,"\n");
    k++;
  }
  if (k%10 != 0) fprintf(fp,"\n");



  fprintf(fp,"Cell_GlobalID 0\n");
  k = 0;
  if (ncells3) {
    idx = 0;
    while ((region = List_Next_Entry(rlist,&idx))) {
      fprintf(fp,"%d ",MEnt_GlobalID((MEntity_ptr)region));
      if ((k+1)%5 == 0) fprintf(fp,"\n");
      k++;
    }
  }
  
  if (ncells2) {
    idx = 0;
    while ((face = List_Next_Entry(flist,&idx))) {
      int fmarked;
#ifdef MSTK_USE_MARKERS
      fmarked = MEnt_IsMarked(face,cellmk);
#else
      MEnt_Get_AttVal(face, cellmkatt, &fmarked, &rval, &pval);
#endif
      if (!fmarked)
        continue;
      else {
        fprintf(fp,"%d ",MEnt_GlobalID((MEntity_ptr)face));
        if ((k+1)%5 == 0) fprintf(fp,"\n");
        k++;
      }
    }
  }

  if (ncells1) {
    idx = 0;
    while ((edge = List_Next_Entry(elist,&idx))) {
      int emarked;
#ifdef MSTK_USE_MARKERS
      emarked = MEnt_IsMarked(edge,cellmk);
#else
      MEnt_Get_AttVal(edge, cellmkatt, &emarked, &rval, &pval);
#endif
      if (!emarked)
        continue;
      else {
        fprintf(fp,"%d ",MEnt_GlobalID((MEntity_ptr)edge));
        if ((k+1)%5 == 0) fprintf(fp,"\n");
        k++;
      }
    }
  }
  if (k%5 != 0) fprintf(fp,"\n");



  fprintf(fp,"Cell_MasterPID 0\n");
  k = 0;
  if (ncells3) {
    idx = 0;
    while ((region = List_Next_Entry(rlist,&idx))) {
      fprintf(fp,"%d ",MEnt_MasterParID((MEntity_ptr)region));
      if ((k+1)%10 == 0) fprintf(fp,"\n");
      k++;
    }
  }
  
  if (ncells2) {    
    idx = 0;
    while ((face = List_Next_Entry(flist,&idx))) {
      int fmarked;
#ifdef MSTK_USE_MARKERS
      fmarked = MEnt_IsMarked(face,cellmk);
#else
      MEnt_Get_AttVal(face, cellmkatt, &fmarked, &rval, &pval);
#endif
      if (!fmarked)
        continue;
      else {
        fprintf(fp,"%d ",MEnt_MasterParID((MEntity_ptr)face));
        if ((k+1)%10 == 0) fprintf(fp,"\n");
        k++;
      }
    }
  }

  if (ncells1) {
    idx = 0;
    while ((edge = List_Next_Entry(elist,&idx))) {
      int emarked;
#ifdef MSTK_USE_MARKERS
      emarked = MEnt_IsMarked(edge,cellmk);
#else
      MEnt_Get_AttVal(edge, cellmkatt, &emarked, &rval, &pval);
#endif
      if (!emarked)
        continue;
      else {
        fprintf(fp,"%d ",MEnt_MasterParID((MEntity_ptr)edge));
        if ((k+1)%10 == 0) fprintf(fp,"\n");
        k++;
      }
    }
  }
  if (k%10 != 0) fprintf(fp,"\n");

#endif /* MSTK_HAVE_MPI */

  fprintf(fp,"endvars \n");
  
  
  /* Finish Export */

  fprintf(fp,"endgmv\n");

  fclose(fp);



  /* Clean up */

  free(gentities);

#ifdef MSTK_USE_MARKERS
  if (ncells2) {
    idx = 0;
    while ((face = List_Next_Entry(flist,&idx))) 
      MEnt_Unmark(face,cellmk);
  }

  if (ncells1) {
    idx = 0;
    while ((edge = List_Next_Entry(elist,&idx)))
      MEnt_Unmark(edge,cellmk);
  }

  MSTK_FreeMarker(cellmk);
#else
  MAttrib_Delete(cellmkatt);
#endif

  if (vidatt) MAttrib_Delete(vidatt);
  if (eidatt) MAttrib_Delete(eidatt);
  if (fidatt) MAttrib_Delete(fidatt);
  if (ridatt) MAttrib_Delete(ridatt);

  if (vlist) List_Delete(vlist);
  if (elist) List_Delete(elist);
  if (flist) List_Delete(flist);
  if (rlist) List_Delete(rlist);

  return 1;
}

#ifdef __cplusplus
  }
#endif

