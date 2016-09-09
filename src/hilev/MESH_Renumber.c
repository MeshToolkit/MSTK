#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MSTK_private.h"
#include "MSTK.h"



void Graph_Build_LevelStructure(int nnodes, int *adj, int *nadj, int **offset, 
                                int start_node, int *newmap, int *level, 
                                int *depth, int *maxwidth);


int Graph_Renumber_GPS(int nnodes, int nstart, int *nadj, int *adj, int *newmap, int *depth, 
                       int *maxwidth);


/* Enforce continuous and possibly optimized numbering for mesh entities */
/* renum_type = 0   --- Sequential renumbering 
              = 1   --- Reverse Cuthill Mckee 
   mtype      = Entity type to reorder (MALLTYPE for all entities)
                MVERTEX, MEDGE, MFACE, MREGION, MALLTYPE
*/



void MESH_Renumber(Mesh_ptr mesh, int renum_type, MType mtype) {
  MVertex_ptr mv, v0=NULL;
  MEdge_ptr me, e0=NULL;
  MFace_ptr mf, f0=NULL;
  MRegion_ptr mr, r0=NULL;
  int idx, idx2, idx3;
  int i, j;
  int done;
  MAttrib_ptr vidatt;
  List_ptr vlist;
  double xyz[3];
  int bandwidth, maxbandwidth1, maxbandwidth2;
  double avebandwidth1, avebandwidth2;

  if (renum_type == 0) {
    if (mtype == MVERTEX || mtype == MALLTYPE) {
      int nv = 0;
      idx = 0; 
      while ((mv = MESH_Next_Vertex(mesh,&idx)))
	MV_Set_ID(mv,++nv);
    }
    
    if (mtype == MEDGE || mtype == MALLTYPE) {
      int ne = 0;
      idx = 0; 
      while ((me = MESH_Next_Edge(mesh,&idx)))
	ME_Set_ID(me,++ne);
    }
    
    if (mtype == MFACE || mtype == MALLTYPE) {
      int nf = 0;
      idx = 0; 
      while ((mf = MESH_Next_Face(mesh,&idx)))
	MF_Set_ID(mf,++nf);
    }
    
    if (mtype == MREGION || mtype == MALLTYPE) {
      int nr = 0;
      idx = 0; 
      while ((mr = MESH_Next_Region(mesh,&idx)))
	MR_Set_ID(mr,++nr);
    }
  }
  else if (renum_type == 1) {
    double minx, miny, minz;
    int mkid = MSTK_GetMarker();
    int minid, maxid;
    int *nadj, *newmap, *adj, *offset, nconn=0;
    int nalloc, depth, maxwidth;
        


    if (mtype == MVERTEX || mtype == MALLTYPE) { 
      int nv = MESH_Num_Vertices(mesh);      

      /* Compute a graph of vertex connections across elements (faces
	 for surface meshes, regions for solid meshes */

      /* Start with the vertex in the lower leftmost corner */
      
      minx = miny = minz = 1.0e+12;
      v0 = NULL;
      idx = 0;
      while ((mv = MESH_Next_Vertex(mesh,&idx))) {
	MV_Coords(mv,xyz);
	if (xyz[0] <= minx && xyz[1] <= miny && xyz[2] <= minz) {
	  minx = xyz[0];
	  miny = xyz[1];
	  minz = xyz[2];
	  v0 = mv;
	}
      }


      nadj = (int *) malloc(nv*sizeof(int));
      nalloc = nv*5;
      adj = (int *) malloc(nalloc*sizeof(int));

      if (!MESH_Num_Regions(mesh)) {
        int nentries = 0;
        i = 0;
        idx = 0;
        while ((mv = MESH_Next_Vertex(mesh,&idx))) {
          List_ptr vfaces, adjvlist;
          MFace_ptr vf;
          MVertex_ptr adjv;

          adjvlist = List_New(0);
          vfaces = MV_Faces(mv);
          idx2 = 0;
          while ((vf = List_Next_Entry(vfaces,&idx2))) {
            List_ptr fverts = MF_Vertices(vf,1,0);
            idx3 = 0;
            while ((adjv = List_Next_Entry(fverts,&idx3))) {
              if (adjv != mv && !MEnt_IsMarked(adjv,mkid)) {
                List_Add(adjvlist,adjv);
                MEnt_Mark(adjv,mkid);
              }
            }
            List_Delete(fverts);
          }
          List_Delete(vfaces);
          List_Unmark(adjvlist,mkid);

          nadj[i] = List_Num_Entries(adjvlist);
          
          if (nentries+nadj[i] > nalloc) {
            nalloc *= 2;
            adj = (int *) realloc(adj,nalloc*sizeof(int));
          }

          idx2 = 0;
          while ((adjv = List_Next_Entry(adjvlist,&idx2)))
            adj[nentries++] = MV_ID(adjv)-1;

          List_Delete(adjvlist);
          i++;
        }
      }
      else {
        int nentries = 0;
        i = 0;
        idx = 0;
        while ((mv = MESH_Next_Vertex(mesh,&idx))) {
          List_ptr vregions, adjvlist;
          MRegion_ptr vr;
          MVertex_ptr adjv;

          adjvlist = List_New(0);
          vregions = MV_Regions(mv);
          idx2 = 0;
          while ((vr = List_Next_Entry(vregions,&idx2))) {
            List_ptr rverts = MR_Vertices(vr);
            idx3 = 0;
            while ((adjv = List_Next_Entry(rverts,&idx3))) {
              if (adjv != mv && !MEnt_IsMarked(adjv,mkid)) {
                List_Add(adjvlist,adjv);
                MEnt_Mark(adjv,mkid);
              }
            }
            List_Delete(rverts);
          }
          List_Delete(vregions);
          List_Unmark(adjvlist,mkid);

          nadj[i] = List_Num_Entries(adjvlist);

          if (nentries+nadj[i] > nalloc) {
            nalloc *= 2;
            adj = (int *) realloc(adj,nalloc*sizeof(int));
          }

          idx2 = 0;
          while ((adjv = List_Next_Entry(adjvlist,&idx2)))
            adj[nentries++] = MV_ID(adjv)-1;

          List_Delete(adjvlist);
          i++;
        }
      }

      /* Compute offsets into adj array */

      offset = (int *) malloc(nv*sizeof(int));
      offset[0] = 0;
      for (i = 1; i < nv; i++)
	offset[i] = offset[i-1] + nadj[i-1];

      /* Compute maximum bandwidth before renumbering */

      maxbandwidth1 = 0;
      avebandwidth1 = 0;
      for (i = 0; i < nv; i++) {
	int off = offset[i];
	int curid = i;
	for (j = 0; j < nadj[i]; j++) {
	  int adjid = adj[off+j];
	  int diff = abs(adjid-curid);
	  maxbandwidth1 = (diff > maxbandwidth1) ? diff : maxbandwidth1;
	  avebandwidth1 += diff;
	  nconn++;
	}
      }
      nconn = offset[nv-1]+nadj[nv-1];
      avebandwidth1 /= nconn;

      fprintf(stderr,
	      "Ave vertex ID difference on elements before renumbering: %-lf\n",
	      avebandwidth1);
      fprintf(stderr,
	      "Max vertex ID difference on elements before renumbering: %-d\n",
	      maxbandwidth1);
    
      fprintf(stderr,"\n");

      newmap = (int *) malloc(nv*sizeof(int));
      Graph_Renumber_GPS(nv, MV_ID(v0)-1, nadj, adj, newmap, &depth, &maxwidth);


      /* Compute bandwidth after renumbering */

      maxbandwidth2 = 0;
      avebandwidth2 = 0;
      for (i = 0; i < nv; i++) {
	int off = offset[i];
	int curid = newmap[i];
	for (j = 0; j < nadj[i]; j++) {
	  int adjid = newmap[adj[off+j]];
	  int diff = abs(adjid-curid);
	  maxbandwidth2 = (diff > maxbandwidth2) ? diff : maxbandwidth2;
	  avebandwidth2 += diff;
	  nconn++;
	}
      }
      nconn = offset[nv-1]+nadj[nv-1];
      avebandwidth2 /= nconn;


      if (maxbandwidth2 < maxbandwidth1 && avebandwidth2 < avebandwidth1) {

        /* Renumber */
      
        idx = 0; i = 0;
        while ((mv = MESH_Next_Vertex(mesh,&idx))) {
          MV_Set_ID(mv,newmap[i]+1);
          i++;
        }

        fprintf(stderr,
                "Ave vertex ID difference on elements after renumbering: %-lf\n",
                avebandwidth2);
        fprintf(stderr,
                "Max vertex ID difference on elements after renumbering: %-d\n",
                maxbandwidth2);    
      }
      else {
        nv = 0;
        idx = 0; 
        while ((mv = MESH_Next_Vertex(mesh,&idx)))
          MV_Set_ID(mv,++nv);

        fprintf(stderr,"Bandwidth did not improve. Keeping old numbering with gaps eliminated\n");
      }
      fprintf(stderr,"\n\n\n");
      
      free(nadj);
      free(adj);
      free(offset);
      free(newmap);

    }
 



    /* Reorder edges according to a breadth first algorithm applied to
       edges (differs from RCM in that it does not add adjacent nodes
       in ascending order of their valence) */

    if (mtype == MEDGE || mtype == MALLTYPE) {
      int ne = MESH_Num_Edges(mesh);
      MEdge_ptr ve;
      List_ptr elist;

      /************************* renumbering code ****************************/

      ne = MESH_Num_Edges(mesh);

      if (mtype == MALLTYPE) { 

	/* RCM algorithm already applied on the vertices. Use an edge
	   connected to the starting vertex as the first edge */

	List_ptr vedges = MV_Edges(v0);
	e0 = List_Entry(vedges,0);
	List_Delete(vedges);
      }
      else {
	/* Find the edge whose mid point is a minimum point */
	minx = miny = minz = 1.0e+12;
	e0 = NULL;
	idx = 0;
	while ((me = MESH_Next_Edge(mesh,&idx))) {
          double exyz[2][3];

	  MV_Coords(ME_Vertex(me,0),exyz[0]);
	  MV_Coords(ME_Vertex(me,1),exyz[1]);
	  xyz[0] = (exyz[0][0]+exyz[1][0])/2.0;
	  xyz[1] = (exyz[0][1]+exyz[1][1])/2.0;
	  xyz[2] = (exyz[0][2]+exyz[1][2])/2.0;
	  if (xyz[0] < minx && xyz[1] < miny && xyz[2] < minz) {
	    minx = xyz[0];
	    miny = xyz[1];
	    minz = xyz[2];
	    e0 = me;
	  }
	}
      }


      nadj = (int *) malloc(ne*sizeof(int));
      nalloc = ne*5;
      adj = (int *) malloc(nalloc*sizeof(int));

      if (!MESH_Num_Regions(mesh)) {
        int nentries = 0;
        i = 0;
        idx = 0;
        while ((me = MESH_Next_Edge(mesh,&idx))) {
          List_ptr efaces, adjelist;
          MFace_ptr ef;
          MEdge_ptr adje;

          adjelist = List_New(0);
          efaces = ME_Faces(me);
          idx2 = 0;
          while ((ef = List_Next_Entry(efaces,&idx2))) {
            List_ptr fedges = MF_Edges(ef,1,0);
            idx3 = 0;
            while ((adje = List_Next_Entry(fedges,&idx3))) {
              if (adje != me && !MEnt_IsMarked(adje,mkid)) {
                List_Add(adjelist,adje);
                MEnt_Mark(adje,mkid);
              }
            }
            List_Delete(fedges);
          }
          List_Delete(efaces);
          List_Unmark(adjelist,mkid);

          nadj[i] = List_Num_Entries(adjelist);

          if (nentries+nadj[i] > nalloc) {
            nalloc *= 2;
            adj = (int *) realloc(adj,nalloc*sizeof(int));
          }

          idx2 = 0;
          while ((adje = List_Next_Entry(adjelist,&idx2)))
            adj[nentries++] = ME_ID(adje)-1;

          List_Delete(adjelist);
          i++;
        }
      }
      else {
        int nentries = 0;
        i = 0;
        idx = 0;
        while ((me = MESH_Next_Edge(mesh,&idx))) {
          List_ptr eregions, adjelist;
          MRegion_ptr er;
          MEdge_ptr adje;

          adjelist = List_New(0);
          eregions = ME_Regions(me);
          idx2 = 0;
          while ((er = List_Next_Entry(eregions,&idx2))) {
            List_ptr redges = MR_Edges(er);
            idx3 = 0;
            while ((adje = List_Next_Entry(redges,&idx3))) {
              if (adje != me && !MEnt_IsMarked(adje,mkid)) {
                List_Add(adjelist,adje);
                MEnt_Mark(adje,mkid);
              }
            }
            List_Delete(redges);
          }
          List_Delete(eregions);
          List_Unmark(adjelist,mkid);

          nadj[i] = List_Num_Entries(adjelist);

          if (nentries+nadj[i] > nalloc) {
            nalloc *= 2;
            adj = (int *) realloc(adj,nalloc*sizeof(int));
          }

          idx2 = 0;
          while ((adje = List_Next_Entry(adjelist,&idx2)))
            adj[nentries++] = ME_ID(adje)-1;

          List_Delete(adjelist);
          i++;
        }
      }

      /* Compute offsets into adj array */

      offset = (int *) malloc(ne*sizeof(int));
      offset[0] = 0;
      for (i = 1; i < ne; i++)
	offset[i] = offset[i-1] + nadj[i-1];


      /* Compute maximum bandwidth before renumbering */

      maxbandwidth1 = 0;
      avebandwidth1 = 0;
      for (i = 0; i < ne; i++) {
	int off = offset[i];
	int curid = i;
	for (j = 0; j < nadj[i]; j++) {
	  int adjid = adj[off+j];
	  int diff = abs(adjid-curid);
	  maxbandwidth1 = (diff > maxbandwidth1) ? diff : maxbandwidth1;
	  avebandwidth1 += diff;
	  nconn++;
	}
      }
      nconn = offset[ne-1]+nadj[ne-1];
      avebandwidth1 /= nconn;

      fprintf(stderr,
	      "Ave edge ID difference on elements before renumbering: %-lf\n",
	      avebandwidth1);
      fprintf(stderr,
	      "Max edge ID difference on elements before renumbering: %-d\n",
	      maxbandwidth1);
    
      fprintf(stderr,"\n");


      /* Call Graph Renumbering algorithm */

      newmap = (int *) malloc(ne*sizeof(int));
      Graph_Renumber_GPS(ne, ME_ID(e0)-1, nadj, adj, newmap, &depth, &maxwidth);


      /* Compute bandwidth after renumbering */

      maxbandwidth2 = 0;
      avebandwidth2 = 0;
      for (i = 0; i < ne; i++) {
	int off = offset[i];
	int curid = newmap[i];
	for (j = 0; j < nadj[i]; j++) {
	  int adjid = newmap[adj[off+j]];
	  int diff = abs(adjid-curid);
	  maxbandwidth2 = (diff > maxbandwidth2) ? diff : maxbandwidth2;
	  avebandwidth2 += diff;
	  nconn++;
	}
      }
      nconn = offset[ne-1]+nadj[ne-1];
      avebandwidth2 /= nconn;

      if (maxbandwidth2 < maxbandwidth1 && avebandwidth2 < avebandwidth1) {
        /* Renumber */
      
        idx = 0; i = 0;
        while ((me = MESH_Next_Edge(mesh,&idx))) {
          ME_Set_ID(me,newmap[i]+1);
          i++;
        }

        fprintf(stderr,
                "Ave edge ID difference on elements after renumbering: %-lf\n",
                avebandwidth2);
        fprintf(stderr,
                "Max edge ID difference on elements after renumbering: %-d\n",
                maxbandwidth2);
        
      }
      else {
        ne = 0;
        idx = 0; 
        while ((me = MESH_Next_Edge(mesh,&idx)))
          ME_Set_ID(me,++ne);

        fprintf(stderr,"Bandwidth did not improve. Keeping old numbering with gaps eliminated\n");
      }
      fprintf(stderr,"\n\n\n");


      free(nadj);
      free(adj);
      free(offset);
      free(newmap);

    }


    /* Reorder faces according to a breadth first algorithm applied to
       edges (differs from RCM in that it does not add adjacent graph nodes
       in ascending order of their valence) */

    if (mtype == MFACE || mtype == MALLTYPE) {      
      int nf = MESH_Num_Faces(mesh);

      if (mtype == MALLTYPE) { 

	/* RCM algorithm already applied on the vertices. Use an edge
	   connected to the starting vertex as the first edge */

	List_ptr vfaces = MV_Faces(v0);
	f0 = List_Entry(vfaces,0);
	List_Delete(vfaces);
      }
      else {
	/* Find the face whose mid point is a minimum point */
	minx = miny = minz = 1.0e+12;
	f0 = NULL;
	idx = 0;
	while ((mf = MESH_Next_Face(mesh,&idx))) {
          double fxyz[MAXPV2][3];
	  int nfv;

	  MF_Coords(mf,&nfv,fxyz);
	  xyz[0] = fxyz[0][0];
	  xyz[1] = fxyz[0][1];
	  xyz[2] = fxyz[0][2];
	  for (i = 1; i < nfv; i++) {
	    xyz[0] += fxyz[i][0];
	    xyz[1] += fxyz[i][1];
	    xyz[2] += fxyz[i][2];
	  }
	  xyz[0] /= nfv; xyz[1] /= nfv; xyz[2] /= nfv;
	  if (xyz[0] < minx && xyz[1] < miny && xyz[2] < minz) {
	    minx = xyz[0];
	    miny = xyz[1];
	    minz = xyz[2];
	    f0 = mf;
	  }
	}
      }


      nadj = (int *) malloc(nf*sizeof(int));
      nalloc = nf*5;
      adj = (int *) malloc(nalloc*sizeof(int));

      if (!MESH_Num_Regions(mesh)) {
        int nentries = 0;
        i = 0;
        idx = 0;
        while ((mf = MESH_Next_Face(mesh,&idx))) {
          List_ptr vfaces, fverts, adjflist;
          MFace_ptr vf, adjf;
          MVertex_ptr fv;

          adjflist = List_New(0);
	  fverts = MF_Vertices(mf,1,0);
          idx2 = 0;
          while ((fv = List_Next_Entry(fverts,&idx2))) {
            List_ptr vfaces = MV_Faces(fv);
            idx3 = 0;
            while ((adjf = List_Next_Entry(vfaces,&idx3))) {
              if (adjf != mf && !MEnt_IsMarked(adjf,mkid)) {
                List_Add(adjflist,adjf);
                MEnt_Mark(adjf,mkid);
              }
            }	    
            List_Delete(vfaces);
          }
          List_Delete(fverts);
          List_Unmark(adjflist,mkid);

          nadj[i] = List_Num_Entries(adjflist);

          if (nentries+nadj[i] > nalloc) {
            nalloc *= 2;
            adj = (int *) realloc(adj,nalloc*sizeof(int));
          }

          idx2 = 0;
          while ((adjf = List_Next_Entry(adjflist,&idx2)))
            adj[nentries++] = MF_ID(adjf)-1;

          List_Delete(adjflist);
          i++;
        }
      }
      else {
        int nentries = 0;
        i = 0;
        idx = 0;
        while ((mf = MESH_Next_Face(mesh,&idx))) {
          List_ptr fregions, adjflist;
          MRegion_ptr fr;
          MFace_ptr adjf;

          adjflist = List_New(0);
          fregions = MF_Regions(mf);
          idx2 = 0;
          while ((fr = List_Next_Entry(fregions,&idx2))) {
            List_ptr rfaces = MR_Faces(fr);
            idx3 = 0;
            while ((adjf = List_Next_Entry(rfaces,&idx3))) {
              if (adjf != mf && !MEnt_IsMarked(adjf,mkid)) {
                List_Add(adjflist,adjf);
                MEnt_Mark(adjf,mkid);
              }
            }
            List_Delete(rfaces);
          }
          List_Delete(fregions);
          List_Unmark(adjflist,mkid);

          nadj[i] = List_Num_Entries(adjflist);

          if (nentries+nadj[i] > nalloc) {
            nalloc *= 2;
            adj = (int *) realloc(adj,nalloc*sizeof(int));
          }

          idx2 = 0;
          while ((adjf = List_Next_Entry(adjflist,&idx2)))
            adj[nentries++] = MF_ID(adjf)-1;

          List_Delete(adjflist);
          i++;
        }
      }

      /* Compute offsets into adj array */

      offset = (int *) malloc(nf*sizeof(int));
      offset[0] = 0;
      for (i = 1; i < nf; i++)
	offset[i] = offset[i-1] + nadj[i-1];


      /* Compute maximum bandwidth before renumbering */

      maxbandwidth1 = 0;
      avebandwidth1 = 0;
      for (i = 0; i < nf; i++) {
	int off = offset[i];
	int curid = i;
	for (j = 0; j < nadj[i]; j++) {
	  int adjid = adj[off+j];
	  int diff = abs(adjid-curid);
	  maxbandwidth1 = (diff > maxbandwidth1) ? diff : maxbandwidth1;
	  avebandwidth1 += diff;
	  nconn++;
	}
      }
      nconn = offset[nf-1]+nadj[nf-1];
      avebandwidth1 /= nconn;

      if (MESH_Num_Regions(mesh)) {
        fprintf(stderr,
                "Ave face ID difference on elements before renumbering: %-lf\n",
                avebandwidth1);
        fprintf(stderr,
                "Max face ID difference on elements before renumbering: %-d\n",
                maxbandwidth1);
      }
      else {
        fprintf(stderr,
                "Ave face ID difference before renumbering: %-lf\n",
                avebandwidth1);
        fprintf(stderr,
                "Max face ID difference before renumbering: %-d\n",
                maxbandwidth1);
      }
    
      fprintf(stderr,"\n");


      /* Call Graph Renumbering algorithm */

      newmap = (int *) malloc(nf*sizeof(int));
      Graph_Renumber_GPS(nf, MF_ID(f0)-1, nadj, adj, newmap, &depth, &maxwidth);


      /* Compute bandwidth after renumbering */

      maxbandwidth2 = 0;
      avebandwidth2 = 0;
      for (i = 0; i < nf; i++) {
	int off = offset[i];
	int curid = newmap[i];
	for (j = 0; j < nadj[i]; j++) {
	  int adjid = newmap[adj[off+j]];
	  int diff = abs(adjid-curid);
	  maxbandwidth2 = (diff > maxbandwidth2) ? diff : maxbandwidth2;
	  avebandwidth2 += diff;
	  nconn++;
	}
      }
      nconn = offset[nf-1]+nadj[nf-1];
      avebandwidth2 /= nconn;


      if (maxbandwidth2 < maxbandwidth1 && avebandwidth2 < avebandwidth1) {
        /* Renumber */
      
        idx = 0; i = 0;
        while ((mf = MESH_Next_Face(mesh,&idx))) {
          MF_Set_ID(mf,newmap[i]+1);
          i++;
        }

        if (MESH_Num_Regions(mesh)) {
          fprintf(stderr,
                  "Ave face ID difference on elements after renumbering: %-lf\n",
                  avebandwidth2);
          fprintf(stderr,
                  "Max face ID difference on elements after renumbering: %-d\n",
                  maxbandwidth2);
        }
        else {
          fprintf(stderr,
                  "Ave face ID difference after renumbering: %-lf\n",
                  avebandwidth2);
          fprintf(stderr,
                  "Max face ID difference after renumbering: %-d\n",
                  maxbandwidth2);
        }
    
      }
      else {
        nf = 0;
        idx = 0; 
        while ((mf = MESH_Next_Face(mesh,&idx)))
          MF_Set_ID(mf,++nf);

        fprintf(stderr,"Bandwidth did not improve. Keeping old numbering with gaps eliminated\n");
      }
      fprintf(stderr,"\n\n\n");

      free(nadj);
      free(adj);
      free(offset);
      free(newmap);
    }


    if (mtype == MREGION || mtype == MALLTYPE) {
      int nr = MESH_Num_Regions(mesh);

      if (nr) {

        if (mtype == MALLTYPE) { 

          /* Renumbering algorithm already applied on the vertices. Use
             a region connected to the starting vertex as the first
             region */

          List_ptr vregions = MV_Regions(v0);
          r0 = List_Entry(vregions,0);
          List_Delete(vregions);
        }
        else {
          /* Find the region whose center point is a minimum point */
          minx = miny = minz = 1.0e+12;
          r0 = NULL;
          idx = 0;
          while ((mr = MESH_Next_Region(mesh,&idx))) {
            double rxyz[MAXPV3][3];
            int nrv;

            MR_Coords(mr,&nrv,rxyz);
            xyz[0] = rxyz[0][0];
            xyz[1] = rxyz[0][1];
            xyz[2] = rxyz[0][2];
            for (i = 1; i < nrv; i++) {
              xyz[0] += rxyz[i][0];
              xyz[1] += rxyz[i][1];
              xyz[2] += rxyz[i][2];
            }
            xyz[0] /= nrv; xyz[1] /= nrv; xyz[2] /= nrv;
            if (xyz[0] < minx && xyz[1] < miny && xyz[2] < minz) {
              minx = xyz[0];
              miny = xyz[1];
              minz = xyz[2];
              r0 = mr;
            }
          }
        }


        nadj = (int *) malloc(nr*sizeof(int));
        nalloc = nr*5;
        adj = (int *) malloc(nalloc*sizeof(int));

        int nentries = 0;
        i = 0;
        idx = 0;
        while ((mr = MESH_Next_Region(mesh,&idx))) {
          List_ptr vregions, rverts, adjrlist;
          MRegion_ptr vr, adjr;
          MVertex_ptr rv;

          adjrlist = List_New(0);
          rverts = MR_Vertices(mr);
          idx2 = 0;
          while ((rv = List_Next_Entry(rverts,&idx2))) {
            List_ptr vregions = MV_Regions(rv);
            idx3 = 0;
            while ((adjr = List_Next_Entry(vregions,&idx3))) {
              if (adjr != mr && !MEnt_IsMarked(adjr,mkid)) {
                List_Add(adjrlist,adjr);
                MEnt_Mark(adjr,mkid);
              }
            }	    
            List_Delete(vregions);
          }
          List_Delete(rverts);
          List_Unmark(adjrlist,mkid);

          nadj[i] = List_Num_Entries(adjrlist);

          if (nentries+nadj[i] > nalloc) {
            nalloc *= 2;
            adj = (int *) realloc(adj,nalloc*sizeof(int));
          }

          idx2 = 0;
          while ((adjr = List_Next_Entry(adjrlist,&idx2)))
            adj[nentries++] = MR_ID(adjr)-1;

          List_Delete(adjrlist);
          i++;
        }

        /* Compute offsets into adj array */

        offset = (int *) malloc(nr*sizeof(int));
        offset[0] = 0;
        for (i = 1; i < nr; i++)
          offset[i] = offset[i-1] + nadj[i-1];


        /* Compute maximum bandwidth before renumbering */

        maxbandwidth1 = 0;
        avebandwidth1 = 0;
        for (i = 0; i < nr; i++) {
          int off = offset[i];
          int curid = i;
          for (j = 0; j < nadj[i]; j++) {
            int adjid = adj[off+j];
            int diff = abs(adjid-curid);
            maxbandwidth1 = (diff > maxbandwidth1) ? diff : maxbandwidth1;
            avebandwidth1 += diff;
            nconn++;
          }
        }
        nconn = offset[nr-1]+nadj[nr-1];
        avebandwidth1 /= nconn;

        fprintf(stderr,
                "Ave region ID difference before renumbering: %-lf\n",
                avebandwidth1);
        fprintf(stderr,
                "Max region ID difference before renumbering: %-d\n",
                maxbandwidth1);
    
        fprintf(stderr,"\n");


        /* Call Graph Renumbering algorithm */

        newmap = (int *) malloc(nr*sizeof(int));
        Graph_Renumber_GPS(nr, MR_ID(r0)-1, nadj, adj, newmap, &depth, &maxwidth);

        /* Compute bandwidth after renumbering */

        maxbandwidth2 = 0;
        avebandwidth2 = 0;
        for (i = 0; i < nr; i++) {
          int off = offset[i];
          int curid = newmap[i];
          for (j = 0; j < nadj[i]; j++) {
            int adjid = newmap[adj[off+j]];
            int diff = abs(adjid-curid);
            maxbandwidth2 = (diff > maxbandwidth2) ? diff : maxbandwidth2;
            avebandwidth2 += diff;
            nconn++;
          }
        }
        nconn = offset[nr-1]+nadj[nr-1];
        avebandwidth2 /= nconn;

        if (maxbandwidth2 < maxbandwidth1 && avebandwidth2 < avebandwidth1) {

          /* Renumber */
      
          idx = 0; i = 0;
          while ((mr = MESH_Next_Region(mesh,&idx))) {
            MR_Set_ID(mr,newmap[i]+1);
            i++;
          }

          fprintf(stderr,
                  "Ave region ID difference after renumbering: %-lf\n",
                  avebandwidth2);
          fprintf(stderr,
                  "Max region ID difference after renumbering: %-d\n",
                  maxbandwidth2);
    
        }
        else {
          nr = 0;
          idx = 0; 
          while ((mr = MESH_Next_Region(mesh,&idx)))
            MR_Set_ID(mr,++nr);

          fprintf(stderr,"Bandwidth did not improve. Keeping old numbering with gaps eliminated\n");
        }
        fprintf(stderr,"\n\n\n");


        free(nadj);
        free(adj);
        free(offset);
        free(newmap);
      }
    }

    MSTK_FreeMarker(mkid);
  }
  

  vidatt = MAttrib_New(mesh,"vidrcm",INT,MVERTEX);
  idx = 0;
  while ((mv = MESH_Next_Vertex(mesh,&idx))) {
    MEnt_Set_AttVal(mv,vidatt,MV_ID(mv),0.0,NULL);
  }
    
  return;
}

