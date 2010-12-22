#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "MSTK.h"
#include "MSTK_util.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MV_Print(MVertex_ptr v, int lvl) {
    double xyz[3];
    MEdge_ptr e;
    MFace_ptr f;
    MRegion_ptr r;
    List_ptr vedges, vfaces, vregions;
    int idx, ne, nf, nr;
    GEntity_ptr gent;
    
    if (!v)
      return;

    if (MEnt_Dim(v) != MVERTEX && MEnt_OrigDim(v) != MVERTEX) {
      fprintf(stderr,"NOT A VERTEX!\n");
      return;
    }

    fprintf(stderr,"\n");
    fprintf(stderr,"Vertex: 0x%-x   ID: %-d   ",(unsigned int)v,MV_ID(v));
    
    
    if (lvl) {
      /* Print more detailed information */
      
      if (MV_GEntDim(v) == -1)
	fprintf(stderr,"Unknown Classification\n");
      else {
	fprintf(stderr,"GEntID: %-d    GEntDim: %-d\n", MV_GEntID(v),
		MV_GEntDim(v));
	if ((gent = MV_GEntity(v)))
	  fprintf(stderr,"Model entity pointer: 0x%-x\n",(unsigned int)gent);	
      }

      MV_Coords(v,xyz);
      fprintf(stderr,"Coordinates: %16.8lf %16.8lf %16.8lf\n",xyz[0],xyz[1],xyz[2]);
    
      
      if (lvl > 1) {
	/* Print upward adjacency information as well */
	
	fprintf(stderr,"\n");
	
	
	vedges = MV_Edges(v);
        if (vedges) {
	  ne = List_Num_Entries(vedges);
	
	  if (ne) {
	    fprintf(stderr,"Edges:%-d\n",ne);
	  
	    fprintf(stderr,"Object        ID      GEntID   GEntDim    Vertex IDs\n");
	    idx = 0;
	    while ((e = List_Next_Entry(vedges,&idx))) {
	      fprintf(stderr,"0x%-8x    %-8d %-8d     %-1d    %-d  %-d\n",
		      (unsigned int)e,ME_ID(e),ME_GEntID(e),ME_GEntDim(e),
		      MV_ID(ME_Vertex(e,0)),MV_ID(ME_Vertex(e,1)));
	    }
	    fprintf(stderr,"\n");
	  }
	  List_Delete(vedges);
        }
	
	
	vfaces = MV_Faces(v);
	if (vfaces) {
	  nf = List_Num_Entries(vfaces);
	  fprintf(stderr,"Faces: %-d\n",nf);
	  
	  fprintf(stderr,"Object        ID      GEntID   GEntDim\n");
	  idx = 0;
	  while ((f = List_Next_Entry(vfaces,&idx)))
	    fprintf(stderr,"0x%-8x    %-8d %-8d     %-1d\n",(unsigned int)f,
		    MF_ID(f),MF_GEntID(f),ME_GEntDim(f));
	  fprintf(stderr,"\n");
	  List_Delete(vfaces);
	}
	
	
	vregions= MV_Regions(v);
	if (vregions) {
	  nr = List_Num_Entries(vregions);
	  fprintf(stderr,"Regions: %-d\n",nr);
	  
	  fprintf(stderr,"Object        ID      GEntID   GEntDim\n");
	  idx = 0;
	  while ((r = List_Next_Entry(vregions,&idx)))
	    fprintf(stderr,"0x%-8x %-8d %-8d     %-1d\n",(unsigned int)r,
		    MR_ID(r),MR_GEntID(r),MR_GEntDim(r));
	  
	  fprintf(stderr,"\n");
	  List_Delete(vregions);      
	}
      }
    }
    else {
      fprintf(stderr,"\n");
      MV_Coords(v,xyz);
      fprintf(stderr,"Coordinates: %16.8lf %16.8lf %16.8lf\n",xyz[0],xyz[1],xyz[2]);
    }

    fprintf(stderr,"\n");
  }
  

  void ME_Print(MEdge_ptr e, int lvl) {
    MVertex_ptr v0, v1; 
    MFace_ptr f;
    MRegion_ptr r;
    List_ptr efaces, eregions;
    int idx, nf, nr;
    GEntity_ptr gent;
    double xyz[3];

    if (!e)
      return;

    if (MEnt_Dim(e) != MEDGE && MEnt_OrigDim(e) != MEDGE) {
      fprintf(stderr,"NOT AN EDGE!\n");
      return;
    }

    fprintf(stderr,"\n");
    fprintf(stderr,"Edge: 0x%-x   ID: %-d   ",(unsigned int)e,ME_ID(e));
    
    
    if (lvl) {
      /* Print more detailed information */
      
      if (ME_GEntDim(e) == -1)
	fprintf(stderr,"Unknown Classification\n");
      else {
	fprintf(stderr,"GEntID: %-d    GEntDim: %-d\n", ME_GEntID(e),
		ME_GEntDim(e));
	if ((gent = ME_GEntity(e)))
	  fprintf(stderr,"Model entity pointer: 0x%-x\n",(unsigned int)gent);
      }

      fprintf(stderr,"\n");

      fprintf(stderr,"Vertices:\n");
      fprintf(stderr,"Object        ID      GEntID   GEntDim\n");
      v0 = ME_Vertex(e,0);
      fprintf(stderr,"0x%-8x    %-8d %-8d     %-1d\n",(unsigned int)v0,
	      MV_ID(v0),MV_GEntID(v0),MV_GEntDim(v0));
      if (lvl > 1) {
	MV_Coords(v0,xyz);
        fprintf(stderr,"Coordinates: %16.8lf %16.8lf %16.8lf\n",xyz[0],xyz[1],xyz[2]);
      }

      v1 = ME_Vertex(e,1);
      fprintf(stderr,"0x%-8x    %-8d %-8d     %-1d\n",(unsigned int)v1,
	      MV_ID(v1),MV_GEntID(v1),MV_GEntDim(v1));
      if (lvl > 1) {
	MV_Coords(v1,xyz);
        fprintf(stderr,"Coordinates: %16.8lf %16.8lf %16.8lf\n",xyz[0],xyz[1],xyz[2]);
      }

      if (lvl > 1) {
	/* Print upward adjacency information as well */
	
	fprintf(stderr,"\n");


	efaces = ME_Faces(e);
	if (efaces) {
	  nf = List_Num_Entries(efaces);
	  fprintf(stderr,"Faces: %-d\n",nf);
	  
	  fprintf(stderr,"Object        ID      GEntID   GEntDim\n");
	  idx = 0;
	  while ((f = List_Next_Entry(efaces,&idx)))
	    fprintf(stderr,"0x%-8x    %-8d %-8d     %-1d\n",(unsigned int)f,
		    MF_ID(f),MF_GEntID(f),ME_GEntDim(f));
	  fprintf(stderr,"\n");
	  List_Delete(efaces);
	}
	
	
	eregions= ME_Regions(e);
	if (eregions) {
	  nr = List_Num_Entries(eregions);
	  fprintf(stderr,"Regions: %-d\n",nr);
	  
	  fprintf(stderr,"Object        ID      GEntID   GEntDim\n");
	  idx = 0;
	  while ((r = List_Next_Entry(eregions,&idx)))
	    fprintf(stderr,"0x%-8x %-8d %-8d     %-1d\n",(unsigned int)r,
		    MR_ID(r),MR_GEntID(r),MR_GEntDim(r));
	  
	  fprintf(stderr,"\n");
	  List_Delete(eregions);      
	}
      }
    }
    else {
      fprintf(stderr,"\n");
      fprintf(stderr,"Vertices: %-d %-d\n",MV_ID(ME_Vertex(e,0)),MV_ID(ME_Vertex(e,1)));
    }

    fprintf(stderr,"\n");
  }

  void MF_Print(MFace_ptr f, int lvl) {
    MVertex_ptr v;
    MEdge_ptr e;
    MRegion_ptr r;
    List_ptr fvertices, fedges, fregions;
    int i, edir, idx, ne, nv, nr;
    GEntity_ptr gent;
    double xyz[3];

    if (!f)
      return;

    if (MEnt_Dim(f) != MFACE && MEnt_OrigDim(f) != MFACE) {
      fprintf(stderr,"NOT A FACE\n");
      return;
    }

    fprintf(stderr,"\n");
    fprintf(stderr,"Face: 0x%-x   ID: %-d   ",(unsigned int)f,MF_ID(f));
    
    
    if (lvl) {
      /* Print more detailed information */
      
      if (MF_GEntDim(f) == -1)
	fprintf(stderr,"Unknown Classification\n");
      else {
	fprintf(stderr,"GEntID: %-d    GEntDim: %-d\n", MF_GEntID(f),
		MF_GEntDim(f));
	if ((gent = MF_GEntity(f)))
	  fprintf(stderr,"Model entity pointer: 0x%-x\n",(unsigned int)gent);
      }

      fprintf(stderr,"\n");

      fedges = MF_Edges(f,1,0);
      ne = List_Num_Entries(fedges);
      fprintf(stderr,"Edges: %-d\n",ne);
      fprintf(stderr,"Object        ID      GEntID   GEntDim    Vertex IDs\n");
      idx = 0; i = 0;
      while ((e = List_Next_Entry(fedges,&idx))) {
	edir = MF_EdgeDir_i(f,i);
	if (edir) 
	  fprintf(stderr,"0x%-8x    %-8d %-8d     %-1d       %-d  %-d\n",
		  (unsigned int)e,ME_ID(e),ME_GEntID(e),ME_GEntDim(e),
		  MV_ID(ME_Vertex(e,0)),MV_ID(ME_Vertex(e,1)));
	else 
	  fprintf(stderr,"0x%-8x    %-8d %-8d     %-1d       %-d  %-d\n",
		  (unsigned int)e,-ME_ID(e),ME_GEntID(e),ME_GEntDim(e),
		  MV_ID(ME_Vertex(e,0)),MV_ID(ME_Vertex(e,1)));
	i++;
      }
      fprintf(stderr,"\n");
      List_Delete(fedges);
      
      if (lvl > 1) {
	/* Print more detailed adjacency information as well */
	
	fprintf(stderr,"\n");

	idx = 0;
	fvertices = MF_Vertices(f,1,0);
	nv = List_Num_Entries(fvertices);
	fprintf(stderr,"Vertices: %-d\n",nv);
	fprintf(stderr,"Object        ID      GEntID   GEntDim\n");
	while ((v = List_Next_Entry(fvertices,&idx))) {
	  fprintf(stderr,"0x%-8x    %-8d %-8d     %-1d\n",(unsigned int)v,
		  MV_ID(v),MV_GEntID(v),MV_GEntDim(v));
	  if (lvl > 2) {
	    MV_Coords(v,xyz);
	    fprintf(stderr,"Coordinates: %16.8lf %16.8lf %16.8lf\n",xyz[0],xyz[1],xyz[2]);
	  }
	}
	fprintf(stderr,"\n");
	List_Delete(fvertices);

	
	fregions= MF_Regions(f);
	if (fregions) {
	  nr = List_Num_Entries(fregions);
	  fprintf(stderr,"Regions: %-d\n",nr);
	  
	  fprintf(stderr,"Object        ID      GEntID   GEntDim\n");
	  idx = 0;
	  while ((r = List_Next_Entry(fregions,&idx)))
	    fprintf(stderr,"0x%-8x %-8d %-8d     %-1d\n",(unsigned int)r,
		    MR_ID(r),MR_GEntID(r),MR_GEntDim(r));
	  
	  fprintf(stderr,"\n");
	  List_Delete(fregions);      
	}
      }
    }
    else {
      fprintf(stderr,"\n");
      fprintf(stderr,"Edges: ");
      fedges = MF_Edges(f,1,0);
      ne = List_Num_Entries(fedges);
      for (i = 0; i < ne; i++) {
	e = List_Entry(fedges,i);
	edir = MF_EdgeDir_i(f,i);
	if (edir)
	  fprintf(stderr," %-d ",ME_ID(e));
	else 
	  fprintf(stderr," %-d ",-ME_ID(e));	
      }
      fprintf(stderr,"\n");
      List_Delete(fedges);
    }

    fprintf(stderr,"\n");
  }

  void MR_Print(MRegion_ptr r, int lvl) {
    MVertex_ptr v;
    MEdge_ptr e;
    MFace_ptr f;
    List_ptr rvertices, redges, rfaces;
    int i, fdir, idx, ne, nv, nf;
    GEntity_ptr gent;
    double xyz[3];

    if (!r)
      return;

    if (MEnt_Dim(r) != MREGION && MEnt_OrigDim(r) != MREGION) {
      fprintf(stderr,"NOT A REGION\n");
      return;
    }

    fprintf(stderr,"\n");
    fprintf(stderr,"Region: 0x%-x   ID: %-d   ",(unsigned int)r,MF_ID(r));
    
    
    if (lvl) {
      /* Print more detailed information */
      
      if (MR_GEntDim(r) == -1)
	fprintf(stderr,"Unknown Classification\n");
      else {
	fprintf(stderr,"GEntID: %-d    GEntDim: %-d\n", MR_GEntID(r),
		MR_GEntDim(r));
	if ((gent = MR_GEntity(r)))
	  fprintf(stderr,"Model entity pointer: 0x%-x\n",(unsigned int)gent);
      }

      fprintf(stderr,"\n");

      fprintf(stderr,"Faces:\n");
      fprintf(stderr,"Object        ID      GEntID   GEntDim\n");
      rfaces = MR_Faces(r);
      nf = List_Num_Entries(rfaces);
      for (i = 0; i < nf; i++) {
	f = List_Entry(rfaces,i);
	fdir = MR_FaceDir_i(r,i);
	if (fdir)
	  fprintf(stderr,"0x%-8x    %-8d %-8d     %-1d\n",(unsigned int)f,
		  MF_ID(f),MF_GEntID(f),ME_GEntDim(f));
	else 
	  fprintf(stderr,"0x%-8x    %-8d %-8d     %-1d\n",(unsigned int)f,
		  -MF_ID(f),MF_GEntID(f),ME_GEntDim(f));
      }
      fprintf(stderr,"\n");
      List_Delete(rfaces);
      

      if (lvl > 1) {
	/* Print more detailed adjacency information as well */
	
	fprintf(stderr,"\n");

	redges = MR_Edges(r);
	ne = List_Num_Entries(redges);
	fprintf(stderr,"Edges:%-d\n",ne);
	fprintf(stderr,"Object        ID      GEntID   GEntDim    Vertex IDs\n");
	idx = 0; i = 0;
	while ((e = List_Next_Entry(redges,&idx))) {
	  fprintf(stderr,"0x%-8x    %-8d %-8d     %-1d       %-d  %-d\n",
		  (unsigned int)e,ME_ID(e),ME_GEntID(e),ME_GEntDim(e),
		  MV_ID(ME_Vertex(e,0)),MV_ID(ME_Vertex(e,1)));
	}
	fprintf(stderr,"\n");
	List_Delete(redges);
      
	idx = 0;
	rvertices = MR_Vertices(r);
	nv = List_Num_Entries(rvertices);
	fprintf(stderr,"Vertices:%-d\n",nv);
	fprintf(stderr,"Object        ID      GEntID   GEntDim\n");
	while ((v = List_Next_Entry(rvertices,&idx))) {
	  fprintf(stderr,"0x%-8x    %-8d %-8d     %-1d\n",(unsigned int)v,
		  MV_ID(v),MV_GEntID(v),MV_GEntDim(v));
	  if (lvl > 2) {
	    MV_Coords(v,xyz);
	    fprintf(stderr,"Coordinates: %16.8lf %16.8lf %16.8lf\n",xyz[0],xyz[1],xyz[2]);
	  }
	}
	fprintf(stderr,"\n");

      }
    }
    else {
      fprintf(stderr,"\n");
      fprintf(stderr,"Faces: ");
      rfaces = MR_Faces(r);
      nf = List_Num_Entries(rfaces);
      for (i = 0; i < nf; i++) {
	f = List_Entry(rfaces,i);
	fdir = MR_FaceDir_i(r,i);
	if (fdir)
	  fprintf(stderr," %-d ",MF_ID(f));
	else 
	  fprintf(stderr," %-d ",-MF_ID(f));	
      }
      fprintf(stderr,"\n");
      List_Delete(rfaces);
    }

    fprintf(stderr,"\n");
  }

#ifdef __cplusplus
}
#endif
