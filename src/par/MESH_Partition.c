#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif

  /* 
     This function partition the mesh into submeshes.

     GlobalID, GEntID, and GEntDim are retained in submesh for
     vertex and element.
     
     MasterParID is assigned for both global and submeshes.
     
     Caller should allocate memory for submeshes: 
     Mesh_ptr *submeshes = (Mesh_ptr*) MSTK_malloc(num*sizeof(Mesh_ptr);

     Author(s): Duo Wang, Rao Garimella

  */



  int MESH_Surf_Partition(Mesh_ptr mesh, int num, int *part, Mesh_ptr *submeshes);
  int MESH_Vol_Partition(Mesh_ptr mesh, int num, int *part, Mesh_ptr *submeshes);


  int MESH_Surf_Partition_FN(Mesh_ptr mesh, int num, int *part, Mesh_ptr *submeshes);
  int MESH_Vol_Partition_FN(Mesh_ptr mesh, int num, int *part, Mesh_ptr *submeshes);
  int MESH_Surf_Partition_R1R2R4(Mesh_ptr mesh, int num, int *part, Mesh_ptr *submeshes);
  int MESH_Vol_Partition_R1R2(Mesh_ptr mesh, int num, int *part, Mesh_ptr *submeshes);
  int MESH_Vol_Partition_R4(Mesh_ptr mesh, int num, int *part, Mesh_ptr *submeshes);


  static int (*MESH_Surf_Partition_jmp[MSTK_MAXREP])(Mesh_ptr mesh, int num,
						     int *part, 
						     Mesh_ptr *submeshes) =
  {MESH_Surf_Partition_FN, MESH_Surf_Partition_FN, MESH_Surf_Partition_R1R2R4,
   MESH_Surf_Partition_R1R2R4, MESH_Surf_Partition_R1R2R4};

  static int (*MESH_Vol_Partition_jmp[MSTK_MAXREP])(Mesh_ptr mesh, int num,
						    int *part, 
						    Mesh_ptr *submeshes) =
  {MESH_Vol_Partition_FN, MESH_Vol_Partition_FN, MESH_Vol_Partition_R1R2,
   MESH_Vol_Partition_R1R2, MESH_Vol_Partition_R4};

  int MESH_Get_Partition(Mesh_ptr mesh, int num, int **part, int method, int rank, MPI_Comm comm) {
    int ok, nf, nr, ncells;
    /* basic mesh information */
    if ( rank == 0 ) {
      nf = MESH_Num_Faces(mesh);
      nr = MESH_Num_Regions(mesh);

      ncells = (nr) ? nr:nf; 
      if (ncells == 0) {
	MSTK_Report("MESH_Get_Partition_Serial",
		    "This is not a valid mstk file for partition",MSTK_ERROR);
	exit(-1);
      }
      *part = (int *) MSTK_malloc(ncells*sizeof(int));
      if ( method == 0)
	ok = MESH_PartitionWithMetis(mesh, num, part);
    }
    if ( method == 1 )
      ok = MESH_PartitionWithZoltan(mesh, num, part,rank,comm);
    return ok;
  }
    
  int MESH_Partition(Mesh_ptr mesh, Mesh_ptr *submeshes, int num, int *part) {
    int i, j, ok, nf, nr, idx, idx2, found, natt, ncomp;
    MVertex_ptr gmv;
    MEdge_ptr gme;
    MFace_ptr gmf;
    MRegion_ptr gmr;
    char attname[256];
    MType mtype;
    MAttType att_type;
    MAttrib_ptr attrib;

    /* basic mesh information */
    nf = MESH_Num_Faces(mesh);
    nr = MESH_Num_Regions(mesh);

    if (!nr && !nf) {
      MSTK_Report("MESH_Partition",
		  "This is not a valid mstk file for partition",MSTK_ERROR);
      exit(-1);
    }

    /* assign global ID equal to ID for global mesh */
    /* and master partition ID of -1                */
    idx = 0;
    while ((gmv = MESH_Next_Vertex(mesh,&idx)))
      MV_Set_GlobalID(gmv,MV_ID(gmv));

    idx = 0;
    while ((gme = MESH_Next_Edge(mesh,&idx)))
      ME_Set_GlobalID(gme,ME_ID(gme));

    idx = 0;
    while ((gmf = MESH_Next_Face(mesh,&idx)))
      MF_Set_GlobalID(gmf,MF_ID(gmf));

    idx = 0;  
    while ((gmr = MESH_Next_Region(mesh,&idx)))
      MR_Set_GlobalID(gmr,MR_ID(gmr));
    
    if (nr) {
      ok = MESH_Vol_Partition(mesh,num,part,submeshes);
    }
    else if (nf) {
      ok = MESH_Surf_Partition(mesh,num,part,submeshes);
    }
    else {
      MSTK_Report("MESH_Partition",
		  "This is not a valid mstk file for partition",MSTK_ERROR);
      exit(-1);
    }


    /* assign attrib list */

    /* We should probably do this outside this routine */

    natt = MESH_Num_Attribs(mesh);
    for (j = 0; j < natt; j++) {
      attrib = MESH_Attrib(mesh,j);
      MAttrib_Get_Name(attrib,attname);
      att_type = MAttrib_Get_Type(attrib);
      ncomp = MAttrib_Get_NumComps(attrib);
      mtype = MAttrib_Get_EntDim(attrib);
      for (i = 0; i < num; i++) {
	if (ncomp == 1)
	  attrib =  MAttrib_New(submeshes[i], attname, att_type, mtype);
	else
	  attrib =  MAttrib_New(submeshes[i], attname, att_type, mtype, ncomp);
      }
    }

    return ok;
  }


  int MESH_Surf_Partition(Mesh_ptr mesh, int num, int *part, Mesh_ptr *submeshes) {
    RepType rtype = MESH_RepType(mesh);
    int i;

    for (i = 0; i < num; i++) submeshes[i] = MESH_New(rtype);

    (*MESH_Surf_Partition_jmp[rtype])(mesh, num, part, submeshes);
  }

  int MESH_Vol_Partition(Mesh_ptr mesh, int num, int *part, Mesh_ptr *submeshes) {
    RepType rtype = MESH_RepType(mesh);
    int i;

    for (i = 0; i < num; i++) submeshes[i] = MESH_New(rtype);

    (*MESH_Vol_Partition_jmp[rtype])(mesh, num, part, submeshes);
  }


  int MESH_Surf_Partition_FN(Mesh_ptr mesh, int num, int *part, Mesh_ptr *submeshes) {
    int i, j, ic, idx, idx2, max_nfv, found;
    int nv, nf, nfv, nfe, mpart_no;
    double xyz[3];
    MVertex_ptr lmv, gmv;
    MEdge_ptr lme, gme;
    MFace_ptr lmf, gmf;
    List_ptr fverts, fedges, lmvlist, lmelist, lmflist;
    MAttrib_ptr g2latt, l2gatt;

    /* basic mesh information */
    RepType rtype = MESH_RepType(mesh);
    nv = MESH_Num_Vertices(mesh);
    nf = MESH_Num_Faces(mesh);


    int part_no;

    /* each face belongs only to one partition, so no need to mark  */
    int mkvid = MSTK_GetMarker();

    g2latt = MESH_AttribByName(mesh,"Global2Local");


    MEdge_ptr *lfedges = (MEdge_ptr *) malloc(MAXPV2*sizeof(MEdge_ptr));
    int *lfedirs = (int *) malloc(MAXPV2*sizeof(int));

    /* for each partition */

    /* In the following we set the master partion ID of global
       entities 1 above the partition number so that we can
       distinguish between no data and part number 0 - Will dial back
       the partition ID by 1 at the end. Don't need to do that for
       local entities because they are not being checked, only being
       set. */

    for (part_no = 0; part_no < num; part_no++) {

      l2gatt = MAttrib_New(submeshes[part_no],"Local2Global",POINTER,MALLTYPE);

      idx = 0; ic = 0;
      while ((gmf = MESH_Next_Face(mesh,&idx))) {
	if (part[ic++] != part_no)
	  continue;

	/* Duplicate the face */

	fverts = MF_Vertices(gmf,1,0);
	nfv = List_Num_Entries(fverts);
	
	for (i = 0; i < nfv; i++) {
	  gmv = List_Entry(fverts,i);
	  MEnt_Get_AttVal(gmv,g2latt,0,0,&lmvlist);
	  if (!lmvlist) {
	    lmvlist = List_New(0);
	    MEnt_Set_AttVal(gmv,g2latt,0,0,lmvlist);
	  }

	  /* Check if we created a copy of this global vertex
	     in _this_ submesh */

	  found = 0;
	  idx2 = 0;
	  while ((lmv = List_Next_Entry(lmvlist,&idx2))) {
	    if (MV_Mesh(lmv) == submeshes[part_no]) {
	      found = 1;
	      break;
	    }
	  }

	  if (!found) {
	    lmv = MV_New(submeshes[part_no]);
	    MV_Coords(gmv,xyz);
	    MV_Set_Coords(lmv,xyz);
	    MV_Set_GEntID(lmv,MV_GEntID(gmv));
	    MV_Set_GEntDim(lmv,MV_GEntDim(gmv));
	    mpart_no = MV_MasterParID(gmv);
	    if (mpart_no == 0) {
	      MV_Set_MasterParID(lmv,part_no);
	      MV_Set_MasterParID(gmv,part_no+1); /* will subract 1 at the end */
	    }
	    else
	      MV_Set_MasterParID(lmv,mpart_no-1);
	    MV_Set_GlobalID(lmv,MV_ID(gmv));	    
	    MEnt_Set_AttVal(lmv,l2gatt,0,0,gmv); 
	    List_Add(lmvlist,lmv);
	  }
	}

	fedges = MF_Edges(gmf,1,0);
	nfe = List_Num_Entries(fedges);
	
	for (i = 0; i < nfe; i++) {
	  gme = List_Entry(fedges,i);
	  MEnt_Get_AttVal(gme,g2latt,0,0,&lmelist);
	  if (!lmelist) {
	    lmelist = List_New(0);
	    MEnt_Set_AttVal(gme,g2latt,0,0,lmelist);
	  }

	  idx2 = 0;
	  found = 0;
	  while ((lme = List_Next_Entry(lmelist,&idx2))) {
	    if (ME_Mesh(lme) == submeshes[part_no]) {
	      found = 1;
	      break;
	    }
	  }

	  if (!found) { 
	    lme = ME_New(submeshes[part_no]);
	    ME_Set_GEntID(lme,ME_GEntID(gme));
	    ME_Set_GEntDim(lme,ME_GEntDim(gme));
	    mpart_no = ME_MasterParID(gme);
	    if (mpart_no == 0) { /* not set */
	      ME_Set_MasterParID(lme,part_no);
	      ME_Set_MasterParID(gme,part_no+1); 
	    }
	    else
	      ME_Set_MasterParID(lme,mpart_no-1);
	    ME_Set_GlobalID(lme,ME_ID(gme));

	    for (j = 0; j < 2; j++) {
	      gmv = ME_Vertex(gme,j);

	      MEnt_Get_AttVal(gmv,g2latt,0,0,&lmvlist);
	      if (!lmvlist)
		MSTK_Report("MESH_Surf_AddGhost_FN","Missing local vertex list",MSTK_FATAL);

	      idx2 = 0;
	      found = 0;
	      while ((lmv = List_Next_Entry(lmvlist,&idx2))) {
		if (MV_Mesh(lmv) == submeshes[part_no]) {
		  found = 1;
		  break;
		}
	      }

	      if (!lmv)
		MSTK_Report("MESH_Surf_AddGhost_FN","Missing local vertex",MSTK_FATAL);
	      ME_Set_Vertex(lme,j,lmv);	      
	    }	    
	    MEnt_Set_AttVal(lme,l2gatt,0,0,gme);
	    List_Add(lmelist,lme);
	  }
	  lfedges[i] = lme;
	  lfedirs[i] = MF_EdgeDir_i(gmf,i);	
	}

	lmf = MF_New(submeshes[part_no]);
	MF_Set_GEntID(lmf,MF_GEntID(gmf));
	MF_Set_GEntDim(lmf,MF_GEntDim(gmf));
	MF_Set_PType(lmf,PINTERIOR);
	MF_Set_MasterParID(lmf,part_no);
	MF_Set_MasterParID(gmf,part_no);
	MF_Set_GlobalID(lmf,MF_ID(gmf));
	MF_Set_Edges(lmf,nfe,lfedges,lfedirs);

	lmflist = List_New(0);
	List_Add(lmflist,lmf);
	MEnt_Set_AttVal(gmf,g2latt,0,0,lmflist);
	MEnt_Set_AttVal(lmf,l2gatt,0,0,gmf);
      }

    } /* for (part_no = 0.... */

    free(lfedges);
    free(lfedirs);


    idx = 0;
    while ((gmv = MESH_Next_Vertex(mesh,&idx)))
      MV_Set_MasterParID(gmv,MV_MasterParID(gmv)-1);

    idx = 0;
    while ((gme = MESH_Next_Edge(mesh,&idx)))
      ME_Set_MasterParID(gme,ME_MasterParID(gme)-1);

    return 1;
  }

    
  int MESH_Vol_Partition_FN(Mesh_ptr mesh, int num, int *part, Mesh_ptr *submeshes) {
    int i, j, ic, idx, idx2, mpart_no, found;
    int nv, nr, nrv, nre, nrf, nfe;
    double xyz[3];
    MVertex_ptr lmv, gmv;
    MEdge_ptr lme, gme;
    MFace_ptr lmf, gmf;
    MRegion_ptr lmr, gmr;
    List_ptr rverts, redges, rfaces, fedges, lmvlist, lmelist, lmflist, lmrlist;
    MAttrib_ptr l2gatt, g2latt;

    /* basic mesh information */
    RepType rtype = MESH_RepType(mesh);
    nv = MESH_Num_Vertices(mesh);
    nr = MESH_Num_Regions(mesh);

    int part_no;

    g2latt = MESH_AttribByName(mesh,"Global2Local");

    MFace_ptr *lrfaces = (MFace_ptr *) malloc(MAXPF3*sizeof(MFace_ptr));
    int *lrfdirs = (int *) malloc(MAXPF3*sizeof(int));
    MEdge_ptr *lfedges = (MEdge_ptr *) malloc(MAXPV2*sizeof(MEdge_ptr));
    int *lfedirs = (int *) malloc(MAXPV2*sizeof(int));

    /* for each partition */
    for (part_no = 0; part_no < num; part_no++) {

      l2gatt = MAttrib_New(submeshes[part_no],"Local2Global",POINTER,MALLTYPE);

      idx = 0; ic = 0;
      while ((gmr = MESH_Next_Region(mesh,&idx))) {
	if (part[ic++] != part_no)
	  continue;

	rverts = MR_Vertices(gmr);
	nrv = List_Num_Entries(rverts);
	
	for (i = 0; i < nrv; i++) {
	  gmv = List_Entry(rverts,i);
	  MEnt_Get_AttVal(gmv,g2latt,0,0,&lmvlist);
	  if (!lmvlist) {
	    lmvlist = List_New(0);
	    MEnt_Set_AttVal(gmv,g2latt,0,0,lmvlist);
	  }

	  /* Check if we created a copy of this global vertex in
	     _this_ submesh */
	  idx2 = 0;
	  found = 0;
	  while ((lmv = List_Next_Entry(lmvlist,&idx2))) {
	    if (MV_Mesh(lmv) == submeshes[part_no]) {
	      found = 1;
	      break;
	    }
	  }

	  if (!found) {
	    lmv = MV_New(submeshes[part_no]);
	    MV_Coords(gmv,xyz);
	    MV_Set_Coords(lmv,xyz);
	    MV_Set_GEntID(lmv,MV_GEntID(gmv));
	    MV_Set_GEntDim(lmv,MV_GEntDim(gmv));
	    mpart_no = MV_MasterParID(gmv);
	    if (mpart_no == 0) { /* Not set */
	      MV_Set_MasterParID(lmv,part_no);
	      MV_Set_MasterParID(gmv,part_no+1);
	    }
	    else
	      MV_Set_MasterParID(lmv,mpart_no-1);
	    MV_Set_GlobalID(lmv,MV_ID(gmv));
	    MEnt_Set_AttVal(lmv,l2gatt,0,0,gmv);
	    List_Add(lmvlist,lmv);
	  }
	}
	List_Delete(rverts);

	redges = MR_Edges(gmr);
	nre = List_Num_Entries(redges);
	
	for (i = 0; i < nre; i++) {
	  gme = List_Entry(redges,i);
	  MEnt_Get_AttVal(gme,g2latt,0,0,&lmelist);

	  if (!lmelist) {
	    lmelist = List_New(0);
	    MEnt_Set_AttVal(gme,g2latt,0,0,lmelist);
	  }

	  idx2 = 0;
	  found = 0;
	  while ((lme = List_Next_Entry(lmelist,&idx2))) {
	    if (ME_Mesh(lme) == submeshes[part_no]) {
	      found = 1;
	      break;
	    }
	  }

	  if (!found) {
	    lme = ME_New(submeshes[part_no]);
	    ME_Set_GEntID(lme,ME_GEntID(gme));
	    ME_Set_GEntDim(lme,ME_GEntDim(gme));
	    mpart_no = ME_MasterParID(gme);
	    if (mpart_no == 0) { /* not set */
	      ME_Set_MasterParID(lme,part_no);
	      ME_Set_MasterParID(gme,part_no+1);
	    }
	    else
	      ME_Set_MasterParID(lme,mpart_no-1);
	    ME_Set_GlobalID(lme,ME_ID(gme));

	    for (j = 0; j < 2; j++) {
	      gmv = ME_Vertex(gme,j);

	      MEnt_Get_AttVal(gmv,g2latt,0,0,&lmvlist);
	      if (!lmvlist)
		MSTK_Report("MESH_Vol_AddGhost_FN","Missing local vertex list",MSTK_FATAL);

	      idx2 = 0;
	      found = 0;
	      while ((lmv = List_Next_Entry(lmvlist,&idx2))) {
		if (MV_Mesh(lmv) == submeshes[part_no]) {
		  found = 1;
		  break;
		}
	      }

	      if (!lmv)
		MSTK_Report("MESH_Vol_AddGhost_FN","Missing submesh vertex",MSTK_FATAL);
	      ME_Set_Vertex(lme,j,lmv);
	    }
	    MEnt_Set_AttVal(lme,l2gatt,0,0,gme);
	    List_Add(lmelist,lme);
	  }
	}
	List_Delete(redges);
	
	rfaces = MR_Faces(gmr);
	nrf = List_Num_Entries(rfaces);
	
	for (i = 0; i < nrf; i++) {
	  gmf = List_Entry(rfaces,i);
	  MEnt_Get_AttVal(gmf,g2latt,0,0,&lmflist);

	  if (!lmflist) {
	    lmflist = List_New(0);
	    MEnt_Set_AttVal(gmf,g2latt,0,0,lmflist);
	  }

	  idx2 = 0;
	  found = 0;
	  while ((lmf = List_Next_Entry(lmflist,&idx2))) {
	    if (MF_Mesh(lmf) == submeshes[part_no]) {
	      found = 1;
	      break;
	    }
	  }

	  if (!found) {
	    lmf = MF_New(submeshes[part_no]);
	    MF_Set_GEntID(lmf,MF_GEntID(gmf));
	    MF_Set_GEntDim(lmf,MF_GEntDim(gmf));
	    mpart_no = MF_MasterParID(gmf);
	    if (mpart_no == 0) { /* not set */
	      MF_Set_MasterParID(lmf,part_no);
	      MF_Set_MasterParID(gmf,part_no+1);
	    }
	    else
	      MF_Set_MasterParID(lmf,mpart_no-1);
	    MF_Set_GlobalID(lmf,MF_ID(gmf));

	    fedges = MF_Edges(gmf,1,0);
	    nfe = List_Num_Entries(fedges);
	    for (j = 0; j < nfe; j++) {
	      gme = List_Entry(fedges,j);

	      MEnt_Get_AttVal(gme,g2latt,0,0,&lmelist);
	      if (!lmelist)
		MSTK_Report("MESH_Vol_Partition_FN","Missing local edge list",MSTK_FATAL);

	      idx2 = 0;
	      found = 0;
	      while ((lme = List_Next_Entry(lmelist,&idx2))) {
		if (ME_Mesh(lme) == submeshes[part_no]) {
		  found = 1;
		  break;
		}
	      }

	      if (!lme)
		MSTK_Report("MESH_Vol_Partition_FN","Missing local edge",MSTK_FATAL);

	      lfedges[j] = lme;
	      lfedirs[j] = MF_EdgeDir_i(gmf,j);
	    }
	    MF_Set_Edges(lmf,nfe,lfedges,lfedirs);
	  }
	  lrfaces[i] = lmf;
	  lrfdirs[i] = MR_FaceDir_i(gmr,i);
	  MEnt_Set_AttVal(lmf,l2gatt,0,0,gmf);
	  List_Add(lmflist,lmf);
	}
	
	
	lmr = MR_New(submeshes[part_no]);
	MR_Set_GEntID(lmr,MR_GEntID(gmr));
	MR_Set_GEntDim(lmr,MR_GEntDim(gmr));
	MR_Set_PType(lmr,PINTERIOR);
	MR_Set_MasterParID(lmr,part_no);
	MR_Set_MasterParID(gmr,part_no);
	MR_Set_GlobalID(lmr,MR_ID(gmr));
	MR_Set_Faces(lmr,nrf,lrfaces,lrfdirs);

	lmrlist = List_New(0);
	List_Add(lmrlist,lmr);
	MEnt_Set_AttVal(gmr,g2latt,0,0,lmrlist);
	MEnt_Set_AttVal(lmr,l2gatt,0,0,gmr);
      }

    } /* for (part_no = 0... */

    free(lrfaces);
    free(lrfdirs);
    free(lfedges);
    free(lfedirs);


    /* Dial back partition ID numbers on global mesh by 1 so that it
       they start at 0 */

    idx = 0;
    while ((gmv = MESH_Next_Vertex(mesh,&idx)))
      MV_Set_MasterParID(gmv,MV_MasterParID(gmv)-1);

    idx = 0;
    while ((gme = MESH_Next_Edge(mesh,&idx)))
      ME_Set_MasterParID(gme,ME_MasterParID(gme)-1);

    idx = 0;
    while ((gmf = MESH_Next_Face(mesh,&idx)))
      MF_Set_MasterParID(gmf,MF_MasterParID(gmf)-1);

    return 1;
  }



  int MESH_Surf_Partition_R1R2R4(Mesh_ptr mesh, int num, int *part, Mesh_ptr *submeshes) {
    MSTK_Report("MESH_Surf_Partition_R1R2R4","Not implemented",MSTK_FATAL);
  }

    
  int MESH_Vol_Partition_R1R2(Mesh_ptr mesh, int num, int *part, Mesh_ptr *submeshes) {
    MSTK_Report("MESH_Vol_Partition_R1R2","Not implemented",MSTK_FATAL);
  }


  int MESH_Vol_Partition_R4(Mesh_ptr mesh, int num, int *part, Mesh_ptr *submeshes) {
    MSTK_Report("MESH_Vol_Partition_R4","Not implemented",MSTK_FATAL);
  }


  
#ifdef __cplusplus
}
#endif

