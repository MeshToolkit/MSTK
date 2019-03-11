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

  /* 
     This function partition the mesh into submeshes.

     GlobalID, GEntID, and GEntDim are retained in submesh for
     vertex and element.
     
     MasterParID is assigned for both global and submeshes.
     
     Caller should allocate memory for submeshes: 
     Mesh_ptr *submeshes = (Mesh_ptr*) malloc(num*sizeof(Mesh_ptr);

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

    
  int MESH_Partition(Mesh_ptr mesh, int num, int *part, Mesh_ptr *submeshes) {
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

    for (i = 0; i < num; i++) {
      submeshes[i] = MESH_New(rtype);
      MESH_Set_Prtn(submeshes[i],i,num);
    }

    return (*MESH_Surf_Partition_jmp[rtype])(mesh, num, part, submeshes);
  }

  int MESH_Vol_Partition(Mesh_ptr mesh, int num, int *part, Mesh_ptr *submeshes) {
    RepType rtype = MESH_RepType(mesh);
    int i;

    for (i = 0; i < num; i++) {
      submeshes[i] = MESH_New(rtype);
      MESH_Set_Prtn(submeshes[i],i,num);
    }

    return (*MESH_Vol_Partition_jmp[rtype])(mesh, num, part, submeshes);
  }


  int MESH_Surf_Partition_FN(Mesh_ptr mesh, int num, int *part, Mesh_ptr *submeshes) {
    int i, j, ic, idx, idx2, max_nfv;
    int nv, nf, nfv, nfe, mpart_no;
    int gfid=1, geid=1, gvid=1;
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

	fedges = MF_Edges(gmf,1,0);
	nfe = List_Num_Entries(fedges);
	
	for (i = 0; i < nfe; i++) {
          int found_lme = 0;

	  gme = List_Entry(fedges,i);
	  MEnt_Get_AttVal(gme,g2latt,0,0,&lmelist);
	  if (!lmelist) {
	    lmelist = List_New(0);
	    MEnt_Set_AttVal(gme,g2latt,0,0,lmelist);
	  }
          else {
            idx2 = 0;
            while ((lme = List_Next_Entry(lmelist,&idx2))) {
              if (ME_Mesh(lme) == submeshes[part_no]) {
                found_lme = 1;
                break;
              }
            }
          }

	  if (!found_lme) { 
	    lme = ME_New(submeshes[part_no]);
	    ME_Set_GEntID(lme,ME_GEntID(gme));
	    ME_Set_GEntDim(lme,ME_GEntDim(gme));
	    mpart_no = ME_MasterParID(gme);
	    if (mpart_no == 0) { /* not set */
	      ME_Set_MasterParID(lme,part_no);
	      ME_Set_MasterParID(gme,part_no+1); 
              if (!ME_GlobalID(gme)) ME_Set_GlobalID(gme,geid++);
	    }
	    else
	      ME_Set_MasterParID(lme,mpart_no-1);
	    ME_Set_GlobalID(lme,ME_GlobalID(gme));

	    for (j = 0; j < 2; j++) {
              int found_lmv = 0;

	      gmv = ME_Vertex(gme,j);

	      MEnt_Get_AttVal(gmv,g2latt,0,0,&lmvlist);
	      if (!lmvlist) {
                lmvlist = List_New(0);
                MEnt_Set_AttVal(gmv,g2latt,0,0,lmvlist);
              }
              else {
                idx2 = 0;
                while ((lmv = List_Next_Entry(lmvlist,&idx2))) {
                  if (MV_Mesh(lmv) == submeshes[part_no]) {
                    found_lmv = 1;
                    break;
                  }
                }
              }

	      if (!found_lmv) {
                lmv = MV_New(submeshes[part_no]);
                MV_Coords(gmv,xyz);
                MV_Set_Coords(lmv,xyz);
                MV_Set_GEntID(lmv,MV_GEntID(gmv));
                MV_Set_GEntDim(lmv,MV_GEntDim(gmv));
                mpart_no = MV_MasterParID(gmv);
                if (mpart_no == 0) {
                  MV_Set_MasterParID(lmv,part_no);
                  MV_Set_MasterParID(gmv,part_no+1); /* will subract 1 at the end */
                  if (!MV_GlobalID(gmv)) MV_Set_GlobalID(gmv,gvid++);	    
                }
                else
                  MV_Set_MasterParID(lmv,mpart_no-1);
                MV_Set_GlobalID(lmv,MV_GlobalID(gmv));	    
                MEnt_Set_AttVal(lmv,l2gatt,0,0,gmv); 
                List_Add(lmvlist,lmv);
              } /* if (!found_lmv) */

	      ME_Set_Vertex(lme,j,lmv);	      
	    } /* for (j = 0; j < 2; j++) */
	    
	    MEnt_Set_AttVal(lme,l2gatt,0,0,gme);
	    List_Add(lmelist,lme);

	  } /* if !found_lme */

	  lfedges[i] = lme;
	  lfedirs[i] = MF_EdgeDir_i(gmf,i);	

	} /* for (i = 0; i < nfe; i++) */
        List_Delete(fedges);

	lmf = MF_New(submeshes[part_no]);
	MF_Set_GEntID(lmf,MF_GEntID(gmf));
	MF_Set_GEntDim(lmf,MF_GEntDim(gmf));
        //	MF_Set_PType(lmf,PINTERIOR);
	MF_Set_MasterParID(lmf,part_no);
	MF_Set_MasterParID(gmf,part_no);
	if (!MF_GlobalID(gmf)) MF_Set_GlobalID(gmf,gfid++);
	MF_Set_GlobalID(lmf,MF_GlobalID(gmf));
	MF_Set_Edges(lmf,nfe,lfedges,lfedirs);

	lmflist = List_New(0);
	List_Add(lmflist,lmf);
	MEnt_Set_AttVal(gmf,g2latt,0,0,lmflist);
	MEnt_Set_AttVal(lmf,l2gatt,0,0,gmf);

      } /* while (gmf = ....) */

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
    int i, j, ic, idx, idx2, mpart_no;
    int nv, nr, nrv, nre, nrf, nfe;
    int grid=1, gfid=1, geid=1, gvid=1;
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
	
	rfaces = MR_Faces(gmr);
	nrf = List_Num_Entries(rfaces);
	
	for (i = 0; i < nrf; i++) {
          int found_lmf = 0;

	  gmf = List_Entry(rfaces,i);
	  MEnt_Get_AttVal(gmf,g2latt,0,0,&lmflist);
          
	  if (!lmflist) {
	    lmflist = List_New(0);
	    MEnt_Set_AttVal(gmf,g2latt,0,0,lmflist);
	  }
          else {
            idx2 = 0;
            while ((lmf = List_Next_Entry(lmflist,&idx2))) {
              if (MF_Mesh(lmf) == submeshes[part_no]) {
                found_lmf = 1;
                break;
              }
            }
          }

	  if (!found_lmf) {
	    lmf = MF_New(submeshes[part_no]);
	    MF_Set_GEntID(lmf,MF_GEntID(gmf));
	    MF_Set_GEntDim(lmf,MF_GEntDim(gmf));
	    mpart_no = MF_MasterParID(gmf);
	    if (mpart_no == 0) { /* not set */
	      MF_Set_MasterParID(lmf,part_no);
	      MF_Set_MasterParID(gmf,part_no+1);
              if (!MF_GlobalID(gmf)) MF_Set_GlobalID(gmf,gfid++);
	    }
	    else
	      MF_Set_MasterParID(lmf,mpart_no-1);
	    MF_Set_GlobalID(lmf,MF_GlobalID(gmf));


	    fedges = MF_Edges(gmf,1,0);
	    nfe = List_Num_Entries(fedges);
	    for (j = 0; j < nfe; j++) {
              int found_lme = 0;

	      gme = List_Entry(fedges,j);
              
	      MEnt_Get_AttVal(gme,g2latt,0,0,&lmelist);
	      if (!lmelist) {
                lmelist = List_New(0);
                MEnt_Set_AttVal(gme,g2latt,0,0,lmelist);
              }
              else {
                idx2 = 0;
                while ((lme = List_Next_Entry(lmelist,&idx2))) {
                  if (ME_Mesh(lme) == submeshes[part_no]) {
                    found_lme = 1;
                    break;
                  }
                }
              }

              if (!found_lme) {
                lme = ME_New(submeshes[part_no]);
                ME_Set_GEntID(lme,ME_GEntID(gme));
                ME_Set_GEntDim(lme,ME_GEntDim(gme));
                mpart_no = ME_MasterParID(gme);
                if (mpart_no == 0) { /* not set */
                  ME_Set_MasterParID(lme,part_no);
                  ME_Set_MasterParID(gme,part_no+1);
                  if (!ME_GlobalID(gme)) ME_Set_GlobalID(gme,geid++);
                }
                else
                  ME_Set_MasterParID(lme,mpart_no-1);
                ME_Set_GlobalID(lme,ME_GlobalID(gme));
                
                int k;
                for (k = 0; k < 2; k++) {
                  int found_lmv = 0;

                  gmv = ME_Vertex(gme,k);
                  
                  MEnt_Get_AttVal(gmv,g2latt,0,0,&lmvlist);
                  if (!lmvlist) {
                    lmvlist = List_New(0);
                    MEnt_Set_AttVal(gmv,g2latt,0,0,lmvlist);
                  }
                  else {
                    /* Check if we created a copy of this global vertex in
                       _this_ submesh */
                    idx2 = 0;
                    while ((lmv = List_Next_Entry(lmvlist,&idx2))) {
                      if (MV_Mesh(lmv) == submeshes[part_no]) {
                        found_lmv = 1;
                        break;
                      }
                    }
                  }

                  if (!found_lmv) {
                    lmv = MV_New(submeshes[part_no]);
                    MV_Coords(gmv,xyz);
                    MV_Set_Coords(lmv,xyz);
                    MV_Set_GEntID(lmv,MV_GEntID(gmv));
                    MV_Set_GEntDim(lmv,MV_GEntDim(gmv));
                    mpart_no = MV_MasterParID(gmv);
                    if (mpart_no == 0) { /* Not set */
                      MV_Set_MasterParID(lmv,part_no);
                      MV_Set_MasterParID(gmv,part_no+1);
                      if (!MV_GlobalID(gmv)) MV_Set_GlobalID(gmv,gvid++);
                    }
                    else
                      MV_Set_MasterParID(lmv,mpart_no-1);
                    MV_Set_GlobalID(lmv,MV_GlobalID(gmv));
                    MEnt_Set_AttVal(lmv,l2gatt,0,0,gmv);
                    List_Add(lmvlist,lmv);
                  } /* if !found_lmv */

                  ME_Set_Vertex(lme,k,lmv);
                }
                MEnt_Set_AttVal(lme,l2gatt,0,0,gme);
                List_Add(lmelist,lme);                
              } /* if !found_lme */

              lfedges[j] = lme;
              lfedirs[j] = MF_EdgeDir_i(gmf,j);

            } /* for (j = 0; j < nfe; j++) */

            List_Delete(fedges);

            MF_Set_Edges(lmf,nfe,lfedges,lfedirs);
            MEnt_Set_AttVal(lmf,l2gatt,0,0,gmf);
            List_Add(lmflist,lmf);
          } /* if !found_lmf */

          lrfaces[i] = lmf;
          lrfdirs[i] = MR_FaceDir_i(gmr,i);

        } /* for (i = 0; i < nrf; i++) */   

        List_Delete(rfaces);
	
	lmr = MR_New(submeshes[part_no]);
	MR_Set_GEntID(lmr,MR_GEntID(gmr));
	MR_Set_GEntDim(lmr,MR_GEntDim(gmr));
        //	MR_Set_PType(lmr,PINTERIOR);
	MR_Set_MasterParID(lmr,part_no);
	MR_Set_MasterParID(gmr,part_no);
        if (!MR_GlobalID(gmr)) MR_Set_GlobalID(gmr,grid++);
	MR_Set_GlobalID(lmr,MR_GlobalID(gmr));
        grid++;
	MR_Set_Faces(lmr,nrf,lrfaces,lrfdirs);

	lmrlist = List_New(0);
	List_Add(lmrlist,lmr);
	MEnt_Set_AttVal(gmr,g2latt,0,0,lmrlist);
	MEnt_Set_AttVal(lmr,l2gatt,0,0,gmr);
      } /* while (gmr ....) */

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
    return 0;
  }

    
  int MESH_Vol_Partition_R1R2(Mesh_ptr mesh, int num, int *part, Mesh_ptr *submeshes) {
    MSTK_Report("MESH_Vol_Partition_R1R2","Not implemented",MSTK_FATAL);
    return 0;
  }


  int MESH_Vol_Partition_R4(Mesh_ptr mesh, int num, int *part, Mesh_ptr *submeshes) {
    MSTK_Report("MESH_Vol_Partition_R4","Not implemented",MSTK_FATAL);
    return 0;
  }


  
#ifdef __cplusplus
}
#endif

