#define _H_MRegion_Private

#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  int MR_Num_Faces_R1R2(MRegion_ptr r) {
    MRegion_DownAdj_R1R2 *downadj;
    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    
    if (downadj->fvtemplate) {
      return downadj->fvtemplate[0][0];
    }
    else {
      int nv = List_Num_Entries(downadj->rvertices);
      MRType regtype = MSTK_nv2rtype[nv];
      return MSTK_nrf_template[regtype];
    }
  }

  List_ptr MR_Faces_R1R2(MRegion_ptr r) {
    MRegion_DownAdj_R1R2 *downadj;
    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    int i, j, k, nf, nfv;
    List_ptr rfaces=NULL;
    MFace_ptr face;
    Mesh_ptr mesh = MEnt_Mesh(r);
    
    if (downadj->fvtemplate) {
      MVertex_ptr fvertices[MAXPV2];

      nf = downadj->fvtemplate[0][0];
      rfaces = List_New(nf);

      for (i = 0; i < nf; i++) {
	nfv = downadj->fvtemplate[1+i][0];
	
	for (j = 0; j < nfv; j++) {
	  k = downadj->fvtemplate[1+i][1+j];
	  fvertices[j] = List_Entry(downadj->rvertices,k);
	}

	face = MF_New(mesh);
	MEnt_Set_Volatile(face);
	MF_Set_Vertices(face,nfv,fvertices);
	MF_Set_GInfo_Auto(face);
	List_Add(rfaces,face);	
      }
    }
    else {
      MVertex_ptr fvertices[4];
      int nv = List_Num_Entries(downadj->rvertices);

      MRType regtype = MSTK_nv2rtype[nv];

      nf = MSTK_nrf_template[regtype];
      rfaces = List_New(nf);

      for (i = 0; i < nf; i++) {
	nfv = MSTK_rfv_template[regtype][i][0];

	for (j = 0; j < nfv; j++) {
	  k = MSTK_rfv_template[regtype][i][1+j];
	  fvertices[j] = List_Entry(downadj->rvertices,k);
	}

	face = MF_New(mesh);
	MEnt_Set_Volatile(face);
	MF_Set_Vertices(face,nfv,fvertices);
	MF_Set_GInfo_Auto(face);
	List_Add(rfaces,face);
      }
    }

    return rfaces;
  }

  int MR_FaceDir_R1R2(MRegion_ptr r, MFace_ptr f) {
    MRegion_DownAdj_R1R2 *downadj;
    downadj = (MRegion_DownAdj_R1R2 *) r->downadj;
    int fvids[MAXPV2], nfv, nfv2, idx, i, j, k, dir=-1, fnd, same, nf;
    List_ptr fverts;
    MVertex_ptr v;

    fverts = MF_Vertices(f,1,0);
    nfv = List_Num_Entries(fverts);

    idx = 0; i = 0;
    while ((v = List_Next_Entry(fverts,&idx)))
      fvids[i++] = MV_ID(v);

    List_Delete(fverts);

    /* Check if a face in the region has the same vertices as the
       given face. Then, assuming that the region-face-vertex
       templates contain face description such that the normal points
       out of the region, we have to check if the vertices of the
       given face are in the same or opposite direction to the
       matching face in region */

    if (downadj->fvtemplate) {

      nf = downadj->fvtemplate[0][0];

      for (i = 0; i < nf; i++) {
	nfv2 = downadj->fvtemplate[1+i][0];
	if (nfv != nfv2)
	  continue;

	/* Are the first pair of vertices of the given face in this
	   face of the region and are they in the same or reverse
	   order */

	for (j = 0, fnd = 0, k = -1; !fnd && (j < nfv); j++)
	  if (downadj->fvtemplate[1+i][1+j] == fvids[0]) {
	    if (downadj->fvtemplate[1+i][1+(j+1)%nfv] == fvids[1]) {
	      fnd = 1;
	      dir = 1;
	      k = j;
	    }
	    else if (downadj->fvtemplate[1+i][1+(j-1+nfv)%nfv] == fvids[1]) {
	      fnd = 1;
	      dir = 0;
	      k = j;
	    }
	  }

	/* Go to next face if even these 2 vertices are not in the face */
	if (!fnd) 
	  continue;

	/* See if the rest of the vertices of given face are there in 
	   the face of the region */

	for (j = 0, same = 1; same && j < nfv; j++) {
	  if (dir) {
	    if (downadj->fvtemplate[1+i][1+(k+j)%nfv] != fvids[j])
	      same = 0;
	  }
	  else {
	    if (downadj->fvtemplate[1+i][1+(k-j+nfv)%nfv] != fvids[j])
	      same = 0;
	  }
	}
	if (!same)
	  continue;

	/* The faces are the same. return the direction */
	return dir;
      }      
    }
    else {
      int nrv = List_Num_Entries(downadj->rvertices);
      MRType rgtyp = MSTK_nv2rtype[nrv];

      nf = MSTK_nrf_template[rgtyp];

      for (i = 0; i < nf; i++) {
	nfv2 = MSTK_rfv_template[rgtyp][i][0];
	if (nfv != nfv2)
	  continue;

	/* Are the first pair of vertices of the given face in this
	   face of the region and are they in the same or reverse
	   order */

	for (j = 0, fnd = 0, k = -1; !fnd && (j < nfv); j++)
	  if (MSTK_rfv_template[rgtyp][i][1+j] == fvids[0]) {
	    if (MSTK_rfv_template[rgtyp][i][1+(j+1)%nfv] == fvids[1]) {
	      fnd = 1;
	      dir = 1;
	      k = j;
	    }
	    else if (MSTK_rfv_template[rgtyp][i][1+(j-1+nfv)%nfv] == fvids[1]){
	      fnd = 1;
	      dir = 0;
	      k = j;
	    }
	  }

	/* Go to next face if even these 2 vertices are not in the face */
	if (!fnd) 
	  continue;

	/* See if the rest of the vertices of given face are there in 
	   the face of the region */

	for (j = 0, same = 1; same && j < nfv; j++) {
	  if (dir) {
	    if (MSTK_rfv_template[rgtyp][i][1+(k+j)%nfv] != fvids[j])
	      same = 0;
	  }
	  else {
	    if (MSTK_rfv_template[rgtyp][i][1+(k-j+nfv)%nfv] != fvids[j])
	      same = 0;
	  }
	}
	if (!same)
	  continue;

	/* The faces are the same. Return the direction */
	return dir;
      }      

    }
    
    return -1; 
 }

  int MR_FaceDir_i_R1R2(MRegion_ptr r, int i) {
    return 1; /* Is this the right thing to do? Yes! */
  }


  int MR_UsesFace_R1R2(MRegion_ptr r, MFace_ptr f) {
    if (MR_FaceDir_R1R2(r,f) == -1)
      return 0; /* Face is not in region */
    else
      return 1;
  }

#ifdef __cplusplus
}
#endif
