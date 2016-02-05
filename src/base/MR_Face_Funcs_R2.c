#define _H_MRegion_Private

#include <stdlib.h>
#include "MRegion.h"
#include "MRegion_jmp.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  int MR_Num_Faces_R2(MRegion_ptr r) {
    MRegion_Adj_R2 *adj;
    adj = (MRegion_Adj_R2 *) r->adj;
    
    if (adj->fvtemplate) {
      return adj->fvtemplate[0][0];
    }
    else {
      int nv = List_Num_Entries(adj->rvertices);
      MRType regtype = MSTK_nv2rtype[nv];
      return MSTK_nrf_template[regtype];
    }
  }

  List_ptr MR_Faces_R2(MRegion_ptr r) {
    MRegion_Adj_R2 *adj;
    adj = (MRegion_Adj_R2 *) r->adj;
    int i, j, k, nf, nfv;
    List_ptr rfaces=NULL;
    MFace_ptr face;
    Mesh_ptr mesh = MEnt_Mesh((MEntity_ptr) r);
    MVertex_ptr fvertices[MAXPV2], fverticestmp[MAXPV2];
    int tmin, ti, tj, tdir;
    
    if (adj->fvtemplate) {
      MVertex_ptr fvertices[MAXPV2];

      nf = adj->fvtemplate[0][0];
      rfaces = List_New(nf);

      for (i = 0; i < nf; i++) {
	nfv = adj->fvtemplate[1+i][0];
	
	for (j = 0; j < nfv; j++) {
	  k = adj->fvtemplate[1+i][1+j];
	  fvertices[j] = List_Entry(adj->rvertices,k);
	}

#ifdef HASHTABLE
	tmin = 0;
	fverticestmp[0] = fvertices[0];
	for (ti = 1; ti < nfv; ti++) {
	  if (fverticestmp[0] > fvertices[ti]) {
	    tmin = ti;
	    fverticestmp[0] = fvertices[ti];
	  }
	}
	tdir = (fvertices[(tmin+1) % nfv] < fvertices[(tmin+nfv-1) % nfv]) ? 1 : nfv-1;

	for (tj = 1, ti=(tmin+tdir) % nfv; tj < nfv; tj++, ti=(ti+tdir) % nfv) {
		fverticestmp[tj] = fvertices[ti];
	}
	face = Hash_Entry(MESH_Hash_Faces(mesh), nfv, fvertices);

	if (face == NULL) {
	  face = MF_New(mesh);
	  MEnt_Set_Volatile(face);
	  MF_Set_Vertices(face,nfv,fvertices);
	  MF_Set_GInfo_Auto(face);

	  Hash_Add(MESH_Hash_Faces(mesh), face, nfv, fvertices);
	}
#else
	face = MF_New(mesh);
	MEnt_Set_Volatile(face);
	MF_Set_Vertices(face,nfv,fvertices);
	MF_Set_GInfo_Auto(face);
#endif
	List_Add(rfaces,face);	
	MF_Lock(face);
      }
    }
    else {
      int nv = List_Num_Entries(adj->rvertices);

      MRType regtype = MSTK_nv2rtype[nv];

      nf = MSTK_nrf_template[regtype];
      rfaces = List_New(nf);

      for (i = 0; i < nf; i++) {
	nfv = MSTK_rfv_template[regtype][i][0];

	for (j = 0; j < nfv; j++) {
	  k = MSTK_rfv_template[regtype][i][1+j];
	  fvertices[j] = List_Entry(adj->rvertices,k);
	}

#ifdef HASHTABLE
	tmin = 0;
	fverticestmp[0] = fvertices[0];
	for (ti = 1; ti < nfv; ti++) {
	  if (fverticestmp[0] > fvertices[ti]) {
	    tmin = ti;
	    fverticestmp[0] = fvertices[ti];
	  }
	}
	tdir = (fvertices[(tmin+1) % nfv] < fvertices[(tmin+nfv-1) % nfv]) ? 1 : nfv-1;

	for (tj = 1, ti=(tmin+tdir) % nfv; tj < nfv; tj++, ti=(ti+tdir) % nfv) {
		fverticestmp[tj] = fvertices[ti];
	}

	face = Hash_Entry(MESH_Hash_Faces(mesh), nfv, fverticestmp);

	if (face == NULL) {
	  face = MF_New(mesh);
	  MEnt_Set_Volatile(face);
	  MF_Set_Vertices(face,nfv,fverticestmp);
	  MF_Set_GInfo_Auto(face);

	  Hash_Add(MESH_Hash_Faces(mesh), face, nfv, fverticestmp);
	}
#else
	face = MF_New(mesh);
	MEnt_Set_Volatile(face);
	MF_Set_Vertices(face,nfv,fvertices);
	MF_Set_GInfo_Auto(face);
#endif
	List_Add(rfaces,face);
	MF_Lock(face);
      }
    }

    if (!MESH_AutoLock(mesh)) {
       i = 0;
       while ((face = List_Next_Entry(rfaces, &i))) {
	 MF_UnLock(face);
       }
    }
    return rfaces;
  }

  void MR_FaceIDs_R2(MRegion_ptr r, int *nrf, int *faceids) {
    MSTK_Report("MR_FaceIDs_R2","Not implemented",MSTK_FATAL);
  }


  int MR_Rev_FaceDir_R2(MRegion_ptr r, MFace_ptr f) {
    return 0;
  }

  int MR_Rev_FaceDir_i_R2(MRegion_ptr r, int i) {
    return 0;
  }

  int MR_Set_FaceDir_R2(MRegion_ptr r, MFace_ptr f, int dir) {
    return 0;
  }

  int MR_Set_FaceDir_i_R2(MRegion_ptr r, int i, int dir) {
    return 0;
  }


  int MR_FaceDir_R2(MRegion_ptr r, MFace_ptr f) {
    MRegion_Adj_R2 *adj;
    adj = (MRegion_Adj_R2 *) r->adj;
    int *rvids, nrv, fvids[MAXPV2], nfv, nfv2, idx, i, j, k, dir=-1, fnd=0, same=0, nf;
    List_ptr fverts;
    MVertex_ptr v;

    fverts = MF_Vertices(f,1,0);
    nfv = List_Num_Entries(fverts);

    idx = 0; i = 0;
    while ((v = List_Next_Entry(fverts,&idx)))
      fvids[i++] = MV_ID(v);

    List_Delete(fverts);

    nrv = List_Num_Entries(adj->rvertices);
    rvids = (int *) malloc(nrv*sizeof(int));
    idx = 0; i = 0;
    while ((v = List_Next_Entry(adj->rvertices, &idx)))
      rvids[i++] = MV_ID(v);

    /* Check if a face in the region has the same vertices as the
       given face. Then, assuming that the region-face-vertex
       templates contain face description such that the normal points
       out of the region, we have to check if the vertices of the
       given face are in the same or opposite direction to the
       matching face in region */

    if (adj->fvtemplate) {

      nf = adj->fvtemplate[0][0];

      for (i = 0; i < nf; i++) {
	nfv2 = adj->fvtemplate[1+i][0];
	if (nfv != nfv2)
	  continue;

	/* Are the first pair of vertices of the given face in this
	   face of the region and are they in the same or reverse
	   order */

	for (j = 0, fnd = 0, k = -1; !fnd && (j < nfv); j++)
	  if (rvids[adj->fvtemplate[1+i][1+j]] == fvids[0]) {
	    if (rvids[adj->fvtemplate[1+i][1+(j+1)%nfv]] == fvids[1]) {
	      fnd = 1;
	      dir = 1;
	      k = j;
	    }
	    else if (rvids[adj->fvtemplate[1+i][1+(j-1+nfv)%nfv]] == fvids[1]) {
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
	    if (rvids[adj->fvtemplate[1+i][1+(k+j)%nfv]] != fvids[j])
	      same = 0;
	  }
	  else {
	    if (rvids[adj->fvtemplate[1+i][1+(k-j+nfv)%nfv]] != fvids[j])
	      same = 0;
	  }
	}
	if (!same)
	  continue;

	/* The faces are the same. return the direction */
	break;
      }      
    }
    else {
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
	  if (rvids[MSTK_rfv_template[rgtyp][i][1+j]] == fvids[0]) {
	    if (rvids[MSTK_rfv_template[rgtyp][i][1+(j+1)%nfv]] == fvids[1]) {
	      fnd = 1;
	      dir = 1;
	      k = j;
	    }
	    else if (rvids[MSTK_rfv_template[rgtyp][i][1+(j-1+nfv)%nfv]] == fvids[1]){
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
	    if (rvids[MSTK_rfv_template[rgtyp][i][1+(k+j)%nfv]] != fvids[j])
	      same = 0;
	  }
	  else {
	    if (rvids[MSTK_rfv_template[rgtyp][i][1+(k-j+nfv)%nfv]] != fvids[j])
	      same = 0;
	  }
	}
	if (!same)
	  continue;

	/* Check if the face is in the correct direction in template */
	if (MSTK_rfdir_template[rgtyp][i]==0)
	  dir = 1-dir;

	/* The faces are the same. Return the direction */
	break;
      }      

    }
    free(rvids);
    if (!same) dir = -1;
    if (!fnd) dir = -1;
    return dir; 
  }

  int MR_FaceDir_i_R2(MRegion_ptr r, int i) {
#ifdef HASHTABLE
    MRegion_Adj_R2 *adj;
    MRType regtype;
    MVertex_ptr fvertices[MAXPV2], vtmp;
    int j, k, nfv, tdir, nv;

    adj = (MRegion_Adj_R2 *) r->adj;
    if (adj->fvtemplate) {

      nfv = adj->fvtemplate[1+i][0];

      for (j = 0; j < nfv; j++) {
	k = adj->fvtemplate[1+i][1+j];
	fvertices[j] = List_Entry(adj->rvertices,k);
      }

      k = 0;
      vtmp = fvertices[0];
      for (j = 1; j < nfv; j++) {
	if (vtmp > fvertices[j]) {
	  k = j;
	  vtmp = fvertices[j];
	}
      }
      tdir = (fvertices[(k+1) % nfv] < fvertices[(k+nfv-1) % nfv]) ? 1 : 0;

    }
    else {
      nv = List_Num_Entries(adj->rvertices);

      regtype = MSTK_nv2rtype[nv];

      nfv = MSTK_rfv_template[regtype][i][0];

      for (j = 0; j < nfv; j++) {
	k = MSTK_rfv_template[regtype][i][1+j];
	fvertices[j] = List_Entry(adj->rvertices,k);
      }

      k = 0;
      vtmp = fvertices[0];
      for (j = 1; j < nfv; j++) {
	if (vtmp > fvertices[j]) {
	  k = j;
	  vtmp = fvertices[j];
	}
      }
      tdir = (fvertices[(k+1) % nfv] < fvertices[(k+nfv-1) % nfv]) ? 1 : 0;

      /* Check if the face is in the correct direction in template */
      if (MSTK_rfdir_template[regtype][i]==0)
	tdir = 1-tdir;

    }
    return tdir;
#else
    MRegion_Adj_R2 *adj;
    MRType rgtyp;
    int nrv;
    adj = (MRegion_Adj_R2 *) r->adj;
    if (adj->fvtemplate) {
      return 1; /* Is this the right thing to do? Yes! */ /* Hmmm, not sure... */
    }
    else {
      nrv = List_Num_Entries(adj->rvertices);
      rgtyp = MSTK_nv2rtype[nrv];
      return MSTK_rfdir_template[rgtyp][i];
    }
#endif
  }

  void MR_Update_ElementType_R2(MRegion_ptr r) {
    MSTK_Report("MR_Update_ElementType_R2","Not yet implemented",MSTK_ERROR);
  }


  int MR_UsesFace_R2(MRegion_ptr r, MFace_ptr f) {
    if (MR_FaceDir_R2(r,f) == -1)
      return 0; /* Face is not in region */
    else
      return 1;
  }

#ifdef __cplusplus
}
#endif
