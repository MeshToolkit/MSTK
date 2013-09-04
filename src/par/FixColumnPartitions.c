#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "MSTK.h"

/* Sometimes partitioners will partition vertical (Z) columns of
   elements into two even though we requested only a horizontal
   partition - Tweak the partition so that columns of elements are on
   the same partition 

   If this procedure is applied to a non-columnar mesh, the results
   are unpredictable
*/


#ifdef __cplusplus
extern "C" {
#endif

int FixColumnPartitions(Mesh_ptr mesh, int *part, MSTK_Comm comm) { 
  int i, j, k, idx, nrf, nfr, nfv, nrv, curid, nxtid, homepid, nxtpid;
  int interior, done, nmove=0, found;
  MRegion_ptr mr, curreg, nxtreg;
  MFace_ptr rf, topf, botf;
  MVertex_ptr fv;
  List_ptr rfaces, fregs, fverts;
  double fxyz[MAXPV2][3], rxyz[MAXPV3][3], rfvec[3], rcen[3], fcen[3], maxdp, norm;
  
  if (MESH_Num_Regions(mesh) == 0) return 1; /* Not implemented for 2D (perhaps not needed either) */

  idx = 0;
  while ((mr = MESH_Next_Region(mesh,&idx))) {
    rfaces = MR_Faces(mr);
    nrf = List_Num_Entries(rfaces);
    
    interior = 1;
    for (i = 0; i < nrf; i++) {
      rf = List_Entry(rfaces,i);      
      fregs = MF_Regions(rf);
      if (List_Num_Entries(fregs) == 1) 
        interior = 0;
      List_Delete(fregs);
      if (!interior) break;
    }

    if (interior) continue;

    /* Found a boundary region - Find the face whose centroid is
       most directly above the region centroid and check if this
       is a boundary face */

    rcen[0] = rcen[1] = rcen[2] = 0.0;
    MR_Coords(mr,&nrv,rxyz);
    for (i = 0; i < nrv; i++) {
      for (k = 0; k < 3; k++)
        rcen[k] += rxyz[i][k];
    }
    for (k = 0; k < 3; k++)
      rcen[k] /= nrv;        

    topf = NULL;
    maxdp = -1.0e+10;
    for (i = 0; i < nrf; i++) {
      rf = List_Entry(rfaces,i);
      MF_Coords(rf,&nfv,fxyz);

      fcen[0] = fcen[1] = fcen[2] = 0.0;
      for (j = 0; j < nfv; j++) {
        for (k = 0; k < 3; k++)
          fcen[k] += fxyz[j][k];
      }
      for (k = 0; k < 3; k++)
        fcen[k] /= nfv;

      for (k = 0; k < 3; k++)
        rfvec[k] = fcen[k]-rcen[k];
      norm = sqrt(rfvec[0]*rfvec[0]+rfvec[1]*rfvec[1]+rfvec[2]*rfvec[2]);
      for (k = 0; k < 3; k++)
        rfvec[k] /= norm;

      /* rfvec[2] is same as rfvec dotted with [0,0,1] */
      if (rfvec[2] > maxdp) {
        maxdp = rfvec[2];
        topf = rf;
      }
    }
    List_Delete(rfaces);

    if (!topf)
      MSTK_Report("FixColumnPartitions","Could not find top face?",MSTK_FATAL);

    /* Check if the top face is also a top boundary face */

    fregs = MF_Regions(topf);
    nfr = List_Num_Entries(fregs);

    if (nfr != 1) {
      List_Delete(fregs);
      continue; /* Not a mesh region at the top of the mesh */
    }
    
    /* Hopefully, found the top face of a column in the mesh */

    curreg = List_Entry(fregs,0);
    curid = MR_ID(curreg);
    homepid = part[curid-1];

    List_Delete(fregs);


    done = 0;
    while (!done) {

      rfaces = MR_Faces(curreg);
      nrf = List_Num_Entries(rfaces);
      botf = 0;
      for (i = 0; i < nrf; i++) {
        rf = List_Entry(rfaces,i);

        fverts = MF_Vertices(rf,1,0);
        nfv = List_Num_Entries(fverts);
        found = 0;
        for (j = 0; j < nfv; j++) {
          fv = List_Entry(fverts,j);
          if (MF_UsesEntity(topf,fv,MVERTEX)) {
            found = 1;  /* found common vertex betweeen two faces */
            break;
          }
        }          
        List_Delete(fverts);

        if (!found) {
          botf = rf;
          break;
        }
      }
      List_Delete(rfaces);

      /* Found bottom face of current mesh region - get
       the next lower mesh region and check if it is on
      the same partition as the top region */

      fregs = MF_Regions(botf);

      if (List_Num_Entries(fregs) == 2) {
        nxtreg = List_Entry(fregs,0);
        if (nxtreg == curreg)
          nxtreg = List_Entry(fregs,1);
        
        nxtid = MR_ID(nxtreg);
        nxtpid = part[nxtid-1];
        if (nxtpid != homepid) {
          part[nxtid-1] = homepid;
          nmove++;
        }

        curreg = nxtreg;
        curid = nxtid;
        topf = botf;
      }
      else 
        done = 1;
      List_Delete(fregs);
    }
    
  }

  if (nmove) {
    char msg[256];
    sprintf(msg,"Redistributed %-d elements to maintain column partitioning",nmove);
    MSTK_Report("FixColumnPartitions",msg,MSTK_MESG);
  }
    
  return 1;
}

#ifdef __cplusplus
  }
#endif

