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

int FixColumnPartitions_IsSideFace(Mesh_ptr mesh, MRegion_ptr mr, MFace_ptr rf) {
  double fxyz[MAXPV2][3];
  double znorm;
  int j, jm, jp, nfv;
  
  MF_Coords(rf,&nfv,fxyz);

  znorm = 0.;
  for (j = 0; j < nfv; j++) {
    /* determine the average z-normal of the face */
    if (j == 0) jm = nfv-1;
    else jm = j-1;

    if (j == nfv-1) jp = 0;
    else jp = j+1;

    znorm += (fxyz[j][1]-fxyz[jm][1]) * (fxyz[jp][0] - fxyz[j][0])
        - (fxyz[jp][1]-fxyz[j][1]) * (fxyz[j][0] - fxyz[jm][0]);
  }
  znorm /= nfv;
  if (fabs(znorm) < 1.e-12) 
    /* z-normal is zero, face is a column side by definition */
    return 1;
  else
    return 0;
}

/* Searches a region mr to find the "upward" and "downward" oriented faces.

   Returns 1 if the relationship is proper, or 0 if the region has
   zero volume and is therefore uncertain.
*/
int FixColumnPartitions_UpDownFaces(Mesh_ptr mesh, MRegion_ptr mr, MFace_ptr* up, MFace_ptr* dn) {
  int found;
  List_ptr rfaces, rfaces2;
  int nrf, i, j, k, nfv, ret;
  MFace_ptr rf; double rcen[3], upcen[3], dncen[3];
  double fxyz[MAXPV2][3];
  MRegion_ptr mr_it;

  *up=NULL; *dn=NULL;
  
  rfaces = MR_Faces(mr);
  nrf = List_Num_Entries(rfaces);

  /* find two faces that are not sides */
  found = 0; i = 0;
  while (found < 2) {
    rf = List_Entry(rfaces,i);
    if (!FixColumnPartitions_IsSideFace(mesh, mr, rf)) {
      if (found < 1) *up = rf;
      else *dn = rf;
      found++;
    }
    i++;
    if (i >= nrf && found < 2) 
      MSTK_Report("FixColumnPartitions","Mesh is not columnar, can't find both up and down faces.",MSTK_FATAL);
  }

  /* figure which is up and which is down */
  MF_Coords(*up,&nfv,fxyz);
  upcen[0] = upcen[1] = upcen[2] = 0.;
  for (j = 0; j < nfv; j++) {
    for (k = 0; k < 3; k++)
      upcen[k] += fxyz[j][k];
  }

  MF_Coords(*dn,&nfv,fxyz);
  dncen[0] = dncen[1] = dncen[2] = 0.;
  for (j = 0; j < nfv; j++) {
    for (k = 0; k < 3; k++)
      dncen[k] += fxyz[j][k];
  }

  /* ASSERT x,y are the same */
  if (fabs(upcen[0] - dncen[0]) > 1.e-12 ||
      fabs(upcen[1] - dncen[1]) > 1.e-12) {
    MSTK_Report("FixColumnPartitions","Mesh is not quite columnar, up and down faces do not align.",MSTK_FATAL);
  }

  if (upcen[2] < dncen[2] - 1.e-12) {
    rf = *up;
    *up = *dn;
    *dn = rf;
    ret = 1;
  } else if (upcen[2] > dncen[2] + 1.e-12) {
    ret = 1;
  } else {
    /* find the first region upward that is proper */
    MSTK_Report("FixColumnPartitions","Partitioning with degenerate cells is not yet supported.",MSTK_FATAL);
    ret = 0;
  }

  List_Delete(rfaces);
  return ret;
}

  
int FixColumnPartitions(Mesh_ptr mesh, int *part, MSTK_Comm comm) { 
  int idx, curid, nxtid, homepid, nxtpid, num_cols, num_cells_in_col;
  int num_touches, num_not_touched;
  int interior, done, nmove=0;
  MRegion_ptr mr, curreg, nxtreg;
  MFace_ptr topf, botf, expected_topf;
  List_ptr fregs;
  int* touched;
  FILE* fid;
  
  if (MESH_Num_Regions(mesh) == 0) return 1; /* Not implemented for 2D (perhaps not needed either) */

  /* ensure every region is touched once */
  touched = calloc(MESH_Num_Regions(mesh),sizeof(int));
  for (idx=0; idx!=MESH_Num_Regions(mesh); ++idx) touched[idx] = 0;
#ifdef DEBUG
  printf("Read %-d regions", MESH_Num_Regions(mesh));
  /* loop over all regions, finding regions whose "up" face is a boundary */
  fid = fopen("col_counts.txt", "w");
#endif
  num_cols = 0;
  num_touches = 0;
  idx = 0;
  //TODO: FixColumnPartitions_UpDownFaces checks for side faces on every invocation, so
  //interior side faces are tested twice
  while ((mr = MESH_Next_Region(mesh,&idx))) {
    FixColumnPartitions_UpDownFaces(mesh, mr, &topf, &botf);
    
    interior = 1;
    fregs = MF_Regions(topf);
    if (List_Num_Entries(fregs) == 1) interior = 0;
    List_Delete(fregs);

    if (interior) continue;

    /* found a top face: iterate to the bottom */
    num_cols++;
    num_cells_in_col = 0;
    done = 0;
    curreg = mr;
    curid = MR_ID(curreg);
    homepid = part[curid-1];

    /* continue to loop through the column */
    expected_topf = topf;
    while (!done) {
      if (touched[curid-1]) 
        MSTK_Report("FixColumnPartitions","Mesh is not columnar, iteration process touches the same cell twice.",MSTK_FATAL);
      touched[curid-1] = 1;
      num_cells_in_col++;
      num_touches++;
      
      /* find the bottom face of that region */
      FixColumnPartitions_UpDownFaces(mesh, curreg, &topf, &botf);
      if (topf != expected_topf)
        MSTK_Report("FixColumnPartitions","Mesh is not columnar, up/down faces aren't consistent.",MSTK_FATAL);

      fregs = MF_Regions(botf);
      if (List_Num_Entries(fregs) == 2) {
        nxtreg = List_Entry(fregs,0);
        if (nxtreg == curreg) nxtreg = List_Entry(fregs,1);
        
        nxtid = MR_ID(nxtreg);
        nxtpid = part[nxtid-1];
        if (nxtpid != homepid) {
          part[nxtid-1] = homepid;
          nmove++;
        }

        curreg = nxtreg;
        curid = nxtid;
        expected_topf = botf;
      } else {
        done = 1;
      }
      List_Delete(fregs);
    }
#ifdef DEBUG
    fprintf(fid, "%-d\n", num_cells_in_col);
#endif
  }

  if (nmove) {
    char msg[256];
    sprintf(msg,"Redistributed %-d elements in %-d columns to maintain column partitioning",nmove, num_cols);
    MSTK_Report("FixColumnPartitions",msg,MSTK_MESG);
  }
#ifdef DEBUG
  fclose(fid);
#endif
  num_not_touched = 0;
  for (idx=0; idx!=MESH_Num_Regions(mesh); ++idx) {
    if (!touched[idx]) {
      num_not_touched++;
      char msg[256];
      sprintf(msg,"  missed %-d", idx);
      MSTK_Report("FixColumnPartitions",msg,MSTK_MESG);
    }
  }
  if (num_not_touched > 0) {
    char msg[256];
    sprintf(msg,"Mesh is not columnar, iteration process did not touch %-d cells, did perform %-d touches.", num_not_touched, num_touches);
    MSTK_Report("FixColumnPartitions",msg,MSTK_MESG);
  }
  free(touched);
  
  return 1;
}

#ifdef __cplusplus
  }
#endif

