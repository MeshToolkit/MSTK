#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MSTK_private.h"
#include "MSTK_globals.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif




  void MSTK_Init(MPI_Comm comm) {

    int initstatus;
    MPI_Initialized(&initstatus);
    if (initstatus == 0) {
      int loc_argc=1;
      char **loc_argv=NULL;
      
      MPI_Init(&loc_argc,&loc_argv);
    }

    MSTK_communicator = comm;
    
    /* Make some dummy calls so that the symbols get defined */
    
    MV_Print(0,0);
    ME_Print(0,0);
    MF_Print(0,0);
    MR_Print(0,0);
    {
      int n = MSTK_rev_template[0][0][0];
    }
    
  }


  MPI_Comm MSTK_Comm() {
    return MSTK_communicator;
  }


  /* MARKERS */

  unsigned int MSTK_marker = 0;
  int MSTK_MAXBITS = 8*sizeof(unsigned int);
  int MSTK_lastbit = -1;
  
  int MSTK_GetMarker() {
    int i;

    if (MSTK_lastbit < MSTK_MAXBITS-1) {
      MSTK_lastbit++;
      MSTK_marker = MSTK_marker | 1<<MSTK_lastbit;
      return (MSTK_lastbit+1); /* lastbit goes from 0 to N-1, markers from 1 to N */
    }
    else {
      /* Search if any other bit is free */
      i = 0;
      while (i <= MSTK_MAXBITS-1) {
	if (MSTK_marker & 1<<i)
	  i++;
	else {
	  MSTK_marker = MSTK_marker | 1<<i;
	  return i+1;
	}
      }

#ifdef DEBUG
      MSTK_Report("MSTK_GetMarker","Out of Markers",MSTK_ERROR);
#endif
      return -1;
    }

#ifdef DEBUG
    MSTK_Report("MSTK_GetMarker","Cannot allocate Marker",MSTK_ERROR);
#endif
    return -1;
  }

  void MSTK_FreeMarker(int mkr) {

#ifdef DEBUG
    if (mkr < 1 || mkr > MSTK_MAXBITS)
      MSTK_Report("MSTK_FreeMarker","Invalid Marker",MSTK_ERROR);
#endif

    MSTK_marker = MSTK_marker & ~(1<<(mkr-1));
    if (mkr-1 == MSTK_lastbit)
      MSTK_lastbit--;
  }


#ifdef __cplusplus
}
#endif
