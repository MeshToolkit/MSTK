#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MSTK_private.h"
#include "Set.h"
#include "MSTK_malloc.h"

#ifdef __cplusplus
extern "C" {
#endif

  unsigned int MSTK_marker = 0;
  int MSTK_MAXBITS = 8*sizeof(unsigned int);
  int MSTK_lastbit = -1;

  void MSTK_Init() {

    MV_Print(0,0);
    ME_Print(0,0);
    MF_Print(0,0);
    MR_Print(0,0);
  }

  int MSTK_GetMarker() {
    int i;

    /*** CHECK CAREFULLY ***/

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
      MSTK_Report("MSTK_GetMarker","Out of Markers",ERROR);
#endif
      return -1;
    }

#ifdef DEBUG
    MSTK_Report("MSTK_GetMarker","Cannot allocate Marker",ERROR);
#endif
    return -1;
  }

  void MSTK_FreeMarker(int mkr) {

#ifdef DEBUG
    if (mkr < 1 || mkr > MSTK_MAXBITS)
      MSTK_Report("MSTK_FreeMarker","Invalid Marker",ERROR);
#endif

    MSTK_marker = MSTK_marker & ~(1<<(mkr-1));
    if (mkr-1 == MSTK_lastbit)
      MSTK_lastbit--;
  }


#ifdef __cplusplus
}
#endif
