#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef MSTK_USE_MARKERS
#include <pthread.h>
#endif

#include "MSTK_private.h"
#include "MSTK_globals.h"

#ifdef MSTK_HAVE_MPI
#include "mpi.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif




  void MSTK_Init() {

#ifdef MSTK_HAVE_MPI

    int initstatus;
    MPI_Initialized(&initstatus);
    if (initstatus == 0) {
      int loc_argc=0;
      char **loc_argv=NULL;
      
      MPI_Init(&loc_argc,&loc_argv);
    }

#endif /* MSTK_HAVE_MPI */

    /* Make some dummy calls so that the symbols get defined */
    
    MV_Print(0,0);
    ME_Print(0,0);
    MF_Print(0,0);
    MR_Print(0,0);
    {
      int n = MSTK_rev_template[0][0][0];
    }
    
  }

#ifdef MSTK_HAVE_MPI
  int MSTK_Comm_size(MSTK_Comm comm) {
    int size;
    MPI_Comm_size(comm,&size);
    return size;
  }

  int MSTK_Comm_rank(MSTK_Comm comm) {
    int rank;
    MPI_Comm_rank(comm,&rank);
    return rank;
  }
#endif /* MSTK_HAVE_MPI */


  /* MARKERS */

#ifdef MSTK_USE_MARKERS
  unsigned int MSTK_marker = 0;
  int MSTK_MAXBITS = 8*sizeof(unsigned int);
  int MSTK_lastbit = -1;
  pthread_mutex_t marker_lock;  /* declared extern in MSTK_private.h */
  
  int MSTK_GetMarker() {
    static int first = 1;
    if (first) {
      first = 0;
      fprintf(stderr,"MSTK using markers\n");
    }
    int i, marker = -1;
    pthread_mutex_lock(&marker_lock);
    
    if (MSTK_lastbit < MSTK_MAXBITS-1) {      
      MSTK_lastbit++;
      MSTK_marker = MSTK_marker | 1<<MSTK_lastbit;
      
      /* lastbit goes from 0 to N-1, markers from 1 to N */
      marker = MSTK_lastbit+1;

    } else {
      /* Search if any other bit is free */
      i = 0;
      while (i <= MSTK_MAXBITS-1) {
	if (MSTK_marker & 1<<i)
	  i++;
	else {
	  MSTK_marker = MSTK_marker | 1<<i;
	  marker = i+1;
          break;
	}
      }

      if (marker == -1)
        MSTK_Report("MSTK_GetMarker","Out of Markers",MSTK_FATAL);
    }

    pthread_mutex_unlock(&marker_lock);
    return marker;
  }
  
  void MSTK_FreeMarker(int mkr) {
    pthread_mutex_lock(&marker_lock);

#ifdef DEBUG
    if (mkr < 1 || mkr > MSTK_MAXBITS)
      MSTK_Report("MSTK_FreeMarker","Invalid Marker",MSTK_ERROR);
#endif

    MSTK_marker = MSTK_marker & ~(1<<(mkr-1));
    if (mkr-1 == MSTK_lastbit)
      MSTK_lastbit--;

    if (MSTK_marker != (int) pow(2,MSTK_lastbit+1)-1) {

      /* Check what the last free bit is - may not be the one that was
         just freed if multiple markers were allocated and freed in
         random order */

      MSTK_lastbit = -1;
      int i = 0;       
      while (i <= MSTK_MAXBITS-1) {
        if (MSTK_marker & 1<<i)
          MSTK_lastbit = i;
        i++;
      }

    }

    pthread_mutex_unlock(&marker_lock);
  }
#endif  /* ifdef MSTK_USE_MARKERS */

#ifdef __cplusplus
}
#endif
