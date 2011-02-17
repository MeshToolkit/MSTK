#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  int compareGlobalID(MEntity_ptr *a, MEntity_ptr *b) {
    if ( MEnt_GlobalID(*a) > MEnt_GlobalID(*b) ) 
      return 1;
    else if ( MEnt_GlobalID(*a) < MEnt_GlobalID(*b) ) 
      return -1;
    else
      return 0;
  }

#ifdef __cplusplus
}
#endif
