#ifndef _H_MSTK_EXTERN
#define _H_MSTK_EXTERN

/* If the compiler complains, change first line to 
   #ifndef _H_MSTK_EXTERN && _H_MSTK_GLOBALVARS
*/

#include "MSTK_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Templates for TETS, PYRAMIDS, PRISMS and HEXS (move to global file?) */

/* Number of faces for standard regions */
extern int const MSTK_nrf_template[5];

/* Number of edges for standard regions */
extern int const MSTK_nre_template[5];

/* Translate number of vertices to standard region types */
extern const MRType MSTK_nv2rtype[9];

/* Face-vertex templates for standard regions */
extern const int MSTK_rfv_template[5][6][5];

/* Face direction templates for standard regions */
extern const int MSTK_rfdir_template[5][6];

/* Face-edge templates for standard regions */
extern const int MSTK_rfe_template[5][6][5]; 

/* Face-edge templates for standard regions */
extern const int MSTK_rfedir_template[5][6][5]; 

/* Edge-vertex templates for standard regions */
extern const int MSTK_rev_template[5][12][2];

extern MAttrib_ptr plnumatt;

#ifdef __cplusplus
}
#endif

#endif
