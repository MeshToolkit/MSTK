#ifndef _H_MSTK_DEFINES
#define _H_MSTK_DEFINES


#define MSTK_UNKNOWN -1

#define MSTK_VER    1.0

#define MSTK_MAXREP 10

/* In reality we can have as many faces and vertices as we want but these
   constants can be used to size some arrays by application programs ?? */

#define MAXPV2 30
#define MAXPV3 50
#define MAXPF3 4*(8*sizeof(unsigned int))  /* Expect this to be 128 */

#define NMAXATT 10




#endif
