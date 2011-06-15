#ifndef _H_MSTK_DEFINES
#define _H_MSTK_DEFINES

#ifdef __cplusplus
extern "C" {
#endif

#define MSTK_UNKNOWN -1

#define MSTK_VER    1.0

#define MSTK_MAXREP 10

/* Number of edges (or vertices) in a polygon and number of faces in a
   polyhedron is determined by how much space we have to store their
   directions */

#define MAXPV2 30
#define MAXPV3 200
#define MAXPF3 4*(8*sizeof(unsigned int))  /* Expect this to be 128 */

#define NMAXATT 10

#ifdef __cplusplus
}
#endif

#endif
