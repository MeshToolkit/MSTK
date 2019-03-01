/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#ifndef _H_MSTK_DEFINES
#define _H_MSTK_DEFINES

#ifdef __cplusplus
extern "C" {
#endif


/* The following lines (MSTK_*_VERSION and if present, #define
 * MSTK_HAVE_MPI) are filled in by CMake during installation */

#define MSTK_VERSION_MAJOR @MSTK_VERSION_MAJOR@
#define MSTK_VERSION_MINOR @MSTK_VERSION_MINOR@
#define MSTK_VERSION_PATCH @MSTK_VERSION_PATCH@

#cmakedefine MSTK_HAVE_MPI
 


#define MSTK_FILE_VER 1.0

#define MSTK_UNKNOWN -1

#define MSTK_MAXREP 10

/* Number of edges (or vertices) in a polygon and number of faces in a
   polyhedron is determined by how much space we have to store their
   directions */

#define MAXPV2 30
#define MAXPV3 200
#define MAXPF3 16*(8*sizeof(unsigned int))  /* Actually it is unlimited */

#define NMAXATT 10

#ifdef __cplusplus
}
#endif

#endif
