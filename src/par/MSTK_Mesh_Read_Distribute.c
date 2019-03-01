/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  

  /* Read a mesh in, partition it and distribute it to 'num' processors 

     Authors: Duo Wang
              Rao Garimella
  */

  int MSTK_Mesh_Read_Distribute(Mesh_ptr *recv_mesh, 
				const char* global_mesh_name, 
				int *dim, int ring, int with_attr,
                                int method, MSTK_Comm comm) {

    int i, ok;

    int rank;
    MPI_Comm_rank(comm,&rank);

    Mesh_ptr mesh=NULL;
    if(rank == 0) {
      mesh = MESH_New(UNKNOWN_REP);
      ok = MESH_InitFromFile(mesh,global_mesh_name,comm);
      if (!ok) {
	fprintf(stderr,"Cannot open input file %s\n\n\n",global_mesh_name);
	exit(-1);
      }
      fprintf(stdout,"mstk mesh %s read in successfully\n",global_mesh_name);

      *dim = MESH_Num_Regions(mesh) ? 3 : 2;
    }

    int del_inmesh = 1;
    MSTK_Mesh_Distribute(mesh, recv_mesh, dim, ring, with_attr, method, 
			 del_inmesh, comm);

    return 1;
  }


#ifdef __cplusplus
}
#endif

