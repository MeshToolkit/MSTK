#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "MSTK.h"


#ifdef __cplusplus
extern "C" {
#endif

  /* 
     Add Global IDs for edges (and faces in 3D meshes)
  */

int MESH_Add_EdgeGIDs(Mesh_ptr mesh, int rank, MPI_Comm comm);
int MESH_Add_FaceGIDs(Mesh_ptr mesh, int rank, MPI_Comm comm);


int MESH_Add_EdgeFaceGIDs(Mesh_ptr mesh, int dim, int rank, MPI_Comm comm) {

  if (MESH_RepType(mesh) <= F4)
    MESH_Add_EdgeGIDs(mesh, rank, comm);
  
  if (dim == 3 && (MESH_RepType(mesh) <= F4 || MESH_RepType(mesh) == R4))
    MESH_Add_FaceGIDs(mesh, rank, comm);

  return 1;
}
  

int MESH_Add_EdgeGIDs(Mesh_ptr mesh, int rank, MPI_Comm comm) {
  int i, ne, idx, id, numprocs, myprocid;
  MEdge_ptr me;
  int *allne;
 

  MPI_Comm_size(comm,&numprocs);
  MPI_Comm_rank(comm,&myprocid);

  /* Count the number of edges on this processor */

  idx = 0; ne = 0;
  while ((me = MESH_Next_Edge(mesh,&idx)))
    if (ME_PType(me) != PGHOST) ne++;


  /* Gather the number of edges on each processor */

  allne = (int *) MSTK_malloc(numprocs*sizeof(int));

  MPI_AllGather(&ne,1,MPI_INT,allne,1,MPI_INT,comm);

  /* Now we can label the edges of this processor */
  
  id = 0;
  for (i = 0; i < myprocid; i++)
    id += allne[i];

  idx = 0; 
  while ((me = MESH_Next_Edge(mesh,&idx)))
    if (ME_PType(me) != PGHOST)
      ME_Set_GlobalID(me,id++);

  MSTK_free(allne);

  /* How do we transmit the IDs of these edges to the ghost copies on
     other processors? */

  return 1;
}

int MESH_Add_FaceGIDs(Mesh_ptr mesh, int rank, MPI_Comm comm) {
  int i, nf, idx, id, numprocs, myprocid;
  MFace_ptr mf;
  int *allnf;

  MPI_Comm_size(comm,&numprocs);
  MPI_Comm_rank(comm,&myprocid);

  /* Count the number of faces on this processor */

  idx = 0; nf = 0;
  while ((mf = MESH_Next_Face(mesh,&idx)))
    if (MF_PType(mf) != PGHOST) nf++;


  /* Gather the number of faces on each processor */

  allnf = (int *) MSTK_malloc(numprocs*sizeof(int));

  MPI_AllGather(&nf,1,MPI_INT,allnf,1,MPI_INT,comm);

  /* Now we can label the faces of this processor */
  
  id = 0;
  for (i = 0; i < myprocid; i++)
    id += allnf[i];

  idx = 0; 
  while ((mf = MESH_Next_Face(mesh,&idx)))
    if (MF_PType(mf) != PGHOST)
      MF_Set_GlobalID(mf,id++);

  MSTK_free(allnf);

  /* How do we transmit the IDs of these faces to the ghost copies on
     other processors? */

  return 1;
}

#ifdef __cplusplus
}
#endif

