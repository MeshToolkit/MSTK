#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "MSTK.h"
#include "MSTK_private.h"
#ifdef __cplusplus
extern "C" {
#endif



  /* 
     Reassign global IDs to mesh entities to make the contiguous on
     each partition and numbered according to some scheme. Assumes
     that the parallel mesh is fully setup with ghost connectivity
     etc.

     For now only method 0 (renumber sequentially) is supported

     Eventually, one could precompute some global IDs and just send it
     to this routine (NOT IMPLEMENTED). Cannot use mtype = MALLTYPE for
     preassigned GIDs

     Author(s): Rao Garimella
  */

  int MESH_Renumber_EntityGlobalIDs(Mesh_ptr mesh, MType mtype, int method,
                                    int *preassigned_gids, MSTK_Comm comm);
  
  int MESH_Renumber_GlobalIDs(Mesh_ptr mesh, MType mtype, int method,
                              int *preassigned_gids, MSTK_Comm comm) {
    if (preassigned_gids) {
      MSTK_Report("MESH_Renumber_GlobalIDs",
                  "Renumbering according to pre-assigned numbering not yet implemented",
                  MSTK_ERROR);
    } else {
      if (MESH_Num_Vertices(mesh) && (mtype == MVERTEX || mtype == MALLTYPE))
        MESH_Renumber_EntityGlobalIDs(mesh, MVERTEX, method, NULL, comm);
      if (MESH_Num_Edges(mesh) && (mtype == MEDGE || mtype == MALLTYPE))
        MESH_Renumber_EntityGlobalIDs(mesh, MEDGE, method, NULL, comm);
      if (MESH_Num_Faces(mesh) && (mtype == MFACE || mtype == MALLTYPE))
        MESH_Renumber_EntityGlobalIDs(mesh, MFACE, method, NULL, comm);
      if (MESH_Num_Regions(mesh) && (mtype == MREGION || mtype == MALLTYPE))
        MESH_Renumber_EntityGlobalIDs(mesh, MREGION, method, NULL, comm);

      // Ghost entity management in the mesh relies on ghost lists
      // being sorted according to their global ID

      MESH_Sort_GhostLists(mesh, compareGlobalID);
    }

    return 1;
  }
  

  int MESH_Renumber_EntityGlobalIDs(Mesh_ptr mesh, MType mtype,
                                    int method, int *preassigned_gids,
                                    MSTK_Comm comm) {
    int i;

    if (method != 0) {
      MSTK_Report("MESH_Renumber_EntityGlobalIDs",
                  "Chosen renumbering scheme not implemented", MSTK_ERROR);
      return 0;
    }
    if (mtype == MALLTYPE) {
      MSTK_Report("MESH_Renumber_EntityGlobalIDs",
                  "Cannot call this routine for MALLTYPE", MSTK_ERROR);
      return 0;
    }

    MAttrib_ptr tmpatt = MAttrib_New(mesh, "tmpatt_renumber", INT, mtype);
 
    int nproc, rank;
    MPI_Comm_size(comm, &nproc);
    MPI_Comm_rank(comm, &rank);

    int idx = 0, nowned = 0;
    MEntity_ptr ment;
    switch (mtype) {
      case MVERTEX:
        while ((ment = MESH_Next_Vertex(mesh, &idx)))
          if (MV_PType(ment) != PGHOST)
            nowned++;
        break;
      case MEDGE:
        while ((ment = MESH_Next_Edge(mesh, &idx)))
          if (ME_PType(ment) != PGHOST)
            nowned++;
        break;
      case MFACE:
        while ((ment = MESH_Next_Face(mesh, &idx)))
          if (MF_PType(ment) != PGHOST)
            nowned++;
        break;
      case MREGION:
        while ((ment = MESH_Next_Region(mesh, &idx)))
          if (MR_PType(ment) != PGHOST)
            nowned++;
        break;
      default: {}
    }

    /* Gather the number of entities on every processor */

    int *nowned_all = NULL;
    if (rank == 0)
      nowned_all = (int *) calloc(nproc, sizeof(int));

    MPI_Gather(&nowned, 1, MPI_INT, nowned_all, 1, MPI_INT, 0, comm);

    int *offset_all = NULL;
    if (rank == 0) {      
      offset_all = (int *) malloc(nproc*sizeof(int));
      offset_all[0] = 0;
      int p;
      for (p = 1; p < nproc; p++)
        offset_all[p] = offset_all[p-1] + nowned_all[p-1];
    }

    int offset = 0;
    MPI_Scatter(offset_all, 1, MPI_INT, &offset, 1, MPI_INT, 0, comm);

    /* At this point if we had different methods for renumbering the
     * global set of entities, we would generate a Global ID map. For
     * now this is just a sequential map (method 0) */

    int *new_index = (int *) malloc(nowned*sizeof(int));
    if (method == 0) {
      for (i = 0; i < nowned; i++)
        new_index[i] = i+1;
    }

    /* Now assign global IDs to the entities */
    int new_gid;
    idx = 0; i = 0;
    switch (mtype) {
      case MVERTEX:
        while ((ment = MESH_Next_Vertex(mesh, &idx)))
          if (MV_PType(ment) != PGHOST) {
            new_gid = offset + new_index[i];
            MEnt_Set_AttVal(ment, tmpatt, new_gid, 0.0, NULL);
            i++;
          }
        break;
      case MEDGE:
        while ((ment = MESH_Next_Edge(mesh, &idx)))
          if (ME_PType(ment) != PGHOST) {
            new_gid = offset + new_index[i];
            MEnt_Set_AttVal(ment, tmpatt, new_gid, 0.0, NULL);
            i++;
          }
        break;
      case MFACE:
        while ((ment = MESH_Next_Face(mesh, &idx)))
          if (MF_PType(ment) != PGHOST) {
            new_gid = offset + new_index[i];
            MEnt_Set_AttVal(ment, tmpatt, new_gid, 0.0, NULL);
            i++;
          }
        break;
      case MREGION:
        while ((ment = MESH_Next_Region(mesh, &idx)))
          if (MR_PType(ment) != PGHOST) {
            new_gid = offset + new_index[i];
            MEnt_Set_AttVal(ment, tmpatt, new_gid, 0.0, NULL);
            i++;
          }
        break;
      default: {}
    }

    /* Now exchange the attribute across processors */
    MESH_Update1Attribute(mesh, tmpatt, comm);

    /* Now assign global IDs of the entities based on the values of
     * the gidatt attribute */
    idx = 0;
    double rval;
    void *pval;
    switch (mtype) {
      case MVERTEX:
        while ((ment = MESH_Next_Vertex(mesh, &idx))) {
          MEnt_Get_AttVal(ment, tmpatt, &new_gid, &rval, &pval);
          MV_Set_GlobalID(ment, new_gid);
        }
        break;
      case MEDGE:
        while ((ment = MESH_Next_Edge(mesh, &idx))) {
          MEnt_Get_AttVal(ment, tmpatt, &new_gid, &rval, &pval);
          ME_Set_GlobalID(ment, new_gid);
        }
        break;
      case MFACE:
        while ((ment = MESH_Next_Face(mesh, &idx))) {
          MEnt_Get_AttVal(ment, tmpatt, &new_gid, &rval, &pval);
          MF_Set_GlobalID(ment, new_gid);
        }
        break;
      case MREGION:
        while ((ment = MESH_Next_Region(mesh, &idx))) {
          MEnt_Get_AttVal(ment, tmpatt, &new_gid, &rval, &pval);
          MR_Set_GlobalID(ment, new_gid);
        }
        break;
      default: {}
    }

    free(nowned_all);
    free(offset_all);

    MAttrib_Delete(tmpatt);

    return 1;
  }

#ifdef __cplusplus
}
#endif

