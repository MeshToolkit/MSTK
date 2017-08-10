#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "MSTK.h"
#include "MSTK_private.h"

#include "exodusII.h"
#ifdef EXODUSII_4
#include "exodusII_ext.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif


  /* Read an Exodus II file into MSTK */
  /* 
     The input arguments parallel_opts indicate if we should build
     parallel adjacencies and some options for controlling that 

     if parallel_opts == NULL, then only processor 0 will read the mesh and
     do nothing else
     
     parallel_opts[0] = 1/0  ---- distribute/Don't distribute the mesh
     parallel_opts[1] = 0-3  ---- 0: read/partition on P0 and distribute
                                  1: read/partition on P0, reread local portion
                                     on each process and weave/connect
                                  2: read on multiple processors Pn where n is 
                                  an arithmetic progression (5,10,15 etc) and
                                  distribute to a subset of processors (like P5
     Opts 1,2,3 NOT ACTIVE        distributes to P0-P4, P10 to P6-P9 etc)
                                  3: read portions of the mesh on each processor
                                  and repartition
     parallel_opts[2] = N    ---- Number of ghost layers around mesh on proc
     parallel_opts[3] = 0/1  ---- partitioning method
                                  0: Metis
                                  1: Zoltan
  */                            

  
  /* Right now we are creating an attribute for material sets, side
     sets and node sets. We are ALSO creating meshsets for each of
     these entity sets. We are also setting mesh geometric entity IDs
     - Should we pick one of the first two? */

  int MESH_ImportFromExodusII(Mesh_ptr mesh, const char *filename, int *parallel_opts, 
                              MSTK_Comm comm) {

  char mesg[256], funcname[32]="MESH_ImportFromExodusII";
  char title[256], elem_type[256], sidesetname[256], nodesetname[256];
  char matsetname[256];
  char **elem_blknames;
  int i, j, k, k1;
  int comp_ws = sizeof(double), io_ws = 0;
  int exoid=0, status;
  int ndim, nnodes, nelems, nelblock, nnodesets, nsidesets;
  int nedges, nedge_blk, nfaces, nface_blk, nelemsets;
  int nedgesets, nfacesets, nnodemaps, nedgemaps, nfacemaps, nelemmaps;
  int *elem_blk_ids, *connect, *node_map, *elem_map, *nnpe;
  int nelnodes, neledges, nelfaces;
  int nelem_i, natts;
  int *sideset_ids, *ss_elem_list, *ss_side_list, *nodeset_ids, *ns_node_list;
  int num_nodes_in_set, num_sides_in_set, num_df_in_set;

  double *xvals, *yvals, *zvals, xyz[3];
  float version;

  int exo_nrf[3] = {4,5,6};
  int exo_nrfverts[3][6] =
    {{3,3,3,3,0,0},{4,4,4,3,3,0},{4,4,4,4,4,4}};
  int exo_rfverts[3][6][4] =
    {{{0,1,3,-1},{1,2,3,-1},{2,0,3,-1},{2,1,0,-1},{-1,-1,-1,-1},{-1,-1,-1,-1}},
     {{0,1,4,3},{1,2,5,4},{2,0,3,5},{2,1,0,-1},{3,4,5,-1},{-1,-1,-1,-1}},
     {{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7},{3,2,1,0},{4,5,6,7}}};

  List_ptr fedges, rfaces;
  MVertex_ptr mv, *fverts, *rverts;
  MEdge_ptr me;
  MFace_ptr mf;
  MRegion_ptr mr;
  MAttrib_ptr nmapatt=NULL, elblockatt=NULL, nodesetatt=NULL, sidesetatt=NULL;
  MSet_ptr faceset=NULL, nodeset=NULL, sideset=NULL, matset=NULL;
  int distributed=0;
  
  ex_init_params exopar;

  FILE *fp;
  char basename[256], modfilename[256], *ext;

  
#ifdef MSTK_HAVE_MPI

  int rank=0, numprocs=1;

  if (comm) {
    MPI_Comm_size(comm,&numprocs);
    MPI_Comm_rank(comm,&rank);
  }

  if (numprocs > 1 && parallel_opts) {

    if (parallel_opts[0] == 1) { /* distribute the mesh */

      int num_ghost_layers = parallel_opts[2];
      
      int have_zoltan = 0, have_metis = 0;
#ifdef _MSTK_HAVE_ZOLTAN
      have_zoltan = 1;
#endif
#ifdef _MSTK_HAVE_METIS
      have_metis = 1;
#endif
      
      PartitionMethod part_method = parallel_opts[3];
      
      if (part_method == METIS) {
        if (!have_metis) {
          MSTK_Report(funcname,"Metis not available. Trying Zoltan with GRAPH",MSTK_WARN);
          part_method = ZOLTAN_GRAPH;
          if (!have_zoltan) 
            MSTK_Report(funcname,"No partitioner defined",MSTK_FATAL);
        }
      }
      else if (part_method == ZOLTAN_GRAPH || part_method == ZOLTAN_RCB) {
        if (!have_zoltan) {
          MSTK_Report(funcname,"Zoltan not available. Trying Metis",MSTK_WARN);
          part_method = METIS;
          if (!have_metis) 
            MSTK_Report(funcname,"No partitioner defined",MSTK_FATAL);
        }
      }
      else if (part_method == METIS_COLUMNAR) {
        if (!have_metis) {
          MSTK_Report(funcname,
            "Metis not available. Trying Zoltan with RCB for columnar meshes",MSTK_WARN);
          part_method = ZOLTAN_RCB_COLUMNAR;
          if (!have_zoltan)
            MSTK_Report(funcname,"No partitioner defined",MSTK_FATAL);
        }
      }
      else if (part_method == ZOLTAN_GRAPH_COLUMNAR ||
               part_method == ZOLTAN_RCB_COLUMNAR) {
        if (!have_zoltan) {
          MSTK_Report(funcname,
            "Zoltan not available. Trying Metis for columnar meshes",MSTK_WARN);
          part_method = METIS_COLUMNAR;
          if (!have_metis)
            MSTK_Report(funcname,"No partitioner defined",MSTK_FATAL);
        }
      }
      
      if (parallel_opts[1] == 0) {
        
        /* Read on processor 0 and distribute to all other processors */
        
        Mesh_ptr globalmesh;
        int topodim;
        
        if (rank == 0) { /* Read only on processor 0 */
          
          globalmesh = MESH_New(MESH_RepType(mesh));
          int read_status = MESH_ReadExodusII_Serial(globalmesh,filename,rank);
          
          if (!read_status) {
            sprintf(mesg,"Could not read Exodus II file %s successfully\n",
                    filename);
            MSTK_Report(funcname,mesg,MSTK_FATAL);
          }
          
          topodim = (MESH_Num_Regions(globalmesh) == 0) ? 2 : 3;
        }
        else
          globalmesh = NULL;
        
        int with_attr = 1;      
        int del_inmesh = 1;
        int dist_status = MSTK_Mesh_Distribute(globalmesh, &mesh, &topodim, 
                                               num_ghost_layers,
                                               with_attr, part_method, 
                                               del_inmesh, comm);
        if (!dist_status)
          MSTK_Report(funcname,
                      "Could not distribute meshes to other processors",
                      MSTK_FATAL);
        
      }
      else if (parallel_opts[1] == 1) {
        
      }
      else if (parallel_opts[1] == 2) {
        
      }
      else if (parallel_opts[1] == 3) {
        
      }

      int parallel_check = MESH_Parallel_Check(mesh,comm);        
      
      if (!parallel_check)
        MSTK_Report(funcname,
                    "Parallel mesh checks failed",
                    MSTK_FATAL);

    }
    else { /* Read the mesh on rank 0 but don't distribute */
      if (rank == 0)
        return (MESH_ReadExodusII_Serial(mesh,filename,0));
    }
  
  } /* if numprocs > 1 */
  else {
    if (rank == 0) 
      return (MESH_ReadExodusII_Serial(mesh,filename,0));
    else
      return 1;
  }

#else /* Serial process */
  return (MESH_ReadExodusII_Serial(mesh,filename,0));
#endif 
  
  return 1;

}


#ifdef __cplusplus
}
#endif
