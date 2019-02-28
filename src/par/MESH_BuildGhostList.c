
#define _H_Mesh_Private
#include "Mesh.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif

  /*
    Build ghost and overlap entity lists and sort them by global id

    If the lists are not empty, delete them first

    Author(s): Duo Wang, Rao Garimella
  */


  int MESH_Surf_Build_GhostLists(Mesh_ptr mesh);
  int MESH_Vol_Build_GhostLists(Mesh_ptr mesh);
  
  int MESH_Build_GhostLists(Mesh_ptr mesh, int topodim) {
    int nf, nr;
    /* basic mesh information */
    nf = MESH_Num_Faces(mesh);
    nr = MESH_Num_Regions(mesh);

    /* first delete the old GHOST or OVERLAP Lists */
    if ( mesh->ghvertex != (List_ptr)NULL ) {
      List_Delete(mesh->ghvertex);
      mesh->ghvertex = (List_ptr)NULL;
    }
    if ( mesh->ghedge != (List_ptr)NULL ) {
      List_Delete(mesh->ghedge);
      mesh->ghedge = (List_ptr)NULL;
    }
    if ( mesh->ghface != (List_ptr)NULL ) {
      List_Delete(mesh->ghface);
      mesh->ghface = (List_ptr)NULL;
    }
    if ( mesh->ghregion != (List_ptr)NULL ) {
      List_Delete(mesh->ghregion);
      mesh->ghregion = (List_ptr)NULL;
    }
    if ( mesh->ovvertex != (List_ptr)NULL ) {
      List_Delete(mesh->ovvertex);
      mesh->ovvertex = (List_ptr)NULL;
    }
    if ( mesh->ovedge != (List_ptr)NULL ) {
      List_Delete(mesh->ovedge);
      mesh->ovedge = (List_ptr)NULL;
    }
    if ( mesh->ovface != (List_ptr)NULL ) {
      List_Delete(mesh->ovface);
      mesh->ovface = (List_ptr)NULL;
    }
    if ( mesh->ovregion != (List_ptr)NULL ) {
      List_Delete(mesh->ovregion);
      mesh->ovregion = (List_ptr)NULL;
    }


    if (topodim == 3) 
      MESH_Vol_Build_GhostLists(mesh);
    else if (topodim == 2)
      MESH_Surf_Build_GhostLists(mesh);
    else
      MSTK_Report("MESH_BuildGhostList",
		  "This is not a valid mesh for building ghost list",
		  MSTK_FATAL);

    /* Ghost and Overlap entity lists must be sorted for efficient
       update of attributes */

    MESH_Sort_GhostLists(mesh,compareGlobalID);

    return 1;
  }


  int MESH_Surf_Build_GhostLists(Mesh_ptr mesh) {
    int count_ghost,count_ov, idx;
    int nv, ne, nf;
    MVertex_ptr mv;
    MEdge_ptr me;
    MFace_ptr mf;

    nv = MESH_Num_Vertices(mesh);
    ne = MESH_Num_Edges(mesh);
    nf = MESH_Num_Faces(mesh);
    
    /* build vertex */
    idx = 0; count_ghost = 0; count_ov = 0;

    while ((mv = MESH_Next_Vertex(mesh,&idx))) {
      if (MV_PType(mv) == PGHOST) {
	count_ghost++;
	MESH_Add_GhostVertex(mesh,mv);
      }
      else if (MV_PType(mv) == POVERLAP) {
	count_ov++;
	MESH_Add_OverlapVertex(mesh,mv);
      }
    }
    
#ifdef DEBUG_MAX
      printf("Mesh 0x%x num of vertex %d , ghost vertex %d, ov vertex %d\n",
	   (unsigned int) mesh, nv, count_ghost, count_ov);  
#endif
    
    /* build ghost edges (if the representation supports it) */
    
    idx = 0; count_ghost = 0; count_ov = 0;
    while ((me = MESH_Next_Edge(mesh,&idx))) {
      if (ME_PType(me) == PGHOST) {
	count_ghost++;
	MESH_Add_GhostEdge(mesh,me);
      }
      else if (ME_PType(me) == POVERLAP) {
	count_ov++;
	MESH_Add_OverlapEdge(mesh,me);
      }
    }

#ifdef DEBUG_MAX
    printf("Mesh 0x%x num of edges %d , ghost edge %d, ov edge %d\n",
	   (unsigned int) mesh, ne, count_ghost, count_ov);  
#endif

    /* build ghost faces */

    idx = 0; count_ghost = 0; count_ov = 0;
    while ((mf = MESH_Next_Face(mesh,&idx))) {
      if (MF_PType(mf) == PGHOST) {
	count_ghost++;
	MESH_Add_GhostFace(mesh,mf);
      }
      else if (MF_PType(mf) == POVERLAP) {
	count_ov++;
	MESH_Add_OverlapFace(mesh,mf);
      }
    }

#ifdef DEBUG_MAX
    printf("Mesh 0x%x num of face %d , ghost face %d, ov face %d\n",
	   (unsigned int) mesh, nf, count_ghost, count_ov);  
#endif

    return 1;
  }

  int MESH_Vol_Build_GhostLists(Mesh_ptr mesh) {
    int count_ghost,count_ov, idx;
    int nv, ne, nf, nr;
    MVertex_ptr mv;
    MEdge_ptr me;
    MFace_ptr mf;
    MRegion_ptr mr;
    

    nv = MESH_Num_Vertices(mesh);
    ne = MESH_Num_Edges(mesh);
    nf = MESH_Num_Faces(mesh);
    nr = MESH_Num_Regions(mesh);

    /* build ghost vertices */
    idx = 0; count_ghost = 0; count_ov = 0;
    while ((mv = MESH_Next_Vertex(mesh,&idx))) {
      if (MV_PType(mv) == PGHOST) {
	count_ghost++;
	MESH_Add_GhostVertex(mesh,mv);
      }
      else if(MV_PType(mv) == POVERLAP) {
	count_ov++;
	MESH_Add_OverlapVertex(mesh,mv);
      }
    }
#ifdef DEBUG_MAX
    printf("Mesh 0x%x  num of vertex %d , ghost vertex %d, ov vertex %d\n",
	   (unsigned int) mesh, nv, count_ghost, count_ov);  
#endif

    /* build ghost edges */
    idx = 0; count_ghost = 0; count_ov = 0;
    while ((me = MESH_Next_Edge(mesh,&idx))) {
      if (ME_PType(me) == PGHOST) {
	count_ghost++;
	MESH_Add_GhostEdge(mesh,me);
      }
      else if (ME_PType(me) == POVERLAP) {
	count_ov++;
	MESH_Add_OverlapEdge(mesh,me);
      }
    }
#ifdef DEBUG_MAX
    printf("Mesh 0x%x num of edges %d , ghost edges %d, ov edges %d\n", 
	   (unsigned int) mesh, ne, count_ghost, count_ov);  
#endif


    /* build ghost faces */
    idx = 0; count_ghost = 0; count_ov = 0;
    while ((mf = MESH_Next_Face(mesh,&idx))) {
      if (MF_PType(mf) == PGHOST) {
	count_ghost++;
	MESH_Add_GhostFace(mesh,mf);
      }
      else if (MF_PType(mf) == POVERLAP) {
	count_ov++;
	MESH_Add_OverlapFace(mesh,mf);
      }
    }

#ifdef DEBUG_MAX
    printf("Mesh 0x%x num of faces %d , ghost face %d, ov face %d\n", 
	   (unsigned int) mesh, nr, count_ghost, count_ov);  
#endif

    /* build ghost regions */
    idx = 0; count_ghost = 0; count_ov = 0;
    while ((mr = MESH_Next_Region(mesh,&idx))) {
      if (MR_PType(mr) == PGHOST) {
	count_ghost++;
	MESH_Add_GhostRegion(mesh,mr);
      }
      else if (MR_PType(mr) == POVERLAP) {
	count_ov++;
	MESH_Add_OverlapRegion(mesh,mr);
      }
    }

#ifdef DEBUG_MAX
    printf("Mesh 0x%x num of region %d , ghost region %d, ov region %d\n", 
	   (unsigned int) mesh, nr, count_ghost, count_ov);  
#endif

    return 1;
  }


#ifdef __cplusplus
}
#endif


