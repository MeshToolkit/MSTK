#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MSTK.h"
#include "metis.h"


/* Parition a MSTK mesh, and write out a parallel X3D file. Also,
   write out the map between the serial and parallel mesh cell
   numbers */


#ifdef __cplusplus
extern "C" {
#endif


int MESH_PartitionWithMetis(Mesh_ptr mesh, int nparts, int **part) {

  MEdge_ptr fedge;
  MFace_ptr mf, oppf, rface;
  MRegion_ptr mr, oppr;
  List_ptr fedges, efaces, rfaces, fregions;
  int  i, ncells;
  int  nv, ne, nf, nr, nfe, nef, nfr, nrf, idx, idx2;
  int  numflag, nedgecut, ipos;
  int  wtflag, metisopts[5] = {0,0,0,0,0};
  idxtype  *xadj, *adjncy, *vwgt, *adjwgt;

  

  /* First build a nodal graph of the mesh in the format required by
     metis */

  nv = MESH_Num_Vertices(mesh);
  ne = MESH_Num_Edges(mesh);
  nf = MESH_Num_Faces(mesh);
  nr = MESH_Num_Regions(mesh);

  ipos = 0;
  
  if (nr == 0) {
    
    if (nf == 0) {
      fprintf(stderr,"Cannot write out FLAG X3D files for wire meshes\n");
      exit(-1);
    }
    
    xadj = (int *) malloc((nf+1)*sizeof(int));
    adjncy = (int *) malloc(2*ne*sizeof(int));
    ncells = nf;

    /* Surface mesh */

    idx = 0; i = 0;
    xadj[i] = ipos;
    while ((mf = MESH_Next_Face(mesh,&idx))) {
      
      fedges = MF_Edges(mf,1,0);
      nfe = List_Num_Entries(fedges);
      
      idx2 = 0;
      while ((fedge = List_Next_Entry(fedges,&idx2))) {
	
	efaces = ME_Faces(fedge);
	nef = List_Num_Entries(efaces);
	
	if (nef > 2) {
	  fprintf(stderr,"Non-manifold surface mesh. Exit!\n");
	  exit(-1);
	}
	else if (nef == 1) {
	  continue;          /* boundary edge; nothing to do */
	}
	else {
	  oppf = List_Entry(efaces,0);
	  if (oppf == mf)
	    oppf = List_Entry(efaces,1);
	  
	  adjncy[ipos] = MF_ID(oppf)-1;
	  ipos++;
	}
	
	List_Delete(efaces);
	
      }
      
      List_Delete(fedges);
      
      i++;
      xadj[i] = ipos;
    }

  }
  else {

    xadj = (int *) malloc((nr+1)*sizeof(int));
    adjncy = (int *) malloc(2*nf*sizeof(int));
    ncells = nr;

    /* Volume mesh */

    idx = 0; i = 0;
    xadj[i] = ipos;
    while ((mr = MESH_Next_Region(mesh,&idx))) {
      
      rfaces = MR_Faces(mr);
      nrf = List_Num_Entries(rfaces);
      
      idx2 = 0;
      while ((rface = List_Next_Entry(rfaces,&idx2))) {
	
	fregions = MF_Regions(rface);
	nfr = List_Num_Entries(fregions);
	
	if (nfr == 1) {
	  continue;          /* boundary face; nothing to do */
	}
	else {
	  oppr = List_Entry(fregions,0);
	  if (oppr == mr)
	    oppr = List_Entry(fregions,1);
	  
	  adjncy[ipos] = MR_ID(oppr)-1;
	  ipos++;
	}
	
	List_Delete(fregions);
	
      }
      
      List_Delete(rfaces);
      
      i++;
      xadj[i] = ipos;
    }

  }
  


  /* Partition the graph */
  
  wtflag = 0;        /* No weights are specified */
  vwgt = adjwgt = NULL;

  numflag = 0;    /* C style numbering of elements (nodes of the dual graph) */
  *part = (int *) malloc(ncells*sizeof(int));

  if (nparts <= 8)
    METIS_PartGraphRecursive(&ncells,xadj,adjncy,vwgt,adjwgt,&wtflag,
			     &numflag,&nparts,metisopts,&nedgecut,*part);
  else
    METIS_PartGraphKway(&ncells,xadj,adjncy,vwgt,adjwgt,&wtflag,&numflag,
			&nparts,metisopts,&nedgecut,*part);
  
  return 1;
}

#ifdef __cplusplus
  }
#endif

