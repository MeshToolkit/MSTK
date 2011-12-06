#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <MSTK.h>


#ifdef __cplusplus
extern "C" {
#endif


int MESH_InitFromGenDesc(Mesh_ptr mesh, int nv, double (*xyz)[3], int nf, 
			 int *nfv, int **fvids, int nr, int *nrv, 
			 int **rvids, int *nrf, int ***rfvtemplate) {

  int i, j, max_nfv, max_nrv;
  MVertex_ptr *verts, *fverts=NULL, *rverts=NULL;
  MFace_ptr mf;
  MRegion_ptr mr;


  verts = (MVertex_ptr *) malloc(nv*sizeof(MVertex_ptr));
  

  /* Create the mesh vertices */
              
  for (i = 0; i < nv; i++) {
    verts[i] = MV_New(mesh);
    MV_Set_Coords(verts[i],xyz[i]);
	
    MV_Set_ID(verts[i],i+1);
  }
      

  if (nf && nr) {
    
    MSTK_Report("MESH_InitFromGenDesc",
		"Cannot specify 'faces' and 'regions' for now",MSTK_WARN);
    return 0;
  }

  if (!nf && !nr) {
    
    MSTK_Report("MESH_InitFromGenDesc",
		"Mesh has only vertices? No faces or regions?",MSTK_WARN);

    return 0;
  }


  max_nfv = 0;
  for (i = 0; i < nf; i++)
    if (nfv[i] > max_nfv)
      max_nfv = nfv[i];

  if (max_nfv)
    fverts = (MVertex_ptr *) malloc(max_nfv*sizeof(MVertex_ptr));

  for (i = 0; i < nf; i++) {
      
    mf = MF_New(mesh);
    
    for (j = 0; j < nfv[i]; j++)
      fverts[j] = verts[fvids[i][j]];
    
    MF_Set_Vertices(mf,nfv[i],fverts);
    
    MF_Set_ID(mf,i+1);
  }

  if (fverts) free(fverts);


  if (nr && rfvtemplate)
    MSTK_Report("MESH_InitFromGenDesc",
		"General mesh import has not been tested yet",MSTK_WARN);

  max_nrv = 0;
  for (i = 0; i < nr; i++)
    if (nrv[i] > max_nrv)
      max_nrv = nrv[i];

  if (max_nrv)
    rverts = (MVertex_ptr *) malloc(max_nrv*sizeof(MVertex_ptr));

  for (i = 0; i < nr; i++) {

    if (!rfvtemplate) {
      /* Assume standard cell types (pyramid, prism, hexahedron) */

      if (nrv[i] <= 6 || nrv[i] == 8) {
	
	mr = MR_New(mesh);
	
	for (j = 0; j < nrv[i]; j++)
	  rverts[j] = verts[rvids[i][j]];
	
	MR_Set_Vertices(mr,nrv[i],rverts,0,NULL);
	
	MR_Set_ID(mr,i+1);
	
      }
      else {
	MSTK_Report("MESH_InitFromGenDesc",
		    "No standard mesh shape for regions with more than 8 vertices",MSTK_WARN);
	
	return 0;
	
      }
    }
    else {

      mr = MR_New(mesh);
      
      for (j = 0; j < nrv[i]; j++)
	rverts[j] = verts[rvids[i][j]];
      
      if (nrv[i] == 4) {
	/* No need to consult template here. It is just a tet */
	
	MR_Set_Vertices(mr,nrv[i],rverts,0,NULL);
	
      }
      else {
	/* Construct region according to given template information */
	
	MR_Set_Vertices(mr,nrv[i],rverts,nrf[i],rfvtemplate[i]);
	  
      }

      MR_Set_ID(mr,i+1);
	  
    }
  }


  if (verts) free(verts);

  return 1;
}

#ifdef __cplusplus
}
#endif
