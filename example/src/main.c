#include <stdio.h>
#include <stdlib.h>
#include "MSTK.h"

#include "test.h"


int main(int argc, char *argv[]) {
  int i, idx, idx2, ok, edir, nv, ne;
  double xyz[3];
  char meshname[256];
  Mesh_ptr mesh;
  MVertex_ptr v;
  MEdge_ptr e;
  MFace_ptr f;
  GEntity_ptr gent;
  Set_ptr fedges;

  
  fprintf(stderr,"\n");


  /* Initialize MSTK - Always do this even if it does seem to
     matter in this version of MSTK */

  MSTK_Init();

  /* Load the mesh */

  strcpy(meshname,argv[1]);
  strcat(meshname,".mstk");

  mesh = MESH_New(UNKNOWN_REP);
  ok = MESH_InitFromFile(mesh,meshname);
  if (!ok) {
    fprintf(stderr,"Cannot file input file %s\n\n\n",meshname);
    exit(-1);
  }


  /* Print some info about the mesh */


  nv = MESH_Num_Vertices(mesh);
  for (i = 0; i < nv; i++) {
    v = MESH_Vertex(mesh,i);

    /* Basic info */
    fprintf(stderr,"\n");
    fprintf(stderr,"Vertex: 0x%-x   ID: %-d   ",v,MV_ID(v));
    
    
    /* Classification w.r.t. geometric model */
      
    if (MV_GEntDim(v) == -1)
      fprintf(stderr,"Unknown Classification\n");
    else {
      fprintf(stderr,"GEntID: %-d    GEntDim: %-d\n", MV_GEntID(v),
	      MV_GEntDim(v));
      if ((gent = MV_GEntity(v)))
	fprintf(stderr,"Model entity pointer: 0x%-x\n",gent);	
    }


    /* Coordinates */
    MV_Coords(v,xyz);
    fprintf(stderr,"Coordinates: %16.8lf %16.8lf %16.8lf\n",xyz[0],xyz[1],xyz[2]);
  }



  idx = 0;
  while (f = MESH_Next_Face(mesh,&idx)) {

    /* Basic info */
    fprintf(stderr,"\n");
    fprintf(stderr,"Face: 0x%-x   ID: %-d   ",f,MF_ID(f));
    
    
    /* Classification w.r.t. geometric model */

    if (MF_GEntDim(f) == -1)
      fprintf(stderr,"Unknown Classification\n");
    else {
      fprintf(stderr,"GEntID: %-d    GEntDim: %-d\n", MF_GEntID(f),
	      MF_GEntDim(f));
      if ((gent = MF_GEntity(f)))
	fprintf(stderr,"Model entity pointer: 0x%-x\n",gent);
    }

    fprintf(stderr,"\n");


    /* Edges of face */
    fedges = MF_Edges(f,1,0);
    ne = Set_Num_Entries(fedges);
    fprintf(stderr,"Edges: %-d\n",ne);
    fprintf(stderr,"Object        ID      GEntID   GEntDim    Vertex IDs\n");
    idx2 = 0; i = 0;
    while (e = Set_Next_Entry(fedges,&idx2)) {
      edir = MF_EdgeDir_i(f,i);
      if (edir) 
	fprintf(stderr,"0x%-8x    %-8d %-8d     %-1d       %-d  %-d\n",
		e,ME_ID(e),ME_GEntID(e),ME_GEntDim(e),
		MV_ID(ME_Vertex(e,0)),MV_ID(ME_Vertex(e,1)));
      else 
	fprintf(stderr,"0x%-8x    %-8d %-8d     %-1d       %-d  %-d\n",
		e,-ME_ID(e),ME_GEntID(e),ME_GEntDim(e),
		MV_ID(ME_Vertex(e,0)),MV_ID(ME_Vertex(e,1)));
      i++;
    }
    fprintf(stderr,"\n");
    Set_Delete(fedges);
  }


  /* Write out a copy of the mesh */

  strcpy(meshname,argv[1]);
  strcat(meshname,"-copy.mstk");
  MESH_WriteToFile(mesh,meshname);

  
  /* Deleting of mesh is not necessary if this the end of the program */

  MESH_Delete(mesh);

  return 1;
}
