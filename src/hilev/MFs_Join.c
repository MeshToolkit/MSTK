#include <stdlib.h>
#include <stdio.h>
#include "MSTK.h"


MFace_ptr MFs_Join(MFace_ptr f1, MFace_ptr f2, MEdge_ptr e) {
  int i, k1, k2, nfe1, nfe2, *fedir1, *fedir2, *fedir3, gdim, gid;
  MEdge_ptr *fe1, *fe2, *fe3;
  MFace_ptr  nuface;
  Mesh_ptr   mesh;
  Set_ptr   fedges;

  mesh = MF_Mesh(f1);
  gid = MF_GEntID(f1);
  gdim = MF_GEntDim(f1);

  if (mesh != MF_Mesh(f2)) {
    MSTK_Report("MFs_Join","Faces not from same mesh",ERROR);
    return 0;
  }
  else if (gid != MF_GEntID(f1) || gdim != MF_GEntDim(f1)) {
    MSTK_Report("MFs_Join","Faces not from same geometric entity",ERROR);
    return 0;
  }

  nfe1 = MF_Num_Edges(f1);
  nfe2 = MF_Num_Edges(f2);

  fe1 = (MEdge_ptr *) malloc(nfe1*sizeof(MEdge_ptr));
  fedir1 = (int *) malloc(nfe1*sizeof(int));
  fe2 = (MEdge_ptr *) malloc(nfe2*sizeof(MEdge_ptr));
  fedir2 = (int *) malloc(nfe2*sizeof(int));
  fe3 = (MEdge_ptr *) malloc((nfe1+nfe2-2)*sizeof(MEdge_ptr));
  fedir3 = (int *) malloc((nfe1+nfe2-2)*sizeof(int));

  fedges = MF_Edges(f1,1,0);
  for (i = 0, k1 = -1; i < nfe1; i++) {
    fe1[i] = Set_Entry(fedges,i);
    fedir1[i] = MF_EdgeDir_i(f1,i);
    if (fe1[i] == e)
      k1 = i;
  }
  Set_Delete(fedges);

  fedges = MF_Edges(f2,1,0);
  for (i = 0, k2 = -1; i < nfe2; i++) {
    fe2[i] = Set_Entry(fedges,i);
    fedir2[i] = MF_EdgeDir_i(f2,i);
    if (fe2[i] == e)
      k2 = i;
  }
  Set_Delete(fedges);

  if (k1 == -1 || k2 == -1) {
    MSTK_Report("MFs_Join","Cannot find edge in face",ERROR);
    return 0;
  }

  for (i = 0; i < nfe1-1; i++) {
    fe3[i] = fe1[(k1+1+i)%nfe1];
    fedir3[i] = fedir1[(k1+1+i)%nfe1];
  }

  for (i = 0; i < nfe2-1; i++) {
    fe3[(nfe1-1)+i] = fe2[(k2+1+i)%nfe2];
    fedir3[(nfe1-1)+i] = fedir2[(k2+1+i)%nfe2];
  }

  
  MF_Delete(f1);
  MF_Delete(f2);
  ME_Delete(e);


  nuface = MF_New(mesh);
  MF_Set_Edges(nuface,(nfe1+nfe2-2),fe3,fedir3);

  MF_Set_GEntDim(nuface,gdim);
  MF_Set_GEntID(nuface,gid);

  free(fe1); free(fedir1);
  free(fe2); free(fedir2);
  free(fe3); free(fedir3);

  return nuface;
}
