#include "MSTK.h"
#include "MSTK_private.h"

#ifdef __cplusplus
extern "C" {
#endif



  int compareINT(const void *a, const void *b) {
    return ( *(int*)a - *(int*)b );
  }

  int compareGlobalID(const void *a, const void *b) {
    return ( MEnt_GlobalID(*(MEntity_ptr*)a) - MEnt_GlobalID(*(MEntity_ptr*)b) );
  }


  int compareCoorDouble(const void * a, const void * b) {
    double tol = 1e-8;
    int i;
    double *coor1 = (double *)a;
    double *coor2 = (double *)b;
    for(i = 0; i < 3; i++) {
      if ( (coor1[i] - coor2[i]) > tol )
	return 1;
      if ( (coor2[i] - coor1[i]) > tol )
	return -1;
    }
    return 0;
  }
  
  int compareVertexCoor(const void *a, const void *b) {
    double coor1[3], coor2[3];
    double tol = 1e-8;
    int i;
    MV_Coords(*(MVertex_ptr*)a,coor1);
    MV_Coords(*(MVertex_ptr*)b,coor2);
    for(i = 0; i < 3; i++) {
      if ( (coor1[i] - coor2[i]) > tol )
	return 1;
      if ( (coor2[i] - coor1[i]) > tol )
	return -1;
    }
    return 0;
  }
  /* compare endpoint global id of an edge */
  int compareEdgeINT(const void *a, const void *b) {
    int *edge1_id = (int*)a; 
    int *edge2_id = (int*)b; 
    int maxgid1 = edge1_id[0];
    int mingid1 = edge1_id[1];
    int maxgid2 = edge2_id[0];
    int mingid2 = edge2_id[1];
    int tmp;
    if(maxgid1 < mingid1) {tmp = maxgid1; maxgid1 = mingid1; mingid1 = tmp;}
    if(maxgid2 < mingid2) {tmp = maxgid2; maxgid2 = mingid2; mingid2 = tmp;}
    if(maxgid1 > maxgid2) return 1;
    if(maxgid1 < maxgid2) return -1;
    if(mingid1 > mingid2) return 1;
    if(mingid1 < mingid2) return -1;

    return 0;
  }

  /* compare larger endpoint global id, then smaller endpoint global id */
  int compareEdgeID(const void *a, const void *b) {
    MEdge_ptr me1 = *(MEdge_ptr*)a;
    MEdge_ptr me2 = *(MEdge_ptr*)b;
    int maxgid1 = MV_GlobalID(ME_Vertex(me1,0));
    int mingid1 = MV_GlobalID(ME_Vertex(me1,1));
    int maxgid2 = MV_GlobalID(ME_Vertex(me2,0));
    int mingid2 = MV_GlobalID(ME_Vertex(me2,1));
    int tmp;
    if(maxgid1 < mingid1) {tmp = maxgid1; maxgid1 = mingid1; mingid1 = tmp;}
    if(maxgid2 < mingid2) {tmp = maxgid2; maxgid2 = mingid2; mingid2 = tmp;}
    if(maxgid1 > maxgid2) return 1;
    if(maxgid1 < maxgid2) return -1;
    if(mingid1 > mingid2) return 1;
    if(mingid1 < mingid2) return -1;

    return 0;
  }
  /* first compare # of vertices, then compare the largest vertex global ID */
  /*
  int compareFaceID(const void *a, const void *b) {
    MFace_ptr af = *(MFace_ptr*)a, bf = *(MFace_ptr*)b;
    List_ptr afv = MF_Vertices(af,0,1),bfv = MF_Vertices(bf,0,1);
    int nafv = List_Num_Entires(afv), nbfv = List_Num_Entries(bfv);
    for(i = 0; i < 3; i++) {
      if ( (coor1[i] - coor2[i]) > tol )
	return 1;
      if ( (coor2[i] - coor1[i]) > tol )
	return -1;
    }
    return 0;
  }
  */
  /* compare two faces, first on number of vertices, then on 
  int compareFaceVids(const void *a, const void *b) {
    double coor1[3], coor2[3];
    double tol = 1e-8;
    int i;
    MV_Coords(*(MVertex_ptr*)a,coor1);
    MV_Coords(*(MVertex_ptr*)b,coor2);
    for(i = 0; i < 3; i++) {
      if ( (coor1[i] - coor2[i]) > tol )
	return 1;
      if ( (coor2[i] - coor1[i]) > tol )
	return -1;
    }
    return 0;
  }
  */
#ifdef __cplusplus
}
#endif
