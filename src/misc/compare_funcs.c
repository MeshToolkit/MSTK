/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdlib.h>
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

  int compareID(const void *a, const void *b) {
    return ( MEnt_ID(*(MEntity_ptr*)a) - MEnt_ID(*(MEntity_ptr*)b) );
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

  int compareFaceINT(const void *a, const void *b) {
    int *faceid1 = (int*)a, *faceid2 = (int*)b;
    int nfv1 = faceid1[0], nfv2 = faceid2[0];
    int i;
    if( nfv1 > nfv2 ) return 1;
    if( nfv1 < nfv2 ) return -1;
    /*
    printf("before sort 1: ");
    for(i = 0; i < nfv1; i++)
      printf("Face ID: %d \t",faceid1[i+1]);
    */
    qsort(faceid1+1,nfv1,sizeof(int),compareINT);
    qsort(faceid2+1,nfv2,sizeof(int),compareINT);
    /*
    printf("after sort 1: ");
    for(i = 0; i < nfv1; i++)
      printf("Face ID: %d \t",faceid1[i+1]);
    */
    for(i = 0; i < nfv1; i++) {
      if ( faceid1[i+1] > faceid2[i+1] )
	return 1;
      if ( faceid1[i+1] < faceid2[i+1] )
	return -1;
    }
    return 0;
  }

  int compareFaceID(const void *a, const void *b) {
    MFace_ptr mf1 = *(MFace_ptr*)a, mf2 = *(MFace_ptr*)b;
    List_ptr mfverts1 = MF_Vertices(mf1,1,0), mfverts2 = MF_Vertices(mf2,1,0);
    int nfv1 = List_Num_Entries(mfverts1), nfv2 = List_Num_Entries(mfverts2);
    int i;
    if( nfv1 > nfv2 ) {
      List_Delete(mfverts1);
      List_Delete(mfverts2);
      return 1;
    }
    if( nfv1 < nfv2 ) {
      List_Delete(mfverts1);
      List_Delete(mfverts2);
      return -1;
    }
    List_Sort(mfverts1,nfv1,sizeof(MFace_ptr),compareGlobalID);
    List_Sort(mfverts2,nfv2,sizeof(MFace_ptr),compareGlobalID);
    for(i = 0; i < nfv1; i++) {
      if ( MV_GlobalID(List_Entry(mfverts1,i)) > MV_GlobalID(List_Entry(mfverts2,i)) ) {
        List_Delete(mfverts1);
        List_Delete(mfverts2);
	return 1;
      }
      if ( MV_GlobalID(List_Entry(mfverts1,i)) < MV_GlobalID(List_Entry(mfverts2,i)) ) {
        List_Delete(mfverts1);
        List_Delete(mfverts2);
	return -1;
      }
    }
    List_Delete(mfverts1);
    List_Delete(mfverts2);
    return 0;
  }


  int compareValence(const void *a, const void *b) {
    return ( MV_Num_Edges(*((MVertex_ptr*)a)) - MV_Num_Edges(*((MVertex_ptr*)b)) );
  }



#ifdef __cplusplus
}
#endif
