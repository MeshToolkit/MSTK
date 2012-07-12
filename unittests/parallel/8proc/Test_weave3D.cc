#include <UnitTest++.h>

#include <algorithm>

#include "../../../include/MSTK.h"

//
//  Build a parallel connected mesh from the following distributed
//  meshes and test for various parallel entities. The meshes on the
//  eight processors have been exploded for clarity. Each edge is of
//  unit length. Lower left corner is at 0,0,0 Pn shows the processor
//  number. Numbers shown are the expected global ID numbers. Vertex
//  numbers in parentheses indicates a master vertex
//                                                                            
//                                                                          
//                       *------*------*. .. . . . . . . . . . . *------*   
//                      /|     /|     /|                        /|     /|   
//                     / | 7  / | 8  / |                       / | 12 / |   
//                   .*------*------*  |                      *-----.*  |   
//                  . |  *---|--*---|--*                      |  *-.-|--*
//                 .  | /    | /    | /                       | / .  | /.   
//                .   |/     |/     |/                        |/ .   |/ .   
//               .    *------*------*                         *-.----*  .   
//              .        .   P4                                . P7     . 
//             .         .                                    .         .
//            .          .                                   .          .
//           .           .                                  .           .
//          .            .                                 .            .
//         .             .                                .             .
//        .              .                               .              .   
//       *------*------* .                       *------*               .   
//      /|     /|     /| .                      /|     /|               .   
//     / | 10 / | 11 / | .                     / | 9  / |               .   
//    *------*------*.........................*------*  |               .   
//    |  *---|--*---|--* .                    |  *---|--*               .   
//    | /    | /    | /  .                    | /    | /                .   
//    |/     |/     |/   .                    |/     |/                 .   
//    *------*------*    .                    *------*                  .   
//    .      P6          .                       P5  .                  .   
//    .                  .                           .                  .
//    .                  .                           .                  .
//    .                  .                           .                  .   
//    .                  *------*------*             .           *------*   
//    .                 /|     /|     /|             .          /|     /|   
//    .                / | 4  / | 5  / |             .         / | 6  / |   
//    .               *------*------*  |             .        *------*  |   
//    .               |  *---|--*---|--*. .. . . . . . . . . .|. *---|--*
//    .               | /    | /    | /              .        | /    | /    
//    .               |/     |/     |/               .        |/     |/     
//    .              .*------*------*                .        *------*      
//    .             .       P2                       .           P3.    
//    .            .                                 .            .      
//    .           .                                  .           .      
//    .          .                                   .          .       
//    .         .                                    .         .        
//    .        .                                     .        .         
//    .       .                                      .       .              
//    .  *------*------*                         *------*   .               
//    . /|  .  /|     /|                        /|   . /|  .                
//    ./ | 1  / |  2 / |                       / | 3 ./ | .                 
//    *------*------*  |                      *------*  |.                  
//    |  *---|--*---|--*                      |  *---|--*                   
//    | /    | /    | /                       | /    | /                    
//    |/     |/     |/                        |/     |/                     
//    *------*------*.........................*------*                      
//           P0                                   P1                      
//                                                                       
//                                                                            
//                                                                            
//                       (28)   (29)   (30)                       30    (36)
//                       *------*------*. .. . . . . . . . . . . *------*   
//                      /|     /|     /|                        /|     /|   
//                  (25) | (26) | (27) |                    27 / |  33/ |   
//                   .*------*------*  |                      *-----.*  |   
//                  . |  *---|--*---|--*                      |  *-.-|--*
//                 .  | /20  | /21  | /22                     | /22  | /24  
//                .   |/     |/     |/                        |/ .   |/ .   
//               .    *------*------*                         *-.----*  .   
//              .    10  .  11     12                        12.    16  . 
//             .         .                                    .         .
//            .          .                                   .          .
//           .           .                                  .           .
//          .            .                                 .            .
//         .             .                                .             .
//        25     26     27                        27    (33)            .   
//       *------*------* .                       *------*               .   
//      /|     /|     /| .                      /|     /|               .   
//  (34) | (35) |  31/ | .                  (31) | (32) |               .   
//    *------*------*.........................*------*  |               .   
//    |  *---|--*---|--* .                    |  *---|--*               .   
//    | /10  | /11  | /12.                    | /12  | /16              .   
//    |/     |/     |/   .                    |/     |/                 .   
//    *------*------*    .                    *------*                  .   
//    7      8      9    .                    9      15                 .   
//    .                  .                           .                  .
//    .                  .                           .                  .
//    .                  (20)   (21)   (22)          .            22    (24)
//    .                  *------*------*             .           *------*   
//    .                 /|     /|     /|             .          /|     /|   
//    .              10/ |  11/ |  12/ |             .       12/ |  16/ |   
//    .               *------*------*  |             .        *------*  |   
//    .               |  *---|--*---|--*. .. . . . . . . . . .|. *---|--*
//    .               | (17) | (18) | (19)           .        | /19  | /(23)
//    .               |/     |/     |/               .        |/     |/     
//    .              .*------*------*                .        *------*      
//    .             . 4      5      6                .        6    . 14 
//    .            .                                 .            .      
//    .           .                                  .           .      
//    .          .                                   .          .       
//    .         .                                    .         .        
//    .        .                                     .        .         
//    .  (10) . (11)   (12)                       12 .  (16) .              
//    .  *------*------*                         *------*   .               
//    . /|  .  /|     /|                        /|   . /|  .                
//  (7)/ | (8)/ | (9)/ |                     9 / | (15) | .                 
//    *------*------*  |                      *------*  |.                  
//    |  *---|--*---|--*                      |  *---|--*                   
//    | (4)  | (5)  | (6)                     | /6   | /(14)                
//    |/     |/     |/                        |/     |/                     
//    *------*------*.........................*------*                      
//   (1)    (2)    (3)                        3       (13)                
//                                                                       
//                                                                            
                                                                              

SUITE(Parallel) {
TEST(Weave3D_from_MSTK) {

  int i, idx;
  int nr, nv, ngr, nor, ngv, nov;
  int nproc, rank, status, dim;
  int *regids, *gregids, *oregids;
  int *vertexids, *gvertexids, *overtexids;
  Mesh_ptr mesh;
  MRegion_ptr mr;
  MVertex_ptr mv;
  char filename[256];


  int expnr[8]={12,8,12,8,12,8,12,8};
  int expngr[8]={10,7,10,7,10,7,10,7};
  int expnor[8]={2,1,2,1,2,1,2,1};

  int expregids[8][12]={
    {1,2,3,4,5,6,7,8,9,10,11,12},
    {3,2,5,6,8,9,11,12,0,0,0,0},
    {4,5,1,2,3,6,7,8,9,10,11,12},
    {6,2,3,5,8,9,11,12,0,0,0,0},
    {7,8,1,2,3,4,5,6,9,10,11,12},
    {9,2,3,5,6,8,11,12,0,0,0,0},
    {10,11,1,2,3,4,5,6,7,8,9,12},
    {12,2,3,5,6,8,9,11,0,0,0,0}};

  int expgregids[8][10]={
    {3,4,5,6,7,8,9,10,11,12},
    {2,5,6,8,9,11,12,0,0,0},
    {1,2,3,6,7,8,9,10,11,12},
    {2,3,5,8,9,11,12,0,0,0},
    {1,2,3,4,5,6,9,10,11,12},
    {2,3,5,6,8,11,12,0,0,0},
    {1,2,3,4,5,6,7,8,9,12},
    {2,3,5,6,8,9,11,0,0,0}};

  int exporegids[8][2]={
    {1,2},
    {3,0},
    {4,5},
    {6,0},
    {7,8},
    {9,0},
    {10,11},
    {12,0}};


  int expnv[8] = {36,27,36,27,36,27,36,27};

  int expvertexids[8][36]={
    {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36},
    {2,3,5,6,8,9,11,12,13,14,15,16,18,19,21,22,23,24,26,27,29,30,31,32,33,35,36,0,0,0,0,0,0,0,0,0},
    {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36},
    {2,3,5,6,8,9,11,12,13,14,15,16,18,19,21,22,23,24,26,27,29,30,31,32,33,35,36,0,0,0,0,0,0,0,0,0},
    {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36},
    {2,3,5,6,8,9,11,12,13,14,15,16,18,19,21,22,23,24,26,27,29,30,31,32,33,35,36,0,0,0,0,0,0,0,0,0},
    {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36},
    {2,3,5,6,8,9,11,12,13,14,15,16,18,19,21,22,23,24,26,27,29,30,31,32,33,35,36,0,0,0,0,0,0,0,0,0}
  };

  int expngv[8] = {24,23,30,25,30,24,34,26};

  int expgvertexids[8][34]={
    {13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,0,0,0,0,0,0,0,0,0,0},
    {2,3,5,6,8,9,11,12,18,19,21,22,23,24,26,27,29,30,31,32,33,35,36,0,0,0,0,0,0,0,0,0,0,0},
    {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,23,24,25,26,27,28,29,30,31,32,33,34,35,36,0,0,0,0},
    {2,3,5,6,8,9,11,12,13,14,15,16,18,19,21,22,26,27,29,30,31,32,33,35,36,0,0,0,0,0,0,0,0,0},
    {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,31,32,33,34,35,36,0,0,0,0},
    {2,3,5,6,8,9,11,12,13,14,15,16,18,19,21,22,23,24,26,27,29,30,35,36,0,0,0,0,0,0,0,0,0,0},
    {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,36},
    {2,3,5,6,8,9,11,12,13,14,15,16,18,19,21,22,23,24,26,27,29,30,31,32,33,35,0,0,0,0,0,0,0,0}
  };

  int expnov[12] = {12,4,6,2,6,3,2,1};

  int expovertexids[12][36]={
    {1,2,3,4,5,6,7,8,9,10,11,12},
    {13,14,15,16,0,0,0,0,0,0,0,0},
    {17,18,19,20,21,22,0,0,0,0,0,0},
    {23,24,0,0,0,0,0,0,0,0,0,0},
    {25,26,27,28,29,30,0,0,0,0,0,0},
    {31,32,33,0,0,0,0,0,0,0,0,0},
    {34,35,0,0,0,0,0,0,0,0,0,0},
    {36,0,0,0,0,0,0,0,0,0,0,0}
  };


  MSTK_Init();
  MSTK_Set_Comm(MPI_COMM_WORLD);

  int debugwait=0;
  while (debugwait);


  nproc = MSTK_Comm_size();
  rank = MSTK_Comm_rank();


  mesh = MESH_New(UNKNOWN_REP);

  sprintf(filename,"parallel/8proc/hex3x2x2.mstk.%-1d",rank);
  status = MESH_InitFromFile(mesh,filename);

  CHECK(status);

  CHECK(MESH_Num_Regions(mesh));


  int input_type = 0;  /* no parallel info present in meshes */
  int num_ghost_layers = 1; /* always */
  status = MSTK_Weave_DistributedMeshes(mesh, num_ghost_layers, input_type);

  CHECK(status);


  nr = MESH_Num_Regions(mesh); /* includes ghost regionss */

  CHECK_EQUAL(expnr[rank],nr);

  regids = (int *) malloc(nr*sizeof(int));

  idx = 0; i = 0;
  while ((mr = MESH_Next_Region(mesh,&idx)))
    regids[i++] = MEnt_GlobalID(mr);

  CHECK_ARRAY_EQUAL(expregids[rank],regids,nr);

  free(regids);


  nv = MESH_Num_Vertices(mesh); /* includes ghost vertices */

  CHECK_EQUAL(expnv[rank],nv);

  vertexids = (int *) malloc(nv*sizeof(int));

  idx = 0; i = 0;
  while ((mv = MESH_Next_Vertex(mesh,&idx))) {
    vertexids[i++] = MEnt_GlobalID(mv);
    double xyz[3];
    MV_Coords(mv,xyz);
  }

  std::sort(vertexids,vertexids+nv);
  CHECK_ARRAY_EQUAL(expvertexids[rank],vertexids,nv);

  free(vertexids);



  ngr = MESH_Num_GhostRegions(mesh);
  CHECK_EQUAL(expngr[rank],ngr);

  gregids = (int *) malloc(ngr*sizeof(int));

  idx = 0; i = 0;
  while ((mr = MESH_Next_GhostRegion(mesh,&idx)))
    gregids[i++] = MEnt_GlobalID(mr);

  CHECK_ARRAY_EQUAL(expgregids[rank],gregids,ngr);

  free(gregids);

  nor = MESH_Num_OverlapRegions(mesh);
  CHECK_EQUAL(expnor[rank],nor);


  oregids = (int *) malloc(nor*sizeof(int));

  idx = 0; i = 0;
  while ((mr = MESH_Next_OverlapRegion(mesh,&idx)))
    oregids[i++] = MEnt_GlobalID(mr);

  CHECK_ARRAY_EQUAL(exporegids[rank],oregids,nor);

  free(oregids);


  ngv = MESH_Num_GhostVertices(mesh);
  CHECK_EQUAL(expngv[rank],ngv);

  gvertexids = (int *) malloc(nv*sizeof(int));

  idx = 0; i = 0;
  while ((mv = MESH_Next_GhostVertex(mesh,&idx)))
    gvertexids[i++] = MEnt_GlobalID(mv);

  std::sort(gvertexids,gvertexids+ngv);
  CHECK_ARRAY_EQUAL(expgvertexids[rank],gvertexids,ngv);

  free(gvertexids);


  nov = MESH_Num_OverlapVertices(mesh);
  CHECK_EQUAL(expnov[rank],nov);

  overtexids = (int *) malloc(nov*sizeof(int));

  idx = 0; i = 0;
  while ((mv = MESH_Next_OverlapVertex(mesh,&idx)))
    overtexids[i++] = MEnt_GlobalID(mv);

  std::sort(overtexids,overtexids+nov);
  CHECK_ARRAY_EQUAL(expovertexids[rank],overtexids,nov);

  free(overtexids);

  CHECK_EQUAL(1,MESH_Parallel_Check(mesh));

  return;
}


}
