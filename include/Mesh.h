#ifndef _H_Mesh
#define _H_Mesh

#include "MSTK_defines.h"
#include "MSTK_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _H_Mesh_Private
  typedef struct Mesh {
    RepType reptype;
    int nv, ne, nf, nr;
    List_ptr mvertex, medge, mface, mregion;
    GModel_ptr geom;
    List_ptr AttribList;
    List_ptr MSetList;
    int max_vid, max_eid, max_fid, max_rid;
    Hash_ptr hedge, hface;
    int autolock;

#ifdef MSTK_HAVE_MPI
    unsigned int mypartn, numpartns;
    unsigned int *par_adj_flags;  /* See note below for explanation */
    unsigned int *par_recv_info;  /* of these two variables */
    List_ptr ghvertex, ghedge, ghface, ghregion;
    int max_ghvid, max_gheid, max_ghfid, max_ghrid;
    List_ptr ovvertex, ovedge, ovface, ovregion;
#endif

  } Mesh, *Mesh_ptr;
#else
   typedef void *Mesh;
#endif



  Mesh_ptr   MESH_New(RepType type);
  int        MESH_InitFromFile(Mesh_ptr mesh, const char *filename);
  int        MESH_WriteToFile(Mesh_ptr mesh, const char *filename, RepType rtype);
  void       MESH_Delete(Mesh_ptr mesh);
  
  GModel_ptr MESH_GModel(Mesh_ptr mesh);
  RepType    MESH_RepType(Mesh_ptr mesh);
  char      *MESH_RepType_Str(Mesh_ptr mesh);

  int         MESH_Num_Attribs(Mesh_ptr mesh);
  MAttrib_ptr MESH_Attrib(Mesh_ptr mesh, int i);
  MAttrib_ptr MESH_Next_Attrib(Mesh_ptr mesh, int *index);
  MAttrib_ptr MESH_AttribByName(Mesh_ptr mesh, const char *name);
  void        MESH_Add_Attrib(Mesh_ptr mesh, MAttrib_ptr attrib);
  void        MESH_Rem_Attrib(Mesh_ptr mesh, MAttrib_ptr attrib);

  int         MESH_Num_MSets(Mesh_ptr mesh);
  MSet_ptr    MESH_MSet(Mesh_ptr mesh, int i);
  MSet_ptr    MESH_Next_MSet(Mesh_ptr mesh, int *index);
  MSet_ptr    MESH_MSetByName(Mesh_ptr mesh, const char *name);
  void        MESH_Add_MSet(Mesh_ptr mesh, MSet_ptr mset);
  void        MESH_Rem_MSet(Mesh_ptr mesh, MSet_ptr mset);


  int        MESH_Num_Vertices(Mesh_ptr mesh);
  int        MESH_Num_Edges(Mesh_ptr mesh);
  int        MESH_Num_Faces(Mesh_ptr mesh);
  int        MESH_Num_Regions(Mesh_ptr mesh);
  
  MVertex_ptr  MESH_Vertex(Mesh_ptr mesh, int i);
  MEdge_ptr    MESH_Edge(Mesh_ptr mesh, int i);
  MFace_ptr    MESH_Face(Mesh_ptr mesh, int i);
  MRegion_ptr  MESH_Region(Mesh_ptr mesh, int i);
  
  MVertex_ptr  MESH_Next_Vertex(Mesh_ptr mesh, int *index);
  MEdge_ptr    MESH_Next_Edge(Mesh_ptr mesh, int *index);
  MFace_ptr    MESH_Next_Face(Mesh_ptr mesh, int *index);
  MRegion_ptr  MESH_Next_Region(Mesh_ptr mesh, int *index);

  void       MESH_Add_Vertex(Mesh_ptr mesh, MVertex_ptr v);
  void       MESH_Add_Edge(Mesh_ptr mesh, MEdge_ptr e);
  void       MESH_Add_Face(Mesh_ptr mesh, MFace_ptr f);
  void       MESH_Add_Region(Mesh_ptr mesh, MRegion_ptr r);

  void       MESH_Rem_Vertex(Mesh_ptr mesh, MVertex_ptr v);
  void       MESH_Rem_Edge(Mesh_ptr mesh, MEdge_ptr e);
  void       MESH_Rem_Face(Mesh_ptr mesh, MFace_ptr f);
  void       MESH_Rem_Region(Mesh_ptr mesh, MRegion_ptr r);

#ifdef MSTK_HAVE_MPI
  void         MESH_Set_Prtn(Mesh_ptr mesh, unsigned int partition, 
                             unsigned int numpartitions);
  unsigned int MESH_Prtn(Mesh_ptr mesh);

  void         MESH_Flag_Has_Ghosts_From_Prtn(Mesh_ptr mesh, unsigned int prtn,
                                              MType mtype);  
  unsigned int MESH_Has_Ghosts_From_Prtn(Mesh_ptr mesh, unsigned int prtn, 
                                         MType mtype);
  void         MESH_Flag_Has_Overlaps_On_Prtn(Mesh_ptr mesh, unsigned int prtn, 
                                              MType mtype);
  unsigned int MESH_Has_Overlaps_On_Prtn(Mesh_ptr mesh, unsigned int prtn, 
                                         MType mtype);

  unsigned int MESH_Num_GhostPrtns(Mesh_ptr mesh);
  void         MESH_GhostPrtns(Mesh_ptr mesh, unsigned int *pnums);

  
  void         MESH_Set_Num_Recv_From_Prtn(Mesh_ptr mesh, unsigned int prtn, 
                                           MType mtype, unsigned int numrecv);
  unsigned int MESH_Num_Recv_From_Prtn(Mesh_ptr mesh, unsigned int prtn, 
                                       MType mtype);

  int        MESH_Num_GhostVertices(Mesh_ptr mesh);
  int        MESH_Num_GhostEdges(Mesh_ptr mesh);
  int        MESH_Num_GhostFaces(Mesh_ptr mesh);
  int        MESH_Num_GhostRegions(Mesh_ptr mesh);


  List_ptr  MESH_GhostVertex_List(Mesh_ptr mesh);
  List_ptr    MESH_GhostEdge_List(Mesh_ptr mesh);
  List_ptr    MESH_GhostFace_List(Mesh_ptr mesh);
  List_ptr  MESH_GhostRegion_List(Mesh_ptr mesh);

  List_ptr  MESH_OverlapVertex_List(Mesh_ptr mesh);
  List_ptr    MESH_OverlapEdge_List(Mesh_ptr mesh);
  List_ptr    MESH_OverlapFace_List(Mesh_ptr mesh);
  List_ptr  MESH_OverlapRegion_List(Mesh_ptr mesh);

  MVertex_ptr  MESH_GhostVertex(Mesh_ptr mesh, int i);
  MEdge_ptr    MESH_GhostEdge(Mesh_ptr mesh, int i);
  MFace_ptr    MESH_GhostFace(Mesh_ptr mesh, int i);
  MRegion_ptr  MESH_GhostRegion(Mesh_ptr mesh, int i);

  MVertex_ptr  MESH_Next_GhostVertex(Mesh_ptr mesh, int *index);
  MEdge_ptr    MESH_Next_GhostEdge(Mesh_ptr mesh, int *index);
  MFace_ptr    MESH_Next_GhostFace(Mesh_ptr mesh, int *index);
  MRegion_ptr  MESH_Next_GhostRegion(Mesh_ptr mesh, int *index);

  void       MESH_Add_GhostVertex(Mesh_ptr mesh, MVertex_ptr v);
  void       MESH_Add_GhostEdge(Mesh_ptr mesh, MEdge_ptr e);
  void       MESH_Add_GhostFace(Mesh_ptr mesh, MFace_ptr f);
  void       MESH_Add_GhostRegion(Mesh_ptr mesh, MRegion_ptr r);

  void       MESH_Rem_GhostVertex(Mesh_ptr mesh, MVertex_ptr v);
  void       MESH_Rem_GhostEdge(Mesh_ptr mesh, MEdge_ptr e);
  void       MESH_Rem_GhostFace(Mesh_ptr mesh, MFace_ptr f);
  void       MESH_Rem_GhostRegion(Mesh_ptr mesh, MRegion_ptr r);

  int        MESH_Num_OverlapVertices(Mesh_ptr mesh);
  int        MESH_Num_OverlapEdges(Mesh_ptr mesh);
  int        MESH_Num_OverlapFaces(Mesh_ptr mesh);
  int        MESH_Num_OverlapRegions(Mesh_ptr mesh);

  MVertex_ptr  MESH_OverlapVertex(Mesh_ptr mesh, int i);
  MEdge_ptr    MESH_OverlapEdge(Mesh_ptr mesh, int i);
  MFace_ptr    MESH_OverlapFace(Mesh_ptr mesh, int i);
  MRegion_ptr  MESH_OverlapRegion(Mesh_ptr mesh, int i);

  MVertex_ptr  MESH_Next_OverlapVertex(Mesh_ptr mesh, int *index);
  MEdge_ptr    MESH_Next_OverlapEdge(Mesh_ptr mesh, int *index);
  MFace_ptr    MESH_Next_OverlapFace(Mesh_ptr mesh, int *index);
  MRegion_ptr  MESH_Next_OverlapRegion(Mesh_ptr mesh, int *index);

  void       MESH_Add_OverlapVertex(Mesh_ptr mesh, MVertex_ptr v);
  void       MESH_Add_OverlapEdge(Mesh_ptr mesh, MEdge_ptr e);
  void       MESH_Add_OverlapFace(Mesh_ptr mesh, MFace_ptr f);
  void       MESH_Add_OverlapRegion(Mesh_ptr mesh, MRegion_ptr r);

#endif

  void       MESH_Set_GModel(Mesh_ptr mesh, GModel_ptr geom);
  int        MESH_Change_RepType(Mesh_ptr mesh, int nurep);
  
  void       MESH_Renumber(Mesh_ptr mesh);

  int MESH_AutoLock(Mesh_ptr mesh);
  void MESH_Set_AutoLock(Mesh_ptr mesh, int autolock);
  
#ifdef __cplusplus
}
#endif


#endif


/*
     par_adj_flags is an integer array of size numprocs
     par_adj_flags[i] on a processor j indicates the adjacency relationship of 
     the mesh on processor j with the mesh on processor i.

     1 and 2 bit from the right indicate relation on vertex
     0(00) no relation,
     1(01) has ghost vertices from processor j,
     2(10) has overlap entities on processor j,
     3(11) both ghosts and overlaps
     3 and 4 bit indicate relation on edge
     5 and 6 bit indicate relation on face
     7 and 8 bit indicate relation on region

     for example: par_adj_flags[3] = 51 (00110011) means mesh on this
     processor has ghost vertices and ghost faces whose master
     processor is 3 and has overlap vertices and overlap faces which
     are ghost entities on processor 3.

     par_recv_info stores some information that tells this processor
     about what sizes of data to expect from an adjacent processor so
     that it can allocate the appropriate sizes

     par_recv_info[0]: number of processors from which this processor
     has ghost entities (called num_recv below)

     par_recv_info[1]-par_recv_info[num_recv+1]: processor ids that
     whose ghosts are on this processor

     Then the subsequent entries in par_recv_info stores the number of
     entities of each type that we exect to receive from another
     processor (4 per processor)

     The numbers of entities that we expect to receive from another
     processor are different from the number of ghost entities from
     that processor because the other processor does not store info
     about the processors corresponding to its overlap entities (it is
     not a one-to-one correspondence). Hence, it sends ALL its overlap
     entities to all processors that it overlaps with even if the
     overlap may be just one or two entities

     Here is an unrealistic but illustrative example

     par_recv_info[0] = 2   --- we have ghosts from two other processors
     par_recv_info[1] = 3   --- we have ghosts from processor 3
     par_recv_info[2] = 99  --- we have ghosts from processor 99
     par_recv_info[3] = 24  --- processor 3 will send us data for 24 vertices
     par_recv_info[4] = 20  --- processor 3 will send us data for 20 edges
     par_recv_info[5] = 12  --- processor 3 will send us data for 12 faces
     par_recv_info[6] =  0  --- processor 3 will send us no data for regions
     par_recv_info[7] = 24  --- processor 99 will send us data for 24 vertices
     par_recv_info[8] =  0  --- processor 99 will send us no data for edges
     par_recv_info[9] =  4  --- processor 99 will send us data for 4 faces
     par_recv_info[10]=  1  --- processor 99 will send us data for 1 region

*/
