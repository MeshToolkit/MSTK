#ifndef _H_Mesh
#define _H_Mesh

#include "MSTK_types.h"
#include "List.h"

#ifdef _cplusplus
extern "C" {
#endif

#ifdef _H_Mesh_Private
  typedef struct Mesh {
    RepType reptype;
    int nv, ne, nf, nr;
    List_ptr mvertex, medge, mface, mregion;
    GModel_ptr geom;
    int max_vid, max_eid, max_fid, max_rid;
  } Mesh, *Mesh_ptr;
#else
  typedef void *Mesh;
#endif



  Mesh_ptr   MESH_New(RepType type);
  int        MESH_InitFromFile(Mesh_ptr mesh, char *filename);
  void       MESH_WriteToFile(Mesh_ptr mesh, char *filename);
  
  GModel_ptr MESH_GModel(Mesh_ptr mesh);
  RepType    MESH_RepType(Mesh_ptr mesh);

  int        MESH_Num_Vertices(Mesh_ptr mesh);
  int        MESH_Num_Edges(Mesh_ptr mesh);
  int        MESH_Num_Faces(Mesh_ptr mesh);
  int        MESH_Num_Regions(Mesh_ptr mesh);
  
  MVertex_ptr  MESH_Vertex(Mesh_ptr mesh, int i);
  MEdge_ptr    MESH_Edge(Mesh_ptr mesh, int i);
  MFace_ptr    MESH_Face(Mesh_ptr mesh, int i);
  MRegion_ptr  MESH_Region(Mesh_ptr mesh, int i);
  
  void       MESH_Add_Vertex(Mesh_ptr mesh, MVertex_ptr v);
  void       MESH_Add_Edge(Mesh_ptr mesh, MEdge_ptr e);
  void       MESH_Add_Face(Mesh_ptr mesh, MFace_ptr f);
  void       MESH_Add_Region(Mesh_ptr mesh, MRegion_ptr r);

  void       MESH_Rem_Vertex(Mesh_ptr mesh, MVertex_ptr v);
  void       MESH_Rem_Edge(Mesh_ptr mesh, MEdge_ptr e);
  void       MESH_Rem_Face(Mesh_ptr mesh, MFace_ptr f);
  void       MESH_Rem_Region(Mesh_ptr mesh, MRegion_ptr r);

  void       MESH_Set_GModel(Mesh_ptr mesh, GModel_ptr geom);
  int        MESH_Change_RepType(Mesh_ptr mesh, int nurep);


#ifdef _cplusplus
}
#endif


#endif
