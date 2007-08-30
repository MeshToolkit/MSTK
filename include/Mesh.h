#ifndef _H_Mesh
#define _H_Mesh

#include "MSTK_defines.h"
#include "MSTK_types.h"
#include "List.h"
#include "Hash.h"

#ifdef _cplusplus
extern "C" {
#endif

#ifdef _H_Mesh_Private
  typedef struct Mesh {
    RepType reptype;
    int nv, ne, nf, nr;
    List_ptr mvertex, medge, mface, mregion;
    GModel_ptr geom;
    List_ptr AttribList;
    int max_vid, max_eid, max_fid, max_rid;
    Hash_ptr hedge, hface;
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

  void       MESH_Set_GModel(Mesh_ptr mesh, GModel_ptr geom);
  int        MESH_Change_RepType(Mesh_ptr mesh, int nurep);
  
  void       MESH_Renumber(Mesh_ptr mesh);

  int        MESH_Init_ParAtts(Mesh_ptr mesh);
  
#ifdef _cplusplus
}
#endif


#endif
