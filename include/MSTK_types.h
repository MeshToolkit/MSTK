#ifndef _H_TYPESMSTK
#define _H_TYPESMSTK

#define MSTK_UNKNOWN -1

#ifdef   __cplusplus
extern "C" {
#endif
#ifndef _H_Mesh_Private
  typedef void     *Mesh_ptr;
#endif
#ifndef _H_MVertex_Private
  typedef void     *MVertex_ptr;
#endif
#ifndef _H_MEdge_Private
  typedef void     *MEdge_ptr;
#endif
#ifndef _H_MFace_Private
  typedef void     *MFace_ptr;
#endif
#ifndef _H_MRegion_Private
  typedef void     *MRegion_ptr;
#endif
#ifndef _H_MEntity_Private
  typedef void     *MEntity_ptr;
#endif
#ifndef _H_GModel_Private
  typedef void     *GModel_ptr;
#endif
#ifndef _H_GEntity_Private
  typedef void     *GEntity_ptr;
#endif

#define MSTK_VER    1.0

#define MSTK_MAXREP 10
typedef enum RepType {UNKNOWN_REP=-1, F1=0, F2, F3, F4, F5, F6, R1, R2, R3, R4} RepType;


typedef enum MFType {FDELETED=-1, FUNKNOWN=0, TRI=3, QUAD, POLYGON} MFType;
typedef enum MRType {RDELETED=-1, RUNKNOWN=0, TET, PYRAMID, PRISM, HEX, POLYHED} MRType;

typedef enum MType {MVERTEX=0, MEDGE, MFACE, MREGION} MType;
typedef enum VType {VIGNORE=21, VPARENT=41, VDELETED=61} VType;


/* typedefs needed for attributes */
typedef enum AttType {INT=0, VINT, VINT2, DOUBLE, VDOUBLE, VDOUBLE2, CHAR, VCHAR} AttType;
typedef enum Rank {SCALAR=4001, VECTOR, TENSOR} Rank;
typedef enum Interp {CONSTANT=5001, COPY, SEQUENCE, LINEAR, LOG, ASINH, MAX,
		     MIN, USER, AND, OR, INCMAX} Interp;

/********
typedef enum ElType {DELETED=-1, UNKNOWN=0, TRI=3, QUAD, TET, PYRAMID, PRISM, HEX, POLYGON, POLYHED} ElType;
*********/


#ifdef __cplusplus
	   }
#endif

#endif
