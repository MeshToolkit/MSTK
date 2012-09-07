#ifndef _H_TYPESMSTK
#define _H_TYPESMSTK

#ifdef MSTK_HAVE_MPI
#include "mpi.h"
typedef MPI_Comm MSTK_Comm;
#else 
/* For serial builds define typedef MSTK_Comm so that codes can send in
   a NULL argument */
typedef void *MSTK_Comm;    
#endif

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
#if !(defined _H_MEntity_Private || defined _H_MVertex_Private ||  \
      defined _H_MEdge_Private || defined _H_MFace_Private || \
      defined _H_MRegion_Private)
  typedef void     *MEntity_ptr;
#endif
#ifndef _H_GModel_Private
  typedef void     *GModel_ptr;
#endif
#ifndef _H_GEntity_Private
  typedef void     *GEntity_ptr;
#endif
#ifndef _H_MAttrib_Private
  typedef void     *MAttrib_ptr;
  typedef void     *MAttIns_ptr;
#endif
#ifndef _H_List_Private
  typedef void     *List_ptr;
#endif
#ifndef _H_MSet_Private
  typedef void     *MSet_ptr;
#endif
#ifndef _H_Hash_Private
  typedef void     *Hash_ptr;
#endif

typedef enum RepType {UNKNOWN_REP=-1, F1=0, F4, R1, R2, R4} RepType;

typedef enum MFType {FDELETED=-1, FUNKNOWN=0, TRI=3, QUAD, POLYGON} MFType;
typedef enum MRType {RDELETED=-1, RUNKNOWN=0, TET, PYRAMID, PRISM, HEX, POLYHED} MRType;

  typedef enum MType {MDELETED=-1, MVERTEX=0, MEDGE=1, MFACE=2, MREGION=3, MUNKNOWNTYPE=4, MALLTYPE=5, MANYTYPE=6} MType;

/* typedefs needed for attributes */
typedef enum MAttType {INT=0, DOUBLE, POINTER, VECTOR, TENSOR} MAttType;


#ifdef MSTK_HAVE_MPI  

/* PINTERIOR is the interior in a submesh, no need to get broadcast 
   POVERLAP also belongs to this submesh, but may be used as ghost entities 
   of another processor
   PBOUNDARY indicates processor boundary
   PGHOST indicates ghost elements;
*/
typedef enum PType {PINTERIOR=0, POVERLAP=1, PBOUNDARY=2, PGHOST=3} PType;

#endif

#ifdef __cplusplus
	   }
#endif

#endif
