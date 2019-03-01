/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

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

  typedef enum {F1=0, F4, R1, R2, R4, UNKNOWN_REP} RepType;

typedef enum {FDELETED=-1, FUNKNOWN=0, TRI=3, QUAD, POLYGON} MFType;
typedef enum {RDELETED=-1, RUNKNOWN=0, TET, PYRAMID, PRISM, HEX, POLYHED} MRType;

  typedef enum {MDELETED=-1, MVERTEX=0, MEDGE=1, MFACE=2, MREGION=3, MUNKNOWNTYPE=4, MALLTYPE=5, MANYTYPE=6} MType;

/* typedefs needed for attributes */
typedef enum {INT=0, DOUBLE, POINTER, VECTOR, TENSOR} MAttType;

/* typedefs to indicate how to aggregate attribute values from different processors */

typedef enum {ATTOP_UNDEF, ATTOP_MAX, ATTOP_MIN, ATTOP_SUM, ATTOP_AVG} MAttOpType;

/* Parallel status of entity 

   PINTERIOR entity is in the interior in a submesh and is independent of any
   other partition. This will be status for serial runs

   POVERLAP entity may be in the partition interior or boundary and is
   the original of a copy on another partition

   PBOUNDARY indicates that an entity is on the processor boundary
   (THIS IS A TEMPORARY FLAG THAT IS USED INTERNALLY AND NO ENTITY
   WILL RETURN PBOUNDARY OUTSIDE OF THE INITIAL CONSTRUCTION AND
   SHOULD NOT BE USED) However, one can ask if an entity is on a
   partition boundary and entities can be marked/unmarked as being on
   partition boundaries.

   PGHOST indicates ghost entity on the boundary of a partition or
   outside the partition in the halo. These types of entities must be
   updated whenever any field variable changes on the original entity.

 */
typedef enum PType {PINTERIOR=0, POVERLAP=1, PBOUNDARY=2, PGHOST=3} PType;


/* Mesh file formats */
typedef enum {MSTK,GMV,EXODUSII,NEMESISI,CGNS,VTK,STL,AVSUCD,DX,X3D} MshFmt;
  
#ifdef __cplusplus
	   }
#endif

#endif
