#define _H_MEntity_Private

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MEntity.h"
#include "MSTK_malloc.h"
#include "MAttrib.h"
#include "MSTK.h"
#include "MSTK_private.h"

/* */

#ifdef __cplusplus
extern "C" {
#endif


  void MEnt_Init_CmnData(MEntity_ptr ent) {
    ent->entdat.mesh = NULL;
    ent->entdat.AttInsList = NULL;

    /* The first bit (from the right) in ent->dim_id contains flag
       indicating if the entity is alive or deleted. Bit 2 indicates
       if the entity is temporary/volatile or permanent (edges in all
       reduced representations are temporary).Bits 3,4,5 combined
       contain the dimension of the entity. If this number is greater
       than or equal to 4, it means the dimension of the entity is
       undefined. The rest of the bits encode the entity ID. Since
       this is an unsigned int (32 bits), the max ID can be 2^(32-5)-1
       = 134217727 (approx. 134.2 million entities of each type) */

    ent->entdat.deleted = 0;
    ent->entdat.temporary = 0;
    ent->entdat.dim = MUNKNOWNTYPE;
    ent->entdat.id = 0;

    ent->entdat.rtype = UNKNOWN_REP;
    ent->entdat.gdim = 4; /* unknown */
    ent->entdat.gid = 0;

     /* Empty marker field */

     ent->entdat.marker = 0;
 
#ifdef MSTK_HAVE_MPI
     ent->entdat.ptype = 0;
     ent->entdat.masterparid = 0;
     ent->entdat.globalid = 0;
#endif
  }

  void MEnt_Free_CmnData(MEntity_ptr ent) {    
    MEnt_Rem_AllAttVals(ent);
  }

  void MEnt_Set_Mesh(MEntity_ptr ent, Mesh_ptr mesh) {
    ent->entdat.mesh = mesh;
  }

  Mesh_ptr MEnt_Mesh(MEntity_ptr ent) {
    return ent->entdat.mesh;
  }


  /* The first bit (from the right) in ent->dim_id contains flag
     indicating if the entity is alive or deleted. Bit 2 indicates if
     the entity is temporary/volatile or permanent (edges in all
     reduced representations are temporary).Bits 3,4,5 combined
     contain the dimension of the entity. If this number is greater
     than or equal to 4, it means the dimension of the entity is
     undefined. The rest of the bits encode the entity ID. Since this
     is an unsigned int (32 bits), the max ID can be 2^(32-5)-1 =
     134217727 (approx. 134.2 million entities of each type) */


  /* MESH ENTITY ID */

  void MEnt_Set_Dim(MEntity_ptr ent, int dim) {
    if (dim == MDELETED)
      ent->entdat.deleted = 1;
    else
      ent->entdat.dim = dim;
  }

  MType MEnt_Dim(MEntity_ptr ent) {
    return (ent->entdat.deleted ? MDELETED : ent->entdat.dim);
  }

  void MEnt_Set_Volatile(MEntity_ptr ent) {
    ent->entdat.temporary = 1;
  }

  int MEnt_IsVolatile(MEntity_ptr ent) {
    return ent->entdat.temporary;
  }

  MType MEnt_OrigDim(MEntity_ptr ent) {
#ifdef DEBUG
    if (!ent->entdat.deleted)
      MSTK_Report("MEnt_OrigDim","This is not a deleted entity",MSTK_WARN);    
#endif

    return ((MType) ent->entdat.dim);
  }	



  /* MESH ENTITY ID */

  void MEnt_Set_ID(MEntity_ptr ent, int id) {
#ifdef DEBUG
    if (id > 134217727 )
      MSTK_Report("MEnt_SetID","ID too large",MSTK_WARN);
#endif

    ent->entdat.id = id;
  }


  int MEnt_ID(MEntity_ptr ent) {
    return (ent->entdat.id);
  }


  /* MESH ENTITY REPRESENTATION TYPE */

  void MEnt_Set_RepType_Data(MEntity_ptr ent, RepType rtype) {
    ent->entdat.rtype = rtype;
  }

  void MEnt_Set_RepType(MEntity_ptr ent, RepType rtype) {
    switch (ent->entdat.dim) {
    case MVERTEX:
      MV_Set_RepType(ent,rtype);
      break;
    case MEDGE:
      ME_Set_RepType(ent,rtype);
      break;
    case MFACE:
      MF_Set_RepType(ent,rtype);
      break;
    case MREGION:
      MR_Set_RepType(ent,rtype);
      break;
    default:
      break;
    }
  }

  RepType MEnt_RepType(MEntity_ptr ent) {
    return ((RepType) ent->entdat.rtype); 
  }



  /* MODEL ENTITY DIMENSION */

  void MEnt_Set_GEntDim(MEntity_ptr ent, int gdim) {
    ent->entdat.gdim = gdim;
  }
    

  int MEnt_GEntDim(MEntity_ptr ent) {
    return ent->entdat.gdim;
  }


  /* MODEL ENTITY ID */

  void MEnt_Set_GEntID(MEntity_ptr ent, int gid) {
    ent->entdat.gid = gid;
  }

  int MEnt_GEntID(MEntity_ptr ent) {
    return ent->entdat.gid;
  }

  /* MODEL ENTITY PTR - We won't store it explicitly */
  
  void MEnt_Set_GEntity(MEntity_ptr ent, GEntity_ptr gent) {
    /* When we have the appropriate functions we should set GEntDim
       and GEntID using this info */
  }

  GEntity_ptr MEnt_GEntity(MEntity_ptr ent) {
    /* When we have the appropriate fucntions we should get GEntity
       from the GEntDim and GEntID info */

    return NULL;
  }



  /* ENTITY DELETION */

  /* Mark an entity as deleted */

  void MEnt_Set_DelFlag(MEntity_ptr ent) {
    ent->entdat.deleted = 1;
  }
    
  /* Mark an entity as undeleted or restored or alive */

  void MEnt_Rem_DelFlag(MEntity_ptr ent) {
    ent->entdat.deleted = 0;
  }
    

  /* Delete an entity */

  void MEnt_Delete(MEntity_ptr ent, int keep) {
    switch (ent->entdat.dim) {
    case MVERTEX:
      MV_Delete(ent,keep);
      break;
    case MEDGE:
      ME_Delete(ent,keep);
      break;
    case MFACE:
      MF_Delete(ent,keep);
      break;
    case MREGION:
      MR_Delete(ent,keep);
      break;
    default:
      break;
    }
  }




  /* ENTITY MARKERS */


  int MEnt_IsMarked(MEntity_ptr ent, int markerID) {
    return (ent->entdat.marker & 1<<(markerID-1));
  }

  void MEnt_Mark(MEntity_ptr ent, int markerID) {
#ifdef DEBUG
    int dim = ent->entdat.dim;
    if (dim < MVERTEX && dim > MREGION)
      MSTK_Report("MEnt_Unmark","Not a valid topological entity",MSTK_ERROR);
#endif

    ent->entdat.marker = ent->entdat.marker | 1<<(markerID-1);
  }

  void MEnt_Unmark(MEntity_ptr ent, int markerID) {
#ifdef DEBUG
    int dim = ent->entdat.dim;
    if (dim < MVERTEX && dim > MREGION)
      MSTK_Report("MEnt_Unmark","Not a valid topological entity",MSTK_ERROR);
#endif

    ent->entdat.marker = ent->entdat.marker & ~(1<<(markerID-1));
  }






  /* ATTRIBUTES */


  /* Set value of attribute */

  void MEnt_Set_AttVal(MEntity_ptr ent, MAttrib_ptr attrib, int ival, 
		       double lval, void *pval) {
    int idx, found;
    MType attentdim, entdim;
    MAttIns_ptr attins;
    List_ptr attinslist;


#ifdef DEBUG
    attentdim = MAttrib_Get_EntDim(attrib);
    entdim = ent->entdat.dim;

    if (attentdim != entdim)
      if (attentdim != MALLTYPE &&
          ((entdim == MDELETED) && (attentdim != MEnt_OrigDim(ent))))
        MSTK_Report("MEnt_Set_AttVal",
                    "Attribute not suitable for this entity type",MSTK_ERROR);
#endif

      attinslist = ent->entdat.AttInsList;
      if (!attinslist)
	ent->entdat.AttInsList = attinslist = List_New(3);
      
      idx = 0; found = 0;
      while ((attins = List_Next_Entry(attinslist,&idx))) {
	if (MAttIns_Attrib(attins) == attrib) {
	  found = 1;
	  break;
	}
      }
      
      if (!found) {
	attins = MAttIns_New(attrib);
	attinslist = List_Add(attinslist,attins);
      }
      
      MAttIns_Set_Value(attins, ival, lval, pval);

  }

  /* Remove an attribute from entity */

  void MEnt_Rem_AttVal(MEntity_ptr ent, MAttrib_ptr attrib) {
    int i, idx, found;
    MType attentdim, entdim;
    MAttIns_ptr attins;
    List_ptr attinslist;

#ifdef DEBUG    
    attentdim = MAttrib_Get_EntDim(attrib);
    entdim = ent->entdat.dim;

    if (attentdim != entdim)
      if (attentdim != MALLTYPE &&
          ((entdim == MDELETED) && (attentdim != MEnt_OrigDim(ent))))
        MSTK_Report("MEnt_Rem_AttVal",
                    "Attribute not suitable for this entity type",MSTK_ERROR);
#endif

      attinslist = ent->entdat.AttInsList;
      if (!attinslist)
	return;
    
      idx = 0; i = 0; found = 0;
      while ((attins = List_Next_Entry(attinslist,&idx))) {
	if (MAttIns_Attrib(attins) == attrib) {
	  found = 1;
	  break;
	}
	else
	  i++;
      }
      
      if (!found)
	return;
      
      List_Remi(attinslist,i);
      MAttIns_Delete(attins);

  }

  void MEnt_Clear_AttVal(MEntity_ptr ent, MAttrib_ptr attrib) {
    int idx, found;
    MType attentdim, entdim;
    MAttIns_ptr attins;
    List_ptr attinslist;

#ifdef DEBUG
    attentdim = MAttrib_Get_EntDim(attrib);
    entdim = ent->entdat.dim;

    if (attentdim != entdim)
      if (attentdim != MALLTYPE &&
          ((entdim == MDELETED) && (attentdim != MEnt_OrigDim(ent))))
        MSTK_Report("MEnt_Clear_AttVal",
                    "Attribute not suitable for this entity type",MSTK_ERROR);
#endif

      attinslist = ent->entdat.AttInsList;
      if (!attinslist) return;
      
      idx = 0; found = 0;
      while ((attins = List_Next_Entry(attinslist,&idx))) {
	if (MAttIns_Attrib(attins) == attrib) {
	  found = 1;
	  break;
	}
      }
      
      if (!found) return;
      
      MAttIns_Set_Value(attins, 0, 0.0, NULL);

  }


  /* Clear value of attribute */

  void MEnt_Rem_AllAttVals(MEntity_ptr ent) {
    int idx;
    MAttIns_ptr attins;
    List_ptr attinslist;
    
    attinslist = ent->entdat.AttInsList;
    if (!attinslist)
      return;
    
    idx = 0;
    while ((attins = List_Next_Entry(attinslist,&idx)))
      MAttIns_Delete(attins);
    List_Delete(attinslist);

    ent->entdat.AttInsList = NULL;
  }


  /* Query the value of the attribute */

  int MEnt_Get_AttVal(MEntity_ptr ent, MAttrib_ptr attrib, int *ival, 
		      double *lval, void **pval) {
    int idx, found;
    MType attentdim, entdim;
    MAttIns_ptr attins;
    List_ptr attinslist;

#ifdef DEBUG    
    attentdim = MAttrib_Get_EntDim(attrib);
    entdim = ent->entdat.dim;
    if (attentdim != entdim)
      if (attentdim != MALLTYPE &&
          ((entdim == MDELETED) && (attentdim != MEnt_OrigDim(ent))))
        MSTK_Report("MEnt_Clear_AttVal",
                    "Attribute not suitable for this entity type",MSTK_ERROR);
#endif

    if (ival) *ival = 0;
    if (lval) *lval = 0;
    if (pval) *pval = NULL;
    
    attinslist = ent->entdat.AttInsList;
    if (!attinslist)
      return 0;
    
    idx = 0; found = 0;
    while ((attins = List_Next_Entry(attinslist,&idx))) {
      if (MAttIns_Attrib(attins) == attrib) {
        found = 1;
        break;
      }
    }
    
    if (!found)
      return 0;
    
    MAttIns_Get_Value(attins, ival, lval, pval);
    
    return 1;
    
  }


  /* Extra functionality for hash-tables */

  MEntity_ptr MEnt_NextInHash(MEntity_ptr ent) {
    switch (ent->entdat.dim) {
    case MVERTEX:
      MSTK_Report("MEnt_NextInHash", "Entity is vertex", MSTK_WARN);
      break;
    case MEDGE:
      return ME_NextInHash(ent);
      break;
    case MFACE:
      return MF_NextInHash(ent);
      break;
    case MREGION:
      MSTK_Report("MEnt_NextInHash", "Entity is region", MSTK_WARN);
      break;
    default:
      MSTK_Report("MEnt_NextInHash", "Unknown entity type", MSTK_WARN);
      break;
    }
    return NULL;
  }

  void MEnt_Set_NextInHash(MEntity_ptr ent, MEntity_ptr next) {
    switch (ent->entdat.dim) {
    case MVERTEX:
      MSTK_Report("MEnt_Set_NextInHash", "Entity is vertex", MSTK_WARN);
      break;
    case MEDGE:
      ME_Set_NextInHash(ent, next);
      break;
    case MFACE:
      MF_Set_NextInHash(ent, next);
      break;
    case MREGION:
      MSTK_Report("MEnt_Set_NextInHash", "Entity is region", MSTK_WARN);
      break;
    default:
      MSTK_Report("MEnt_Set_NextInHash", "Unknown entity type", MSTK_WARN);
      break;
    }
  }
  void MEnt_HashKey(MEntity_ptr ent, unsigned int *pn, void* **pp) {
    switch (ent->entdat.dim) {
    case MVERTEX:
      MSTK_Report("MEnt_HashKey", "Entity is vertex", MSTK_WARN);
      break;
    case MEDGE:
      ME_HashKey(ent, pn, pp);
      break;
    case MFACE:
      MF_HashKey(ent, pn, pp);
      break;
    case MREGION:
      MSTK_Report("MEnt_HashKey", "Entity is region", MSTK_WARN);
      break;
    default:
      MSTK_Report("MEnt_HashKey", "Unknown entity type", MSTK_WARN);
      break;
    }
  }

  int MEnt_IsLocked(MEntity_ptr ent) {
    switch (ent->entdat.dim) {
    case MEDGE:
      return ME_IsLocked(ent);
    case MFACE:
      return MF_IsLocked(ent);
    case MVERTEX:
    case MREGION:
      return 0;
    default:
#ifdef DEBUG
      MSTK_Report("MEnt_HashKey", "Unknown entity type", MSTK_WARN);
#endif
      return 0;
    }
  }


#ifdef MSTK_HAVE_MPI

  PType MEnt_PType(MEntity_ptr ent) {
    return ent->entdat.ptype;
  }

  void  MEnt_Set_PType(MEntity_ptr ent, PType ptype) {
    ent->entdat.ptype = ptype;
  }

  int   MEnt_MasterParID(MEntity_ptr ent) {
    return ent->entdat.masterparid;
  }

  void  MEnt_Set_MasterParID(MEntity_ptr ent, int masterparid) {
    ent->entdat.masterparid = masterparid;
  }

  int   MEnt_GlobalID(MEntity_ptr ent) {
    return ent->entdat.globalid;
  }

  void  MEnt_Set_GlobalID(MEntity_ptr ent, int globalid) {
    ent->entdat.globalid = globalid;
  }

  /*
  int  MEnt_MasterLocalID(MEntity_ptr ent) {
    if (ent->entdat.globalid < 0)
      return -(ent->entdat.globalid);
    else
      return -1;
  }


  void  MEnt_Set_MasterLocalID(MEntity_ptr ent, int localid) {

#ifdef DEBUG
    if (ent->entdat.globalid > 0)
      MSTK_Report("MEnt_Set_MasterLocalID","Global ID already set for entity but overwriting with Local ID on master processor",MSTK_WARN);
#endif

    ent->entdat.globalid = -abs(localid);
  }
  */

#endif /* MSTK_HAVE_MPI */

#ifdef __cplusplus
}
#endif
