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
    ent->entdat.dim_id = 0;
    ent->entdat.rtype_gdim_gid = 0;
    ent->entdat.marker = 0;
    ent->entdat.AttInsList = 0;

#ifdef MSTK_HAVE_MPI
    ent->entdat.ptype_masterparid = (0<<2 | PINTERIOR); /* ?? */
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
    unsigned int data = ent->entdat.dim_id;

    if (dim == MDELETED)
      data = data | 1; /* Set the first bit to 1 */
    else {
      /* zero out bits 3,4,5 by bitwise AND with inverse of 0...011100 (28) */
      data = (data & ~28); 
      data = (data | dim<<2); /* set bits 3,4,5 */
    }

    ent->entdat.dim_id = data;
  }

  MType MEnt_Dim(MEntity_ptr ent) {
    int del_flag=0, dim=0;
    unsigned int data = ent->entdat.dim_id;

    del_flag = data & 1; /* check the first bit */
    if (del_flag)
      return MDELETED;
    
    dim = (data >> 2) & 7; /* check bits 2,3 and 4 */
    if (dim <= 3)
      return dim;
    else
      return MUNKNOWNTYPE;
  }

  void MEnt_Set_Volatile(MEntity_ptr ent) {
    unsigned int data = ent->entdat.dim_id;

    data = data | 2; /* 2 is 0000.....00010, i.e. second bit is 1 */
  }

  int MEnt_IsVolatile(MEntity_ptr ent) {
    unsigned int data = ent->entdat.dim_id;
    
    return (data & 2);
  }

  MType MEnt_OrigDim(MEntity_ptr ent) {
    int del_flag=0, dim=0;
    unsigned int data = ent->entdat.dim_id;

    del_flag = data & 1;
    dim = (data>>2) & 7;

    if (!del_flag)
      MSTK_Report("MEnt_OrigDim","This is not a deleted entity",WARN);    

    if (dim <= 3)
      return dim;
    else 
      return MUNKNOWNTYPE;
  }	



  /* MESH ENTITY ID */

  void MEnt_Set_ID(MEntity_ptr ent, int id) {
    unsigned int data = ent->entdat.dim_id;

    if (id > 134217727 )
      MSTK_Report("MEnt_SetID","ID too large",WARN);

    /* Zero out all but the first 5 bits by bitwise AND with 0..11111 (31) */
    data = (data & 31);
    
    /* Set the ID */
    data = (data | (id<<5));

    ent->entdat.dim_id = data;
  }


  int MEnt_ID(MEntity_ptr ent) {
    unsigned int data = ent->entdat.dim_id;
    return (data>>5);
  }




  /* The first 4 bits (from the right) in ent->reptype_gdim_gid encode
     the representation type (maximum of 14; 15, i.e., 1111 is
     reserved to indicate that it is an unknown representation). The
     next 3 bits from the right contain the dimension of the geometric
     model entity that the mesh entity is classified on. If this
     number is greater than 3 (specifically 7 or binary 111), then the
     dimension of the geometric entity is unknown. The rest of the
     bits are reserved for the ID of the geometric model entity */

  /* MESH ENTITY REPRESENTATION TYPE */

  void MEnt_Set_RepType_Data(MEntity_ptr ent, RepType rtype) {
    unsigned int data = ent->entdat.rtype_gdim_gid;

    /* first zero out the first four bits */
    data = data & ~(15);

    if (rtype == UNKNOWN_REP) 
      data = data | 15; /* Set the first 4 bits to 1 */
    else
      data = data | rtype;

    ent->entdat.rtype_gdim_gid = data;
  }

  void MEnt_Set_RepType(MEntity_ptr ent, RepType rtype) {
    unsigned int dim, data = ent->entdat.dim_id;

    /* Get the dimension of the entity */

    dim = (data>>1) & 7;

    switch (dim) {
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
    int rtype=0;
    unsigned int data = ent->entdat.rtype_gdim_gid;

    rtype = data & 15;
    
    if (rtype == 15)
      return UNKNOWN_REP;
    else
      return rtype;
  }



  /* MODEL ENTITY DIMENSION */

  void MEnt_Set_GEntDim(MEntity_ptr ent, int gdim) {
    unsigned int data = ent->entdat.rtype_gdim_gid;

    /* zero out bits 5,6,7 by doing a bitwise AND with inverse of
       0...01110000 (7<<4) */

    data = data & ~(7<<4);  

    /* now set these bits to the model entity dimension */

    if (gdim <= 3)
      data = data | (gdim<<4);
    else
      data = data | (7<<4);  /* set all 3 bits to 1 */

    ent->entdat.rtype_gdim_gid = data;
  }
    

  int MEnt_GEntDim(MEntity_ptr ent) {
    int gdim=0;
    unsigned int data = ent->entdat.rtype_gdim_gid;

    gdim = (data>>4) & 7;

    if (gdim <= 3)
      return gdim;
    else
      return 4; /* Code for unknown type */
  }


  /* MODEL ENTITY ID */

  void MEnt_Set_GEntID(MEntity_ptr ent, int gid) {
    unsigned int data = ent->entdat.rtype_gdim_gid;

    /* first zero out all bits higher than the 7th bit by performing a
       bitwise AND with 127 (0...01111111) */
    data = data & 127;

    if (gid == -1)
      data = data | ~(127); /* all bits higher than 7 are 1 */
    else
      data = data | (gid<<7);
    
    ent->entdat.rtype_gdim_gid = data;
  }

  int MEnt_GEntID(MEntity_ptr ent) {
    unsigned int data = ent->entdat.rtype_gdim_gid;
    return (data>>7);
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

    /* Set first bit to 1 to indicate that this is a deleted entity */    

    ent->entdat.dim_id = ent->entdat.dim_id | 1;
  }
    
  /* Mark an entity as undeleted or restored or alive */

  void MEnt_Rem_DelFlag(MEntity_ptr ent) {

    /* Set first bit to 0 to indicate that this is a valid entity */    

    ent->entdat.dim_id = ent->entdat.dim_id & 0;
  }
    

  /* Delete an entity */

  void MEnt_Delete(MEntity_ptr ent, int keep) {
    switch ( MEnt_Dim(ent) ) {
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
    int dim = MEnt_Dim(ent);
    if (dim < MVERTEX && dim > MREGION)
      MSTK_Report("MEnt_Unmark","Not a valid topological entity",ERROR);
#endif

    ent->entdat.marker = ent->entdat.marker | 1<<(markerID-1);
  }

  void MEnt_Unmark(MEntity_ptr ent, int markerID) {
#ifdef DEBUG
    int dim = MEnt_Dim(ent);
    if (dim < MVERTEX && dim > MREGION)
      MSTK_Report("MEnt_Unmark","Not a valid topological entity",ERROR);
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

    attentdim = MAttrib_Get_EntDim(attrib);
    entdim = MEnt_Dim(ent);

    if ((attentdim == MALLTYPE) || (attentdim == entdim) || 
	(entdim == MDELETED && attentdim == MEnt_OrigDim(ent))) {

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
    else
      MSTK_Report("MEnt_Set_AttVal",
		  "Attribute not suitable for this entity type",ERROR);
  }

  /* Clear value of attribute */

  void MEnt_Rem_AttVal(MEntity_ptr ent, MAttrib_ptr attrib) {
    int i, idx, found;
    MType attentdim, entdim;
    MAttIns_ptr attins;
    List_ptr attinslist;
    
    attentdim = MAttrib_Get_EntDim(attrib);
    entdim = MEnt_Dim(ent);

    if ((attentdim == MALLTYPE) || (attentdim == entdim) || 
	(entdim == MDELETED && attentdim == MEnt_OrigDim(ent))) {

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
    else
      return;
 
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
    
    if (ival) *ival = 0;
    if (lval) *lval = 0;
    if (pval) *pval = NULL;
    
    attentdim = MAttrib_Get_EntDim(attrib);
    entdim = MEnt_Dim(ent);

    if ((attentdim == MALLTYPE) || (attentdim == entdim) || 
	(entdim == MDELETED && attentdim == MEnt_OrigDim(ent))) {

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
    else
      return 0;

  }


  /********************************************************************/
  /* Some functions to emulate parallel functionality for now. We will
     do this right later on */
  /********************************************************************/

  /* TIME TO GET RID OF THIS???? */



  /* Get the number of processors connected to entity */

  int MEnt_NumProcs(MEntity_ptr ent) {
    int ival;
    double rval;
    void *pval;
    int *plnumdat;

    MEnt_Get_AttVal(ent,plnumatt,&ival,&rval,&pval);
    plnumdat = (void *) pval;

    if (!plnumdat)
      return 0;
    else
      return plnumdat[0];

  }
    



  /* Add processor to list of processors connected to entity (if not
     already there) */


  void MEnt_Set_ProcIDs(MEntity_ptr ent, int np, int *procids) {
    int ival, i;
    double rval;
    void *pval;
    int *plnumdat;

    MEnt_Get_AttVal(ent,plnumatt,&ival,&rval,&pval);
    plnumdat = (int *) pval;

    if (!plnumdat) {
      plnumdat = (int *) malloc((1+2*np)*sizeof(int));

      MEnt_Set_AttVal(ent,plnumatt,0,0.0,(void *)plnumdat);      
    }

    /* Number of processors connected to this entity */

    plnumdat[0] = np;

    /* Initialize find this processor ID */

    for (i = 0; i < np; i++) {
      plnumdat[1+2*i] = procids[i];
      plnumdat[1+2*i+1] = 0;
    }
  }




  /* IDs of processors connected to this entity (could be just 1) */


  int MEnt_ProcIDs(MEntity_ptr ent, int *np, int *procids) {
    int ival, i;
    double rval;
    void *pval;
    int *plnumdat;

    MEnt_Get_AttVal(ent,plnumatt,&ival,&rval,&pval);
    plnumdat = (int *) pval;

    if (!plnumdat) {
      MSTK_Report("MEnt_LocalID","No processor specific data. Call MEnt_Set_Num{Procs first",WARN);
      return 0;
    }

    /* Number of processors connected to this entity */

    *np = plnumdat[0];

    for (i = 0; i < *np; i++)
      procids[i] = plnumdat[1+2*i];

    return 1;
  }




  /* Set the local number of an entity on a processor. (Add processor
     to list of processors connected to entity if it is not already
     there */


  void MEnt_Set_LocalID(MEntity_ptr ent, int procid, int lnum) {
    int ival, np, i;
    double rval;
    void *pval;
    int *plnumdat;

    MEnt_Get_AttVal(ent,plnumatt,&ival,&rval,&pval);
    plnumdat = (int *) pval;

    if (!plnumdat) {
      MSTK_Report("MEnt_Set_LocalID",
		  "No processor specific data. Call MEnt_Set_ProcIDs first",
		  WARN);
      return;
    }

    /* Number of processors connected to this entity */

    np = plnumdat[0];

    /* Search to find this processor ID */

    for (i = 0; i < np; i++) {

      if (plnumdat[1+2*i] == procid) {
	plnumdat[1+2*i+1] = lnum;
	return;
      }      
      
    }

    /* This processor ID is not found */

    MSTK_Report("MEnt_Set_LocalID","Processor not connected to entity",WARN);
  }




  /* Local number of the entity on this processor */


  int MEnt_LocalID(MEntity_ptr ent, int procid) {
    int ival, i, np, lnum;
    double rval;
    void *pval;
    int *plnumdat;

    MEnt_Get_AttVal(ent,plnumatt,&ival,&rval,&pval);
    plnumdat = (int *) pval;

    if (!plnumdat) {
      MSTK_Report("MEnt_LocalID",
		  "No processor specific data. Call MEnt_Set_ProcIDs first",
		  WARN);
      return MEnt_ID(ent);
    }

    /* Number of processors connected to this entity */

    np = plnumdat[0];

    /* Search to find this processor ID */

    for (i = 0; i < np; i++) {

      if (plnumdat[1+2*i] == procid) {
	lnum = plnumdat[1+2*i+1];
	return lnum;
      }
      
    }

    /* Processor ID was not found */

    MSTK_Report("MEnt_LocalID","Processor not connected to entity",WARN);
    return -1;
  }



  /* Extra functionality for hash-tables */

  MEntity_ptr MEnt_NextInHash(MEntity_ptr ent) {
    switch ( MEnt_Dim(ent) ) {
    case MVERTEX:
      MSTK_Report("MEnt_NextInHash", "Entity is vertex", WARN);
      break;
    case MEDGE:
      return ME_NextInHash(ent);
      break;
    case MFACE:
      return MF_NextInHash(ent);
      break;
    case MREGION:
      MSTK_Report("MEnt_NextInHash", "Entity is region", WARN);
      break;
    default:
      MSTK_Report("MEnt_NextInHash", "Unknown entity type", WARN);
      break;
    }
    return NULL;
  }

  void MEnt_Set_NextInHash(MEntity_ptr ent, MEntity_ptr next) {
    switch ( MEnt_Dim(ent) ) {
    case MVERTEX:
      MSTK_Report("MEnt_Set_NextInHash", "Entity is vertex", WARN);
      break;
    case MEDGE:
      ME_Set_NextInHash(ent, next);
      break;
    case MFACE:
      MF_Set_NextInHash(ent, next);
      break;
    case MREGION:
      MSTK_Report("MEnt_Set_NextInHash", "Entity is region", WARN);
      break;
    default:
      MSTK_Report("MEnt_Set_NextInHash", "Unknown entity type", WARN);
      break;
    }
  }
  void MEnt_HashKey(MEntity_ptr ent, unsigned int *pn, void* **pp) {
    switch ( MEnt_Dim(ent) ) {
    case MVERTEX:
      MSTK_Report("MEnt_HashKey", "Entity is vertex", WARN);
      break;
    case MEDGE:
      ME_HashKey(ent, pn, pp);
      break;
    case MFACE:
      MF_HashKey(ent, pn, pp);
      break;
    case MREGION:
      MSTK_Report("MEnt_HashKey", "Entity is region", WARN);
      break;
    default:
      MSTK_Report("MEnt_HashKey", "Unknown entity type", WARN);
      break;
    }
  }

  int MEnt_IsLocked(MEntity_ptr ent) {
    switch ( MEnt_Dim(ent) ) {
    case MEDGE:
      return ME_IsLocked(ent);
    case MFACE:
      return MF_IsLocked(ent);
    case MVERTEX:
    case MREGION:
      return 0;
    default:
#ifdef DEBUG
      MSTK_Report("MEnt_HashKey", "Unknown entity type", WARN);
#endif
      return 0;
    }
  }


#ifdef MSTK_HAVE_MPI

  PType MEnt_PType(MEntity_ptr ent) {
    unsigned int data = ent->entdat.ptype_masterparid;
    return data & 3 /*00...00011*/;
  }

  void  MEnt_Set_PType(MEntity_ptr ent, PType ptype) {
    unsigned int data = ent->entdat.ptype_masterparid;
    data = (data & ~3); /* zero out last 2 digits*/
    data = (data | ptype);
    ent->entdat.ptype_masterparid = data;
  }

  int   MEnt_MasterParID(MEntity_ptr ent) {
    unsigned int data = ent->entdat.ptype_masterparid;
    return data>>2;
  }

  void  MEnt_Set_MasterParID(MEntity_ptr ent, int masterparid) {
    unsigned int data = ent->entdat.ptype_masterparid;
    data = (data & 3); /* zero out high bits*/
    data = (data | masterparid<<2);
    ent->entdat.ptype_masterparid = data;
  }

  int   MEnt_GlobalID(MEntity_ptr ent) {
    return ent->entdat.globalid;
  }

  void  MEnt_Set_GlobalID(MEntity_ptr ent, int globalid) {
    ent->entdat.globalid = globalid;
  }
#endif /* MSTK_HAVE_MPI */


#ifdef __cplusplus
}
#endif
