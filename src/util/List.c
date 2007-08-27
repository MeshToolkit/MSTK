#define _H_List_Private

#include <stdio.h>
#include <stdlib.h>
#include "List.h"

#ifdef __cplusplus
extern "C" {
#endif

  /*  typedef struct List {
    int nentdat;
    int *remdat;
    void *entry;
  } List, *List_ptr;
  */

  /* NOTES: The List structure contains a linear array of pointers to
     the entities in the List. If a entry is added to the List it is
     added after the last entry in the List. 

     The field 'nentdat' keeps track of the number of valid entries
     and the number of allocated entries. At any time the number of
     allocated entries is 2^N-1; So N is stored in the lowest (or
     rightmost) 5 bits and can have a maximum value of 31. The higher 27
     bits are used for storing the actual number of valid entries in
     the list. That means we can have a list with a maximum of 2^27-1
     = 134,217,727 or just above 134 million elements. That should be
     plenty since we will move to a 64 bit architecture before we hit
     that limit

     If an entry is removed, the removed entry is replaced by NULL.
     If any entries are removed, they are kept track of in the field
     remdat, which is a pointer to two integers. Normally, this is
     initialized to NULL and when the first element is removed from
     the list, the integers are allocated. 'remdat[0]' will keep track
     of how many entries are removed and 'remdat[1]' will keep track
     of the location of the first removed entry. We assume that most
     lists will not have entries removed from them and for the few
     that are, we have the remdat
     
  */

  /* nent - number of valid entries 
     p - parameter to say how many space are allocated i.e. 2^p-1
     nrem - number of removed entities
     rem1 - index of first removed entity 
  */

  void pvtList_Get_Pars(List_ptr l, int *nent, int *p, int *nrem, int *rem1){

    *nent = (l->nentdat)>>5;    
    if (l->remdat) {
      *nrem = l->remdat[0];
      *rem1 = l->remdat[1];
    }
    else {
      *nrem = 0;
      *rem1 = -1;
    }
    *p = l->nentdat & 31;  /* 31 is 0000...00011111 i.e. 1 in first 5 spaces */
  }

  int pvtList_Get_Nent(List_ptr l) {
    return (l->nentdat>>5);
  }


  void pvtList_Set_Pars(List_ptr l, int nent, int p, int nrem, int rem1) { 

#ifdef DEBUG
      if (p > 27) 
	MSTK_Report("List_Add","Maximum list length exceeded",FATAL);
#endif

    l->nentdat = ((nent<<5) | p);
    if (nrem) {
      if (!l->remdat)
	l->remdat = (unsigned int *) MSTK_malloc(2*sizeof(void *));
      l->remdat[0] = nrem;
      l->remdat[1] = rem1;
    }
    else {
      if (l->remdat)
	MSTK_free(l->remdat);
      l->remdat = NULL;
    }
  }

  List_ptr List_New(int inisize) {
    int p, nalloc;
    List_ptr newl;

    if (inisize) { /* Find power of 2 greater than inisize */
      p = 0;
      while (1) {
	if ((1<<p) > inisize) break;
	p++;
      }
    }
    else 
      p = 3;

    nalloc = (1<<p)-1; /* This is equivalent to 2^np-1 */

    newl = (List_ptr) MSTK_malloc(sizeof(List));
    newl->entry = (void **) MSTK_malloc(nalloc*sizeof(void *));
    newl->nentdat = 0;
    newl->remdat = NULL;

    pvtList_Set_Pars(newl,0,p,0,-1);

    return newl;
  }
    
  
  void List_Delete(List_ptr l) {
    pvtList_Set_Pars(l,0,0,0,-1);
    free(l->entry);
    free(l);
  }

  List_ptr List_Copy(List_ptr oldl) {
    int nalloc, nent, p, nrem, rem1, ntot;
    List_ptr newl;

    pvtList_Get_Pars(oldl,&nent,&p,&nrem,&rem1);
    nalloc = (1<<p)-1;
    ntot = nent+nrem;

    newl = (List_ptr) MSTK_malloc(sizeof(List));
    newl->entry = (void **) MSTK_malloc(nalloc*sizeof(void *));
    newl->nentdat = 0;
    newl->remdat = NULL;

    memcpy(newl->entry,oldl->entry,ntot*sizeof(void *));
    pvtList_Set_Pars(newl,nent,p,nrem,rem1);

    return newl;
  }


  int List_Num_Entries(List_ptr l) {
    return pvtList_Get_Nent(l);
  }

  /* Return i'th valid entry (not necessarily the i'th element of the
     array if the array has gaps) */
  /* returns NULL pointer if index is greater than number of valid entries */
  /* index goes from 0 to (NumberOfEntries-1)                              */

  void *List_Entry(List_ptr l, int i) {
    void *entry=NULL;
    int ic=0, ip=0, p, ntot, nent, nrem, rem1;

    pvtList_Get_Pars(l,&nent,&p,&nrem,&rem1);

    if (i >= nent)
      return (void *) NULL;

    ntot = nent + nrem;   /* total number of list entries */

    ic = -1; /* (valid entry count)-1 (1st entry may be invalid)*/
    ip = 0;  /* position */

    if (!nrem || rem1 > i) {
      /* The array has no removed entities (gaps) or has no gaps upto
	 the i'th position. We can just return the i'th entry as is */
      ic = i;
      ip = i;
    }
    else {
      /* The array has gaps before position i. We have to step through
	 the set and count entities so as to return the i'th 'valid'
	 entry */

      for (ip = 0; ip < ntot; ip++) {
	entry = l->entry[ip];
    
	if (entry != (void *) NULL) {
	  ic++;
	  if (ic == i)
	    break;
	}
      }
    }
    if (ic != i) { 
      MSTK_Report("List_Entry","Could not find the i'th entry in set",ERROR);
      return 0;
    }

    return l->entry[ip];
  }

  void *List_Next_Entry(List_ptr l, int *index) {
    int ip=0, nent, nrem, p, rem1, ntot=0;
    void *lentry=NULL;
    
    ip = *index;

    pvtList_Get_Pars(l,&nent,&p,&nrem,&rem1);

    ntot = nent + nrem;

    if (ip < 0 || ip >= ntot)
      return (void *) NULL;

    lentry = l->entry[ip];

    if (lentry != (void *) NULL) {
      (*index)++;
      return lentry;
    }
    else {
      for (ip = *index; ip < ntot; ip++) {
	lentry = l->entry[ip];
	if (lentry != (void *) NULL) {
	  *index = ip+1;
	  return lentry;
	}
      }

      return (void *) NULL;
    }
  }

  int List_Contains(List_ptr l, void *entry) {
    int ip, nent, nrem, p, rem1, ntot=0;
    void *lentry;

    pvtList_Get_Pars(l,&nent,&p,&nrem,&rem1);

    ntot = nent + nrem;

    for (ip = 0; ip < ntot; ip++) {
      lentry = l->entry[ip];

      if (lentry == entry)
	return 1;
    }

    return 0;
  }


  int List_Locate(List_ptr l, void *entry) {
    int ip, ic=-1, nent, nrem, ntot, p, rem1;
    void *lentry;
    
    pvtList_Get_Pars(l,&nent,&p,&nrem,&rem1);

    ntot = nent + nrem;

    for (ip = 0; ip < ntot; ip++) {
      lentry = l->entry[ip];

      if (lentry != (void *) NULL)
	ic++;
      if (lentry == entry)
	return ic;
    }

    return -1;
  }
      
      


  List_ptr List_Add(List_ptr l, void *entry) {    
    int nent, nrem, rem1, ntot, p, nalloc;

    pvtList_Get_Pars(l,&nent,&p,&nrem,&rem1);

    ntot = nent + nrem;

    nalloc = (1<<p)-1;

    if (ntot == nalloc) {
      p++;
      nalloc = (1<<p)-1;
      l->entry = (void **) MSTK_realloc(l->entry,nalloc*sizeof(void *));
    }

    l->entry[ntot] = entry;
    nent++;

    pvtList_Set_Pars(l,nent,p,nrem,rem1);

    return l;
  }


  List_ptr List_ChknAdd(List_ptr l, void *entry) {
    int i, nent, nrem, p, ntot, rem1, nalloc;

    pvtList_Get_Pars(l,&nent,&p,&nrem,&rem1);

    ntot = nent + nrem;

    for (i = 0; i < ntot; i++)
      if (l->entry[i] == entry)
	return l;

    /* not found in set - add it */

    nalloc = (1<<p)-1;

    if (ntot == nalloc) {
      p++;
      nalloc = (1<<p)-1;
      l->entry = (void **) MSTK_realloc(l->entry,nalloc*sizeof(void *));
    }

    l->entry[ntot] = entry;
    nent++;    

    pvtList_Set_Pars(l,nent,p,nrem,rem1);

    return l;
  }


  List_ptr List_Inserti(List_ptr l, void *nuentry, int i) {
    void *entry;
    int ic, ip, nent, p, ntot, nrem=0, rem1=-1, nalloc;

    pvtList_Get_Pars(l,&nent,&p,&nrem,&rem1);

    ntot = nent + nrem;

    nalloc = (1<<p)-1;
    
    if (ntot == nalloc) {
      p++;
      nalloc = (1<<p)-1;
      l->entry = (void **) MSTK_realloc(l->entry,nalloc*sizeof(void *));
    }

    if (i == nent) { /* Just append */
      l->entry[ntot] = nuentry;
    }
    else {
      ic = -1; /* (valid entry count)-1 (1st entry may be invalid)*/
      ip = 0;  /* position */
      
      if (!nrem || (rem1 > i)) {
	/* The array has no removed entities (gaps) or has no gaps upto
	   the i'th position. We can just return the i'th entry as is */
	ic = i;
	ip = i;
      }
      else {
	/* The array has gaps before position i. We have to step through
	   the set and count entities so as to return the i'th 'valid'
	   entry */
	
	for (ip = 0; ip < ntot; ip++) {
	  entry = l->entry[ip];
	  
	  if (entry != (void *) NULL) {
	    ic++;
	    if (ic == i)
	      break;
	  }
	}
      }
      if (ic != i) { 
	MSTK_Report("List_Remi","Could not find the i'th entry in set",ERROR);
	return 0;
      }
      
      /* Use efficient block memory movement */
      /* Move memory block starting from the start of the last valid
	 block (or end of deleted block + 1) to the end of the
	 set */
      memmove(&(l->entry[ip+1]),&(l->entry[ip]),(ntot-ip)*sizeof(void *));
      
      l->entry[ip] = nuentry;
    }
    nent++;
    
    pvtList_Set_Pars(l,nent,p,nrem,rem1);

    return l;
  }

  List_ptr List_Insert(List_ptr l, void *nuentry, void *b4entry) {
    void *entry;
    int ip, ntot, nent, p, nrem, rem1, nalloc;

    pvtList_Get_Pars(l,&nent,&p,&nrem,&rem1);

    ntot = nent + nrem;

    nalloc = (1<<p)-1;

    if (ntot == nalloc) {
      p++;
      nalloc = (1<<p)-1;
      l->entry = (void **) MSTK_realloc(l->entry,nalloc*sizeof(void *));
    }

    if (!b4entry) { /* add to the end of the list */
      l->entry[ntot] = nuentry;
    }
    else {
      for (ip = 0; ip < ntot; ip++) {
	entry = l->entry[ip];
	
	if (entry != b4entry)
	  break;
      }

      if (ip == ntot) {
	MSTK_Report("List_Insert","Could not find 'b4entry'",WARN);
	return l;
      }
      
      /* Use efficient block memory movement */
      /* Move memory block starting from the insert point to the end
	 back by 1 */
      memmove(&(l->entry[ip+1]),&(l->entry[ip]),(ntot-ip)*sizeof(void *));
    
      l->entry[ip] = nuentry;
    }
    nent++;

    pvtList_Set_Pars(l,nent,p,nrem,rem1);

    return l;
  }



  int List_Rem(List_ptr l, void *entry) {
    int  ip, ntot, nent, p, nrem, rem1;

    pvtList_Get_Pars(l,&nent,&p,&nrem,&rem1);

    ntot = nent + nrem;

    /* search for this entry */    
    for (ip = 0; ip < ntot; ip++)
      if (l->entry[ip] == entry) {
	l->entry[ip] = (void *) NULL;
	nrem++;
	nent--;
	if (rem1 == -1 || rem1 > ip)
	  rem1 = ip;
	pvtList_Set_Pars(l,nent,p,nrem,rem1);
        return 1;
      }

    /* We will not use List_Compress here since we cannot predict
       whether some of the entries yet to be seen in List iterators
       (List_Entry and List_Next_Entry) will get pushed to before the
       current location. This could be particularly problematic if
       multiple processes are iterating through the same set */
    /* HOWEVER IT WOULD BE GOOD TO FIGURE OUT SOME WAY TO USE THIS SAFELY!!! */
    
    /*
    if (l->nrem > 0.5*(l->nent))
      l = List_Compress(l);
    */

    return 0;
  }



  int List_Replace(List_ptr l, void *entry, void *nuentry) {
    int k, ntot, nent, p, nrem, rem1;

    pvtList_Get_Pars(l,&nent,&p,&nrem,&rem1);

    ntot = nent + nrem;

    /* search for this entry */    
    for (k = 0; k < ntot; k++)
      if (l->entry[k] == entry) {
	l->entry[k] = nuentry;
	return 1;
      }

    return 0;
  }


  /* Remove the i'th valid entry (not necessarily the i'th entry in
     the set which could be a pointer to some other location in the
     set) */

  int List_Remi(List_ptr l, int i) {
    int ic, ip, ntot, nent, p, nrem, rem1;
    void *entry;

    pvtList_Get_Pars(l,&nent,&p,&nrem,&rem1);

    if (i >= nent) {
#ifdef DEBUG
      MSTK_Report("List_Remi","Not that many entries in set",ERROR);
#endif
      return 0;
    }

    ntot = nent + nrem;

    ic = -1; /* (valid entry count)-1 (1st entry may be invalid)*/
    ip = 0;  /* position */

    /* Find the i'th entry */
    if (!nrem || (rem1 > i)) {
      ic = i;
      ip = i;
    }
    else {    
      /* The array has gaps before position i. We have to step through
	 the set and count entities so as to get the i'th 'valid'
	 entry */

      for (ip = 0; ip < ntot; ip++) {
	entry = l->entry[ip];

	if (entry != (void *) NULL) {
	  ic++;
	  if (ic == i)
	    break;
	}
	/* ip++; */ /* I think this is erroneous line */
      }
    }
    if (ic != i) { 
      MSTK_Report("List_Remi","Could not find the i'th entry in set",ERROR);
      return 0;
    }


    l->entry[ip] = (void *) NULL;
    nent--;
    nrem++;
    if (rem1 == -1 || rem1 > ip)
      rem1 = ip;

    pvtList_Set_Pars(l,nent,p,nrem,rem1);

    return 1;
  }


  /* Replace the i'th valid entry (not necessarily the i'th entry in
     the set which could be a pointer to some other location in the
     set) */

  int List_Replacei(List_ptr l, int i, void *nuentry) {
    int ic, ip, ntot, nent, p, nrem, rem1;
    void *entry;

    pvtList_Get_Pars(l,&nent,&p,&nrem,&rem1);

    ntot = nent + nrem;

    ic = -1; /* (valid entry count)-1 (1st entry may be invalid)*/
    ip = 0; /* Position */
    if (!nrem || (rem1 > i)) {
      /* The array has no removed entities (gaps) or has no gaps upto
	 the i'th position. We can just replace the i'th entry as is */
      ic = i;
      ip = i;
    }
    else {
      /* The array has gaps before position i. We have to step through
	 the set and count entities so as to return the i'th 'valid'
	 entry */

      for (ip = 0; ip < ntot; ip++) {
	entry = l->entry[ip];
    
	if (entry != (void *) NULL) {
	  ic++;
	  if (ic == i)
	    break;
	}
      }
    }
    if (ic != i) { 
      MSTK_Report("List_Remi","Could not find the i'th entry in set",ERROR);
      return 0;
    }

    l->entry[ip] = nuentry;
    return 1;
  }


  List_ptr List_Compress(List_ptr l) {
    int i, kb, ke, ntot, ndel, found, nent, p, nrem, rem1;

    pvtList_Get_Pars(l,&nent,&p,&nrem,&rem1);

    if (nrem == 0)
      return l;

    ntot = nent + nrem;

    /* Compress the set starting from the tail */

    kb = ke = ntot;
    while (1) {

      /* Find the end of a block of deleted entries */
      found = 0;
      ke = kb-1;
      while (ke >= 0 && !found) {
	if (l->entry[ke] == (void *) NULL) {
	  /* found an NULL or deleted entry */
	  found = 1;
	}
	else
	  ke--;
      }
      if (!found) {
	/* we are done */
	break;
      }
      
      /* Find the beginning of the block - worst case, it will go all
         the way to the beginning of the array */
      i = ke;
      while (i >= 0) {
	if (l->entry[i] == (void *) NULL) {
	  kb = i;
	  i--;
	}
	else {
	  /* found a valid entry - next entry is start of deleted block */
	  break;
	}
      }

      ndel = ke - kb + 1;

      
      /* Check if the deleted block is at the end of the set. Cut the
	 set off at the beginning of the deleted block */

      if (ke == ntot-1) {
	ntot = ntot - ndel;
	continue;
      }


      /* Use efficient block memory movement */
      /* Move memory block starting from the start of the last valid
	 block (or end of deleted block + 1) to the end of the
	 set */
      memmove(&(l->entry[kb]),&(l->entry[ke+1]),(ntot-(ke+1))*sizeof(void *));
      ntot = ntot - ndel;
    }
       
    nrem = 0;
    rem1 = -1;
    
    /* Should we realloc here so that p is just sufficiently big? */
    pvtList_Set_Pars(l,nent,p,nrem,rem1);

    return l;
  }

  List_ptr List_Cat(List_ptr dest, List_ptr src) {
    int ntot, ntot_d, ntot_s, nent_d, nent_s, nrem_d, nrem_s, rem1_d, rem1_s;
    int p, p_d, p_s, nalloc;

    pvtList_Get_Pars(dest,&nent_d,&p_d,&nrem_d,&rem1_d);
    pvtList_Get_Pars(src,&nent_s,&p_s,&nrem_s,&rem1_s);

    ntot_d = nent_d + nrem_d;
    ntot_s  = nent_s + nrem_s;
    ntot = ntot_d + ntot_s;

    p = (p_d > p_s) ? p_d : p_s;
    while (1) {
      if ((1<<p) > ntot) break;
      p++;
    }

    if (p > p_d) {
      p_d = p;
      nalloc = (1<<p)-1;    
      dest->entry = (void **) MSTK_realloc(dest->entry,nalloc*sizeof(void *));
    }

    memcpy(dest->entry+ntot_d,src->entry,ntot_s*sizeof(void *));
    nent_d += nent_s;
    nrem_d += nrem_s;
    if (!nrem_d && nrem_s)
      rem1_d = rem1_s;

    pvtList_Set_Pars(dest,nent_d,p_d,nrem_d,rem1_d);

    /* List_Compress(dest); */

    return dest;    
  }


#ifdef DEBUG

  void List_Print(List_ptr l) {
    int ip, k, ntot, nent, p, nrem, rem1;
    void *lentry;
    
    pvtList_Get_Pars(l,&nent,&p,&nrem,&rem1);

    ntot = nent + nrem;


    fprintf(stderr,"List 0x%x:\n",(unsigned int) l);
    for (ip = 0, k = 0; ip < ntot; ip++) {
      lentry = l->entry[ip];

      if (lentry != (void *) NULL) {
	fprintf(stderr,"0x%x ",(unsigned int)lentry);
	k++;
	if (k%10 == 0)
	  fprintf(stderr,"\n");
      }
    }
    fprintf(stderr,"\n");
  }

#endif

  
#ifdef __cplusplus
}
#endif
