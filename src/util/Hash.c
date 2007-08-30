#define _H_Hash_Private

#include <stdio.h>
#include <stdlib.h>
#include "Hash.h"

#ifdef __cplusplus
extern "C" {
#endif


	/* External functions from MSTK */
	void *MEnt_NextInHash(void *ent);
	void  MEnt_Set_NextInHash(void *ent, void *next);
	void  MEnt_HashKey(void *ent, unsigned int *pn, void* **pp);

	/*  typedef struct Hash {
	    unsigned int nentdat;
	    void **entry;
	    } Hash, *Hash_ptr;
	    */

	/* The field 'nentdat' keeps track of the number of valid entries
	   and the number of allocated entries. At any time the number of
	   allocated entries is 2^N-1; So N is stored in the lowest (or
	   rightmost) 5 bits and can have a maximum value of 31. The higher 27
	   bits are used for storing the actual number of valid entries in
	   the list. That means we can have a list with a maximum of 2^27-1
	   = 134,217,727 or just above 134 million elements. That should be
	   plenty since we will move to a 64 bit architecture before we hit
	   that limit

*/

	void pvtHash_Get_Pars(Hash_ptr h, unsigned int *nent, unsigned int *p, int *t){
		int d = h->nentdat;

		*p = d & 0x1f;
		d >>= 5;
		*t = d & 0x01;
		d >>= 1;
		*nent = d;
	}

	unsigned int pvtHash_Get_Nent(Hash_ptr h) {
		return h->nentdat >> 6;
	}

	unsigned int pvtHash_Get_Size(Hash_ptr h) {
		return 1 << (h->nentdat & 0x1f);
	}

	int pvtHash_Get_Type(Hash_ptr h) {
		return (h->nentdat >> 5) & 0x01;
	}


	void pvtHash_Set_Pars(Hash_ptr h, unsigned int nent, unsigned int p, int t) { 

#ifdef DEBUG
		if (p > 26) 
			MSTK_Report("Hash_Add","Maximum hash size exceeded",FATAL);
#endif

		h->nentdat = ((nent << 6) | (t << 5) | p);
	}


	Hash_ptr Hash_New(unsigned int inisize, int type) {
		unsigned int p, nalloc, i;
		Hash_ptr newh;

		if (inisize) { /* Find power of 2 greater than inisize */
			for (p=0; 1<<p < inisize; p++);
		}
		else {
			p = 10;
		}

		nalloc = 1<<p;

		newh = (Hash_ptr) MSTK_malloc(sizeof(Hash));
		newh->entry = (void **) MSTK_malloc(nalloc*sizeof(void *));
		pvtHash_Set_Pars(newh, 0, p, type);
		for (i=0; i<nalloc; i++) {
			newh->entry[i] = NULL;
		}

		return newh;
	}


	void Hash_Delete(Hash_ptr h) {
		free(h->entry);
		free(h);
	}

	int Hash_Num_Entries(Hash_ptr h) {
		return pvtHash_Get_Nent(h);
	}

	unsigned int pvtHash_Function(unsigned int np, void* *p) { /* very simple */
		unsigned int i, j;
		unsigned int hash=0;
		unsigned long ptr;

		for (i=0; i<np; i++) {
			ptr = (unsigned long) p[i];
			for (j=0; j<sizeof(unsigned long); j++) {
				hash += ptr & 0xff;
				hash += (hash << 10);
				hash ^= (hash >> 6);
				ptr >>= 8;
			}
		}
		hash += (hash << 3);
		hash ^= (hash >> 11);
		hash += (hash << 15);
		return hash;
	}

	int pvtHash_CheckKeys(unsigned int np1, void* *p1, unsigned int np2, void* *p2) {
		unsigned int i;
		if (np1 != np2) return 0;
		for (i=0; i<np1; i++) {
			if (p1[i] != p2[i]) return 0;
		}
		return 1;
	}


	void *Hash_Entry(Hash_ptr h, unsigned int np, void* *p) {
		unsigned int enp;
		void* *ep;
		unsigned int hash, mask, i;
		void *ent;

		hash = pvtHash_Function(np, p);
		mask = pvtHash_Get_Size(h) - 1;

		i = hash & mask;
		ent = h->entry[i];
		while (ent) {
			MEnt_HashKey(ent, &enp, &ep);
			if (pvtHash_CheckKeys(np, p, enp, ep)) return ent;
			ent = MEnt_NextInHash(ent);
		}
		return NULL;
	}

	unsigned int pvtHash_Enlarge(Hash_ptr h) {
		unsigned int nent, pwr, nalloc, i, hash, nold, mask;
		int t;
		void *ent, *prev, *next;
		unsigned int enp;
		void* *ep;

		pvtHash_Get_Pars(h, &nent, &pwr, &t);

		nold = 1 << pwr;
		pwr++;
		nalloc = 1 << pwr;
		mask = nalloc - 1;
		h->entry = (void **) MSTK_realloc(h->entry,nalloc*sizeof(void *));

/*		printf("New hash size = %d (%d)\n", nalloc, pwr); */

		for (i=0; i<nold; i++) {
			ent = h->entry[i];
			prev = NULL;
			h->entry[i+nold] = NULL;
			while (ent) {
				MEnt_HashKey(ent, &enp, &ep);
				hash = pvtHash_Function(enp, ep);
				if (i != (hash & mask)) {
					if (i + nold != (hash & mask)) {
						MSTK_Report("Hash_Enlarge","New hash differs",FATAL);
					}
					next = MEnt_NextInHash(ent);
					MEnt_Set_NextInHash(ent, h->entry[i+nold]);
					h->entry[i+nold] = ent;
					if (prev) {
						MEnt_Set_NextInHash(prev, next);
					}
					else {
						h->entry[i] = next;
					}
					ent = next;
				}
				else {
					prev = ent;
					ent = MEnt_NextInHash(ent);
				}
			}
		}

		pvtHash_Set_Pars(h, nent, pwr, t);
		return nalloc;
	}

	Hash_ptr Hash_Add(Hash_ptr h, void *entry, unsigned int np, void* *p) {    
		unsigned int nent, pwr, nalloc, i, hash, mask;
		int t;
		void *ent;

		hash = pvtHash_Function(np, p);

		pvtHash_Get_Pars(h, &nent, &pwr, &t);

		nalloc = 1 << pwr;


		mask = nalloc - 1;

		i = hash & mask;
		ent = h->entry[i];
		MEnt_Set_NextInHash(entry, ent);
		h->entry[i] = entry;
		nent++;

		pvtHash_Set_Pars(h, nent, pwr, t);

		if (1.0*nent >= 2.0*nalloc) {
			pvtHash_Enlarge(h);
		}

		return h;
	}


	Hash_ptr Hash_ChknAdd(Hash_ptr h, void *entry, unsigned int np, void* *p) {
		unsigned int enp;
		void* *ep;
		unsigned int nent, pwr, nalloc, i, hash, mask;
		int t;
		void *ent;

		hash = pvtHash_Function(np, p);

		pvtHash_Get_Pars(h, &nent, &pwr, &t);

		nalloc = (1<<pwr);
		mask = nalloc - 1;

		i = hash & mask;
		ent = h->entry[i];
		while (ent) {
			MEnt_HashKey(ent, &enp, &ep);
			if (pvtHash_CheckKeys(np, p, enp, ep)) return ent;
			ent = MEnt_NextInHash(ent);
		}

		MEnt_Set_NextInHash(entry, h->entry[i]);
		h->entry[i] = entry;
		nent++;

		pvtHash_Set_Pars(h, nent, pwr, t);

		if (1.0*nent >= 2.0*nalloc) {
			pvtHash_Enlarge(h);
		}

		return h;
	}


	int Hash_Rem(Hash_ptr h, void *entry, unsigned int np, void* *p) {
		(void) h, (void) entry, (void) np, (void) p;
		MSTK_Report("Hash_Rem","Not implemented",FATAL);

		return 0;
	}

	List_ptr Hash_Entries(Hash_ptr h) {
		unsigned int ne, i, p, size;
		int t;
		void *ent;
		List_ptr l;

		pvtHash_Get_Pars(h, &ne, &p, &t);
		size = 1 << p;

		l = List_New(ne);
		for (i=0; i<size; i++) {
			ent = h->entry[i];
			while (ent) {
				List_Add(l, ent);
				ent = MEnt_NextInHash(ent);
			}
		}

		return l;
	}

#ifdef DEBUG

	void Hash_Print(Hash_ptr h) {
		unsigned int i, nent, p, size, hits, maxhits=0, zerohits=0, nzhits=0, tothits=0, ghits=0;;
		int type;
		void *ent;

		pvtHash_Get_Pars(h, &nent, &p, &type);

		size = 1 << p;

		fprintf(stderr,"Hash 0x%lx:\n", (unsigned long) h);
		for (i = 0; i < size; i++) {
			ent = h->entry[i];
			hits=0;
			while (ent) {
				hits++;
				ent = MEnt_NextInHash(ent);
			}
			if (hits) {
				nzhits++;
				tothits += hits;
				ghits += hits*hits;
				if (maxhits < hits) maxhits = hits;
			}
			else {
				zerohits++;
			}
		}
		fprintf(stderr, "nent=%d\n", nent);
		if (nent != tothits) fprintf(stderr, "Total hits: %d\n", tothits);
		fprintf(stderr, "Max hits: %d\n", maxhits);
		fprintf(stderr, "Avg hits:  %.2lf\n", ((double)tothits)/nzhits);
		fprintf(stderr, "GAvg hits: %.2lf\n", ((double)ghits)/tothits);
		fprintf(stderr, "   Empty / NonEmpty /    Total\n");
		fprintf(stderr, "%8d   %8d   %8d\n", zerohits, nzhits, size);

		fprintf(stderr,"\n");
	}

#endif


#ifdef __cplusplus
}
#endif

