#include <stdio.h>
#include <stdlib.h>
#include "MSTK.h"

#ifdef DEBUG

  void List_PrintID(List_ptr l) {
    int i, k, num;
    MEntity_ptr entry;

    num = List_Num_Entries(l);

    fprintf(stderr,"Set 0x%x:\n",(unsigned int) l);
    for (i = 0, k = 0; i < num; i++) {
      entry = List_Entry(l,i);

      fprintf(stderr,"%d ",MEnt_ID(entry));
      k++;
      if (k%10 == 0)
	fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n");
  }

#endif
