#include <stdio.h>
#include <stdlib.h>
#include "MSTK.h"

#ifdef DEBUG

  void Set_PrintID(Set_ptr l) {
    int i, k, num;
    MEntity_ptr entry;

    num = Set_Num_Entries(l);

    fprintf(stderr,"Set 0x%x:\n",l);
    for (i = 0, k = 0; i < num; i++) {
      entry = Set_Entry(l,i);

      fprintf(stderr,"%d ",MEnt_ID(entry));
      k++;
      if (k%10 == 0)
	fprintf(stderr,"\n");
    }
    fprintf(stderr,"\n");
  }

#endif
