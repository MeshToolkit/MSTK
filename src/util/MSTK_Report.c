/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "MSTK_util.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MSTK_Report(const char *funcname, const char *message, ErrType type) {
    static char last_message[1024];
    int len;

    if (strcmp(last_message,message) == 0)
      return;

    
    len = strlen(message);
    memcpy(last_message,message,len*sizeof(char));
    last_message[len] = '\0';

    switch (type) {
    case MSTK_MESG:
      fprintf(stdout, "\n%s: %s\n",funcname,message);
      break;
    case MSTK_WARN:
      fprintf(stderr, "\nWarning!! in %s: %s\n",funcname,message);
      break;
    case MSTK_ERROR:
      fprintf(stderr, "\nERROR!! in %s: %s\n",funcname,message);
      break;
    case MSTK_FATAL:
      fprintf(stderr, "\nFATAL ERROR!! in %s: %s\n",funcname,message);
      exit(-1);
    }
  }

#ifdef __cplusplus
}
#endif
