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
    case MESG:
      fprintf(stdout, "\n%s: %s\n",funcname,message);
      break;
    case WARN:
      fprintf(stderr, "\nWarning!! in %s: %s\n",funcname,message);
      break;
    case ERROR:
      fprintf(stderr, "\nERROR!! in %s: %s\n",funcname,message);
      break;
    case FATAL:
      fprintf(stderr, "\nFATAL ERROR!! in %s: %s\n",funcname,message);
      exit(-1);
    }
  }

#ifdef __cplusplus
}
#endif
