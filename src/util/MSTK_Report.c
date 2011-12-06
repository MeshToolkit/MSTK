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
      fprintf(stdout, "%s: %s\n",funcname,message);
      break;
    case MSTK_WARN:
      fprintf(stderr, "Warning!! in %s: %s\n",funcname,message);
      break;
    case MSTK_ERROR:
      fprintf(stderr, "MSTK_ERROR!! in %s: %s\n",funcname,message);
      break;
    case MSTK_FATAL:
      fprintf(stderr, "MSTK_FATAL MSTK_ERROR!! in %s: %s\n",funcname,message);
      exit(-1);
    }
  }

#ifdef __cplusplus
}
#endif
