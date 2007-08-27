#include <stdio.h>
#include <stdlib.h>


#ifdef __cplusplus
extern "C" {
#endif

void *MSTK_malloc(size_t size) {
  return malloc(size);
}

void *MSTK_calloc(size_t num, size_t size) {
  return calloc(num,size);
}

void *MSTK_realloc(void *ptr, size_t size) {
  return realloc(ptr,size);
}

void MSTK_free(void *ptr) {
  free(ptr);
}

#ifdef __cplusplus
}
#endif
