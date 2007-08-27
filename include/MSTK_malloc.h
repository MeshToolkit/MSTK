#ifndef _H_MSTK_ALLOC
#define _H_MSTK_ALLOC

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

void *MSTK_malloc(size_t size);
void *MSTK_calloc(size_t num, size_t size);
void *MSTK_realloc(void *ptr, size_t size);
void MSTK_free(void *ptr);

#ifdef __cplusplus
}
#endif

#endif
