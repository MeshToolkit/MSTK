#ifndef _H_MSTKUTIL
#define _H_MSTKUTIL

#ifdef __cplusplus
extern "C" {
#endif 

  typedef enum ErrType {MSTK_MESG=0, MSTK_WARN, MSTK_ERROR, MSTK_FATAL} ErrType;
  
  void MSTK_Report(const char *, const char *, ErrType);
  
#ifdef __cplusplus
}
#endif

#endif
