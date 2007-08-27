#ifndef _H_MSTKUTIL
#define _H_MSTKUTIL

#ifdef __cplusplus
extern "C" {
#endif 

  typedef enum ErrType {MESG=0, WARN, ERROR, FATAL} ErrType;
  
  void MSTK_Report(const char *, const char *, ErrType);
  
#ifdef __cplusplus
}
#endif

#endif
