/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

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
