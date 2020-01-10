/* 
Copyright 2019 Triad National Security, LLC. All rights reserved.

This file is part of the MSTK project. Please see the license file at
the root of this repository or at
https://github.com/MeshToolkit/MSTK/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "MSTK.h"
#include "MSTK_util.h"

#ifdef __cplusplus
extern "C" {
#endif

  void MSTK_Version(int *major_version, int *minor_version, int *patch_version,
                  char **version_string) {
    *major_version = MSTK_VERSION_MAJOR;
    *minor_version = MSTK_VERSION_MINOR;
    *patch_version = MSTK_VERSION_PATCH;
    sprintf(*version_string,"%s.%s.%s",*major_version, *minor_version, *patch_version);
  }

#ifdef __cplusplus
}
#endif
