# Copyright 2019 Triad National Security, LLC. All rights reserved.
# 
# This file is part of the MSTK project. Please see the license file at
# the root of this repository or at
# https://github.com/MeshToolkit/MSTK/blob/master/LICENSE

# -*- mode: cmake -*-
#
# Mesh Utilities
#
# Author: Rao Garimella (rao@lanl.gov)
#



#############################################################################
# Mesh conversion utility
###############################################################################

add_subdirectory(meshconvert)


#############################################################################
# Mesh checking utility
###############################################################################

add_subdirectory(chkmesh)


###############################################################################
# Utility to attach auxilliary attributes to exodus mesh
###############################################################################

if (ENABLE_ExodusII)
  add_subdirectory(exoatt)
endif (ENABLE_ExodusII)

#############################################################################
# Partition visualization
###############################################################################

if (ENABLE_PARALLEL)
  add_subdirectory(vizpart)
endif (ENABLE_PARALLEL)


##############################################################################
# Rectangular mesh creation
##############################################################################

add_subdirectory(mkstruc)

##############################################################################
# Tet 2 Hex mesh conversion
##############################################################################

add_subdirectory(tet2hex)

##############################################################################
# Tri 2 Quad mesh conversion
##############################################################################

add_subdirectory(tri2quad)
