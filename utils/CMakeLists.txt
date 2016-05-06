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

if (ENABLE_EXODUSII)
  add_subdirectory(exoatt)
endif (ENABLE_EXODUSII)

#############################################################################
# Partition visualization
###############################################################################

if (ENABLE_PARALLEL)
  add_subdirectory(vizpart)
endif (ENABLE_PARALLEL)
