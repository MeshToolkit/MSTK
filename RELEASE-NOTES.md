04/06/2016
MSTK v2.26rc1

* Fixed a "bug" in the FLAG X3D writer. Apparently the matnames field in the file must only be string version of integers. Misunderstanding by some legacy code a long time ago and the tradition has stuck.
* Minor fixup in MESH_Update1Attribute
'chkmesh' now reads meshes in any of the importable formats

02/07/2015
MSTK v2.25

* Fixed up the parallel import/export of FLAG X3D files
* Incorporated utility to generate structured meshes (mkstruc)
* Added functionality to exchange attributes on entities that have 1-1 mapping on parallel boundaries (faces in 3D, edges in 2D)
* Added functionality to do gather-scatter operation on attributes on entities that have many-to-one mapping on parallel boundaries (vertices in 2D, vertices and edges in 3D)
* Added new operator to query if entity is on parallel boundary
* Removed use of MSTK_malloc, MSTK_calloc, MSTK_realloc and MSTK_free
* Added call that will invert knowledge about which processors will send a processor data to figure out what processors it needs to send data to
* Use CMAKE_INSTALL_PREFIX as installlation directory if INSTALL_DIR is not defined
* Other bug fixes

11/23/2015
MSTK 2.25 rc1

* Added capability to write out vector fields (component by component) to Exodus II meshes
* Added capability to read scalar and vector fields (component by component and aggregated) from Exodus II meshes

10/23/2015
MSTK 2.24

* Fixed some minor bugs in interpretation of natt argument in MESH_ExportToGMV
* Fixed tests to specify the Exodus II format as exodusii instead of exo"
* Bug fix in MFs_Join
* Added unit test for getting the dihedral angle between two mesh faces

6/26/2015
MSTK 2.24rc1

* Added routine MRs_Join to joint two mesh regions along a common face
* MESH_ImportToFile and MESH_ExportToFile enhanced to figure out format from file extension if format is given as NULL
* MESH_ExportToGMV will not write out special MSTK attributes vidatt, ridatt etc.

5/8/2015
MSTK 2.23

* MAJOR: Enhanced the exoatt utility to read text file with element and side set info and introduce them into the mesh
* MINOR: Some tweaks to Exodus II export routine including an option to not renumber vertices to ensure contiguous IDs

5/5/2015
MSTK v2.23rc2

* MAJOR: Added capability to read/write elements sets from/to Exodus II files
* MINOR: Eliminated hash table while reading Exodus II files, freed some memory during mesh partitioning, not checking topological consistency automatically during meshconvert

3/30/2015
MSTK v2.23rc1

* Fixed an input parsing error and added option to recognize .par for output extensions

1/13/2015
MSTK v2.22

* Install utilities such as meshconvert in INSTALL_DIR/bin if specified or in the bin directory under the project source
* Restore MSTK_UpdateAttr call temporarily for backward compatibility but with a warning in DEBUG mode

12/8/2014
MSTK v2.22 rc4

* Major changes to mesh partitioning code to interleave data collection with MPI communication for greater efficiency
* On some machines (mustang at LANL), we see dramatic decreases in partitioning time (4 times faster) but on others (edison.nersc.gov) it is showing a slight increase. This may be due to differences in MPI latency.

11/26/2014
MSTK v2.22 rc3

* Fixed a bug in the restructured code for MESH_ExportToExodusII

11/25/2014
MSTK v2.22 rc2

* Fixed a SERIOUS bug in parallel partitioning that showed up only on machines with high MPI latency
* More specifically, fixed a serious flaw in MESH_SendMesh where I was using a statically declared array inside a routine to do a buffered send but was waiting for all buffered sends to complete outside the routine.
**IF YOU CARE ABOUT PARALLEL PARTITIONING, DO NOT USE ANY VERSION PRIOR TO THIS - THEY MAY WORK, BUT IT IS NOT GUARANTEED**

11/21/2014
MSTK v2.22 rc1

* Added a new utility 'exoatt' to import an Exodus II mesh, augment it with attributes read from a file and export it back out, possibly partitioning it in the process.
* As part of the above, added ability to export real valued attributes as properties to Exodus II files
* Fixed some serious bugs lurking in MESH_CopyAttr (I don't think it ever worked before)

09/08/2014
MSTK v2.21

* Improve efficiency of MESH_Delete
* Improve efficiency of MESH_CopySet in parallel distribution
* Both result in significant improvement in performance of parallel distribution
* Fixed bug in creation of element sets for surface meshes imported from Exodus II files

09/05/2014
MSTK v 2.21 rc4

* A minor change to MESH_RecvMesh to used non-blocking receives to the extent possible
* Small improvement in performance

09/03/2014
MSTK v 2.21 rc3

* Another enhancement in mesh Send/Recv? routines to reduce memory usage at minimal computational cost. Can reduce memory usage on large 3D meshes by 10-15%
* Some minor fixes to how an empty communicator is checked

09/02/2014
MSTK v 2.21 rc2

* A few more tweaks to reduce some MPI communications during mesh distribution. See note below about backward incompatibility

08/27/2014
MSTK v 2.21 rc1

**NOTE: This release is backward INCOMPATIBLE with v 2.20 rc1 and older**
* Changes to make partitioning and distribution of meshes from processor 0 to other processor more efficient
* MSTK_Mesh_Distribute, MESH_Partition function arguments changed to allow caller to ask for deletion of input mesh shortly after partitioning
* MSTK_Send* and MESH_Send* function arguments modified to send back pending MPI requests as well as pointers to memory allocated. This allows code to defer the blocking call after sending mesh information.

06/07/2014
MSTK v 2.20 rc1

**NOTE: This is a backward INCOMPATIBLE release if you use mesh modification functions in MSTK**
* Added argument to ME_Collapse routines to return list of deleted entities (API change)
* Added option to Merge routines to ignore conformity with topology of the underlying domain if application requests it (API change). This is required for ME_Collapse to be able to ignore conformity with the topology of the underlying domain.
* Added code to prefer collapsing from the interior to the exterior and preserving exterior entities

05/07/2014
MSTK v 2.12

* Make the default library name libmstk.a for both debug and optimized builds
* Can request '-d' suffix for debug library by setting option INSTALL_ADD_DEBUG_SUFFIX to true
**WARNING: This release could break your application build if you link to the MSTK debug library**
* Support upto 512 faces per region
* Enhance MESH_Renumber to enable Reverse Cuthill Mckee and Gibbs-Poole-Stockmeyer algorithms
* Add routine ME_MultiSplit to split edge at multiple points at once
* Avoid renumbering in GMV export if IDs are already contiguous
* Minor suppression of error messages in parallel partitioning
* Read MSTK file only on processor 0 in parallel runs
* quicksort, graph renumbering utilities
* New routines to return lists of entity IDs instead of entity pointers
* Fix some memory leaks
* Made MEntity functions easier to understand by using bitfields
* Change CMake policy CMP0017 to new in order to avoid warnings

03/05/2014
MSTK v 2.12 rc1

* Modified to allow faces to have upto 64 edges

02/26/2014
MSTK v 2.11

* Fixed a bug in which MR_Set_Faces called on an existing region was not updating the element type. Refs Ticket #85.

02/21/2014
MSTK v 2.11 rc5

* A few more minor leak fixes to discovered by valgrind tool 'massif'

02/18/2014
MSTK v 2.11 rc4

* Fixed some old memory leaks and new one introduced when adapting the code to use Metis 5.x.x

02/10/2014
MSTK v 2.11 rc3

* Fixed import of polygonal meshes which was failing
* Some other minor fixes

02/06/2014
MSTK v 2.11 rc2

* Reinstated checking for mesh set name 'sideset_*' before counting it as a sideset
* Added option of using INSTALL_PREFIX_ARCHOS as way of making the install go to a machine specific subdir instead of the generic INSTALL_DIR/lib location

01/27/2014
MSTK v 2.11 rc1

* Added code to support Metis 5.x.y (Metis interface has changed)
* If you use Metis 5.x.x as your partitioner then you must add the line -D METIS_MAJOR_VER=5 to your CMake configure line.* Overlapping sidesets correctly handled when partitioning Exodus II files
* Added functions for splitting faces and edges along with connected elements in simplex meshes.
* Changed definition of MF_Split. Previous MF_Split function is now renamed MF_Split_with_EdgeLoop

11/18/2013
MSTK v 2.10

* Fixed some CMake issues with HDF5 and CMake 2.8.12
* Fixed a Mac compilation issue by removing some obsolete functions

11/14/2013
MSTK v 2.10 rc5

* Minor commit to avoid normalization of vector in MSTK_VNormalize3 in case of zero length.

11/04/2013
MSTK v 2.10 rc4

* Fixed bug in writing out polyhedral cells
* Fixed top level CMakeLists.txt to look for local headers before going to other sources
* Minor clean ups

09/19/2013
MSTK v2.10 rc3

* Minor edit to suppress debug messages from Zoltan. Also, tweak of unit test which was failing after last commit.

09/04/2013
MSTK v2.10 rc2

* Changed parameters of Zoltan RCB partitioning in XY to RCB_RECTILINEAR_BOX=1 to ensure mesh columns are on the same processor
* Added function FixColumnPartitions? to enforce above condition if Zoltan failed to meet this requirement

08/21/2013
MSTK v 2.10 rc1

* Had to bump up minor revision number because of change in arguments of one public function MESH_BuildClassification
* MESH_BuildClassification now takes use_geometry argument indicating
if it should refine the mesh classification based on expensive geometric checks
* 'meshconvert' now accepts --classify=2 in addition to 1 to indicate
that it should do these more expensive classification checks. If the option is --classify=1, mesh classification will be computed based only on topological information
* Import from Nemesis I files now discards the global IDs in nemesis files
which can be discontinuous on a processor. Instead it recomputes them.
* ParallelCheck checks for contiguous global IDs on a processor

08/19/2013
MSTK v 2.03 rc1

* Fixed bug related to reading parallel Exodus II/Nemesis I files
* Fixed bug related to exporting parallel Exodus II/Nemesis I files
* Fixed bug in meshconvert related to partitioning an imported Exodus II file

08/13/2013
MSTK v 2.02

* Added partitioning method 2 to mesh partitioning and distribution This uses Zoltan's Recursive Coordinate Bisection to partition the mesh only in XY plane by suppressing the third coordinate.
* Fixed bug in writing out vertex IDs in Exodus II file in parallel meshes

07/17/2013
MSTK v 2.02 release candidate 3

* Fixed bugs in reading and writing of polyhedral elements to Exodus II files.

05/20/2013
MSTK v 2.02 release candidate 2

Added the file MSTKConfig.cmake.in that is essential for writing out the options and TPLs that MSTK is configured with into the MSTK install dir. It was claimed that this option was available in MSTK 2.01 but since this file was omitted the configuration file was not being written out.

If you want to use the MSTK configuration file to simplify building of MSTK based applications, use this version

05/14/2013
MSTK v 2.02 release candidate 1

* Added function MSet_Rename
* Fixed some memory leaks
* CMake searches for libexodus.a before libexoIIv2c
* Fixed error in exporting X3D files in meshconvert utility
* If your code does not complain about MSet_Rename being undefined while linking with MSTK, you can continue using MSTK v 2.01

NOTE 5/21/2013: There was an error in this push. I forgot to add MSTKConfig.cmake.in As a result, the configuration files were not being written

03/22/2013
MSTK v 2.01 (final)

* make install now installs a file called MSTKConfig.cmake containing important information about the MSTK build. This can be used to simplify the process of building an application using MSTK
* added example do-configure and CMakeLists.txt file to help build MSTK based application
NOTE 5/21/2013: There was an error in this push. I forgot to add MSTKConfig.cmake.in As a result, the configuration files were not being written

03/19/2013
MSTK v 2.01 release candidate 2

* Minor fix to MESH_ExportToFile to not call MPI_Comm_* functions if MSTK_comm is NULL
* Some fixes of unitialized variable warnings

03/06/2013
MSTK v 2.01 release candidate 1

* Update Ghost and Overlap lists if necessary when an entity is removed from the mesh
* Added a parallel adjacency status flag to Mesh to indicate if the info of the
number of entities to be received from other processors is current or not
* Set parallel adjacency status flag to stale when entities are removed from the mesh
* Set parallel adjacency status flag to current after calling MESH_Update_ParallelAdj
* Check Parallel Ajdacency flag before updating attributes and if it is stale, call MESH_Update_ParallelAdj

11/26/2012
MSTK v 2.0 stable

* MAJOR RELEASE: Note the version number change
* MSTK 2.0 is NOT backward compatible
* Bug fixes to MESH_SendRecvMSet that was causing failure with MPICH but not OpenMPI
* Minor fixes to suppress some warnings outside of debug mode

10/22/2012
MSTK v 2.0 release candidate 4

* MAJOR RELEASE: Note the version number change
* MSTK 2.0 is NOT backward compatible
* Bug fixes to rc4


MSTK v 2.0 release candidate 3
* Converted mesh send/receive to non-blocking calls
* Efficiency improved for large meshes distributed over a large number of processors

10/02/2012
MSTK v 2.0 release candidate 2

* MAJOR RELEASE: Note the version number change
* MSTK 2.0 is NOT backward compatible
* Improvements to importing of meshes from Exodus II files
* Ability to improve Nemesis I files (these are augmented Exodus II files for parallel applications)
* Efficiency improvements to parallel distribution
* Eliminated some major leaks

09/07/2012
MSTK v 2.0 release candidate 1

* MAJOR RELEASE: Note the version number change
* MSTK 2.0 is NOT backward compatible
* This version has parallel read capability built into Exodus II and FLAG X3D formats. Depending on the files on disk, the code will read in a serial file and automatically distribute it to the number of processors or read in parallel files and weave them together to make parallel connections
* Can also write parallel Exodus meshes
* meshconvert can now be used to read in a mesh in some format, partition the mesh and write it out in the same or another format - For example, this can be used to partition Exodus files (instead of nemspread)
* Zoltan can now be used as partitioner for MSTK meshes
* Arguments of functions like MESH_ImportFromExodusII have changed - some of these need to be supplied with the communicator.
* Some minor bug fixes

07/26/2012
MSTK v 1.86 release candidate 2

* Major change in the distribution of meshes to multiple processors. Anyone reading and distributing Exodus II files with sidesets will see a huge performance improvement for large meshes
* Minor fixes to keep compilers happy

07/16/2012
MSTK v 1.86 release candidate 1

* New parallel functionality to weave distributed meshes
* Many related changes under the hood
* Developers see only the addition of the functions MESH_Weave_DistributedMeshes
and MESH_UpdVtxCoords

07/16/2012
MSTK v 1.85 final release

* Bug fix in MESH_ImportFromExodusII.c
* Bug fix in MESH_AddGhost.c
* New unit tests for parallel partitioning in 3D

06/26/2012
MSTK v 1.85 release candidate 4 released

* Bug fix to sidesets in Exodus II output

06/21/2012
MSTK v 1.85 release candidate 3 released.

* Global IDs of entities in meshes partitioned on processor 0 are now sequential and without holes
* Changed the standard numbering of regions' vertices for standard element types
* New Find modules for CMake
* New unit test for removing edge from mesh face and face from mesh region
* Number of bug fixes

02/08/2012
MSTK v 1.85 release candidate 2 released.

* Removes use of MPI_COMM_WORLD in PMSTK.c

02/02/2012
MSTK v 1.85 release candidate 1 released.

* This release contains an important fix for parallel mesh partitioning!! Other changes include fixes for memory leaks.

12/06/2011
MSTK v 1.84
* Minor release to prefix internal symbols ERROR, WARN, MESG and FATAL with MSTK_ due to name conflict with another package

11/21/2011
MSTK v 1.83
* Has improvements and bug fixes especially for the 2D parallel mesh distribution

06/17/2011
MSTK 1.83rc3

* Build system revamped to use Find*.cmake modules for package discovery.

Version 1.83rc2
* CMakeLists.txt file now takes user hints about location of HDF5 libraries

06/15/2011
MSTK 1.83rc1

* Fixed CMakeLists.txt so that it can recognize HDF5 hints Fixed bugs in recognizing tetrahedral and hexahedral elements Put in code to properly answer queries for prisms and pyramids Introduced unit tests

* Changed chkmesh utility so that it distinguishes between truly inverted elements and those that are just not star shaped (some tets in the decomposition of non-simplex elements are inside out but the total volume of the element is posit ive)

06/02/2011
MSTK 1.82

* Bug fix in ME_Faces_R1R2 Can read R1 type representation into F1 type mesh Added utils with meshconvert and chkmesh Removed use of unistd.h from MESH_ExportToDXBin to enable Windows compile Fixed GMV export of Prisms and Pyramids

05/11/2011
MSTK 1.81 

* Bug fixes and improvement to parallel functionality Support for edge collapse

02/18/2011
MSTK 1.8

* MSTK with support for distributed meshes released.


