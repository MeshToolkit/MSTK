What is MSTK
============

MSTK is a mesh framework that allows users to represent, manipulate
and query unstructured 3D arbitrary topology meshes in a general
manner without the need to code their own data structures. MSTK is a
flexible framework in that it allows a variety of underlying
representations for the mesh while maintaining a common interface. It
allows users to choose from different mesh representations either at
initialization (implemented) or during the program execution (not
yet implemented) so that the optimal data structures are used for the
particular algorithm. The interaction of users and applications with
MSTK is through a functional interface that acts as though the mesh
always contains vertices, edges, faces and regions and maintains
connectivity between all these entities.

MSTK allows for the simultaneous existence of an arbitrary number of
meshes. However, any entity in MSTK can belong to only one mesh at a
time. MSTK also allows applications to attach application or field
data to entities. This data may be integers, reals and pointers.

MSTK supports distributed meshes for parallel computing. In the future
it will allow for parallel mesh modification.

MSTK is not a mesh generator but the infrastructure it provides can be
used to develop sophisticated mesh generators and other mesh based
applications.

**MSTK is not related in anyway to STK or STKmesh in the Trilinos suite
of software products**

I have tested MSTK only on Linux systems, on HPC Unix systems and on
MacOS with GNU, PGI and Intel compilers. If it doesn't work on a
particular system, please let me know. If you have a patch, feel free
to submit a pull request.

Install MSTK
============

Starting from MSTK v1.8, we are including CMakeLists.txt files to allow
users to build and install MSTK on a wide variety of platforms with various 
options enable/disabled easily. If you are unfamiliar with CMake, it is an 
open-source project configuration system developed by Kitware, Inc. (See 
[http://www.cmake.org](http://www.cmake.org)). 

The following are the steps for build and installing MSTK using cmake

1. Verify that you have cmake or install it on your system (system-wide or
   locally)

2. Examine and Modify the file config/do-configure-mstk. This is the driver 
   file for telling the build system what options are needed in the MSTK 
   build. DO NOT HACK the CMakeLists.txt files. Try to control the build 
   system behavior through the configure script. Using this script, you can 
   turn on/off Exodus support, parallel mesh support, specify where to find
   third-party libraries, where to install the mstk libraries etc. The most
   up-to-date options for the configuration are given in the example configure
   script supplied with the distribution (config/do-configure-mstk).

3. Choose a build directory. This could be a subdirectory named build in
   the mstk source tree or some other directory. If you want to build both
   debug and optimized targets you have to run cmake and build in separate
   subdirectories or one build will clobber the other. I like to create 
   subdirectories build/Debug and build/Release in the mstk directory and
   build the different versions there.

4. cd to the build directory and run the do-configure-mstk script there.
   This will create the appropriate Makefiles.

5. Then run 'make' followed by 'make install' or just run 'make install'
   right away. 'make VERBOSE=1' will display details about what compile 
   commands make is using.

6. 'make install' not only installs the library and the include files
   but also a CMake configuration file called MSTKConfig.cmake. This
   file contains many important CMake variables used in building MSTK
   and can simplify the process of building an executable based on
   MSTK. See section C for more detail.


Learn MSTK
==========

1. Read the manual ! It is included with the distribution. While it
may not be the most up-to-date, it has the essentials


Use MSTK for developing applications
====================================

1. Include "MSTK.h" in any code that uses MSTK

2. Link with the MSTK libraries and any third party libraries that 
   were used to build MSTK.

   To simplify the process of building an end application using MSTK,
   a CMake configuration file called MSTKConfig.cmake is written out
   when you 'make install' MSTK.

   The directory config/application_cmake contains a CMakeLists.txt
   file and a do-configure script to use as a "rough" template for
   developing a CMake based build system for an application using MSTK

Contact:
========

Feel free to contact me if you need any assistance at rao@lanl.gov. If
you are successfully using MSTK in your work, drop me a note.

--------------------------------------------------------------------
Rao V Garimella, T-5, MS B284, Los Alamos National Laboratory  
Los Alamos, NM 87544 USA  
Tel: (505) 665-2928  
Email: rao@lanl.gov  
http://math.lanl.gov/~rao  	

--------------------------------------------------------------------

	


