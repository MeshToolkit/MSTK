
------------------------------------------------------------------------------
A.  -------- Install MSTK ----------
------------------------------------------------------------------------------

Starting from MSTK v1.8, we are including CMakeLists.txt files to allow
users to build and install MSTK on a wide variety of platforms with various 
options enable/disabled easily. If you are unfamiliar with CMake, it is an 
open-source project configuration system developed by Kitware, Inc. (See 
http://www.cmake.org). 

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

------------------------------------------------------------------------
B.      --------- Learn MSTK  --------------
------------------------------------------------------------------------

1. Read the manual !!!!! It is in .../mstk-N.n/docs/


------------------------------------------------------------------------
B.      --------- Use MSTK for developing applications ----------------
------------------------------------------------------------------------

1. Include "MSTK.h" in any code that uses MSTK

2. Link with the MSTK libraries and any third party libraries that 
   were used to build MSTK.

   To simplify the process of building an end application using MSTK,
   a CMake configuration file called MSTKConfig.cmake is written out
   when you 'make install' MSTK.

   The directory config/application_cmake contains a CMakeLists.txt
   file and a do-configure script to use as a "rough" template for
   developing a CMake based build system for an application using MSTK



NOTES:

I have tested MSTK only on a Linux systems, on HPC Unix systems and on
MacOS with GNU, PGI and Intel compilers. If it doesn't work on a
particular system, please let me know. 

Feel free to contact me if you need any assistance. Also, see the MSTK
website https://software.lanl.gov/MeshTools/trac/

Regards
Rao


-- 
--------------------------------------------------------------------
Rao V Garimella  			Tel: (505) 665-2928
T-5, MS B284 	 	   '\   '\ 	FAX: (505) 665-5757
Los Alamos National Lab	     ( )-( )   	Email: rao@lanl.gov
Los Alamos, NM 87545	                http://math.lanl.gov/~rao		
--------------------------------------------------------------------

	


