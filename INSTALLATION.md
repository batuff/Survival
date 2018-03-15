Linux
=====

Prerequisites
-------------

To correctly compile the software, please be sure to have:

- The GNU Compiler Collection (GCC) (http://gcc.gnu.org) compatible with the C++11 version (information about C++11 compatible compilers at: http://en.cppreference.com/w/cpp/compiler_support)
- The GNU Scientific Libraries (GSL) installed (Free download from the gsl website: http://www.gnu.org/software/gsl/)  
- The OpenMP libraries (information about OpenMP at: http://openmp.org)

    
### Installation of prerequired libraries

Usually, in the Linux distributions the GCC is already installed.

#### Debian/Ububtu systems

You can install the GSL libraries using the following command:

    sudo apt-get install libgsl-dev

Youe can install the OpenMP libraries using the following command:

    sudo apt-get install libpomp-dev
    
#### RedHat/CentOS systems

You can install the GSL libraries using the following command:

    sudo yum install gsl-devel

You can install the OpenMP libraries using the following command:

    sudo yum install libgomp-devel


Compilation
-----------

1. Download the whole *Survival* folder
2. In the command line change the active directory to the folder and simply type:

        make


Running
-------

1. Modify the *setenv.sh* file with the correct path
2. In the command line, type:

        source setenv.sh

3. In the command line, type:

        survival --help

Specific usage details are given in the README.md file.

 
MacOS X
=======

Prerequisites
-------------

To correctly compile the software, please be sure to have:

- XCode and the Command Line Tools installed on the system
- The GNU Scientific Libraries (GSL) installed (Free download from the gsl website: http://www.gnu.org/software/gsl/)  
- The OpenMP libraries (information about OpenMP at: http://openmp.org)
    
### Installation of prerequired libraries

Install XCode from the App Store.
Once installed, from the Terminal type:

    xcode-select --install

to install the Command Line Tools.

Note that the GSL and OpenMP pre-compiled libraries for MacOS High Sierra are already present in the "Survival" GitHub repository.

Compilation
-----------

1. Download the whole *Survival* folder
2. In the command line change the active directory to the folder and simply type:

        make -f Makefile_clang


Running
-------

1. Modify the *setenv.sh* file with the correct path
2. In the command line, type:

        source setenv.sh
        
2. In the command line, type:

        survival --help

Specific usage details are given in the README.md file.


Windows
=================

The software has never been tested on *Windows*. Nevertheless, it technically works (under the previous requirements).  
The test is currently in progress...


**NOTES:**
The correct compilation and execution of "Survival" has been tested only on Ubuntu 14.04, Ubuntu 16.04 and macOS High Sierra (v10.13.3)