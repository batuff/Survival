Requirements
=================

To correctly compile the software, please be sure to have:
 - The GNU Scientific Libraries (GSL) installed (Free download from the gsl website: http://www.gnu.org/software/gsl/)
 - A compiler supporting the OpenMP and compatible with the C++11 version. (Information about OpenMP compilers at: http://openmp.org/wp/openmp-compilers/; and about C++11 compatible compilers at: http://en.cppreference.com/w/cpp/compiler_support)

Linux and MacOS X
=================

The current version supports the MacOS X up to the 10.10 version (Yosemite).
Compatibility for MacOS X 10.11 (El Capitan) is not guaranteed (under development).

Please:
 1. Download the whole *Survival* folder
 2. Modify the *setenv.sh* file with the correct path
 3. In the *command line*, type:  
 `> source setenv.sh` to load the setenv.sh settings  
 `> make` to create the executable file named *survival*

Windows
=================

The software has never been tested on *Windows*. Nevertheless, it technically works (under the previous requirements).
The test is currently in progress...
