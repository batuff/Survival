Licence
=================

Before starting the installation process, please read the LICENCE of the software (./LICENCE).

Unix and MacOS X
=================

To correctly compile the software, please be sure to have:
 - The GNU Scientific Libraries (GSL) installed (Free download from the gsl website: http://www.gnu.org/software/gsl/)
 - A compiler supporting the OpenMP and compatible with the C++11 version. (Information about OpenMP compilers at: http://openmp.org/wp/openmp-compilers/)

Then:
 1. Download the whole *Survival* folder
 2. Modify the *setenv.sh* file with the correct path
 3. In the *command line*, type:  
 `> source setenv.sh` to load the setenv.sh settings  
 `> make` to create the executable file named *survival*

Windows
=================

Not supported yet. Work in progress...
