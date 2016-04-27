Requirements
=================

To correctly compile the software, please be sure to have:
 - The GNU Scientific Libraries (GSL) installed (Free download from the gsl website: http://www.gnu.org/software/gsl/)  
 - A version of the gcc compiler installed that:
    1. Supports the OpenMP libraries (information about OpenMP compilers at: http://openmp.org/wp/openmp-compilers/)
    2. Is compatible with the C++11 version (information about C++11 compatible compilers at: http://en.cppreference.com/w/cpp/compiler_support)  

**NOTES:**  
- To **check your gcc version** type from the command line: `gcc --version` (recommended 4.9 or newer).  
- To **test if your gcc version is compatible with the OpenMP libraries**:
   1. cd to the *test_OpenMP* directory
   2. From the command line type: `./compile_test_OpenMP.sh`: Compilation should proceed with no errors or warnings.
   3. Test the result typing `./test_OpenMP`  
- To **install gcc** (full guide at https://gcc.gnu.org/install/):  
   1. Download gcc-4.xx.yy-bin.tar.gz download or newer from http://hpc.sourceforge.net/
   2. cd to your downloads folder and un-gzip the archive: `gunzip gcc-xx.yy-bin.tar.gz`
   3. In the same folder run `sudo tar -xvf gcc-xx.yy-bin.tar -C /` - this will place new executable to /usr/local/bin
   4. Add the following to your ~/.bash_profile (*Linux Users*) or ~/.profile (*MacOS X Users*)  
      `export PATH=/usr/local/bin:$PATH`  
      `export C_INCLUDE_PATH=/usr/local/include:$C_INCLUDE_PATH`  
      `export CPLUS_INCLUDE_PATH=/usr/local/include:$CPLUS_INCLUDE_PATH`  
      `export LIBRARY_PATH=/usr/local/lib:$LIBRARY_PATH`  
      `export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH`\*\*  
      \*\* **MacOS X Users** must replace LD\_LIBRARY\_PATH with DYLD\_LIBRARY\_PATH
   5. Open new terminal and run `which gcc`. This should point to /usr/local/bin/gcc.
- For **MacOS X 10.11 (El Capitan) Users** (or higher): In the wake of the introduction of the SIP (System Integrity Protection) security feature, the user who wants to install the GSL libraries from the command line has to do it as super user (1. download the GSL; 2. `./configure`; 3. `sudo make`; 4. `sudo make install`). More information about the SIP at https://support.apple.com/en-us/HT204899

Linux and MacOS X
=================

Please:
 1. Download the whole *Survival* folder
 2. Modify the *setenv.sh* file with the correct path
 3. In the command line, type:  
 `> source setenv.sh` to load the setenv.sh settings  
 `> make` to create the executable file named *survival*

Windows
=================

The software has never been tested on *Windows*. Nevertheless, it technically works (under the previous requirements).  
The test is currently in progress...
