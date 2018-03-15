#
# Macro to compile the OpenMP test (Mac)
#
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$install_folder/ext_lib

clang++ -Xpreprocessor -fopenmp -std=c++11 -I../ext_include/omp/ -L../ext_lib/ -lomp test_OpenMP.cpp -o test_OpenMP