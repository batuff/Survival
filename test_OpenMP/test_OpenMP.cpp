/* Test to check if the compiler support the installed version od the OpenMP libraries */

#include <omp.h>
#include <iostream>

int main() {
    #pragma omp parallel
    {
        printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
    }
    return 0;
}