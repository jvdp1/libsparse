# Fortran library to handle sparse matrices  


## Overview  
The Fortran 2003 library __libsparse__ is a library that provides objects to create and handle rectangular and square sparse matrices using different formats:  

 * Linked List (LL);  


 * COOrdinate storage (COO) (with elements stored using a hashing function);  


 * Compressed Row Storage (CRS).   


The library is written following an object-oriented approach. It has been tested mainly on small datasets.  



## Compilation  
To build the `libsparse` you need (at least):

 * at least a Fortran 2008 compliant compiler (GCC Fortran 11 and Intel Fortran
   classic compilers have been tested successfully);
 * Intel MKL library;
 * Make or CMake or fpm.

The library relies on different libraries, such as BLAS/LAPACK libraries (currently on Intel MKL library), and optionally on PARDISO (at this stage, Intel MKL PARDISO), and on [METIS 5](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview).  
The library can be built with the compilers `gfortran` and `ifort`.


See the brief [documentation](doc/documentation.md) for more details regarding the compilation.  


## Documentation  
The brief documentation is available in the directory *doc* (see *mainpage.md*). An extended documentation can be generated with *Doxygen*.  


## Usage examples  

To __evaluate the quadratic form__ `x' A x` for a sparse matrix `A` and a dense vector `x`:  

````  
use modsparse, only: crssparse  
type(crssparse) :: A  
real(8) :: x(100), q  

!... build A and set x ...  
q = A%xAx(x)  
````  

To __evaluate the bilinear form__ `x' A y` for a sparse matrix `A` and two dense vectors `x` and `y` (`xAx` is the special case `x = y`):  

````  
use modsparse, only: crssparse  
type(crssparse) :: A  
real(8) :: x(100), y(100), q  

!... build A and set x, y ...  
q = A%xAy(x, y)  
````  

To __compute the trace of the product of a square sub-block of a sparse matrix with a dense matrix__, e.g. `trace(A(3:6, 3:6) * B)`:  

````  
use modsparse, only: crssparse  
type(crssparse) :: A  
real(8) :: B(4,4), t  

!... build A and set B ...  
t = A%traceproduct(3, 6, 3, 6, B)  
````  

*xAx*, *xAy*, and *traceproduct* work on matrices in COO and CRS format and run in O(nnz). `xAy` also works on non-square matrices (length of `x` = rows, length of `y` = columns). Indices in *traceproduct* are 1-based and inclusive.  
 

## Acknowledgements  

The code in modspainv.f90 is based on code originally written by Karin Meyer
(didgeridoo.une.edu.au/womwiki/doku.php?id=fortran:fortran). The incorporated and
modified code is re-licensed under the MIT license with express permission from
Karin Meyer. We are grateful for her permission to use and adapt this work.

This library was inspired by several sources:  


 * http://burtleburtle.net/bob/hash/index.html#lookup  


 * https://didgeridoo.une.edu.au/km/homepage.php  


 * https://genomeek.wordpress.com/  


 * https://gist.github.com/n-s-k/522f2669979ed6d0582b8e80cf6c95fd  


 * https://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=sparse90.pdf  


 * https://www.netlib.org/lapack/explore-html/index.html  


 * https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/755612  


 * and by many courses related to object-oriented programming and Fortran 2003/2008.  

## To be implemented  

 * Check for symmetric matrix  

 * Allow the option spainv + single precision

 * Full support of 8-byte integers
