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
 * Make or fpm.

The library relies on different libraries, such as BLAS/LAPACK libraries (currently on Intel MKL library), and optionally on PARDISO (at this stage, Intel MKL PARDISO), and on [METIS 5](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview).  
The library can be built with the compilers `gfortran` and `ifort`


See the brief [documentation](doc/documentation.md) for more details regarding the compilation.  


## Documentation  
The brief documentation is available in the directory *doc* (see *mainpage.md*). An extended documentation can be generated with *Doxygen*.  


## Acknowledgements  
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
