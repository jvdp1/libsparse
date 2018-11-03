# Fortran library to handle sparse matrices  


## Overview  
The Fortran 2003 library __libsparse__ is a library that provides objects to create and handle rectangular and square sparse matrices using different formats:  

 * Linked List (LL);  


 * COOrdinate storage (COO) (with elements stored using a hashing function);  


 * Compressed Row Storage (CRS).   


The library is written following an object-oriented approach. It has been tested mainly on small datasets.  



## Compilation  
The library relies on different libraries, such as BLAS/LAPACK libraries, PARDISO (at this stage, Intel MKL PARDISO), and [METIS 5] (http://glaros.dtc.umn.edu/gkhome/metis/metis/overview).  

For compilation, go in the directory *libsparse/src/lib* and type *make* for the default options. By default, it will not be compiled against METIS 5. To compile with METIS 5, *make METISENABLE=1*  


Compilation with debug options (`-g -check all -traceback`) is possible by adding the argument `DEBUGENABLE=1`.  

Compilation for single precision is possible by adding the argument `DPENABLE=0`.  


## Documentation  
The documentation is available in the directory *doc* (see *mainpage.md*). An extended documentation can be generated with *Doxygen*.  


## Acknowledgements  
This library was inspired by several sources:  


 * https://didgeridoo.une.edu.au/km/homepage.php  


 * https://genomeek.wordpress.com/  


 * https://gist.github.com/n-s-k/522f2669979ed6d0582b8e80cf6c95fd  


 * https://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=sparse90.pdf  


 * https://www.netlib.org/lapack/explore-html/index.html  


 * https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/755612  


 * https://stackoverflow.com/questions/466204/rounding-up-to-next-power-of-2   


 * and by many courses related to object-oriented programming and Fortran 2003/2008.  



