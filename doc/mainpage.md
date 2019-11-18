## Libsparse Documentation  


## Overview  
The Fortran 2003 library __libsparse__ is a library that provides objects to create and handle rectangular and square sparse matrices using different formats:  

 * Linked List (LL);  


 * COOrdinate storage (COO) (with elements stored using a hashing function);  


 * Compressed Row Storage (CRS).   


The library is written following an object-oriented approach. It has been tested mainly on small datasets.  



## Compilation  
The library relies on different libraries, such as BLAS/LAPACK libraries, PARDISO (at this stage, Intel MKL PARDISO), and [METIS 5] (http://glaros.dtc.umn.edu/gkhome/metis/metis/overview).  


See the brief [documentation] (documentation.md) for more details regarding the compilation.  


## Documentation  
The brief [documentation] (documentation.md) is available in the directory *doc*. An extended documentation can be generated with *Doxygen*.  


## Acknowledgements  
This library was inspired by several sources:  


 * http://burtleburtle.net/bob/hash/index.html#lookup  


 * https://didgeridoo.une.edu.au/km/homepage.php  


 * https://genomeek.wordpress.com/  


 * https://gist.github.com/n-s-k/522f2669979ed6d0582b8e80cf6c95fd  


 * https://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=sparse90.pdf  


 * https://www.netlib.org/lapack/explore-html/index.html  


 * https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/755612  


 * https://stackoverflow.com/questions/466204/rounding-up-to-next-power-of-2   


 * and by many courses related to object-oriented programming and Fortran 2003/2008.  



