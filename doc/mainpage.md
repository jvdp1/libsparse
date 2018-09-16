## Short documentation for the Fortran 2003 library libsparse  


## Overview  
The Fortran 2003 library __libsparse__ is a library that provides objects to create and handle rectangular and square sparse matrices using different formats:  

 * Linked List (LL);  


 * COOrdinate storage (COO);  


 * Compressed Row Storage (CRS).   


The library is written following an object-oriented approach.  



This library was inspired by several sources:  


 * https://didgeridoo.une.edu.au/km/homepage.php  


 * https://genomeek.wordpress.com/  


 * https://gist.github.com/n-s-k/522f2669979ed6d0582b8e80cf6c95fd  


 * https://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=sparse90.pdf  


 * https://www.netlib.org/lapack/explore-html/index.html  


 * https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/755612  


 * https://stackoverflow.com/questions/466204/rounding-up-to-next-power-of-2 


 * and by many courses related to object-oriented programming and Fortran 2003/2008.  


## Use  

Below the term *mat* refers to one of the three available sparse objects (i.e., llsparse, coosparse, or crssparse), except if it is mentioned otherwise.  

To __construct__ (__initiate__) a sparse matrix, the constructor of the same name as the object, must be used, e.g.:  

````   
integer::dim1=100     !number of rows     
type(coosparse)::mat  
  
mat=coosparse(dim1)  
````  

Other options are possible for the constructor (see details in the module pages). For sparse matrices with upper storage (default: fulled storage), it must be mentioned as:  

````
mat=coosparse(dim1,lupper=.true.)  
````

To __add__ a value *val* at the position (*row*,*col*) of a sparse matrix, the method *add* must be used, e.g.:  

````   
call mat%add(row,col,val)
````  

To __get__ a value *val* at the position (*row*,*col*) of a sparse matrix, the method *get* must be used, e.g.:  

````   
val=mat%get(row,col)
````  

To __convert__  a sparse matrix from a format to another format, the assignment *=* can be used. E.g., to convert from COO to CSR:  

````   
type(coosparse)::coomat  
type(crssparse)::crsmat  
  
csrmat=coomat
````  

To __print__ a sparse matrix as stored __to default output__ (screen), the method *print* must be used, e.g.:  


````
call mat%print()  
````

To __print__ a sparse matrix as stored __to a file__ called *file.dat*, the method *printtofile* must be used, e.g.:  


````
call mat%printtofile('file.dat')  
````


To __set__ an entry to a specified value (even __0__), the method *set* can be used:  

````
call mat%set(row,col,val)
````

To __sort__ (in an ascending order elements) within a row of a CRS sparse matrix, the method *sort* must be used:  

````
call crsmat%sort()  
````


To __deallocate__  a sparse matrix, the method *reset* can be used, e.g.:  

````   
call mat%reset()  
````  


## Compilation  
For the compilation, go in the directory *libsparse/src/lib* and type *make* for the default options.  

Compilation with debug options (`-g -check all -traceback`) is possible by adding the argument `DEBUGENABLE=1`.  










