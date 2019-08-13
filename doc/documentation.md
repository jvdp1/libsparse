## Brief documentation  


## Overview  
The Fortran 2003 library __libsparse__ is a library that provides objects to create and handle rectangular and square sparse matrices using different formats:  

 * Linked List (LL);  


 * COOrdinate storage (COO) (with elements stored using a hashing function);  


 * Compressed Row Storage (CRS).   


The library is written following an object-oriented approach.  


## Compilation  
The library relies on different libraries, such as BLAS/LAPACK libraries, PARDISO (at this stage, Intel MKL PARDISO), and [METIS 5] (http://glaros.dtc.umn.edu/gkhome/metis/metis/overview).  

For compilation, go in the directory *libsparse/src/lib* and type *make* for the default options. By default, it will not be compiled against METIS 5. To compile with METIS 5, *make METISENABLE=1*  


Compilation with debug options (`-g -check all -traceback`) is possible by adding the argument `DEBUGENABLE=1`.  

Compilation for single precision is possible by adding the argument `DPENABLE=0`.  


## Methods  

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

This option can also be used to __load__ a matrix from a stream file.  


To __add__ a value *val* at the position (*row*,*col*) of a sparse matrix, the method *add* must be used, e.g.:  

````   
call mat%add(row,col,val)  
````  

To __convert__  a sparse matrix from a format to another format, the assignment = can be used. E.g., to convert from COO to CSR:  

````   
type(coosparse)::coomat  
type(crssparse)::crsmat  
  
csrmat=coomat  
````  

To __copy__ a sparse matrix, the assignment = can be used. E.g.,:  

````   
type(crssparse)::crsmat  
type(crssparse)::crsmatcopy  
  
csrmatcopy=crsmat  
````  

To __deallocate__  a sparse matrix, the method *destroy* can be used, e.g.:  

````   
call mat%destroy()  
````  

To extract the __diagonal elements__ of a sparse matrix into a array, the method *diag* must be used:  

````  
array=mat%diag()  
````  

To extract the __diagonal elements__ and __off-diagonals__ of a sparse matrix into a array, the method *diag* must be used:  

````  
extractedmat=mat%diag(x)  
````  


where *x* is the number of off-diagonals (next to the diagonal) that must be extraced. If *x* is equal to 0, only the diagonal will be extracted and stored into a sparse matrix.  

To __get a value *val*__ at the position (*row*,*col*) of a sparse matrix, the method *get* must be used, e.g.:  

````   
val=mat%get(row,col)  
````  

To __get one of the dimensions__ of a sparse matrix, the method *getdim* must be used:  

````   
val=mat%getdim(x)  
````  
where *x* is 1 (=*number of rows*) or 2 (=*number of columns*).  

To __get a permutation vector__ from the METIS 5 fill-reducing ordering approach, the method *getordering* can be used:  

````  
permarray=mat%getordering()  
````  

Options for METIS 5 can be changed through optional arguments of this method.  


To __create__ a CRS matrix from __existing arrays__, the method *external* must be used:  

````  
call crsmat%external(ia,ja,a)  
````  

where the arrays *ia*, *ja*, and *a* must have the same size as the ones of the sparse matrix *crsmat*.  

To __multiply__ a sparse matrix by a vector as, *y = alpha \* mat(trans) \* v + val \* y* , the method *multbyv* must be used:  

````  
call mat%multbyv(alpha,trans,v,val,y)  
````  
where *alpha* and *val* are double-precision real values, *v* and *y* are vectors, and *trans* (= 'n' or 't') relates to the transposition of the matrix.   
The method *multbyv* is based on the MKL Sparse BLAS library.  

To get the number of __non-zero elements__ of a sparse matrix, the method *nonzero* must be used, e.g.:  

````  
nonzeros=mat%nonzero()  
````  


To __print__ a sparse matrix as stored __to default output__ (screen), the method *print* must be used, e.g.:  


````  
call mat%print()  
````  

To __print__ a sparse matrix as stored __to a file__ called *file.dat*, the method *printtofile* must be used, e.g.:  


````  
call mat%printtofile('file.dat')  
````  


To __save__ a matrix in the internal format, the method *save* can be used:  

````  
call mat%save('file.stream')  
````  

To __set an entry__ to a specified value (even __0__), the method *set* can be used:  

````  
call mat%set(row,col,val)  
````  

To __set a permutation vector__,the method *setpermutation* can be used:  

````  
 call mat%setpermutation(array)  
````  

It is possible to get and set a permutation vector in one go as follows:  

````  
 call mat%setpermutation(mat%getordering())  
````  


To __solve__ a linear system of equations of the form *mat \* x = y*, the method *solve* must be used:  

````  
call mat%solve(x,y)  
````  
The method *solve* is based on Intel MKL Pardiso. If a permutation vector was provided with the method `setpermutation`, this permutation vector will be used by Intel MKL Pardiso (instead of determining it internally).  


To __sort__ a column (in an ascending order) within a row of a CRS sparse matrix, the method *sort* must be used:  

````  
call crsmat%sort()  
````  

To check if the sparse matrix is a *square* matrix, the method __lsquare__ must be used:  

````  
square=mat%lsquare()  
````  

where the variable __square__ is a logial.  


To extract a __submatrix__ from a sparse matrix, the method *submatrix* must be used, e.g.:  

````  
submatrix=mat%submatrix(startrow,endrow,startrow,endrow,lupper=log)  
````  

where *log* is a logical to extract the full matrix (lupper=.false.) or the upper triangular matrix (lupper=.true.).  
