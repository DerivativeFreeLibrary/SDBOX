-----------------------------------------------------------
 How to use the derivative-free optimizer SDBOX for 
 bound constrained optimization problems
-----------------------------------------------------------
 The package provides a C and FORTRAN90 version of the code.

0- Gunzip and untar the archive in a folder on your computer by
   issuing in a directory of your choice (ex. curdir) the command

   $> tar -xvf SDBOX.tar.gz

2- FORTRAN VERSION OF THE CODE:

   Edit file curdir/F90_src/problem.f90 to define your own objective function.
   In particular, modify the subroutines 
   setdim    : which sets problem dimension
   setbounds : which sets upper and lower bounds on the variables
   startp    : which sets the starting point
   funct     : which defines the objective function

   PYTHON VERSION OF THE CODE:

   Edit file curdir/PYTHON_src/problem.py to define your own objective function.
   In particular, modify procedures 
   setbounds : which sets upper and lower bounds on the variables
   startp    : which sets the starting point
   funct     : which defines the objective function

   C VERSION OF THE CODE:

   Edit file curdir/C_src/problem_c.c to define your own objective function.
   In particular, modify the subroutines 
   setdim    : which sets problem dimension
   setbounds : which sets upper and lower bounds on the variables
   startp    : which sets the starting point
   funct     : which defines the objective function

   JULIA INTERFACE
   Edit file curdir/Julia_interface/example.jl to define your own objective function
   along with appropriate number of variables, upper and lower bounds on the variables,
   starting point, maximum number of iterations, tolerance.

   
2- At command prompt in curdir execute 

     $> make
 
   which will create the executables 'sdbox_c' (using C sources),
   'sdbox_f' (using FORTRAN90 sources) and will create 'libsdbox.a' and
   place it in the directory 'Julia_interface'

4- execute

     $> ./sdbox_f
     $> ./sdbox_c
	 $> python PYTHON_src/sdbox.py

5- execution within Julia. From within the curdir/Julia_interface directory let Julia start.
   At the Julia prompt type

     julia> include("example.jl")

6- A Matlab version of the code is also provided. File sdbox.m can be found in the folder

     curdir/Matlab

   Also, file mani.m contains example usage of sdbox to optimize the powell function
