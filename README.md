# numerical_cpp
#1
Poisson Equation Solver with many new features of C++11 and Matrix Template Library 4

#2
Solved a system of linear equations without using the Standard Template Library
Basically my own vector class (although named it array), with some aspects like 
resizing, {} intialization, fill, etc that have analogues in STL.
The test problem I used was from my Numerical Analysis book by Burden & Faires
but the code is written such that a larger size matrix can be solved too. 
When I get some more time I plan to implement better solver algorithms like
Gauss-Sidel with Successive Over Relaxation (SOR). After that I would like to 
write my own finite element solver for Navier Stokes equations, but still
need to learn more fluid dynamics first.
