# Bachelor's thesis
This repository contains the code I have written for my Bachelor's thesis in Mathematics at the University of Trieste.
The goal of the (small) programming part of the dissertation was to demonstrate the efficiency of Krylov's method when used to compute the action of the matrix exponential.
All the code is written in MATLAB.

The file `Arnoldi.m` contains the implementation of the Arnoldi method.
Additionally, `driver_Krylov_laplacian.m` and `driver_Krylov_random.m` provide examples of the method's application to the matrix resulting from the 1D discretization of the Laplacian and to a random matrix, respectively.

The file `funm2.m` contains the code for MATLAB's built-in function used to compute functions of matrices.
