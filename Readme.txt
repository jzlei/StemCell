Readme:

This fold includes the codes (MATLAB and C++) used for the general mathematical framework model of heterogeneous stem cell regeneration. For details of the model, refer to bioRxiv 592139 (https://doi.org/10.1101/592139).

1. The matlab code (MATLAB2014b) include three part:

    MATLAB-----1_StemCellgeneration
           |
	   |---2_Tissuegrowth
           |
           |---3_Cancerdevelopment


1_StemCellgeneration include the codes to simulate the general process of stem cell regeneration

2_Tissuegrowth includes the codes to simulate tissue growth with time dependent differentiation rate, and the terminally differentiated cells.

3_Cancerdevelopment includes the codes to simulation cancer development with mutations to the pathways of proliferation and differentiation.

The matlab code is mainly designed to solve the differential-integral model numerically. Here, only one dimensional epigenetic state x (0<x<1) is considered.

Each folder includes following m-files:

--main.m		The main function

--control.m		The file to control the simulation
--parameter.m		The file to define the parameters
--Initialization.m	The file to initialize the system

--beta0.m, beta.m	The files to define the function beta(Q,x)
--eta.m			The function eta(x)
--fnu.m			The function nu(x)
--kappa.m		The function kappa(t,x)
--fkappa0.m		The function for the time dependent parameter kappa0(t)
--mu.m			The function mu(x)
--p.m			The inheritance function p(x,y)
--phi.m			The function phi(x)
--tau.m			The delay tau(x)
--xi.m			The function xi(x) 

Please edit the control, initialization, and parameter files in according to your own applications.

2. The C++ code to perform single-cell-based stochastic simulation. In this code, we do not solve the differential-integral equation directly, but apply a stochastic simulation to model the growth process of a multiple cell system. A multiple cell system is represented as a collation of epigenetic states, and the single-cell-based stochastic simulation tracks the behaviors of each cell according to their own epigenetic states. 

	C++ --------StemCell    include all source codes. 
               |
               |----Run         the script to run the program


The program has been test in MacOS (10.11.6)

Usage:

In the folder StemCell, use the command sh compile.sh to generate the execute file bct_StemCell.

In the folder Run, edit the input files md.in and par.dat to change the control parameter and the parameter values, and use the command ./bct_StemCell md.in to run the program. 
	
