Readme:

This fold include matlab codes used for the StemCell project. 
There are three main functions:
1. stemcellreg.m:   The based function used to general simulation.
2. cancerdev.m:     The code to simulate cancer development with mutations to the pathways of proliferation and differentiation.
3. tissuegrowth.m   The code to simulate tissue growth with time dependent differentiation rate, and the terminally differentiated cells are also included.

Other files include the control, initialization, and parameter setting functions used for simulation:
--control.m		The file to control the simulation
--parameter.m		The file to define the parameters
--Initialization.m	The file to initiazlize the system

--beta0.m, beta.m	The files to define the function beta(Q,x)
--eta.m			The function eta(x)
--fnu.m			The function nu(x)
--kappa.m		The function kappa(x)
--mu.m			The function mu(x)
--p.m			The inheritance function p(x,y)
--phi.m			The function phi(x)
--tau.m			The delay tau(x)
--xi.m			The function xi(x) 

Please edit the control, initialization, and parameter files in according to your own simulation.
