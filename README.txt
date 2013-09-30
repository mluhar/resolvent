This series of MATLAB functions and scripts performs the resolvent analysis proposed by Dr. Beverley McKeon and Dr. Ati Sharma (see McKeon, B.J. and Sharma, A.S., 2010, A critical-layer framework for turbulent pipe flow. Journal of Fluid Mechanics, 658, 336-382), and offers a few related utilities.  

McKeon and Sharma (2010) employed a projection onto a series of divergence-free basis functions to satisfy mass continuity and eliminate the pressure term in the Navier-Stokes equations, following Meseguer and Trefethen (2003, Linearized pipe flow to Reynolds number 10^7, Journal of Computational Physics, 186 (1): 178-197).  The software provided here does not employ this projection.  The Navier-Stokes resolvent is formulated in terms of primitive variables -- pressure and mass continuity are retained explicitly.  This provides some information on the pressure field, and also allows for more flexibility on the boundary conditions imposed.

In general, the code has been designed to be modular, so that changes in geometry, boundary conditions etc., can be incorporated without too much trouble.  All the functions should be relatively well commented, with lists of inputs and outputs.  In case of any confusion, do not hesitate to contact Mitul (mluhar@caltech.edu).

-----
List of Functions 

1) resolventSVD - constructs the discretized resolvent operator (using a collocation method) and performs a singular value decomposition of the resolvent operator, returning a (user-defined) number of forcing and response modes.  It calls functions 2-5

2) pipeCoords: generates the coordinate system and differentiation matrices

3) pipeVel: generates the mean velocity profile used in the resolvent operator.  This needs the log-interpolated velocity profiles in allProfilesLogInterp.txt  (data for Re > 75000 from Princeton superpipe experiments, and for Re = 5000-44000 from Wu and Moin DNS).

4) pipeOperators: generates the linear Navier-Stokes operator and the mass matrix

5) pipeBC: imposes the appropriate no-slip and no-through-flow boundary conditions at the wall.  Optionally, it can also be used to impose boundary conditions corresponding to amplitude- and phase-varying opposition control

6) chebdif - Written by J.A.C. Weideman and S.C. Reddy, generates the appropriate differencing matrices (called in pipeCoords).

7) clencurt - Written by Lloyd Trefethen, generates integration weights (called in pipeCoords)

8) fourier2physical - utility that generates the physical velocity field from the resolvent modes

9) pipeGreenFn - utility that generates the Green's function for the pressure Poisson equation in cylindrical coordinates

10) pipeUGT - utility that calculates the fourier-transformed velocity gradient tensor for the pipes (i.e., not in physical space)

11) pipePacketSwirl - utility that calculates the swirl field for a packet of resolvent modes

-----
Developed by Mitul Luhar*, Beverley McKeon*, and Ati Sharma+. 
Written and maintained by Mitul Luhar (mluhar@caltech.edu, mluhar@cantab.net)

* Graduate Aerospace Laboratories, California Institute of Technology
+ Engineering and the Environment, University of Southampton

SEE ATTACHED LICENSING AGREEMENT
