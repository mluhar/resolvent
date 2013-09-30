function [ R, RP, G ] = pipeGreenFn(k,n,r,rp)
% Function to calculate the Green Function for the the pressure poisson
%equation in a pipe
% Written by Mitul Luhar, 10/22/2012
% For derivation, see Luhar, Sharma & McKeon (2013), Wall pressure 
% fluctuations induced by coherent structures in turbulent pipe flow,
% in Eigth International Symposium on Turbulence and Shear Flow Phenomena.

% k: streamwise wave number
% n: azimuthal wave number
% r: radial coordinate
% rp: location of source

% Note that this Green function corresponds to the Neumann condition, dG/dr
% = 0 at r = 1.  i.e., it does not include the Stokes component of pressure
% (see e.g., Kim 1989, J. Fluid Mech.)

% One of the terms that crops up frequently
A = (-besselk(n-1,k)-besselk(n+1,k))/(besseli(n-1,k)+besseli(n+1,k));

% Set up grid and initialize
[R,RP] = meshgrid(r,rp);
G = zeros(size(R));

G(R>RP) = A*besseli(n,k*RP(R>RP)).*besseli(n,k*R(R>RP)) - besseli(n,k*RP(R>RP)).*besselk(n,k*R(R>RP));
G(R<=RP) = (A*besseli(n,k*RP(R<=RP)) - besselk(n,k*RP(R<=RP))).*besseli(n,k*R(R<=RP));
end

