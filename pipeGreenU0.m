function [ R, RP, G ] = pipeGreenU0(r,rp)
% Function to calculate the Green Function for Laplace's equation in polar
% coordinates
% Written by Mitul Luhar, 10/22/2012

% r: radial coordinate
% rp: location of source

% Note that this Green function corresponds to G = 0 at r = 1, and bounded
% G at r =0

% Set up grid and initialize
[R,RP] = meshgrid(r,rp);
G = zeros(size(R));

% Calculate Green's function
G(R>RP) = log(R(R>RP));
G(R<=RP)= log(RP(R<=RP));

end

