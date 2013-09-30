function [L,M] = pipeOperators(Re,k,n,r,D1O,D1E,D2O,D2E,U0,dU0)
% Function to calculate the linear operator for the Navier-Stokes Eqn
% M(dx/dt) = Lx + Mf
% x = [u;v;w;p]
% f = nonlinear forcing
% Re: Reynolds number based on bulk-avg velocity and diameter
% k,n: axial/azimuthal wave number
% D1O,D1E,D2O,D2E: finite difference matrices

% Written by Mitul Luhar on 02/11/2013

% A few basic matrices to improve readability
I   = eye(length(r));
Z   = zeros(length(r));
Rm1 = diag(r.^-1);
Rm2 = diag(r.^-2);
ikU0 = 1i*k*diag(U0);

% Laplacian components
AE   = -(k^2)*I - (n^2)*Rm2 + Rm1*D1E + D2E;
AO   = -(k^2)*I - (n^2)*Rm2 + Rm1*D1O + D2O;
B    = -2*1i*n*Rm2;

% Change between centerline and bulk scaling
Ret = 1*Re; 
% NOTE: There may be a factor of 2 involved here

% The last column is the pressure gradient, last row is continuity
L1 = [-ikU0+AE/Ret, -dU0              , Z                 ,   -1i*k*I]; %u (ax.)
L2 = [Z           , -ikU0+(AO-Rm2)/Ret, B/Ret             ,      -D1E]; %v (rad.)
L3 = [Z           , -B/Ret            , -ikU0+(AO-Rm2)/Ret, -1i*n*Rm1]; %w (az.)
L4 = [1i*k*I      , Rm1+D1O           , 1i*n*Rm1          ,         Z]; %contin.
L  = [L1; L2; L3; L4];

% Mass matrix M with zeros representing continuity
M = [I Z Z Z; Z I Z Z; Z Z I Z; Z Z Z Z];

end

