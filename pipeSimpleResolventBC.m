%% Collocation code based on full Navier-Stokes Equations
% Written by Mitul Luhar
% Re: Reynolds number
% n : azimuthal wave number (has to be an integer!)
% k : axial wave number (k > 0)
% om: frequency
% Fourier decomposition: u(r)exp(i*[k*x+n*theta-om*t])
% N : number of grid points in r:(0,1]
function [r,su,ss,sv,sp,L,H,U0,yP,UP,dU0,dr] = pipeSimpleResolventBC(Re,k,n,om,N,yPc)
% The Navier-Stokes equations can be written as: M(dx/dt) = Lx + Mf
% Here x = [u;v;w;p], L is a linear operator and f is the forcing
% comprising the nonlinear terms. M is a modified mass matrix

% Eventually, each cell will be modularized.

%% Pseudospectral Chebyshev differentiation matrices D1 and D2 at points x
[x,DM] = chebdif(2*N,2); % x:[-1 1]
half   = 1:N;            % indices corresponding to x:(0,1]
r  = x(half);            % r:(0 1]
D1 = DM(:,:,1);          % First differential
D2 = DM(:,:,2);          % Second differential

%Split D1 and D2 into odd and even components (cf. Meseguer Trefethen 2003)
%For even (odd) n, u is even (odd) over r:[-1 1] 
%For even (odd) n, v and w are odd (even) over r:[-1 1]
s = (-1)^mod(n,2);
D1E = D1(half,half) + s*D1(half,2*N+1-half);
D1O = D1(half,half) - s*D1(half,2*N+1-half);
D2E = D2(half,half) + s*D2(half,2*N+1-half);
D2O = D2(half,half) - s*D2(half,2*N+1-half);

% A few basic matrices
I   = eye(length(half));
Z   = zeros(length(half));
Rm1 = diag(r.^-1);
Rm2 = diag(r.^-2);

%% Velocity profile. Eventually create module with switch for laminar/turb.
% % Laminar
% U0  = (I-R.^2); dU0 = -2*R;
% Turbulent
[U0,yP,UP]  = pipeVel(Re,1-r);
dU0 = D1E*U0;
dU0 = diag(smooth(dU0));


%% Create important block components
ikU0 = 1i*k*diag(U0);
AE   = -(k^2)*I - (n^2)*Rm2 + Rm1*D1E + D2E;
AO   = -(k^2)*I - (n^2)*Rm2 + Rm1*D1O + D2O;
B    = -2*1i*n*Rm2;

% Block matrix L representing linearized NS equations
% The last column is the pressure gradient, last row is continuity
Ret = 2*Re; % Change between centerline and bulk scaling
L1 = [-ikU0+AE/Ret, -dU0              , Z                 ,   -1i*k*I]; %u (ax.)
L2 = [Z           , -ikU0+(AO-Rm2)/Ret, B/Ret             ,      -D1E]; %v (rad.)
L3 = [Z           , -B/Ret            , -ikU0+(AO-Rm2)/Ret, -1i*n*Rm1]; %w (az.)
L4 = [1i*k*I      , Rm1+D1O           , 1i*n*Rm1          ,         Z]; %contin.
L  = [L1; L2; L3; L4];

% Mass matrix M with zeros representing continuity
M = [I Z Z Z; Z I Z Z; Z Z I Z; Z Z Z Z];

%% Impose Boundary Conditions 
% Dirichlet at the wall (r = 1)
% Get rid of all the columns corresponding to u,v,w at r = 1;

% Scale omega with centerline velocity
omt = om*max(U0);

% Basic matrices for Resolvent
LHS = -1i*omt*M-L;
RHS = M; 

% Implement boundary conditions
LHS([1 N+1 2*N+1],:) = 0;
% u - no slip
LHS(1,1) = 1;
% v - opposition control
ci = find(yP>yPc,1);

LHS(N+1,N+1) =-1;
LHS(N+1,N+ci)= 1;
% w - no slip
LHS(2*N+1,2*N+1)=1;

% Make sure the forcing modes have the same behavior
RHS(1:N:end,:) = 0;

% Create modified resolvent matrix H
H = LHS\RHS;

%% Choose appropriate norm for the the SVD
% Integral weights for norm based on the Clenshaw-Curtis quadrature
[~,dr] = clencurt(2*N-1);
dr   = (dr(half))';
IW   = sqrtm(diag(r.*dr));
iIW  = eye(length(half))/IW;
% Do not include pressure in the norm
sqW  = [IW Z Z Z; Z IW Z Z; Z Z IW Z; Z Z Z Z];
isqW = [iIW Z Z Z; Z iIW Z Z; Z Z iIW Z; Z Z Z Z];

% Weighted resolvent
HW = sqW*H*isqW;

%% Singular value decomposition
% nsvd = 4*N-3;
nsvd = 5;
[suW ssW svW] = svds(HW,nsvd);
su = isqW*suW; 
sv = isqW*svW;
ss = diag(ssW);

% set phase at first r-location to be 0
phase_shift = -1i*angle(su(1,:));
su = su*diag(exp(phase_shift));
sv = sv*diag(exp(phase_shift));

% calculate pressure
sp = H*sv;
sp = sp*diag(1./ss);
sp = sp(end-N+1:end,:);