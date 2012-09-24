%% Collocation code based on full Navier-Stokes Equations
% Written by Mitul Luhar
% Re: Reynolds number
% n : azimuthal wave number (has to be an integer!)
% k : axial wave number (k > 0)
% om: frequency
% Fourier decomposition: u(r)exp(i*[om*t-k*x-n*theta])
% N : number of grid points in r:(0,1]
function [r,su,ss,sv,DIV,L1,L2,L3] = pipeSimpleResolvent(Re,k,n,om,N)
% The Navier-Stokes equations can be written as: M(dx/dt) = Lx + Mf
% Here x = [u;v;w;p], L is a linear operator and f is the forcing
% comprising the nonlinear terms.

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
R   = diag(r);
Rm1 = diag(r.^-1);
Rm2 = diag(r.^-2);

%% Velocity profile. Eventually create module with switch for laminar/turb.
% % Laminar
% U0  = (I-R.^2); dU0 = -2*R;
% Turbulent
U0  = pipeVel(Re,1-r);
dU0 = D1E*U0;
U0  = diag(U0);
dU0 = diag(smooth(dU0));


%% Create important block components
ikU0 = 1i*k*U0;
AE   = -(k^2)*I - (n^2)*Rm2 + Rm1*D1E + D2E;
AO   = -(k^2)*I - (n^2)*Rm2 + Rm1*D1O + D2O;
B    = 2*1i*n*Rm2;

% Block matrix L representing linearized NS equations
% The last column is the pressure gradient, last row is continuity
Ret = 2*Re; % Change between centerline and bulk scaling
L1 = [ikU0+AE/Ret, -dU0             , Z                ,   1i*k*I]; %u (ax.)
L2 = [Z          , ikU0+(AO-Rm2)/Ret, B/Ret            ,     -D1E]; %v (rad.)
L3 = [Z          , -B/Ret           , ikU0+(AO-Rm2)/Ret, 1i*n*Rm1]; %w (az.)
L4 = [-1i*k*I    , Rm1+D1O          , -1i*n*Rm1        ,        Z]; %contin.
L  = [L1; L2; L3; L4];
DIV = L4;

% Mass matrix M with zeros representing continuity
M = [I Z Z Z; Z I Z Z; Z Z I Z; Z Z Z Z];

%% Impose Boundary Conditions 
% Dirichlet at the wall (r = 1)
% Get rid of all the columns corresponding to u,v,w at r = 1;
LHS = 1i*om*M-L;
LHS(:,1:N:3*N) = [];

RHS = M; 
RHS(:,1:N:3*N) = []; 

% At this point the problem is overspecificed (4*N eqns for 4*N-3 vars)
% The velocities at the wall are known, but we have retained the momentum-
% and continuity equations at the wall. 
% In theory, we should only need one of the three momentum equations at the
% wall to specify pressure. See e.g. Gresho and Sani, Int J Numer Methods
% Fluids, 7:1111 (1987). 
% Alternatively, we can retain all the conditions.  This will lead to
% MATLAB doing the inversion below (LHS\RHS) in a least-squares sense. 

% % Keep only axial stress boundary condition
% LHS([N+1 2*N+1 3*N+1],:) = []; 
% RHS([N+1 2*N+1 3*N+1],:) = []; 
% 
% % Keep only radial stress boundary condition
% LHS([1 2*N+1 3*N+1],:) = []; 
% RHS([1 2*N+1 3*N+1],:) = []; 
% 
% % Keep only azimuthal stress boundary condition
% LHS([1 N+1 3*N+1],:) = []; 
% RHS([1 N+1 3*N+1],:) = []; 
% 
% Keep only continuity at the wall
LHS([1 N+1 2*N+1],:) = []; 
RHS([1 N+1 2*N+1],:) = []; 

% Create modified resolvent matrix H
H = LHS\RHS;

%% Choose appropriate norm for the the SVD
% Integral weights for norm based on the Clenshaw-Curtis quadrature
[~,dr] = clencurt(2*N-1);
dr   = (dr(half))';
IW   = sqrtm(diag(r.*dr));
iIW  = eye(length(half))/IW;
% Include pressure in the norm
sqW  = [IW Z Z Z; Z IW Z Z; Z Z IW Z; Z Z Z 1e0*IW];
isqW = [iIW Z Z Z; Z iIW Z Z; Z Z iIW Z; Z Z Z 1e-0*iIW];

% Get rid of rows and columns corresponding to u,v,w at r = 1
sqW(:,1:N:3*N)=[];
sqW(1:N:3*N,:)=[];
isqW(:,1:N:3*N)=[];
isqW(1:N:3*N,:)=[];

% Weighted resolvent
HW = sqW*H*isqW;

%% Singular value decomposition
nsvd = 5;
[suW ssW svW] = svds(HW,nsvd);
su = isqW*suW; 
sv = isqW*svW;
ss = diag(ssW);

% Add the zero velocity and forcing at r = 1 back in
su = [zeros(1,nsvd); su(1:N-1,:); zeros(1,nsvd); su(N:2*N-2,:); zeros(1,nsvd); su(2*N-1:end,:)];
sv = [zeros(1,nsvd); sv(1:N-1,:); zeros(1,nsvd); sv(N:2*N-2,:); zeros(1,nsvd); sv(2*N-1:end,:)];