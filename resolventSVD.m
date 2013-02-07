function [r,su,ss,sv,U0,yP,UP,dU0,dr] = resolventSVD(Re,k,n,om,N,nsvd,varargin)
% This code computes the singular value decomposition of  the Navier-Stokes
% resolvent for turbulent pipe flow
% Written by Mitul Luhar on 02/06/2013

% After Fourier decomposition, u(r)exp(i*[k*x+n*theta-om*t]), the 
% Navier-Stokes equations for turbulent pipe flow are: 

% (-1i*om*M)x = Lx + Mf  --> x = (-1i*om*M - L)\M f

% Here x = [u;v;w;p], denotes the velocities and pressure
% L is a linear operator and f are the nonlinear 'forcing' terms. 
% M is a modified mass matrix.  The resolvent is: (-1i*om*M - L)\M

% Re: Reynolds number based on pipe radius and 2x bulk-avg. velocity
% n : azimuthal wave number (has to be an integer!)
% k : axial wave number (k > 0)
% om: radian frequency scaled such that phase speed c = om/k is normalized 
% based on centerline velocity
% N : number of grid points in r:(0,1]
% nsvd: number of singular modes to compute

% varargin specifies the boundary condition
% varargin = {} - no slip
% varargin = {yPD,AD}: opposition control with detection at yPD from the
% wall (in plus units), with amplitude AD, such that v(yPD) = -AD*v(0)

%% Coordinate system based on Chebyshev collocation
[r,dr,D1E,D1O,D2E,D2O] = pipeCoords(n,N);
% r: radial coordinate (0,1]
% dr: integration weights
% D1E,D1O:  Even and odd first difference matrices 
% D2E,D2O:  Even and odd second difference matrices

% NOTE: definition of even and odd changes with n   
% D1E/D2E correspond to the behavior of the axial velocity
% D1O/D2O correspond to the behavior of the radial/azimuthal velocity

%% Load velocity profile. 
[U0,yP,UP]  = pipeVel(Re,1-r);
% U0: profile rescaled to match laminar (i.e. bulk average = 0.5)
% y = (1-r): distance from the wall, normalized by pipe radius
% UP, yP: velocity and y in plus units

% Calculate mean shear
D1R = mod(n+1,2)*D1E + mod(n,2)*D1O; % This is the 'true' even matrix
dU0 = diag(smooth(D1R*U0)); % Shear must be smoothed as data can be noisy
% NOTE: Need to create a better solution than 'smooth'

%% Calculate linear operator, L, and mass matrix, M
% M(dx/dt) = Lx + Mf 
[L, M] = pipeOperators(Re,k,n,r,D1O,D1E,D2O,D2E,U0,dU0);

%% Computer Resolvent
% Scale omega such that c = om/k scales with centerline velocity
omt = om*max(U0);

% The governing equation reads: (-i*om*M - L)x = Mf
LHS = -1i*omt*M-L;
RHS = M; 

% Impose boundary conditions (BC) and calculate resolvent, H = (LHS/RHS)
if(isempty(varargin))
    % No slip
    H = pipeBC(LHS,RHS,yP,N);
else
    % Opposition Control
    H = pipeBC(LHS,RHS,yP,N,varargin{1},varargin{2});
end

%% Scale Resolvent
% Currently, the resolvent is scaled to yield unit l2 norm for the singular
% forcing and response modes.  Alternative norms may be considered.
IW   = sqrtm(diag(r.*dr));
iIW  = eye(length(r))/IW;
Z = zeros(length(r));
sqW  = [ IW Z Z Z; Z  IW Z Z; Z Z  IW Z; Z Z Z Z];
isqW = [iIW Z Z Z; Z iIW Z Z; Z Z iIW Z; Z Z Z Z];

% Weighted resolvent
HW = sqW*H*isqW;

%% Singular value decomposition
[suW ssW svW] = svds(HW,nsvd); 
su = isqW*suW; 
sv = isqW*svW;
ss = diag(ssW);
% ss: singular values
% su: singular response (velocity) modes
% sv: singular forcing modes

% set phase near critical layer to be zero
phase_shift = -1i*angle(su(find((omt/k)>U0,1),:));
sv = sv*diag(exp(phase_shift));

% Because of the l2 norm used to scale the resolvent, we do not have any
% pressure data.  Calculate pressure modes using the un-scaled resolvent.
su = H*sv;
su = su*diag(1./ss);
