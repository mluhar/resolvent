% Function to impose boundary conditions for Navier-Stokes Resolvent
% for a compliant surface
% Ati Sharma 2013-02-15
%
% effDamping is the damping coefficient divided by the mass;
% effModulus is the spring constant divided by mass;
% u_tau is the friction velocity
% nu is the kinematic viscosity
%
% TODO modify pipeBC*.m so separate file for each BC case
% TODO modify pipeBC*.m so that LHS, RHS formed internally (i.e. take L,M as args)
% TODO or even call BC modifications from within pipeOperators since is really a modification of these

function [H] = pipeBCCompliantSurface(L,M,omt,u_tau,Re_tau,effModulus,effDamping)

N=length(L)/4;

% r-deriv at wall of mean vel
% dUdr0 = dU0(1); % if doing from known mean vel profile, but careful as can be noisy.
% Or, analytically,
% U+ = y+ very near wall, =>
% U/u_tau = (1-r) u_tau / nu, =>
% dU_dr = - u_tau^2 / nu = -u_tau * Re_tau, where Re_tau = u_tau R / nu
dUdr0 = -u_tau * Re_tau;

% since we need to jimmy around with L and M both, for LHS and RHS here
% and return modified H.

% After Fourier decomposition, the Navier-Stokes eqns are:
% (-1i*omega*M - L)x = Mf , or LHS x = RHS f
% where x = [u;v;w;p] and f = [fu;fv;fw;-]

%------ equations for surface -----%
% h(theta,z,t) is vertical displacement of surface (Fourier coeff of);
% vh is dh/dt;
% p is pressure at surface;
% note that effModulus can be complex to allow for hysteretic damping
% (positive imaginary part for damping);
% and, in principle, negative mass may be allowed.
%
% d^2/dt^2 h(t) + effDamping d/dt h(t) + effModulus h(t) = p(r=1,t)

%----- equations for small deformations at wall -----%
% roughly after Gaster, Grosch and Jackson 1994 JFM and others
% see also Woodcock et al 2012
% this finds the velocity field at r=1 *equivalent* to a small displacement h
% with no-slip at y=h.
%
% u(r=1,t) + h(t) u'(r=1,t) + h(t) U0'(r=1) = 0         + O(h^2)
% v(r=1,t) + h(t) v'(r=1,t)                 = d/dt h(t) + O(h^2)
% w(r=1,t) + h(t) w'(r=1,t)                 = 0         + O(h^2)
% 
% where ' indicates r-derivative.
%
% The terms bilinear in h and u,v,w are NOT negligible,
% and contribute to streaming, change the mean etc.
% But we will neglect them to keep in a linear setting.
% This gives:
%
% u(r=1,t) = -U0'(r=1) h(t)
% v(r=1,t) =      d/dt h(t)

%----- coupling of fluid with surface -----%
% substituting together, noting that v(r=1,t) = d/dt h(t), gives
%
% -i*om [-U0'(r=1)   0   - ] [ u(r=1) ]   =   [0                         1       ] [ u(r=1) ]  + [0] p(r=1)
%       [    0       1   - ] [ v(r=1) ]       [-effModulus/U0'(r=1)  -effDamping ] [ v(r=1) ]  + [1]
%
% which relates p, u and v at the boundary.
%

%------ boundary conditions implementation -----%
% Implement boundary conditions by changing rows 1, N+1, 2*N+1 in LHS/RHS. 
% These correspond to the x,r,theta momentum balance at r = 1 (y = 0). 

% set BC by leaving mass matrix entry and setting L components to 0
%
%--- u(r=1) mass matrix
M(1,1) = -dUdr0;
% v(r=1) mass matrix element remains the same.
%
%--- u(r=1) dynamic matrix
L(1,1) = 0;
L(1,N+1) = 1;
%
%--- v(r=1) dynamic matrix
L(N+1,1) = effModulus/dUdy0; % coupling to u due to oscillations of h
L(N+1,N+1) = -effDamping; % self-coupling due to damping
% forcing from pressure BC
L(N+1,3*N+1) = 1;
%
%--- w(r=1) dynamic matrix
L(:,2*N+1) = 0;

%----- form resolvent -----%
% The governing equation reads: (-i*om*M - L)x = Mf
LHS = -1i*omt*M-L;
RHS = M; 
H = LHS\RHS;
