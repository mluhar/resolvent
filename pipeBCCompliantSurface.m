% Function to impose boundary conditions for Navier-Stokes Resolvent
% for a compliant surface
% Ati Sharma 2013-02-15
%
% effDamping is the damping coefficient divided by the mass;
% effModulus is the spring constant divided by mass;
%
% TODO modify pipeBC*.m so separate file for each BC case
% TODO modify pipeBC*.m so that LHS, RHS formed internally (i.e. take L,M as args)
% TODO or even call BC modifications from within pipeOperators since is really a modification of these

function [H] = pipeBCCompliantSurface(L,M,Re,omt,effModulus,effDamping)

N=length(L)/4;

% r-deriv at wall of mean vel
% dUdr0 = dU0(1); % if doing from known mean vel profile
% or
% dU+/dy+=1, dU = (dU/dy) (U+/U0)/R+ = -U+/(R+.U0) dU/dr
[U0,RP,U0P]  = pipeVel(Re,1); % y+ = R+ at CL.
dUdr0 = -RP*U0/U0P;

% After Fourier decomposition, the Navier-Stokes eqns are:
% (-1i*omega*M - L)x = Mf , or LHS x = RHS f
% where x = [u;v;w;p] and f = [fu;fv;fw;-]
% we need to jimmy around with L here
% and return modified H.

%------ equations for surface -----%
% h(theta,z,t) is vertical displacement of surface (Fourier coeff of);
% vh is dh/dt;
% p is pressure at surface;
% note that effModulus can be complex to allow for hysteretic damping
% (positive imaginary part for damping);
% and, in principle, negative mass may be allowed.
%
% d^2/dt^2 h(t) + effDamping d/dt h(t) + effModulus h(t) = p(r=1,t)		[1]

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
% u(r=1,t) = -U0'(r=1) h(t)												[2]
% v(r=1,t) =      d/dt h(t)												[3]

%----- coupling of fluid with surface -----%
% [3] -> d/dt [2] gives first line
% [2], [3], d/dt[3] -> [1] gives second line
%
% -i*omt [ 1  0 ] [ u(r=1) ] = [0                         -U0'    ] [ u(r=1) ]  + [0] p(r=1)
%        [ 0  1 ] [ v(r=1) ]   [effModulus/U0'(r=1)   -effDamping ] [ v(r=1) ]  + [1]
%
% which relates p, u and v at the boundary.
%

%------ boundary conditions implementation -----%
% Implement boundary conditions by changing rows 1, N+1, 2*N+1 in LHS/RHS. 
% These correspond to the x,r,theta momentum balance at r = 1 (y = 0). 
%
%--- u(r=1) mass matrix and v(r=1) mass matrix elements remain the same.
%
%--- u(r=1) dynamic matrix
L(1,1) = 0;
L(1,N+1) = -dUdr0;
%--- v(r=1) dynamic matrix
L(N+1,1) = effModulus/dUdr0; % coupling to u due to oscillations of h
L(N+1,N+1) = -effDamping; % self-coupling due to damping
% forcing on v from pressure BC
L(N+1,3*N+1) = 1;
%--- w(r=1) dynamic matrix no slip
L(:,2*N+1) = 0;

%----- form resolvent -----%
% The governing equation reads: (-i*om*M - L)x = Mf
LHS = -i*omt*M - L;
RHS = M; 
H = LHS\RHS;

