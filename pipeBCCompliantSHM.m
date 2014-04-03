function [H] = pipeBCCompliantSHM(LHS,RHS,dU0,omt,freqRatio,dampRatio,massRatio)
% Function to impose boundary conditions for Navier-Stokes Resolvent
% This accounts for a simple compliant surface model, which essentially
% leads to a harmonic-motion relationship between pressure and velocity
% Original: Ati Sharma 02/15/2013
% Modified: Mitul Luhar 03/02/2014

% After Fourier decomposition, the Navier-Stokes eqns are:
% (-1i*omt*M - L)x = Mf , or LHS x = RHS f
% where x = [u;v;w;p] and f = [fu;fv;fw;-]
% dU0 is the mean shear at the wall-required in the linearized BC
% omt is the frequency of oscillation

%----- Dynamic equations for wall-deformation ------%
% A simple spring-dashpot model leads to the following dynamic relation 
% between pressure, p, and wall deflection, h
% h'' + 2 [om0 dampRatio] h' + [om0^2] h = [1/massRatio] p(r=1)         [1]

% massRatio = [rho_w d_w] / [rho R] is the ratio of wall to fluid mass
% freqRatio = om0/omt is the ratio of natural and oscillatory frequency
% dampRatio: damping Ratio
% NOTE: om0 can be complex to account for hysteritic damping.

%----- Wall BC -----%
% roughly after Gaster, Grosch and Jackson 1994 JFM and others
% see also Woodcock et al 2012
% Essentially a Taylor expansion around r=1 for small h.
%
% u(r=1,t) + h(t) du/dr(r=1,t) + h(t) dU0(r=1) = 0         + O(h^2)
% v(r=1,t) + h(t) dv/dr(r=1,t)                 = d/dt h(t) + O(h^2)
% w(r=1,t) + h(t) dw/dr(r=1,t)                 = 0         + O(h^2)
% 
% The terms bilinear in h and u,v,w are NOT negligible
% But we will neglect them to keep in a linear setting.  This gives:
% u(r=1,t) = -dU0(r=1) h(t)												[2]
% v(r=1,t) = h'(t)                                                      [3]

%----- After a harmonic analysis, we get the following -----% 
% Combining [1] and [3]
% -i omt massRatio [1 - freqRatio^2 + 2 i freqRatio dampRatio] v = p    [4]
% Combining [2] and [3]
% -i omt u = -v dUO                                                     [5]

% Implement boundary conditions by changing rows 1, N+1, 2*N+1 in LHS/RHS. 
% These correspond to the x,r,theta momentum balance at r = 1 (y = 0). 

N = length(dU0);
% Set rows [1,N+1,2*N+1] in LHS = 0.  To be modified below to reflect BC
LHS(1:N:3*N,:) = 0;
% Set rows [1,N+1,2*N+1] in RHS = 0.
RHS(1:N:3*N,:) = 0;

% Row 1 of mom. balance: Equation [5]
LHS(1,1) = -1i*omt;
LHS(1,N+1) = dU0(1,1);
% Row N+1 of mom. balance: Equation [4]
LHS(N+1,N+1) = -1i*omt*massRatio*(1-freqRatio^2+2*1i*freqRatio*dampRatio);
LHS(N+1,3*N+1) = -1;
% Row 2*N+1 of mom. balance: w(r=1) = 0
LHS(2*N+1,2*N+1)=1; 

% Compute resolvent which reflects these BCs
H = LHS\RHS;
end