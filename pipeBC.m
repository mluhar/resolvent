function [H] = pipeBC(LHS,RHS,yP,N,varargin)
% Function to impose boundary conditions for Navier-Stokes Resolvent
% Written by Mitul Luhar on 02/06/2013

% After Fourier decomposition, the Navier-Stokes eqns are:
% (-1i*om*M - L)x = Mf , or LHS x = RHS f
% where x = [u;v;w;p] and f = [fu;fv;fw;-]

% varargin specifies the boundary condition
% varargin = {} - no slip
% varargin = {yPD,AD}: opposition control, detection at yPD (plus units)
% and amplitude AD, such that v(yPD) = -AD*v(0)

% yP: coordinates (in plus units)

% Implement boundary conditions by changing rows 1, N+1, 2*N+1 in LHS/RHS. 
% These correspond to the x,r,theta momentum balance at r = 1 (y = 0). 

% Set rows [1,N+1,2*N+1] in LHS = 0.  To be modified below to reflect BC
LHS(1:N:3*N,:) = 0;
% Set rows [1,N+1,2*N+1] in RHS = 0.
RHS(1:N:3*N,:) = 0;

% No slip for u and w
LHS(1,1) = 1;
% Row 1 of mom. balance now represents: u(y=0) = 0
LHS(2*N+1,2*N+1)=1;
% Row 2*N+1 of mom. balance now represents: w(y=0) = 0

if(size(varargin,2)>1)
    % Opposition control
    yPD = varargin{1};
    AD  = varargin{2};
    ci = find(yP>yPD,1);
    LHS(N+1,N+1) = 1;
    LHS(N+1,N+ci)= AD*1;
    % Row N+1 of mom. balance now represents: v(y=0) + AD*v(y=yD) = 0
else
    % No slip
    LHS(N+1,N+1) = 1;
    % Row N+1 of mom. balance now represents: v(y=0) = 0
end

% Compute resolvent which reflects these BCs
H = LHS\RHS;
end