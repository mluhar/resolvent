function [H] = pipeBC(LHS,RHS,yP,N,varargin)
% Function to impose boundary conditions for Navier-Stokes Resolvent
% Written by Mitul Luhar on 02/06/2013

% After Fourier decomposition, the Navier-Stokes eqns are:
% (-1i*om*M - L)x = Mf , or LHS x = RHS f

% varargin specifies the boundary condition
% varargin = {} - no slip
% varargin = {yPD,AD}: opposition control, detection at yPD (plus units)
% and amplitude AD, such that v(yPD) = -AD*v(0)

% yP: coordinates (in plus units)

% Implement boundary conditions
LHS(1:N:3*N,:) = 0;
RHS(1:N:3*N,:) = 0;

% No slip for u and w
LHS(1,1) = 1;
LHS(2*N+1,2*N+1)=1;

if(size(varargin,2)>1)
    % Opposition control
    yPD = varargin{1};
    AD  = varargin{2};
    ci = find(yP>yPD,1);
    LHS(N+1,N+1) = 1;
    LHS(N+1,N+ci)= AD*1;
else
    % No slip for v
    LHS(N+1,N+1) = 1;
end

% Compute resolvent
H = LHS\RHS;
end