% Calculate the velocity gradient tensor
% Note that the tensor is still Fourier transformed. 
function [UGT]=pipeUGT(u,k,n,r,D1E,D1O)
% u(r) is the fourier-transforced velocity field, u = [u;v;w];
% k is the axial wave number
% n is the azimuthal wave number
% D1E (D1O) are the even (odd) difference matrices corresponding to u (v,w)

% UGT is the velocity gradient tensor

% Written by Mitul Luhar on 03/29/2013

%Indices corresponding to velocity component
N  = length(r);
ax = 1:N; 
ra = (N+1):(2*N);
az = (2*N+1):(3*N);

%Velocities 
uu = u(ax); %axial velocity
vv = u(ra); %radial velocity
ww = u(az); %azimuthal velocity

%The last N columns of L represent the gradient operator on pressure. We
%will use them here to calculate the velocity gradient tensor
%Gradient
DX = diag(1i*k*ones(N,1));
DT = diag(1i*n./r);

%3x3 Velocity gradient at N points corresponding to ri
UGT = zeros(N,3,3);
UGT(:,1,1) = DX*uu; UGT(:,1,2) = D1E*uu; UGT(:,1,3) = DT*uu;
UGT(:,2,1) = DX*vv; UGT(:,2,2) = D1O*vv; UGT(:,2,3) = DT*vv;
UGT(:,3,1) = DX*ww; UGT(:,3,2) = D1O*ww; UGT(:,3,3) = DT*ww;
return