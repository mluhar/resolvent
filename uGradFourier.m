% Calculate the velocity gradient tensor
% Note that the tensor is still Fourier transformed. 

function [UGT]=uGradFourier(U,L,N,Ni)
% U(r) is the fourier-transforced velocity field, U = [u;v;w];
% N is the number of points in R
% Ni is the number of points in R at which the velocity gradient is needed
% L is the linear operator from the Resolvent analysis.

%UGT is the velocity gradient tensor of size [3,3,Ni] at the Ni points in R

%Indices corresponding to velocity component
ax = 1:N; 
ra = (N+1):(2*N);
az = (2*N+1):(3*N);
pr = (3*N+1):(4*N);

%Velocities 
u = U(ax); %axial velocity
v = U(ra); %radial velocity
w = U(az); %azimuthal velocity

%The last N columns of L represent the gradient operator on pressure. We
%will use them here to calculate the velocity gradient tensor
%Gradient
ddx = -L(ax,pr);
ddr = -L(ra,pr);
ddt = -L(az,pr);

%3x3 Velocity gradient at N points corresponding to ri
UGT = zeros(3,3,Ni);
for ri = 1:Ni
    UGT(:,:,ri) = [ddx(ri,:)*u ddr(ri,:)*u ddt(ri,:)*u;
        ddx(ri,:)*v ddr(ri,:)*v ddt(ri,:)*v;
        ddx(ri,:)*w ddr(ri,:)*w ddt(ri,:)*w];
end

return
