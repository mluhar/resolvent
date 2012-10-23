function [P,U,V,W] = multiModePU(K,x,theta,N)
%This function calculates the velocity, pressure, and reynolds stress
%fields for multiple modes superposed together

%Ni is the number of points in R the fields are required at
%K(j).a  : amplitude of mode j
%K(j).u  : Fourier-transformed velocities corresponding to mode j
%K(j).p  : Fourier-transformed pressure for mode j
%K(j).k  : (k,n) are the axial and azimuthal wave numbers

% Indices
ax = 1:N;
ra = N+1:2*N;
az = 2*N+1:3*N;

%Calculate velocity gradient tensors
U = zeros(N,length(theta),length(x));
V = zeros(N,length(theta),length(x));
W = zeros(N,length(theta),length(x));
P = zeros(N,length(theta),length(x));

for ji = 1:length(K)
    U = U + fourier2physical(K(ji).a*K(ji).u(ax), K(ji).k(1) , x, K(ji).k(2), theta);
    V = V + fourier2physical(K(ji).a*K(ji).u(ra), K(ji).k(1) , x, K(ji).k(2), theta);
    W = W + fourier2physical(K(ji).a*K(ji).u(az), K(ji).k(1) , x, K(ji).k(2), theta);
    P = P + fourier2physical(K(ji).a*K(ji).p,     K(ji).k(1) , x, K(ji).k(2), theta);
end

end