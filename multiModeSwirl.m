function [swirl,azvor] = multiModeSwirl(a0,UGT0,K,x,theta,N,Ni)
%This function calculates the swirl (lambda_ci) of the superposed modes
%defined by the structure K, at coordinates x and theta. 
%Also calculates the azimuthal vorticity, azvor

%a0 and UGT0 are the amplitude and velocity gradient tensor for the mean
%flow
%Ni is the number of points in R required for the swirl
%K(j).a  : amplitude of mode j
%K(j).u  : Fourier-transformed velocities corresponding to mode j
%K(j).L  : Linear NS operator
%K(j).k  : (k,n) are the axial and azimuthal wave numbers

%Calculate velocity gradient tensors in Fourier space
for ji = 1:length(K)
    tmp(ji).UGT = uGradFourier(K(ji).u, K(ji).L, N, Ni);
end

swirl = zeros(Ni,length(theta),length(x));
azvor = zeros(Ni,length(theta),length(x));

for xi = 1:length(x)
    for ti = 1:length(theta)
        for ri = 1:Ni
            totUGT = a0*UGT0(:,:,ri);
            for ji = 1:length(K)
                totUGT(:,1) = totUGT(:,1) + fourier2physical(K(ji).a*tmp(ji).UGT(:,1,ri),K(ji).k(1), x(xi), K(ji).k(2), theta(ti));
                totUGT(:,2) = totUGT(:,2) + fourier2physical(K(ji).a*tmp(ji).UGT(:,2,ri),K(ji).k(1), x(xi), K(ji).k(2), theta(ti));
                totUGT(:,3) = totUGT(:,3) + fourier2physical(K(ji).a*tmp(ji).UGT(:,3,ri),K(ji).k(1), x(xi), K(ji).k(2), theta(ti));         
            end
            
            swirl(ri,ti,xi) = max(imag(eig(totUGT)));
            azvor(ri,ti,xi) = totUGT(2,1)-totUGT(1,2);
        end
    end
end


end