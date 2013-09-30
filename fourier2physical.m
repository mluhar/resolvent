% Function to convert Fourier modes to physical coordinates
% function [U]=fourier2physical(u(r),k,z,n,theta,[omega,t]);
function [U]=fourier2physical(u,k,x,n,theta,varargin)
% U is the fourier-transformed velocity field
% u is the input (complex) mode shape
% k,n : streamwise, azimuthal wavenumber
% x,theta: streamwise, azimuthal coordinates

% Optional can also include frequency (omega) and time (t) as inputs

if length(varargin)
	exp_omegat=exp(-1i*varargin{1}*varargin{2});
else
	exp_omegat=1;
end

for ix=1:length(x)
	for it=1:length(theta)
		U(:,it,ix)=u*exp(1i*k*x(ix)+1i*n*theta(it))*exp_omegat;
	end
end

U = squeeze(real(U));

return
