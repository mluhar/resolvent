function [P] = pipePacketSwirlU(a0,UGT0,K,P,x,th,N)
%This function calculates the swirl, azimuthal vorticiity (azvor) and
%velocity fields (U) for a packet (P) of modes (K).

% x and th are the axial and azimuthal locations for the fourier-transform

% a0 and UGT0 are the mean flow amplitude and velocity gradient tensor
%K(ik).mode(1)  : axial wave number for mode ik
%K(ik).mode(2)  : azimuthal wave number for mode ik
%K(ik).su       : Fourier-transformed velocities corresponding to mode ik
%K(ik).UGT      : Velocity gradient tensor for mode ik

%P(ip).a        : amplitudes, a, of modes K in packet P.

% Written by Mitul Luhar on 03/29/2013

for ip = 1:length(P)
    % Initialize velocity field
    P(ip).U = zeros(4*N,length(th),length(x));
    % Initialize swirl and vorticity fields
    P(ip).sw = zeros(N,length(th),length(x));
    P(ip).az = zeros(N,length(th),length(x));
    % Include mean shear
    P(ip).sw0 = zeros(N,length(th),length(x));
    P(ip).az0 = zeros(N,length(th),length(x));
    % Initialize temporary array to store velocity gradient tensor
    tempUGT  = zeros(N,length(th),length(x),3,3);
    
    for ik = 1:length(K)
        % Only run calculations if amplitude > 0
        if(abs(P(ip).a(ik))>0)
            % Calculate fourier transformed velocity and pressure fields
            P(ip).U = P(ip).U + fourier2physical(P(ip).a(ik)*K(ik).su,K(ik).mode(1),x,K(ik).mode(2),th);
            % Store velocity gradient tensor at each location
            if(P(ip).calcswirl)
                for iu = 1:3
                    for ix = 1:3
                        tempUGT(:,:,:,iu,ix) = tempUGT(:,:,:,iu,ix) + fourier2physical(P(ip).a(ik)*K(ik).UGT(:,iu,ix),K(ik).mode(1),x,K(ik).mode(2),th);
                    end
                end
            end
        end
    end
    
    % Now calculate swirl and vorticity fields
    if(P(ip).calcswirl)
        for ix = 1:length(x)
            for it = 1:length(th)
                for ir = 1:N
                    %local velocity gradient tensor
                    lUGT = squeeze(tempUGT(ir,it,ix,:,:));
                    lUGT0 = squeeze(tempUGT(ir,it,ix,:,:))+squeeze(a0*UGT0(ir,:,:));
                    %local swirl and vorticity fields
                    P(ip).sw(ir,it,ix)  = max(imag(eig(lUGT)));
                    P(ip).sw0(ir,it,ix) = max(imag(eig(lUGT0)));
                    P(ip).az(ir,it,ix)  = lUGT(2,1)-lUGT(1,2);
                    P(ip).az0(ir,it,ix) = lUGT0(2,1)-lUGT0(1,2);
                end
            end
        end
    end
    
end

return