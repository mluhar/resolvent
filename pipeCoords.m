function [ r, dr, D1E, D1O, D2E, D2O, IW ] = pipeCoords( n, N )
% Function to set up the coordinate system, differentiation and integration
% matrices for the pipe geometry.

% Created by Mitul Luhar, 10/22/2012

% r: radial coordinate (0,1]
% dr: integration weights
% D1E,D1O:  Even and odd first difference matrices
% D2E,D2O:  Even and odd second difference matrices
% IW is the integration weights matrix

[x,DM] = chebdif(2*N,2); % x:[-1 1]
half   = 1:N;            % indices corresponding to x:(0,1]
r  = x(half);            % r:(0 1]
D1 = DM(:,:,1);          % First differential
D2 = DM(:,:,2);          % Second differential

%Split D1 and D2 into odd and even components (cf. Meseguer Trefethen 2003)
%For even (odd) n, u is even (odd) over r:[-1 1] 
%For even (odd) n, v and w are odd (even) over r:[-1 1]
s = (-1)^mod(n,2);
D1E = D1(half,half) + s*D1(half,2*N+1-half);
D1O = D1(half,half) - s*D1(half,2*N+1-half);
D2E = D2(half,half) + s*D2(half,2*N+1-half);
D2O = D2(half,half) - s*D2(half,2*N+1-half);

% Find integration weights
[~,dr] = clencurt(2*N-1);
dr   = (dr(half))';
IW   = sqrtm(diag(r.*dr));

end

