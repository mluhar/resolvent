% Data from Princeton Superpipe, McKeon et al 2004,
% http://gasdyn.princeton.edu/data/e248/mckeon_data.html
% Interpolation by A Sharma, 2009

function [U0] = pipeVel(Re, y)
% linear interpolation on y, yp; log interpolation on Re
% U0: profile rescaled to match laminar (area under curve)
% y = (1-r): distance from the wall, normalized by pipe radius

vels = load('pipeTurbProfiles.txt');
% 1st column: log10(Re) 
% 2nd column: y
% 6th column: U0

y = [ones(length(y),1)*log10(Re) y];
U0 = griddatan(vels(:,1:2), vels(:,6), y);

return