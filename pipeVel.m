% Code to calculate mean velocity profile for given Reynolds number, Re
% Re defined based on bulk-averaged velocity, and pipe diameter.
% Written by Mitul Luhar, 02/06/2013

% Re>= 75000 Data from Princeton Superpipe, McKeon et al 2004,
% http://gasdyn.princeton.edu/data/e248/mckeon_data.html
% Re = 5000 - 4400 Data from DNS, Wu and Moin 2008, JFM

% Interpolation by A Sharma, 2009

function [U0, yP, UP] = pipeVel(Re, y)
% linear interpolation on y, yP
% log interpolation on Re
% U0: profile rescaled to match laminar (i.e. bulk average = 0.5)
% y = (1-r): distance from the wall, normalized by pipe radius
% UP, yP: velocity and y in plus units

vels = load('allProfilesLogInterp.txt');
% 1st column: log10(Re) 
% 2nd column: y
% 3rd column: yP
% 4th column: UP
% 5th column: U0 (normalized based on 2x bulk-averaged velocity)

%Calculate all gridded Re values
availableRe = 10.^(vels(vels(:,2)==0,1));
%Calculate minimum % difference
minDif = min(abs(Re-availableRe)/Re);
%If difference > 10% display warning about potential interpolation
if(minDif>0.1)
    warning('Warning: No DNS or Superpipe data available within +/-10% of requested Re value.')
    warning('Interpolating based on available data, which could result in noisy profiles.') 
end

y = [ones(length(y),1)*log10(Re) y];
yP = griddatan(vels(:,1:2), vels(:,3), y);
UP = griddatan(vels(:,1:2), vels(:,4), y);
U0 = griddatan(vels(:,1:2), vels(:,5), y);
return
