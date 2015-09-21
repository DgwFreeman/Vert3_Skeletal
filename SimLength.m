function [maxt] = SimLength(Ca)

% Do some figuring about how long to run the simulation
% based on the Ca2+ value
if Ca<=10^-6.2
    maxt=6;
 elseif Ca<=10^-5.51
     maxt=6; %(Mids to 10?)
elseif Ca<=10^-5.5
    maxt=6;
elseif Ca<=10^-5
    maxt=4;
else
    maxt=2;
end