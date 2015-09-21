function [TotalRuns] = SimRuns(DataPoints, DataSpan, NSTEPS, ~)

% The ~ is maxt.

% Now we need to calculate how many runs we want.
% Say that we are going to want to average so many data points say:
%DataPoints=40000; Were used in the two filament model.
% Say that our tessellation now has 24 times as many nodes interacting
% this number would then be 40k/24= 1.6k.  But we think we need at least
% 5 runs.so make DataPoints:
%TotalRuns=ceil(DataPoints./(DataSpan*(NSTEPS/2))); %% Used for longer
TotalRuns=ceil(DataPoints./(DataSpan*(NSTEPS/1)));
% For testing on the longer duration runs at the mid-pCa levels
% make them run 10 times.
%      if maxt==7
%          TotalRuns=15;
%      end