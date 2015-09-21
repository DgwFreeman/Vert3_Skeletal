function [OutTotal]=GatherData_v1(OutTotal, DataDir, pCa)

%% Gather the data that was written to text files, into the OutTotal cell
%% array, then pass the cell array back to the function.

%% Our sizing of OutTotal, should be correct, it just might not contain the
%% proper information from the distribution to the VACC cluster and back
%% again.  This should bypass the problem.

% Initialize OutTotal as a cell array to hold the
% Col 1) pCa Value
% Col 2) Average Time Series Output at a pCa value
% Col 3) SS Individual Trial Output
% Col 4) Time Half from Individual Trials

% The number of rows for OutTotal will be equal to the number of pCa values
% run for the simulation.
for i=1:length(pCa)
    OutTotal{i, 1}=pCa(i,1);

    %% First bring in the time series data; strip the file name exactly as
    %% used to write the text file, at the end of Run_v1.m
    InFile=sprintf('%sTimeSeriesAvg_pCa_%s.txt', DataDir, num2str(pCa(i,1), '%3.2f'));
    %% Import the time series data, strip the numeric data and load it into
    %% the cell array (col 2)
    InData=importdata(InFile);
    OutTotal{i, 2}=InData.data;

    % Now repeat the same for the SSData at col 3 and the HalfTimeData at
    % col 4
    InFile=sprintf('%sSSData_pCa_%s.txt', DataDir, num2str(pCa(i,1), '%3.2f'));
    InData=importdata(InFile);
    OutTotal{i,3}=InData.data;

    InFile=sprintf('%sHalfTimeData_pCa_%s.txt', DataDir, num2str(pCa(i,1), '%3.2f'));
    InData=importdata(InFile);
    OutTotal{i,4}=InData.data;
end



