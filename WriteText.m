function [] = WriteText(OutDir, pCa, dt, Binder, Steps, Means, Vars, IndexThalf)

TotalRuns = size(Means, 2);
NSTEPS = size(Steps, 2);

t = (dt*(1:NSTEPS))-dt;

% Package Up the Averaged (Index) information into a vector.
% Also, flip to column vectors, so we can write it to a outfile
OutIndex=zeros(TotalRuns, 20);
%% 24 rows for the first half results
OutIndex(:, 1)=(1:TotalRuns)';
OutIndex(:, 2)=pCa*ones(TotalRuns, 1);
%% Concat the SS mean data, followed by Variance
OutIndex(:, 3)=Means(1,:)';
OutIndex(:, 4)=Vars(1,:)';
OutIndex(:, 5)=Means(2,:)';
OutIndex(:, 6)=Vars(2,:)';
OutIndex(:, 7)=Means(3,:)';
OutIndex(:, 8)=Vars(3,:)';
OutIndex(:, 9)=Means(4,:)';
OutIndex(:, 10)=Vars(4,:)';
OutIndex(:, 11)=Means(5,:)';
OutIndex(:, 12)=Vars(5,:)';
OutIndex(:, 13)=Means(6,:)';
OutIndex(:, 14)=Vars(6,:)';
OutIndex(:, 15)=Means(7,:)';
OutIndex(:, 16)=Vars(7,:)';
OutIndex(:, 17)=Means(8,:)';
OutIndex(:, 18)=Vars(8,:)';
OutIndex(:, 19)=Means(9,:)';
OutIndex(:, 20)=Vars(9,:)';


%% NOW Concat the half times
OutHalfTimes=zeros(TotalRuns, 11);
OutHalfTimes(:, 1)=(1:TotalRuns)';
OutHalfTimes(:, 2)=pCa*ones(TotalRuns, 1);
OutHalfTimes(:, 3)=IndexThalf(1,:)';
OutHalfTimes(:, 4)=IndexThalf(2,:)';
OutHalfTimes(:, 5)=IndexThalf(3,:)';
OutHalfTimes(:, 6)=IndexThalf(4,:)';
OutHalfTimes(:, 7)=IndexThalf(5,:)';
OutHalfTimes(:, 8)=IndexThalf(6,:)';
OutHalfTimes(:, 9)=IndexThalf(7,:)';
OutHalfTimes(:, 10)=IndexThalf(8,:)';
OutHalfTimes(:, 11)=IndexThalf(9,:)';
%%%%%%%%%%%%%%%%%%%% END INDEX (All Runs) OUTPUT PARAMS %%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%BEGIN TIME SERIES OUPUT TO FILE %%%%%%%%%%%%%%%%%%
% Where every row here is a data input.

TimeSeries_DataOut=[t; Steps(1,:); Steps(2,:); Steps(8,:); Steps(3,:)+Steps(4,:); Steps(3,:); Steps(4,:); Steps(5,:); Steps(6,:); Steps(7,:)]';
% Build output filename
OutFile=sprintf('%sTimeSeriesAvg_pCa_%s.txt', OutDir, num2str(pCa, '%3.2f'));
% Open File
fid=fopen(OutFile, 'wt' );	%open outfile--tab delimited text
%Create the output file column header:
fprintf(fid, 'Time (sec)\tThick F pN\tThin F (pN\tATP per dt\tFrac Bound\tFract. XB1\tFract. XB2\tActins Ca0\tActins Ca1\tActins Ca2\n');	%header for outfile
%Create the output file format string, based on the size and data type.
%Force the default precision is 10 places wide, with 6 decimal places,
% a floating point number, with a tab delimiter.
% This comes out as %10.6f\t, in the Format String.
FormatString=[];
[~, ColOut]=size(TimeSeries_DataOut);
for i=1:ColOut-1 %for all but last
    FormatString=[FormatString, '%10.6f\t'];
end
FormatString=[FormatString, '%10.6f\n'];
% Write Data File and close file
fprintf(fid, FormatString, TimeSeries_DataOut');
fclose(fid);     %close the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%% TIME SERIES OUPUT TO FILE %%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%BEGIN INDEXED SS OUTPUT TO FILE %%%%%%%%%%%%%%%%%%
% Build output filename, open file, write header, write data, close file
% Just as above:
OutFile=sprintf('%sSSData_pCa_%s.txt', OutDir, num2str(pCa, '%3.2f'));
fid=fopen(OutFile,'w');	%open outfile--tab delimited text
%Create the output file column header:
fprintf(fid, 'Run Index\tpCa Value\tThickF(pN)\tVARThick F\tThinF (pN)\tVAR Thin F\tFract. XB1\tVAR F XB1\tFract. XB2\tVAR F XB2\tActins Ca0\tVAR ActCa0\tActins Ca1\tVAR ActCa1\tActins Ca2\tVAR ActCa2\tATP per dt\tVAR ATP dt\tFrct Bound\tVAR XB Bnd\n');	%header for outfile
FormatString=[];
[~, ColOut]=size(OutIndex);
for i=1:ColOut-1 %for all but last
    FormatString=[FormatString, '%10.6f\t'];
end
FormatString=[FormatString, '%10.6f\n'];
% Write Data File and close file
fprintf(fid, FormatString, OutIndex');
fclose(fid);     %close the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%% INDEXED SS OUTPUT TO FILE %%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%BEGIN INDEXED Time Half OUTPUT TO FILE %%%%%%%%%%%%%%%%%%
% Build output filename, open file, write header, write data, close file
% Just as above:
OutFile=sprintf('%sHalfTimeData_pCa_%s.txt', OutDir, num2str(pCa, '%3.2f'));
fid=fopen(OutFile,'w');	%open outfile--tab delimited text
%Create the output file column header:
fprintf(fid, 'Run Index\tpCa Value\tThickF(pN)\tThinF (pN)\tFract. XB1\tFract. XB2\tActins Ca0\tActins Ca1\tActins Ca2\tATP per dt\tFrct Bound\n');	%header for outfile
FormatString=[];
[~, ColOut]=size(OutHalfTimes);
for i=1:ColOut-1 %for all but last
    FormatString=[FormatString, '%10.6f\t'];
end
FormatString=[FormatString, '%10.6f\n'];
% Write Data File and close file
fprintf(fid, FormatString, OutHalfTimes');
fclose(fid);     %close the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%% INDEXED Time Half OUTPUT TO FILE %%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%BEGIN Binding Events OUTPUT TO FILE %%%%%%%%%%%%%%%%%%
% Build output filename, open file, write header, write data, close file
% Just as above:
% Only do so if we had binding events
if not(isempty(Binder))
    OutFile=sprintf('%sXBBindingData_pCa_%s.txt', OutDir, num2str(pCa, '%3.2f'));
    fid=fopen(OutFile,'w');	%open outfile--tab delimited text
    %Create the output file column header:
    fprintf(fid, 'Run Index\tStartNSTEP\tTot_dt_Bnd\tXB1_dt_Bnd\tXB2_dt_Bnd\tXB1--->XB2\tXB2--->XB1\tXB2--->XB0\tAVG Xthick\tAVG X_thin\tAvgXBstrain\tAvgXB1strn\tAvgXB2strn\tThick Node\tThin FNode\tkxscaler  \treach(xb0)\tCompleted \tVAR Xthick\tVAR X_thin\tVARXBstrain\tVARXB1strn\tVARXB2strn\n');	%header for outfile
    FormatString=[];
    ColOut=size(Binder, 2); % not sure if valid
    for i=1:ColOut-1 %for all but last
        FormatString=[FormatString, '%10.6f\t'];
    end
    FormatString=[FormatString, '%10.6f\n'];
    % Write Data File and close file
    fprintf(fid, FormatString, Binder');
    fclose(fid);     %close the file
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Binding Events OUTPUT TO FILE %%%%%%%%%%%%%%%%%%