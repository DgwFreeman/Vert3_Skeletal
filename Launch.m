%% BCWT, August 17, 2010
%% Mod from Skeltal--to Drosophila 
%% Port distributed code back to a single machine, and directory.
%% Run isometric force-pCa curves, and automtically save the data.

%% User inputs, toggles:
% Set number
% Cooperativity Scaler
% Functional Tn density
% Functional XB density
% Tm Span
% pCa range with maximum time wanted to run the simulation at the pCa
% Desired number of runs per pCa

%% Simulation Outputs:
% Code and output will be written to an associated
% Data directory:
% <CurrentDir\DataFiles\SetXXXX\>

%% User inputs, toggles:
% Set number
% Cooperativity Scaler
% Functional Tn density
% Functional XB density
% Tm Span
% pCa range with maximum time wanted to run the simulation at the pCa
% Desired number of runs per pCa

%% Simulation Outputs:
% Code and output will be written to an associated
% Data directory:
% <CurrentDir\DataFiles\SetXXXX\>

function Output_Location = Launch(Dir, SetNumber, RequestedRuns, kxscaler, TnKOType, TnFraction, START_LENGTH, L_TITIN)

PathToGit = '"C:\Program Files (x86)\Git\bin\git.exe"';

close all

disp(sprintf('Set %d\nkxscaler=%f\nTnFraction=%f\nSTART_LENGTH=%d\nL_TITIN=%d\n', SetNumber, kxscaler, TnFraction, START_LENGTH, L_TITIN))

%% Make working directory:
% cd('J:\TannerBertGroup\Sims\Vert_3Simple')
% cd('J:\TannerBertGroup\UndergradStudents\AlexWood\Vert_3Simple')
% Just stay in pwd

% SetNumber=1000; % Set number


pCaRange=[
   4.0000
    4.5
    5.0
%    5.25
     5.5 
     5.7 
    5.8
    5.9
    6.0
    6.1
    6.25
    6.5
    7.0
    7.5
    8.0
%    9.0000
    ];

% Pass the user inputs to an intiaization function, to set up the output
% directories:
% Current_Dir=pwd;
[Output_Location] = Initialize_v1(Dir, SetNumber);

% Initialize OutTotal as a cell array to hold the 
% Col 1) pCa Value
% Col 2) Average Time Series Output at a pCa value
% Col 3) SS Individual Trial Output
% Col 4) Time Half from Individual Trials
tic
OutTotal=cell(length(pCaRange), 4);
for i=1:length(pCaRange)
    OutTotal{i, 1}=pCaRange(i,1);
    [OutTotal{i,2}, OutTotal{i,3}, OutTotal{i,4}] = Run_v4_Alex(Output_Location, pCaRange(i,1), RequestedRuns, kxscaler, TnKOType, TnFraction, START_LENGTH, L_TITIN);
end

toc

OutMatFile=sprintf('%sSimulationData.mat', Output_Location);
clear('i'); % don't need to save this
save(OutMatFile)

%% Store the commit hash.
% What this does: The source should be under version control, which means
% every version of the code has a "hash" descriptor attached to it.
% We store this descriptor with simulation data so that later we know what
% code ran a given simulation.
% Our version control is git, so we can get the hash with "git rev-parse"
% on the last commit before the working tree, which is called "HEAD". See
% git's docs for more info.
% Unfortunately, msys-git's stuff is not in PATH, so we have to give a full
% path (set above). Good luck, other users. I apologize.
[success, hash] = system(sprintf('%s rev-parse HEAD', PathToGit));
if (success == 0)
    OutVC = sprintf('%scommit.txt', Output_Location);
    fd = fopen(OutVC, 'w');
    fprintf(fd, '%s', hash);
    fclose(fd);
else
    disp(sprintf('Warning: Could not store commit hash: git rev-parse HEAD returned unsuccessfully\n%s', out));
end
    
%% Plot
clf(figure(1))
plot(OutTotal{1,2}(:,1), OutTotal{1,2}(:,2), 'k-'), hold on %%%Plot Thick filament force vs. Time
plot(OutTotal{1,2}(:,1), OutTotal{1,2}(:,3), 'b-') %%%Plot Thin filament force vs. Time
legend('Thick', 'Thin', 'location', 'southeast')
ylabel('Force (pN)')
xlabel('Time (s)')

% Pass to Process_HillCurves_v1
% 1) pCa Range
% 2) Output directory
% 3) Finished Data cell array
% 4) Figure Number for plotting
% Process_Vert_3SimpleHillCurves_v2(pCaRange, Output_Location, OutTotal, 1, SetNumber, Tn_den, XB_den)

% % Following few lines are beta testing on tracking the attachment time.
% % Comment out for standard running, unless using BCWT's machine.
% pCa_string='4.00';
% data=importdata(['J:\TannerBertGroup\Sims\Vert_3Simple\DataFiles\Set' num2str(SetNumber) '\XBBindingData_pCa_' pCa_string '.txt']);
% Plot_ton_ExpDist_2013_7_31
