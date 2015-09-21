%% Initialize Output Variables:


%TmSpan=1; % Number of Symetrical nodes outside of the 'super actin/Tn' to activate.

%nStates=2; % Number of States the Actin BS can hold

% Load the random seed as a function of time
%time=clock;
%RandStream.setDefaultStream(RandStream('mcg16807','Seed',sum(time(3:end))));

% Define the thin filament binding state, where set to zero initially.
% The first col. will be the sate of the thin filament (Ca binding)
% The second col. will keep track of each actin node's XB binding (XB state).
% The third col. will keep track of which XB each actin node is binding.
TFstate=zeros(ACTNodes, 3);
% Load in some to keep track of the final tallies.
%% It turns out that without the LaCie drive, we cannot accomodate
%% this large set:
%BigTFState=zeros(ACTNodes, InitCol, length(Ca));
%BigBoundState=BigTFState;
%BigBoundPair=BigTFState;

%%%%%%%%%%%%%%%%%%%%%WILL BE TIME SERIES OUPUT TO FILE %%%%%%%%%%%%%%%%%%
%% However, we must take care of these variables:
AFvec=zeros(length(Ca),InitCol);  % Will become the average force at time for Actin
MFvec=AFvec;% Will become the average force at time for Myosin
%X=zeros(TOTNodes, 1); %Because we need to know possition.
ATPuse=AFvec;% To store the average ATP useage.
FractCa0=AFvec; % To store the total number in this state over the run, to ensure that
% nothing is numerically increasing.  It shouldn't we are a discrete
% system.  Anyway, NCa0 is the number in Ca State 0, for thin filament
% activation.  There is also NCa1, NCa2, NXB0, NXB1, and NXB2 for the rest
% of the Ca thin filament states and then the 3 (0-2) cross-bridge states.
FractCa1=AFvec;
FractCa2=AFvec;% What used to be FractAvailAct
% FractXB0=AFvec;  We must take this one out, because we cannot explicitly
% account for it given the code we have.  Outlined in the paper notes.
% This can however be calculated from 1-(FractXB1+FractXB2)
FractXB1=AFvec;
FractXB2=AFvec;% What used to be FractBoundXB
%BigX=zeros(TOTNodes, InitCol, length(Ca));
% Cannot cut the BigX, just too big for right now
% Horse walks into a bar.  Bartender says; Why the long face?
t=zeros(length(Ca), InitCol);
for iCa=1:length(Ca)
    t(iCa, 1:NSTEPS(iCa)) = (dt*(1:NSTEPS(iCa)))-dt;
end
%%%%%%%%%%%%%%%%%%%%%%%%% END TIME SERIES OUPUT TO FILE %%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% WILL BE INDEX (All Runs) OUTPUT PARAMS %%%%%%%%%%%%%%%%%
%%% Now Outline the intermediate indexes:
%%%%%%%%%%%%%%%%%%%%%%%%% All for the 1st Half of the Simulation:
%%%%%%%%%%%%%%%%%%%%%%%%% All for the 1st Half of the Simulation:
%%%%%%%%%%%%%%%%%%%%%%%%% All for the 1st Half of the Simulation:
%% First, Intermediate Indexing for the Averages:
IndexSSFractBound=zeros(1, TotalRuns);
IndexVARFractBound=zeros(1, TotalRuns);

IndexSSMFvec=zeros(1, TotalRuns);
IndexSSAFvec=zeros(1, TotalRuns);
IndexSSFractXB1=zeros(1, TotalRuns);
IndexSSFractXB2=zeros(1, TotalRuns);
IndexSSFractCa0=zeros(1, TotalRuns);
IndexSSFractCa1=zeros(1, TotalRuns);
IndexSSFractCa2=zeros(1, TotalRuns);
IndexSSATPuse=zeros(1, TotalRuns);
% Continues with Intermediate Indexing for the
% standard deviations of the index forming the Averages:
IndexVARMFvec=zeros(1, TotalRuns);
IndexVARAFvec=zeros(1, TotalRuns);
IndexVARFractXB1=zeros(1, TotalRuns);
IndexVARFractXB2=zeros(1, TotalRuns);
IndexVARFractCa0=zeros(1, TotalRuns);
IndexVARFractCa1=zeros(1, TotalRuns);
IndexVARFractCa2=zeros(1, TotalRuns);
IndexVARATPuse=zeros(1, TotalRuns);

%%%%%%%%%%%%%%%%%%%% END INDEX (All Runs) OUTPUT PARAMS %%%%%%%%%%%%%%%%%

    %% Also calculate the time to reach 1/2 of the SS value for each of the
    %% variables.
        IndexThalfMFvec=zeros(1, TotalRuns);
        IndexThalfAFvec=zeros(1, TotalRuns);
        IndexThalfFractXB1=zeros(1, TotalRuns);
        IndexThalfFractXB2=zeros(1, TotalRuns);
        IndexThalfFractCa0=zeros(1, TotalRuns);
        IndexThalfFractCa1=zeros(1, TotalRuns);
        IndexThalfFractCa2=zeros(1, TotalRuns);
        IndexThalfATPuse=zeros(1, TotalRuns);
        IndexThalfFractBound=zeros(1, TotalRuns);

