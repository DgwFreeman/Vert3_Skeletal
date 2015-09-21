%% AMW, move individual run to separate function
%% Not sure about inputs and outputs yet

function [Binder_3, FractCa0, FractCa1, FractCa2, FractXB1, FractXB2, MFvec, AFvec, ATPuse] = OneRun(Tm_Type, dt, NSTEPS, StepRecord, LACT, KACT, NACT, LMYO, KMYO, NMYO, NTn, Tn, TnFraction, HARDY, LAUREL, XB_Fraction, ACTNodes, MYONodes, NumSites, NumBridges, Mcore, EndVeccore, N_Thick_Start, CoopPass, TFRatePass, Ca, iCa, kxscaler, thermochem, TnKOType, XBKOType)
    
%% First we do the knockouts.
% We coerce the random arrays to doubles (they start out boolean) because
% for some reason I don't want to think too hard about, that makes things
% faster.

if (TnKOType == 0)
    % Random knockout
    TnKO = double(rand(NACT, NTn) <= TnFraction);
else
    % "Uniform" knockout from the free end of the thin thilament
    TnKO=ones(NACT, NTn); %% Fill as though no KOs
    [~, jTest]=find(TnKO==1); %% Gather row and col indicies to
    %% filter on the KO
    TnKO=zeros(NACT, NTn);
    for iii=1:length(jTest) %% Loop down each row, looking at first half
        %% then the second half to fill with zero if NTn<(1-TnKO)
        %% Now the matrix is zeros, update with the Tn that will actually
        %% exist.
        if(iii<=(NTn/2)) %% First Half
            %if( (iii/(NTn/2))<(1-TnFraction) ) %% ERROR--9/3/2013
            if( (iii/(NTn/2))<=(TnFraction) )
                TnKO(:, iii)=1;
            end
        end
        if(iii>(NTn/2)) %% 2nd Half
            %if( ((iii-NTn/2)/(NTn/2))<(1-TnFraction) ) %% ERROR--9/3/2013
            if( ((iii-NTn/2)/(NTn/2))<=(TnFraction) )
                TnKO(:, iii)=1;
            end
        end
    end
end


if (XBKOType == 0)
    % Random knockout
    XBKO = double(rand(size(LAUREL)) <= XB_Fraction);
else
    % "Uniform" knockout from the thick filament's free end
    XBKO=ones(size(LAUREL)); %% Fill as though no KOs
    [iTest, jTest]=find(XBKO==1); %% Gather row and col indicies to filter the KO on
    for iii=1:length(iTest)
        % If we are XB_Fraction near the free end of the filament
        % Set XBKO to 0.  This will also set XBK) to zeros
        % for the 61st XB Node (handle) at the M-line
        if(mod(iii,NumBridges)<=NumBridges*(1-XB_Fraction))
            XBKO(iTest(iii,1), jTest(iii,1) )=0;
        end
    end
end

%% Now all the variables, before the loop starts proper

CountOff=0; % To count the total number of reverse transitions from
% unbound to strongly bound in the XB cycle.

TFstate=zeros(ACTNodes, 3);  % Refresh initial value
TFstate(1:NumSites:ACTNodes, :)=-1; % To set the handle nodes to -1, so as
% to not interfere with thinking they are counted as an unboud state.
Prior_TFstate=[];
Binder_3=[];

% Set up our temperary variables/zero them out.
% NOTE: These names used to be for all-run averages, but that has been
% moved out of this function.
FractCa0=zeros(1, NSTEPS);
FractCa1=zeros(1, NSTEPS);
FractCa2=zeros(1, NSTEPS);
FractXB1=zeros(1, NSTEPS);
FractXB2=zeros(1, NSTEPS);
MFvec=zeros(1, NSTEPS);
AFvec=zeros(1, NSTEPS);
ATPuse=zeros(1, NSTEPS);

M = Mcore;

X = M\EndVeccore;

%% Do every step.

for iStep = 1:NSTEPS
    TFstate = DispatchCaRegCoop(Tm_Type, NACT, NTn, TFstate, Tn, dt, CoopPass, TFRatePass, Ca, iCa, TnKO);
    if iStep>=StepRecord
        %% If we want to keep track of time-on - only at end of sim, where we also will average the steady-state information
        [Prior_TFstate, Binder_3] = Time_On_Calculator(X, TFstate, Prior_TFstate, Binder_3, StepRecord, iStep, thermochem.reach, kxscaler);
    end
    AFvec(iStep) = ComputeActinForce(X, NumSites, LACT, KACT);
    MFvec(iStep) = ComputeMyosinForce(X, NumBridges, ACTNodes, LMYO, KMYO);
    [TempCount, TFstate, ATPuse(iStep)] = OneStep(dt, HARDY, LAUREL, X, TFstate, N_Thick_Start, MYONodes, ACTNodes, XBKO, thermochem);
    CountOff = CountOff + TempCount;
    M = alterM(Mcore, TFstate, kxscaler);
    X = M\EndVeccore;
    [FractXB1(iStep), FractXB2(iStep)] = ComputeXBFracts(TFstate, N_Thick_Start, MYONodes, NMYO);
    [FractCa0(iStep), FractCa1(iStep), FractCa2(iStep)] = ComputeCaFracts(TFstate, ACTNodes, NACT);
end % end of one run (iStep)

% Output from a given run.

% Update the fitting in the matrix.
% Collapse the Temp. vectors into the big vector that will be used
% for final averages:
% FractXB1(1, 1:NSTEPS)=FractXB1(1, 1:NSTEPS)+TempFractXB1/TotalRuns;
% FractXB2(1, 1:NSTEPS)=FractXB2(1, 1:NSTEPS)+TempFractXB2/TotalRuns;
% FractCa0(1, 1:NSTEPS)=FractCa0(1, 1:NSTEPS)+TempFractCa0/TotalRuns;
% FractCa1(1, 1:NSTEPS)=FractCa1(1, 1:NSTEPS)+TempFractCa1/TotalRuns;
% FractCa2(1, 1:NSTEPS)=FractCa2(1, 1:NSTEPS)+TempFractCa2/TotalRuns;
% MFvec(1, 1:NSTEPS)=MFvec(1, 1:NSTEPS)+TempMFvec/TotalRuns;
% AFvec(1, 1:NSTEPS)=AFvec(1, 1:NSTEPS)+TempAFvec/TotalRuns;


