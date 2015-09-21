%% BCWT, August 17, 2010-- Hack the Isometric Skeltal code:
%% Include new Geomtry--should be in the CaLookup File (4-start)
%% Prepare Istometric code for David to run some isometric force-pCa
%% simulations.  Make this more of a function use more initialization as a
%% script too.
%% Changes to this code:
%% 1) Force isometric, and only one time-half
%% 2) Functionalize
%% 3) Comment Out/Delete 2nd half averages
%% 4) Add in variances as output
%% 5) For ktr rate; remove the data fitting to an Exp and simply do the
%%    t1/2 from the SS isometric rate
%% 6) Bring into simulation runs a toggle for XBKO fraction.

%% Function Inputs:

%% Function Outputs:
% OutIndex%%(24, TotalRuns)
% OutExpFitInfo%%(28,1)
% OutSimpAve%%(20, 1)

%% Dependent scripts:


function [TimeSeries_DataOut, OutIndex, OutHalfTimes] = Run_v4_Alex(OutDir, pCa, RequestedRuns, kxscaler, TnKOType, TnFraction, START_LENGTH, L_TITIN)

% TnKOType = 0; % 0 for random, 1 for uniform
XBKOType = 0; % ditto

Psi=100; % Cooperativity Scaler
%TnFraction=1; % Functional Tn density
XB_Fraction=1; % Functional XB density
% Tm Span; 4 options; User Input 0, 1, 2, 3, 4
%% August 18, 2010--Work out the below details a bit more tidy
% 0 ) 1/2 Regulatory Unit , Spanning ~ 3-4 actins and ~ 17 nm
% 1 ) 1 Regulatory Unit , Spanning ~ 7 actins and ~ 37 nm
% 2 ) 1.33 Regulatory Unit Activation, Spanning ~ 9 actins and ~ 50 nm
% 3 ) 1.66 Regulatory Unit Activation, Spanning ~ 11 actins and ~ 62 nm
% 4 ) 2 Regulatory Units Activation, Spanning ~ 14 actins and ~ 74a nm
Tm_Type=2;

%% ADDITIONAL USER DEFINED, BUT NOT MESS WITH JUST YET
Coop_Type=7;
ScaleK1=1/100;
ScaleK2=1;
ScaleK3=1;
tic


% % Build our lattice, stiffness values
% Lattice_v1
% % Load all the constants, coefficients for XB cycling
% % paramters and Free energies and such.  Via script:
% XB_pCa_Initialize_v1
% Output_Variables_v1
% CaLookup

%% Scripts don't work with parfor's lexical analyzing, so the above have to be copied into this very file.
% TFstate is actually only used within runs, so don't need to init here.

% Reseed the SWB PRNG based on the clock.
RandStream.setDefaultStream(RandStream('swb2712','Seed',sum(100*clock)));

% Begin with Scalers of Stiffness Constants (spring constants):
kascaler = 1.0;
kmscaler = 1.0;%*XB_Fraction;
%kxscaler = 3.0; %3.0


%START_LENGTH = 1150;  %1200 in nm
%START_LENGTH = 1300;  %1200 in nm

LACT = 37.3/3; % L IS LENGTH, IN nm; rest length for thin fil lengths
LMYO = 42.9/3; % Rest length for thick fil lengths

KACT_0=1743; %thin fil stiff
KMYO_0=2020; %Thick fil stiff
KACT = kascaler*3*KACT_0;%1743; % K IS FOR SPRING CONSTANT, ALL IN pN / nm
KMYO = kmscaler*3*KMYO_0;%2020;

NACT = 8; % N IS 'NUMBER OF' filaments actin
NMYO = 4; % N IS 'NUMBER OF' filaments myosin
NumBridges=60 + 1; % 60 bridges per filament + handel
NumSites= 90 + 1; % 90 per filament + handel
ACTNodes=NACT*NumSites;
MYONodes=NMYO*NumBridges;
TOTNodes=ACTNodes+MYONodes;
HeadIndex=12; % The total number of counted indicies on the heads (12 rows)
% AMW note: no idea where HeadIndex is supposed to be used

%% Add some variable to keep track of our new arraingment of heads/thin fil
%% rows (Oct 20, 2010):
N_Thick_Start=3; %% 3-start thick filaments (3 heads per crown)
Angle_Thick_Start=120; %% Each of the heads will be separated by 120 degrees within a crown.
Angle_Crowns=40; % Each of the adjacent crowns will rotate by 40 degrees.  Pitch.
ThinRow_Index=24;  % Here we are modifying which row to look at.  This will begin to replace
% our venerable LAUREL and  HARDY.
% AMW note: unused, potentially used with HARDY

%% Play with titin--thick fil length =858 nm
%% START_LENGTH-858: 1200-858=342 nm
%% Say we want 20 pN from titin; and stiffness of 0.1 pN/nm
%% Then we must have the titin extended 200 nm
%% this would be a rest length of 342-200=142nm
K_TITIN=0.1344; % pN/nm
%L_TITIN=342; % 342 nm for 50 pN passive contribution at 2.6 um SL
%L_TITIN=247; % 247 nm for 24 pN passive contribution at 2.3 um SL

Mcore = MakeMcore(NACT, KACT, ACTNodes, NMYO, KMYO, NumSites, NumBridges, K_TITIN, TOTNodes);
EndVeccore = MakeEndVeccore(NACT, NumSites, NMYO, NumBridges, KACT, LACT, KMYO, LMYO, START_LENGTH, K_TITIN, L_TITIN, ACTNodes, TOTNodes);

NTn = 46;
Tn = MakeTn(NTn, NACT, NumSites);

HARDY = MakeHARDY(); % more moved out for space than anything.
LAUREL = MakeLAUREL(Angle_Crowns, NumBridges, Angle_Thick_Start, MYONodes, N_Thick_Start, NMYO);

dt=(1/1)*1e-3; % simulation timestep, in seconds
Ca = 10^(-pCa); % calcium level
maxt = SimLength(Ca); % simulation length, in seconds
NSTEPS = ceil(maxt/dt); % number of simulation steps to do.
DataSpan = 0.1; % what fraction of end of sim we'll average over.

if (RequestedRuns == 0)
    TotalRuns = SimRuns(6400, DataSpan, NSTEPS, maxt);
else
    TotalRuns = RequestedRuns;
end

t=zeros(length(Ca), NSTEPS);
for iCa=1:length(Ca)
    t(iCa, 1:NSTEPS(iCa)) = (dt*(1:NSTEPS(iCa)))-dt;
end

% outline the fraction multipliers to be used when creating the free energy
% drops in the the free energy diagrams.  This will then be used as
% efficienies and fraction of free energy differences.

%%                        New Rate Constants

% Rate constants are used in XBRates3/DooItRates/a similar function.
% They are packed into a structure together.
% Only the ones actually used are put into the thermochem struct; the rest
% are intermediates for computing the values in that struct.

thermochem.SXB1=1;
thermochem.SXB2=1;
thermochem.SXB3=1;

thermochem.eta = 0.685; % The mechanical efficiency fraction for the total amount of work performed.
thermochem.f1 = 0.28; % The fraction drop from the first energy to the second minimum
%thermo.f1 = 0.28*(1/XB_Fraction);
% f2 = eta - f1; % To complete the balance of drop from second E min. to the 3rd E min.

%% Parameters for k20
y0=20; % Offset
SlopeL=-100; %Slope of line left of 0
SlopeR=20; % Slope of line right of 0
xL=-1; % xposition left of 0;
xR=1; % xposition right of 0;
yL=SlopeL*xL + y0; % yposition left of 0;
yR=SlopeR*xR + y0; % yposition right of 0;
%So, we will have two parameters, m and A;
%These have been developed from Mathmatica for the
%solution of two systems of equations one with xL and yL
% and one with xR and yR, where x and y appear below:
%sqrt(A*(x^2)) + m*x + y0 = y
m=(-y0 + yR - (0.5*sqrt( ((xL*(y0-yR)-xR*(y0-yL))/xL)^2 )))/xR;
A=((xL*(y0-yR)-xR*(y0-yL))^2)/(4*(xL^2)*(xR^2));

thermochem.AA=1;   %% k12 offset
thermochem.BB=A;   %% k20 ???
thermochem.CC=m;   %% k20 ???
thermochem.DD=y0;  %% k20 offset--right?

ATP=5e-3; 		%  Intracellular ATP in Molar
ADP=30e-6;		%  Intracellular ADP in Molar
phos=3e-3;		%  Intracellular Pi in Molar


J2pNnm=1e21;
kCal2Joule=4.1868e3;
N_Avo=6.022e23;
% Gnot = 13*4.1868e3*1e21/6.022e23; % 13 kCal/mol (E of hyrolysis of ATP) * kCal2Joule * J2pNnm / N_Avagadro's (molecules/mol)
% % Gnot then becomes free energy of hydrolysis in pN*nm from the above line
Gnot=7.8*kCal2Joule; % Gnot of ATP hydrolysis in Joule/mol
% Gnot=13*kCal2Joule; % Gnot of ATP hydrolysis in Joule/mol
T=273.15; % Temperature in K
J2RT=1/(8.314*T); %To convert J/mol to units of RT [R=8.314 J/(mol K) ]
GnotRT=Gnot*J2RT; %Gnot in units of RT ~ 21
thermochem.dGtot=-GnotRT-log(ATP/(ADP*phos)); %Total free energy drop from the
% hydrolysis and concentration changes in units of RT
BridgeNRG=abs(thermochem.eta*thermochem.dGtot); % The total amount of energy allowed for the XB displacement
% givent he efficiency term on the total free energy available.  In units
% of RT.

RTnm2 = 8.314*T*J2pNnm/N_Avo; %  to convert RT/nm^2 into pN/nm, such that 1 RT/nm^2 =~ 4 pN/nm @ 310K
thermochem.kRT = kxscaler/RTnm2; % cross brige spring constant in RT/nm^2
thermochem.reach = sqrt(BridgeNRG./thermochem.kRT); % The cross-brige reach, in nm


%% Now deal with the thin filament rates and the cooperativity:
%% Arange to pass into function.  FYI, the rates influenced
% CoopPass=zeros(13, 1);
% CoopPass(1, 1)=TF3TF12;
% CoopPass(2, 1)=XB2TF12;
% CoopPass(3, 1)=XB3TF12;
% CoopPass(4, 1)=XB2TF23;
% CoopPass(5, 1)=XB3TF23;
% CoopPass(6, 1)=TF3TF23;
% CoopPass(7, 1)=XB2TF21;
% CoopPass(8, 1)=XB3TF21;
% CoopPass(9, 1)=TF3TF21;
% CoopPass(10, 1)=XB2TF32;
% CoopPass(11, 1)=XB3TF32;
% CoopPass(12, 1)=TF3TF32;
% CoopPass(13, 1)=CooperativeType;
if Coop_Type==7
    CoopPass=ones(13,1);
    CoopPass(1:6, 1)=Psi^(1/3);
else
    CoopPass=ones(13,1);
end
%% Seven types of cooperativity with numerical values 1-7
%% 1: RURU Only
%% 2: XB3RU Only
%% 3: Both RURU or XB3RU
%% 4: XB2RU Only
%% 5: Both RURU or XB2RU
%% 6: Both XB3RU or XB2RU
%% 7: All three possible; RURU, XB3RU, XB2RU
CoopPass(13, 1)=Coop_Type;


%% Call the initial RU Rates
%% Where 1=ScaleK1=ScaleK2=ScaleK3
%% into the function.
TFRatePass=zeros(6,1);
%[koff, kon, RuOff, RuOn, CaOff, CaOn] = ScaleThinFilRates( 1,1,1 );
% Unpacked later as:
%%% Unpack Thin Filament Rates
% koff=TFRatePass(1,1);
% kon=TFRatePass(2,1);
% RuOff=TFRatePass(3,1);
% RuOn=TFRatePass(4,1);
% CaOff=TFRatePass(5,1);
% CaOn=TFRatePass(6,1);
[TFRatePass(1,1), TFRatePass(2,1), TFRatePass(3,1), TFRatePass(4,1), TFRatePass(5,1), TFRatePass(6,1)] = ScaleThinFilRates( ScaleK1, ScaleK2, ScaleK3 );

%% Ca should have a single value
iCa=1;

%% Big Binder
Big_Binder_3=[];

%% Steps we do more recording on
EndStepStart = (1-DataSpan)*NSTEPS;
EndSteps = EndStepStart:NSTEPS;

%% Preinitialized output variables
%% Time series (i.e. per step), averaged across runs
% originally the dimensions were "length(Ca),InitCol"
% but I'm not sure how Ca could be anything but scalar
% and InitCol is just a synonym for NSTEPS
MFvec = zeros(1,NSTEPS);
AFvec = zeros(1,NSTEPS);
FractXB1 = zeros(1,NSTEPS);
FractXB2 = zeros(1,NSTEPS);
FractCa0 = zeros(1,NSTEPS);
FractCa1 = zeros(1,NSTEPS);
FractCa2 = zeros(1,NSTEPS);
ATPuse = zeros(1,NSTEPS);

%% Means
MFMean = zeros(1,TotalRuns);
AFMean = zeros(1,TotalRuns);
FractXB1Mean = zeros(1,TotalRuns);
FractXB2Mean = zeros(1,TotalRuns);
FractCa0Mean = zeros(1,TotalRuns);
FractCa1Mean = zeros(1,TotalRuns);
FractCa2Mean = zeros(1,TotalRuns);
ATPuseMean = zeros(1,TotalRuns);
FractBoundMean = zeros(1,TotalRuns);

%% Variances
MFVar = zeros(1,TotalRuns);
AFVar = zeros(1,TotalRuns);
FractXB1Var = zeros(1,TotalRuns);
FractXB2Var = zeros(1,TotalRuns);
FractCa0Var = zeros(1,TotalRuns);
FractCa1Var = zeros(1,TotalRuns);
FractCa2Var = zeros(1,TotalRuns);
ATPuseVar = zeros(1,TotalRuns);
FractBoundVar = zeros(1,TotalRuns);

%% Time to reach half the mean (I think)
IndexThalfMFvec = zeros(1,TotalRuns);
IndexThalfAFvec = zeros(1,TotalRuns);
IndexThalfFractXB1 = zeros(1,TotalRuns);
IndexThalfFractXB2 = zeros(1,TotalRuns);
IndexThalfFractCa0 = zeros(1,TotalRuns);
IndexThalfFractCa1 = zeros(1,TotalRuns);
IndexThalfFractCa2 = zeros(1,TotalRuns);
IndexThalfATPuse = zeros(1,TotalRuns);
IndexThalfFractBound = zeros(1,TotalRuns);

parfor iTotRun=1:TotalRuns
    % too many variables
    % ATPuse (Output_Variables.m) looks like an output but is averaged over all runs
    % TotalRuns should really not be necessary
    % DataSpan is also an Output_Variables thing
    % the way Tm_Type is used is problematic
    % Tn is a big ol' matrix in CaLookup, "the equivalent of HARDY for the Troponin Nodes"
    % dt is the timestep in seconds and is tied to NSTEPS
    % no real idea about CoopPass. or TFRatePass
    % Ca is just the calcium concentration (in molar i think)
    % iCa is defined upfile here. no idea.
    [Binder_3, TempFractCa0, TempFractCa1, TempFractCa2, TempFractXB1, TempFractXB2, TempMFvec, TempAFvec, TempATPuse] = OneRun(Tm_Type, dt, NSTEPS, EndStepStart, LACT, KACT, NACT, LMYO, KMYO, NMYO, NTn, Tn, TnFraction, HARDY, LAUREL, XB_Fraction, ACTNodes, MYONodes, NumSites, NumBridges, Mcore, EndVeccore, N_Thick_Start, CoopPass, TFRatePass, Ca, iCa, kxscaler, thermochem, TnKOType, XBKOType);
    Big_Binder_3 = [Big_Binder_3; Binder_3];
    
    MFvec=MFvec+TempMFvec/TotalRuns;
    AFvec=AFvec+TempAFvec/TotalRuns;
    FractXB1=FractXB1+TempFractXB1/TotalRuns;
    FractXB2=FractXB2+TempFractXB2/TotalRuns;
    FractCa0=FractCa0+TempFractCa0/TotalRuns;
    FractCa1=FractCa1+TempFractCa1/TotalRuns;
    FractCa2=FractCa2+TempFractCa2/TotalRuns;
    ATPuse=ATPuse+TempATPuse/TotalRuns;
    
    MFMean(iTotRun) = mean(TempMFvec(EndSteps));
    AFMean(iTotRun) = mean(TempAFvec(EndSteps));
    FractXB1Mean(iTotRun) = mean(TempFractXB1(EndSteps));
    FractXB2Mean(iTotRun) = mean(TempFractXB2(EndSteps));
    FractCa0Mean(iTotRun) = mean(TempFractCa0(EndSteps));
    FractCa1Mean(iTotRun) = mean(TempFractCa1(EndSteps));
    FractCa2Mean(iTotRun) = mean(TempFractCa2(EndSteps));
    ATPuseMean(iTotRun) = mean(TempATPuse(EndSteps));
    FractBoundMean(iTotRun) = FractXB1Mean(iTotRun) + FractXB2Mean(iTotRun); % sums: convenient
    
    MFVar(iTotRun)=var(TempMFvec(EndSteps), 1);
    AFVar(iTotRun)=var(TempAFvec(EndSteps), 1);
    FractXB1Var(iTotRun)=var(TempFractXB1(EndSteps), 1);
    FractXB2Var(iTotRun)=var(TempFractXB2(EndSteps), 1);
    FractCa0Var(iTotRun)=var(TempFractCa0(EndSteps), 1);
    FractCa1Var(iTotRun)=var(TempFractCa1(EndSteps), 1);
    FractCa2Var(iTotRun)=var(TempFractCa2(EndSteps), 1);
    ATPuseVar(iTotRun)=var(TempATPuse(EndSteps), 1);
    FractBoundVar(iTotRun)=var(TempFractXB1(EndSteps) + TempFractXB2(EndSteps), 1); % these are correlated probably, so we can't just sum variances
    
    %% Also calculate the time to reach 1/2 of the SS value for each of the
    %% variables.
    %% Load them with time value, closest to half of the max avg value
    %% This is irritatingly verbose but I don't know how to fix it.
    [~,Index] = min(abs(TempMFvec-(0.5*MFMean(iTotRun))));
    IndexThalfMFvec(iTotRun)=dt*Index;
    [~,Index] = min(abs(TempAFvec-(0.5*AFMean(iTotRun))));
    IndexThalfAFvec(iTotRun)=dt*Index;
    [~,Index] = min(abs(TempFractXB1-(0.5*FractXB1Mean(iTotRun))));
    IndexThalfFractXB1(iTotRun)=dt*Index;
    [~,Index] = min(abs(TempFractXB2-(0.5*FractXB2Mean(iTotRun))));
    IndexThalfFractXB2(iTotRun)=dt*Index;
    [~,Index] = min(abs(TempFractCa0-(0.5*FractCa0Mean(iTotRun))));
    IndexThalfFractCa0(iTotRun)=dt*Index;
    [~,Index] = min(abs(TempFractCa1-(0.5*FractCa1Mean(iTotRun))));
    IndexThalfFractCa1(iTotRun)=dt*Index;
    [~,Index] = min(abs(TempFractCa2-(0.5*FractCa2Mean(iTotRun))));
    IndexThalfFractCa2(iTotRun)=dt*Index;
    [~,Index] = min(abs(TempATPuse-(0.5*ATPuseMean(iTotRun))));
    IndexThalfATPuse(iTotRun)=dt*Index;
    [~,Index] = min(abs((TempFractXB1(EndSteps) + TempFractXB2(EndSteps))- (0.5*FractBoundMean(iTotRun))));
    IndexThalfFractBound(iTotRun)=dt*Index;

end

%%%%%%%%%%%%%%%%%%%% BEGIN INDEX (All Runs) OUTPUT PARAMS %%%%%%%%%%%%%%%%%
% Package Up the Averaged (Index) information into a vector.
% Also, flip to column vectors, so we can write it to a outfile
OutIndex=zeros(TotalRuns, 20);
%% 24 rows for the first half results
OutIndex(:, 1)=(1:TotalRuns)';
OutIndex(:, 2)=pCa*ones(TotalRuns, 1);
%% Concat the SS mean data, followed by Variance
OutIndex(:, 3)=MFMean';
OutIndex(:, 4)=MFVar';
OutIndex(:, 5)=AFMean';
OutIndex(:, 6)=AFVar';
OutIndex(:, 7)=FractXB1Mean';
OutIndex(:, 8)=FractXB1Var';
OutIndex(:, 9)=FractXB2Mean';
OutIndex(:, 10)=FractXB2Var';
OutIndex(:, 11)=FractCa0Mean';
OutIndex(:, 12)=FractCa0Var';
OutIndex(:, 13)=FractCa1Mean';
OutIndex(:, 14)=FractCa1Var';
OutIndex(:, 15)=FractCa2Mean';
OutIndex(:, 16)=FractCa2Var';
OutIndex(:, 17)=ATPuseMean';
OutIndex(:, 18)=ATPuseVar';
OutIndex(:, 19)=FractBoundMean';
OutIndex(:, 20)=FractBoundVar';


%% NOW Concat the half times
OutHalfTimes=zeros(TotalRuns, 11);
OutHalfTimes(:, 1)=(1:TotalRuns)';
OutHalfTimes(:, 2)=pCa*ones(TotalRuns, 1);
OutHalfTimes(:, 3)=IndexThalfMFvec';
OutHalfTimes(:, 4)=IndexThalfAFvec';
OutHalfTimes(:, 5)=IndexThalfFractXB1';
OutHalfTimes(:, 6)=IndexThalfFractXB2';
OutHalfTimes(:, 7)=IndexThalfFractCa0';
OutHalfTimes(:, 8)=IndexThalfFractCa1';
OutHalfTimes(:, 9)=IndexThalfFractCa2';
OutHalfTimes(:, 10)=IndexThalfATPuse';
OutHalfTimes(:, 11)=IndexThalfFractBound';
%%%%%%%%%%%%%%%%%%%% END INDEX (All Runs) OUTPUT PARAMS %%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%BEGIN TIME SERIES OUPUT TO FILE %%%%%%%%%%%%%%%%%%
% Where every row here is a data input.

TimeSeries_DataOut=[...
    t(1, 1:NSTEPS); MFvec(1, 1:NSTEPS); AFvec(1, 1:NSTEPS);...
    ATPuse(1, 1:NSTEPS); FractXB1(1, 1:NSTEPS)+FractXB2(1, 1:NSTEPS); ...
    FractXB1(1, 1:NSTEPS); FractXB2(1, 1:NSTEPS); ...
    FractCa0(1, 1:NSTEPS); FractCa1(1, 1:NSTEPS); FractCa2(1, 1:NSTEPS)]';
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
if not(isempty(Big_Binder_3))
    OutFile=sprintf('%sXBBindingData_pCa_%s.txt', OutDir, num2str(pCa, '%3.2f'));
    fid=fopen(OutFile,'w');	%open outfile--tab delimited text
    %Create the output file column header:
    fprintf(fid, 'Run Index\tStartNSTEP\tTot_dt_Bnd\tXB1_dt_Bnd\tXB2_dt_Bnd\tXB1--->XB2\tXB2--->XB1\tXB2--->XB0\tAVG Xthick\tAVG X_thin\tAvgXBstrain\tAvgXB1strn\tAvgXB2strn\tThick Node\tThin FNode\tkxscaler  \treach(xb0)\tCompleted \tVAR Xthick\tVAR X_thin\tVARXBstrain\tVARXB1strn\tVARXB2strn\n');	%header for outfile
    FormatString=[];
    ColOut=size(Big_Binder_3, 2); % not sure if valid
    for i=1:ColOut-1 %for all but last
        FormatString=[FormatString, '%10.6f\t'];
    end
    FormatString=[FormatString, '%10.6f\n'];
    % Write Data File and close file
    fprintf(fid, FormatString, Big_Binder_3');
    fclose(fid);     %close the file
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Binding Events OUTPUT TO FILE %%%%%%%%%%%%%%%%%%
