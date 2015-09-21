function [Steps, Means, Vars, IndexThalf, Binder] = RunSeveral(RequestedRuns, DataParams, StartLength, pCa, StiffScale, filaments, knockout, coop, TFRateScale, tcparam)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Ca = 10^(-pCa); % calcium level

%% Data/sim level stuff

DataSpan = DataParams.DataSpan;
dt=DataParams.dt; % simulation timestep, in seconds
maxt = SimLength(Ca); % simulation length, in seconds
NSTEPS = ceil(maxt/dt); % number of simulation steps to do.
EndStepStart = (1-DataSpan)*NSTEPS;
EndSteps = EndStepStart:NSTEPS; % the steps we record more closely because they are steady-state

if (RequestedRuns == 0)
    TotalRuns = SimRuns(6400, DataSpan, NSTEPS, maxt);
else
    TotalRuns = RequestedRuns;
end

%% Knockout!

TnKOType = knockout.TnKOType; % 0 for random, 1 for uniform
XBKOType = knockout.XBKOType; % ditto

TnFraction=knockout.TnFraction; % Functional Tn density
XB_Fraction=knockout.XB_Fraction; % Functional XB density

%% Filament and stiffness info

kxscaler = StiffScale.kxscaler;

Tm_Type = filaments.Tm_Type;

LACT = filaments.LACT;
LMYO = filaments.LMYO;
L_TITIN = filaments.L_TITIN;
NACT = filaments.NACT; % N IS 'NUMBER OF' filaments actin
NMYO = filaments.NMYO; % N IS 'NUMBER OF' filaments myosin
NumBridges = filaments.NumBridges; % bridges per filament + handel
NumSites = filaments.NumSites; % sites per filament + handel
N_Thick_Start = filaments.N_Thick_Start; %% 3-start thick filaments (3 heads per crown)

NTn = filaments.NTn;

ACTNodes=NACT*NumSites;
MYONodes=NMYO*NumBridges;

TOTNodes=ACTNodes+MYONodes;

KACT = StiffScale.kascaler*3*filaments.KACT_0;
KMYO = StiffScale.kmscaler*3*filaments.KMYO_0;
K_TITIN = filaments.K_TITIN; % internal to this file, but L_TITIN should rely on it

%% Cooperativity

CoopPass=zeros(13, 1);
CoopPass(1, 1)=coop.TF3TF12;
CoopPass(2, 1)=coop.XB2TF12;
CoopPass(3, 1)=coop.XB3TF12;
CoopPass(4, 1)=coop.XB2TF23;
CoopPass(5, 1)=coop.XB3TF23;
CoopPass(6, 1)=coop.TF3TF23;
CoopPass(7, 1)=coop.XB2TF21;
CoopPass(8, 1)=coop.XB3TF21;
CoopPass(9, 1)=coop.TF3TF21;
CoopPass(10, 1)=coop.XB2TF32;
CoopPass(11, 1)=coop.XB3TF32;
CoopPass(12, 1)=coop.TF3TF32;
CoopPass(13, 1)=coop.Coop_Type;

%% thermochemical stuff

thermochem.SXB1=tcparam.SXB1;
thermochem.SXB2=tcparam.SXB2;
thermochem.SXB3=tcparam.SXB3;

thermochem.eta = tcparam.eta; % The mechanical efficiency fraction for the total amount of work performed.
thermochem.f1 = tcparam.f1; % The fraction drop from the first energy to the second minimum
%thermo.f1 = 0.28*(1/XB_Fraction);
% f2 = eta - f1; % To complete the balance of drop from second E min. to the 3rd E min.

y0=tcparam.y0; % Offset
xL=tcparam.xL; % xposition left of 0;
xR=tcparam.xR; % xposition right of 0;
yL=tcparam.SlopeL*xL + y0; % yposition left of 0;
yR=tcparam.SlopeR*xR + y0; % yposition right of 0;
m=(-y0 + yR - (0.5*sqrt( ((xL*(y0-yR)-xR*(y0-yL))/xL)^2 )))/xR;
A=((xL*(y0-yR)-xR*(y0-yL))^2)/(4*(xL^2)*(xR^2));

thermochem.AA=1;   %% k12 offset, not sure if should be parametrized
thermochem.BB=A;   %% k20 ???
thermochem.CC=m;   %% k20 ???
thermochem.DD=y0;  %% k20 offset--right?

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
thermochem.dGtot=-GnotRT-log(tcparam.ATP/(tcparam.ADP*tcparam.phos)); %Total free energy drop from the
% hydrolysis and concentration changes in units of RT
BridgeNRG=abs(thermochem.eta*thermochem.dGtot); % The total amount of energy allowed for the XB displacement
% givent he efficiency term on the total free energy available.  In units
% of RT.

RTnm2 = 8.314*T*J2pNnm/N_Avo; %  to convert RT/nm^2 into pN/nm, such that 1 RT/nm^2 =~ 4 pN/nm @ 310K
thermochem.kRT = kxscaler/RTnm2; % cross brige spring constant in RT/nm^2
thermochem.reach = sqrt(BridgeNRG./thermochem.kRT); % The cross-brige reach, in nm

%% Construct some big stuff
Tn = MakeTn(NTn, NACT, NumSites);

Mcore = MakeMcore(NACT, KACT, ACTNodes, NMYO, KMYO, NumSites, NumBridges, K_TITIN, TOTNodes);
EndVeccore = MakeEndVeccore(NACT, NumSites, NMYO, NumBridges, KACT, LACT, KMYO, LMYO, StartLength, K_TITIN, L_TITIN, ACTNodes, TOTNodes);

HARDY = MakeHARDY(); % more moved out for space than anything.
LAUREL = MakeLAUREL(filaments.Angle_Crowns, NumBridges, filaments.Angle_Thick_Start, MYONodes, N_Thick_Start, NMYO);

[TFRatePass(1,1), TFRatePass(2,1), TFRatePass(3,1), TFRatePass(4,1), TFRatePass(5,1), TFRatePass(6,1)] = ScaleThinFilRates(TFRateScale.ScaleK1, TFRateScale.ScaleK2, TFRateScale.ScaleK3 );

%% XB binding stuff
Big_Binder_3 = [];

%% Data per time step, averaged across runs
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
    % iCa should be removed
    [Binder_3, TempFractCa0, TempFractCa1, TempFractCa2, TempFractXB1, TempFractXB2, TempMFvec, TempAFvec, TempATPuse] = OneRun(Tm_Type, dt, NSTEPS, EndStepStart, LACT, KACT, NACT, LMYO, KMYO, NMYO, NTn, Tn, TnFraction, HARDY, LAUREL, XB_Fraction, ACTNodes, MYONodes, NumSites, NumBridges, Mcore, EndVeccore, N_Thick_Start, CoopPass, TFRatePass, Ca, 1, kxscaler, thermochem, TnKOType, XBKOType);
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

Binder = Big_Binder_3;
Steps = [MFvec; AFvec; FractXB1; FractXB2; FractCa0; FractCa1; FractCa2; ATPuse];
Means = [MFMean; AFMean; FractXB1Mean; FractXB2Mean; FractCa0Mean; FractCa1Mean; FractCa2Mean; ATPuseMean; FractBoundMean];
Vars = [MFVar; AFVar; FractXB1Var; FractXB2Var; FractCa0Var; FractCa1Var; FractCa2Var; ATPuseVar; FractBoundVar];
IndexThalf = [IndexThalfMFvec; IndexThalfAFvec; IndexThalfFractXB1; IndexThalfFractXB2; IndexThalfFractCa0; IndexThalfFractCa1; IndexThalfFractCa2; IndexThalfATPuse; IndexThalfFractBound];

end

