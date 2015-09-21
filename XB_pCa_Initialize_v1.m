%% BCWT, May 29, 2009
%% XB_pCa_Intialize

% Load all the constants, coefficients, and stiffness values for the
% lattice.  A good degree of this code will also load the XB cycling
% paramters and Free energies and such.  Also this will initialize some of
% our output/saving variables.

% Set time step for the simulation
dt=(1/1)*1e-3; % Milli second time steps.

% Reseed the SWB PRNG based on the clock.
RandStream.setDefaultStream(RandStream('swb2712','Seed',sum(100*clock)));

Ca=10^(-pCa);

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


%We need to initialize everything for the maximum sizes of NSTEPS:
NSTEPS=ceil(maxt/dt);
InitCol=NSTEPS;  % This will be the initail maximum columns.
% Saw we are going to average over the last tenth of the simulation
DataSpan=0.1;

%% We we don't want just 1 run
if Run_bool==1
    TotalRuns=1;
else
    % Now we need to calculate how many runs we want.
    % Say that we are going to want to average so many data points say:
    %DataPoints=40000; Were used in the two filament model.
    % Say that our tessellation now has 24 times as many nodes interacting
    % this number would then be 40k/24= 1.6k.  But we think we need at least
    % 5 runs.so make DataPoints:
%    DataPoints=3200;
    DataPoints=6400;
    %TotalRuns=ceil(DataPoints./(DataSpan*(NSTEPS/2))); %% Used for longer
    TotalRuns=ceil(DataPoints./(DataSpan*(NSTEPS/1)));
    % For testing on the longer duration runs at the mid-pCa levels
    % make them run 10 times.
%      if maxt==7
%          TotalRuns=15;
%      end
end



% outline the fraction multipliers to be used when creating the free energy
% drops in the the free energy diagrams.  This will then be used as
% efficienies and fraction of free energy differences.

%%___________________________New Rate Constants _________________________
eta = 0.685; % The mechanical efficiency fraction for the total amount of work performed.

f1 = 0.28; % The fraction drop from the first energy to the second minimum
%f1 = 0.28*(1/XB_Fraction);
% f2 = eta - f1; % To complete the balance of drop from second E min. to the 3rd E min.
%%
%%%% % These will be used by the DooItRates.m
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

AA=1;   %% k12 offset
BB=A;   %% k20 ???
CC=m;   %% k20 ???
DD=y0;  %% k20 offset--right?

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
dGtot=-GnotRT-log(ATP/(ADP*phos)); %Total free energy drop from the
% hydrolysis and concentration changes in units of RT
BridgeNRG=abs(eta*dGtot); % The total amount of energy allowed for the XB displacement
% givent he efficiency term on the total free energy available.  In units
% of RT.

RTnm2 = 8.314*T*J2pNnm/N_Avo; %  to convert RT/nm^2 into pN/nm, such that 1 RT/nm^2 =~ 4 pN/nm @ 310K
kRT = kxscaler/RTnm2; % cross brige spring constant in RT/nm^2
reach = sqrt(BridgeNRG./kRT); % The cross-brige reach, in nm


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

