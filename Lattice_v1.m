%% Alex's Notes

%% Inout: script

%% Dependencies: none (apparently?)

%% Defined variables:

% kascaler, kmscaler, kxscaler
% START_LENGTH
% LACT, LMYO
% KACT_0, KMYO_0, KACT, KMYO
% NACT, NMYO, NumBridges, NumSites, ACTNodes, MYONodes, TOTNodes, HeadIndex
% N_Thick_Start, Angle_Thick_Start, Angle_Crowns, ThinRow_Index
% Mcore [IMPORTANT: made sparse later. why not now?]
% EndVeccore
% K_TITIN, L_TITIN

%% File local variables

% kascaler, kmscaler
% START_LENGTH, KACT_0, KMYO_0
% HeadIndex, K_TITIN, L_TITIN

%% Defined here but used in only CaLookup
% Angle_Crowns, ThinRow_Index, Angle_Thick_Start
% N_Thick_Start

%% Mcore and EndVeccore, the premier results of this file, are only used read-only (in Run)

%% Purpose
% Seems to be setting up Mcore and EndVeccore for later use
% Other variables also used outside, but I'm not sure which yet

%% Code

% Begin with Scalers of Stiffness Constants (spring constants):
kascaler = 1.0;
kmscaler = 1.0;%*XB_Fraction;
kxscaler = 3.0; %3.0


START_LENGTH = 1150;  %1200 in nm
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

%% Add some variable to keep track of our new arraingment of heads/thin fil
%% rows (Oct 20, 2010):
N_Thick_Start=3; %% 3-start thick filaments (3 heads per crown)
Angle_Thick_Start=120; %% Each of the heads will be separated by 120 degrees within a crown.
Angle_Crowns=40; % Each of the adjacent crowns will rotate by 40 degrees.  Pitch.
ThinRow_Index=24;  % Here we are modifying which row to look at.  This will begin to replace
% our venerable LAUREL and  HARDY.

Mcore = zeros(TOTNodes);% the part of the spring constant matrix that doesn't change
EndVeccore = zeros(TOTNodes,1); % column vector (for now) of end conditions

% initializes Mcore

%% Play with titin--thick fil length =858 nm
%% START_LENGTH-858: 1200-858=342 nm
%% Say we want 20 pN from titin; and stiffness of 0.1 pN/nm
%% Then we must have the titin extended 200 nm
%% this would be a rest length of 342-200=142nm
K_TITIN=0.1344; % pN/nm
%L_TITIN=342; % 342 nm for 50 pN passive contribution at 2.6 um SL
L_TITIN=247; % 247 nm for 24 pN passive contribution at 2.3 um SL

% Actin routine
for cntr = 1:NACT

    for trac = 1:NumSites

        switch trac

            case 1 % end; set to 0
                Mcore((cntr-1)*NumSites + 1, (cntr-1)*NumSites + 1) = 1;
        
            case NumSites
                Mcore(cntr*NumSites,cntr*NumSites - 1) = -KACT;
                Mcore(cntr*NumSites,cntr*NumSites) = KACT;
                EndVeccore(cntr*NumSites) = KACT*LACT;

            otherwise
                Mcore((cntr-1)*NumSites + trac,(cntr-1)*NumSites + trac - 1) = -KACT;
                Mcore((cntr-1)*NumSites + trac,(cntr-1)*NumSites + trac) = 2*KACT;
                Mcore((cntr-1)*NumSites + trac,(cntr-1)*NumSites + trac + 1) = -KACT;

        end

    end


end

% Myosin routine
for cntr = 1:NMYO

    for trac = 1:NumBridges

        switch trac

            case 1
                %Mcore(ACTNodes + (cntr-1)*NumBridges + 1, ACTNodes + (cntr-1)*NumBridges + 1) = KMYO;
                %Mcore(ACTNodes + (cntr-1)*NumBridges + 1, ACTNodes + (cntr-1)*NumBridges + 2) = -KMYO;
                %EndVeccore(ACTNodes + (cntr-1)*NumBridges + 1) = - KMYO * LMYO;
                Mcore(ACTNodes + (cntr-1)*NumBridges + 1, ACTNodes + (cntr-1)*NumBridges + 1) = KMYO + K_TITIN;
                Mcore(ACTNodes + (cntr-1)*NumBridges + 1, ACTNodes + (cntr-1)*NumBridges + 2) = -KMYO;
                EndVeccore(ACTNodes + (cntr-1)*NumBridges + 1) = (- KMYO * LMYO) + (K_TITIN*L_TITIN);

            case NumBridges
                Mcore(ACTNodes + cntr*NumBridges,ACTNodes + cntr*NumBridges) = 1;
                EndVeccore(ACTNodes + cntr*NumBridges) = START_LENGTH;

            otherwise
                Mcore(ACTNodes + (cntr-1)*NumBridges + trac,ACTNodes + (cntr-1)*NumBridges + trac - 1) = -KMYO;
                Mcore(ACTNodes + (cntr-1)*NumBridges + trac,ACTNodes + (cntr-1)*NumBridges + trac) = 2*KMYO;
                Mcore(ACTNodes + (cntr-1)*NumBridges + trac,ACTNodes + (cntr-1)*NumBridges + trac + 1) = -KMYO;

        end

    end


end


