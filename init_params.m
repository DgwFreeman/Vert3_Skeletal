RequestedRuns = 4;
StartLength = 1150;

pCa = 4;

Psi = 100; % not actually passed!

DataParams.DataSpan = 0.1;
DataParams.dt = 10^-3;

knockout.TnKOType = 0;
knockout.XBKOType = 0;
knockout.TnFraction = 1;
knockout.XB_Fraction = 1;

StiffScale.kxscaler = 3;
StiffScale.kascaler = 1;
StiffScale.kmscaler = 1;

filaments.Tm_Type = 2;
filaments.LACT = 37.3/3;
filaments.LMYO = 42.9/3;
filaments.L_TITIN = 247;
filaments.NACT = 8;
filaments.NMYO = 4;
filaments.NumBridges = 60 + 1; % 60 bridges per filament + handel
filaments.NumSites = 90 + 1; % 90 per filament + handel
filaments.N_Thick_Start = 3;
filaments.NTn = 46;
filaments.KACT_0 = 1743;
filaments.KMYO_0 = 2020;
filaments.K_TITIN = 0.1344;
filaments.Angle_Crowns = 40;
filaments.Angle_Thick_Start = 120;

tcparam.SXB1=1;
tcparam.SXB2=1;
tcparam.SXB3=1;
tcparam.eta = 0.685; % The mechanical efficiency fraction for the total amount of work performed.
tcparam.f1 = 0.28; % The fraction drop from the first energy to the second minimum
tcparam.y0 = 20; % Offset
tcparam.SlopeL=-100; %Slope of line left of 0
tcparam.SlopeR=20; % Slope of line right of 0
tcparam.xL = -1; % xposition left of 0;
tcparam.xR = 1; % xposition right of 0;
tcparam.ATP=5e-3; 		%  Intracellular ATP in Molar
tcparam.ADP=30e-6;		%  Intracellular ADP in Molar
tcparam.phos=3e-3;		%  Intracellular Pi in Molar

coop.TF3TF12 = Psi^(1/3);
coop.XB2TF12 = Psi^(1/3);
coop.XB3TF12 = Psi^(1/3);
coop.XB2TF23 = Psi^(1/3);
coop.XB3TF23 = Psi^(1/3);
coop.TF3TF23 = Psi^(1/3);
coop.XB2TF21 = 1;
coop.XB3TF21 = 1;
coop.TF3TF21 = 1;
coop.XB2TF32 = 1;
coop.XB3TF32 = 1;
coop.TF3TF32 = 1;
coop.Coop_Type = 7;

TFRateScale.ScaleK1 = 1/100;
TFRateScale.ScaleK2 = 1;
TFRateScale.ScaleK3 = 1;