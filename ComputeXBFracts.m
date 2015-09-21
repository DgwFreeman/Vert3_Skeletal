function [FractXB1, FractXB2] = ComputeXBFracts(TFstate, N_Thick_Start, MYONodes, NMYO)

FractXB1=length(find(TFstate(:,2)==1))/(N_Thick_Start*(MYONodes-NMYO)); % Fract XB Low force bearing
FractXB2=length(find(TFstate(:,2)==2))/(N_Thick_Start*(MYONodes-NMYO)); % Fract XB High force bearing