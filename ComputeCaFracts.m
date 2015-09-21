function [FractCa0, FractCa1, FractCa2] = ComputeCaFracts(TFstate, ACTNodes, NACT)

FractCa0=length(find(TFstate(:,1)==0))/(ACTNodes-NACT);
FractCa1=length(find(TFstate(:,1)==1))/(ACTNodes-NACT);
FractCa2=length(find(TFstate(:,1)==2))/(ACTNodes-NACT); % Fract Avail BS.