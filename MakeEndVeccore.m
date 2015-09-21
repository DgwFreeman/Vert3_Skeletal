function [EndVeccore] = MakeEndVeccore(NACT, NumSites, NMYO, NumBridges, KACT, LACT, KMYO, LMYO, START_LENGTH, K_TITIN, L_TITIN, ACTNodes, TOTNodes)

EndVeccore = zeros(TOTNodes,1);

% Actin routine
for cntr = 1:NACT
    EndVeccore(cntr*NumSites)=KACT*LACT;
end

% Myosin routine
for cntr = 1:NMYO
    EndVeccore(ACTNodes + (cntr-1)*NumBridges + 1) = (- KMYO * LMYO) + (K_TITIN*L_TITIN);
    EndVeccore(ACTNodes + cntr*NumBridges) = START_LENGTH;
end