function [Mcore] = MakeMcore(NACT, KACT, ACTNodes, NMYO, KMYO, NumSites, NumBridges, K_TITIN, TOTNodes)

Mcore = zeros(TOTNodes);
% Actin routine
for cntr = 1:NACT
    Mcore((cntr-1)*NumSites + 1, (cntr-1)*NumSites + 1) = 1;
    for trac = 2:NumSites-1
        Mcore((cntr-1)*NumSites + trac,(cntr-1)*NumSites + trac - 1) = -KACT;
        Mcore((cntr-1)*NumSites + trac,(cntr-1)*NumSites + trac) = 2*KACT;
        Mcore((cntr-1)*NumSites + trac,(cntr-1)*NumSites + trac + 1) = -KACT;
    end
    Mcore(cntr*NumSites,cntr*NumSites - 1) = -KACT;
    Mcore(cntr*NumSites,cntr*NumSites) = KACT;
end

% Myosin routine
for cntr = 1:NMYO
    Mcore(ACTNodes + (cntr-1)*NumBridges + 1, ACTNodes + (cntr-1)*NumBridges + 1) = KMYO + K_TITIN;
    Mcore(ACTNodes + (cntr-1)*NumBridges + 1, ACTNodes + (cntr-1)*NumBridges + 2) = -KMYO;
    for trac = 2:NumBridges-1
        Mcore(ACTNodes + (cntr-1)*NumBridges + trac,ACTNodes + (cntr-1)*NumBridges + trac - 1) = -KMYO;
        Mcore(ACTNodes + (cntr-1)*NumBridges + trac,ACTNodes + (cntr-1)*NumBridges + trac) = 2*KMYO;
        Mcore(ACTNodes + (cntr-1)*NumBridges + trac,ACTNodes + (cntr-1)*NumBridges + trac + 1) = -KMYO;
    end
    Mcore(ACTNodes + cntr*NumBridges,ACTNodes + cntr*NumBridges) = 1;
end

Mcore = sparse(Mcore);

% %% Gonna see if I can do a proper sparse initialization.
% % The actin routine assigns NACT*(1+(NumSites-2)*3+2) indices.
% % The myosin routine assigns NMYO*(2+(NumSites-2)*3+1).
% 
% indices = NACT*(1+(NumSites-2)*3+2) + NMYO*(2+(NumSites-2)*3+1);
% 
% i = zeros(1,indices);
% j = zeros(1,indices);
% s = zeros(1,indices);
