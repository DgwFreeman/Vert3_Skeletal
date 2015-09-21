function Tn = MakeTn(NTn, NACT, NumSites)

Tn = zeros(NACT, NTn);

for i=1:NACT
    Tn(i, 1:NTn)=[(((i-1)*NumSites)+2):4:(i*NumSites), (((i-1)*NumSites)+3):4:(i*NumSites)];
end