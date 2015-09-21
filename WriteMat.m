function [] = WriteMat(OutDir, pCa, Binder, Steps, Means, Vars, IndexThalf)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

MFvec = Steps(1,:);
AFvec = Steps(2,:);
FractXB1 = Steps(3,:);
FractXB2 = Steps(4,:);
FractCa0 = Steps(5,:);
FractCa1 = Steps(6,:);
FractCa2 = Steps(7,:);
ATPuse = Steps(8,:);

MFMean = Means(1,:);
AFMean = Means(2,:);
FractXB1Mean = Means(3,:);
FractXB2Mean = Means(4,:);
FractCa0Mean = Means(5,:);
FractCa1Mean = Means(6,:);
FractCa2Mean = Means(7,:);
ATPuseMean = Means(8,:);
FractBoundMean = Means(9,:);

MFVar = Vars(1,:);
AFVar = Vars(2,:);
FractXB1Var = Vars(3,:);
FractXB2Var = Vars(4,:);
FractCa0Var = Vars(5,:);
FractCa1Var = Vars(6,:);
FractCa2Var = Vars(7,:);
ATPuseVar = Vars(8,:);
FractBoundVar = Vars(9,:);

IndexThalfMFvec = IndexThalf(1,:);
IndexThalfAFvec = IndexThalf(2,:);
IndexThalfFractXB1 = IndexThalf(3,:);
IndexThalfFractXB2 = IndexThalf(4,:);
IndexThalfFractCa0 = IndexThalf(5,:);
IndexThalfFractCa1 = IndexThalf(6,:);
IndexThalfFractCa2 = IndexThalf(7,:);
IndexThalfATPuse = IndexThalf(8,:);
IndexThalfFractBound = IndexThalf(9,:);

save([OutDir filesep 'SimulationData_pCa_' num2str(pCa, '%3.2f') '.mat'], 'Binder', 'MF*', 'AF*', 'Fract*', 'ATP*', 'Index*');

end

