function TFstate = DispatchCaRegCoop(Tm_Type, NACT, NTn, TFstate, Tn, dt, CoopPass, TFRatePass, Ca, iCa, TnKO)

%%%*****************************************************************
%This is the loop over the actin nodes to check and see their
%state, and assing them Ca binding or not.  Regulation loop.

switch Tm_Type
    case 0 %% 18 nm
        TFstate=CaRegCoop_0a_TnKO(NACT, NTn, TFstate, Tn, dt, CoopPass, TFRatePass, Ca, iCa, TnKO);
    case 1
        TFstate=CaRegCoop_1a_TnKO(NACT, NTn, TFstate, Tn, dt, CoopPass, TFRatePass, Ca, iCa, TnKO);
    case 2
        TFstate=CaRegCoop_2a_TnKO(NACT, NTn, TFstate, Tn, dt, CoopPass, TFRatePass, Ca, iCa, TnKO);
    case 3
        TFstate=CaRegCoop_3a_TnKO(NACT, NTn, TFstate, Tn, dt, CoopPass, TFRatePass, Ca, iCa, TnKO);
    case 4 %% 74 nm
        TFstate=CaRegCoop_3b_TnKO(NACT, NTn, TFstate, Tn, dt, CoopPass, TFRatePass, Ca, iCa, TnKO);
end
%%%*****************************************************************