%% This function will take care of cooperativity for the TmSpan=56 nm
%% which is 1.5*37 nm.

function [TFstate]=CaRegCoop_2a_TnKO(NACT, NTn, TFstate, Tn, dt, CoopPass, TFRatePass, Ca, TnKO)

%%% Unpack coopeartativity parameters
% CooperativeType=CoopPass(13, 1);

%%% Unpack Thin Filament Rates
koff=TFRatePass(1,1);
kon=TFRatePass(2,1);
RuOff=TFRatePass(3,1);
RuOn=TFRatePass(4,1);
CaOff=TFRatePass(5,1);
CaOn=TFRatePass(6,1);

BaseTransitions = zeros(3,3);

BaseTransitions(1,2) = dt*kon*Ca; % 0->1
BaseTransitions(1,3) = dt*CaOn*Ca; % 0->2

BaseTransitions(2,1) = dt*koff; % 1->0
BaseTransitions(2,3) = dt*RuOn; % 1->2

BaseTransitions(3,1) = dt*CaOff; % 2->0
BaseTransitions(3,2) = dt*RuOff; % 2->1

%% Actual loop
% Look over each BS
% Look over each BS
% We will house the Ca binding in the first col. of the
% TFstate matrix.  The binding is
% 0==CaFree + TnC + TnI.A
% 1==CaTnC + TnI.A
% 2==CaTnC.TnI + A --> This is the state needed to make
% actin nodes available to bind.
for iRow=1:NACT  % Loop over each filament
    for iTn=1:NTn % Loop down the Tn index of each filament
        %% Now test whether this Tn can be modified.  TnKO will be zero if
        %% there are not Tn avail, and 1 if Tn is avail
        if TnKO(iRow, iTn) == 1
            ii=Tn(iRow, iTn); % This will dump a given Actin node into ii
            orientation = mod(iTn, (NTn/2));
            %% decide which nodes can cooperate
            switch orientation
                case 1 % if iTn==1 or 24
                    neighbors = [Tn(iRow, iTn+1) Tn(iRow, iTn+1)+2];
                    XB3nearby = (TFstate(ii, 2)==2)||(TFstate(ii+2, 2)==2);
                    affected = [ii+2, Tn(iRow, iTn + 1)];
                case 0  %if iTn=23 or 46
                    neighbors = [Tn(iRow, iTn-1) Tn(iRow, iTn-1)+2];
                    XB3nearby = (TFstate(ii, 2)==2)||(TFstate(ii-2, 2)==2);
                    affected = [ii-2, Tn(iRow, iTn - 1)];
                case 22 %if iTn=22 or 45
                    neighbors = [Tn(iRow, iTn-1) Tn(iRow, iTn-1)+2 Tn(iRow, iTn+1)];
                    XB3nearby = (TFstate(ii-2, 2)==2)||(TFstate(ii, 2)==2)||(TFstate(ii+2, 2)==2);
                    affected = [ii+2, Tn(iRow, iTn+1); ii-2, Tn(iRow, iTn - 1)];
                otherwise % If we are not at the either end, we will have Tm
                    neighbors = [Tn(iRow, iTn-1) Tn(iRow, iTn-1)+2 Tn(iRow, iTn+1) Tn(iRow, iTn+1)+2];
                    XB3nearby = (TFstate(ii-2, 2)==2)||(TFstate(ii, 2)==2)||(TFstate(ii+2, 2)==2);
                    affected = [ii+2, Tn(iRow, iTn + 1); ii-2, Tn(iRow, iTn - 1)];
            end
            %% do the sample to get new state
            transitions = CaRegTransitionProbabilities(BaseTransitions, CoopPass, neighbors, TFstate);
            OldState = TFstate(ii,1);
            NewState = OldState;
            sample = rand;
            if ~((OldState == 2) && XB3nearby) % MOD June 28, 2007  Here we will modify the condition here.  Therefore, if a XB is bound within the RU (TmSpan) then we cannot have the RU turn off.
                if (sample < transitions(OldState+1,1))
                    NewState = 0;
                elseif (sample < transitions(OldState+1,2))
                    NewState = 1;
                else
                    NewState = 2;
                end
            end
            % here we alter the states of the Tn and its collateral damage,
            % so if there's no actual change just move on instead
            if (NewState == OldState)
                continue
            end
            % do the actual state update
            TFstate(ii,1) = NewState;
            % now we update nearby affected Tns - which is more complex
            for nearby = 1:size(affected, 1); % # of nearby to alter
                % if it's pulled to the 1 state, we can't go through with
                % it if the current state is 2, which Coop_2aTo1 handles.
                % If it's pulled to the zero state, we copy it.
                % if it's pulled to the 2 state we just put in 2.
                switch NewState
                    case 0
                        TFstate(affected(nearby, 1), 1) = TFstate(affected(nearby, 2), 1);
                    case 1
                        TFstate(affected(nearby, 1), 1) = Coop_2aTo1(TFstate(affected(nearby, 2), 1));
                    case 2
                        TFstate(affected(nearby, 1), 1) = NewState;
                end
                % if we transitioned out of the 2 state, we have to change
                % the other TFstate fields.
                % The old state we check is ours, not the nearby. I don't
                % know why this is the case.
                if (OldState == 2 && TFstate(affected(nearby, 1), 1) ~= 2)
                    TFstate(affected(nearby, 1), 2) = 0;
                    TFstate(affected(nearby, 1), 3) = 0;
                end
            end
        end % End the If on whether we can go into this guy or not based on TnKO
    end % end loop over the Tn columns.  Where the columns of each row is given in Tn
    % by the Tn on 1 helix, then followed by the Tn on the
    % other helix.  Together this makes up all the Tn on one of
    % the 8 filaments. (Note: there are 46 col. in Tn matrix)
end % end loop over the Tn rows.  Each row of the Tn matrix satisfys one of the
% Thin filaments.  There are eight rows.

