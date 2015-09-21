%% June 26, 2007 BCWT
%% function [shared] = Coop_2aTo0(inState);
%% This function will work on the nearest-neighbor thin filament
%% cooperativity.  Query the neighboring Tn to see what state it is in,
%% then return where the shared actin should transition (for 56 nm TmSpan).
%%
%% This function works on the cooperativity when a particular Tn is
%% switching to the 1 (mid state, Ca bound but unavailable to bind myosin)
%% thin filament state for the RU.
%%
%% Dependencies:
%% State of the neighboring Tn Tn in line

function [shared] = Coop_2aTo0(inState)

switch inState
  %% If the other Tn State is 0 (not Ca bound, inactivated)
  %% then the shared one will become inactivated.
  case 0
    shared=0;
    %% If the other Tn State is 1 ( Ca bound, inactivated)
    %% then the shared one will become 1 as well
  case 1
    shared=1;
    %% If the other Tn State is 2 (activated)
    %% then the shared one will remain activated.
  case 2
    shared=2;
end