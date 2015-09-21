%% June 26, 2007 BCWT
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

function [shared] = Coop_2aTo1(inState)

switch inState
  %% If the other Tn State is 2 (activated)
  %% then the shared one will remain activated.
  case 2
    shared=2;
    %% Otherwise, the shared node will go to 1, as the Tn in question here
    %% is either at 1 or 0.
  otherwise
    shared=1;
end