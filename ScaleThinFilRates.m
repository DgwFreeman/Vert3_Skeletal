%% BCWT July 15, 2007
%% Modify the thin filament rates, to keep the detailed balance as things
%% scale.

function [koff, kon, RuOff, RuOn, CaOff, CaOn] = ScaleThinFilRates( ScaleK1, ScaleK2, ScaleK3 )

%%%%% Thin Filament Rates
%% Thin filament transition State 0 to 1
%% Ca binding dirctly to TnC
% Ca binding, fast equilibrium, K=1e5 M^-1*sec^-1
Keq=1e5*ScaleK1;
koff=5*1e1; %sec-1
kon = Keq*koff;  % (M*sec)^-1,
% Second Step, the TnI onto CaTnC
% Transition State 1 to 2
KaTnI=1e1*ScaleK2;
RuOff=1e0;
RuOn=KaTnI*RuOff;
% The loop where Ca comes off and the Regulatory unit
% of the thin filament goes back to the CaFree state.
% Thin filament State transition 2 to 0
Kd=(1/(KaTnI*Keq))*ScaleK3;
% The Ca being torn off and back to ground state
CaOff=5;%2*1e1;
CaOn=CaOff/Kd; % The reverse, which essentially cannot happen