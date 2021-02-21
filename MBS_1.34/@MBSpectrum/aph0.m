% aph0 - Autophase zero-order
%
%   sp = aph0(sp). Performs the default phase calculation and applies the
%   phase shift.

% Created: 8/7/2002 Patrick Bolan, with revisions by Curt Corum
% MBS - Minnesota Breast Spectroscopy package
function [sp, phs] = aph0(sp)

% In previous versions I had many different variations on aph0. Instead
% there is now one quick-and-dirty aph0. For more explicit phase
% corrections, keep the estimation and correction steps separate.
phs = calcZeroOrderPhaseTD(sp); 

% My original autophase (Bolan MRM 2003) takes the mean of these two 
%   approaches:
%phsTD = calcZeroOrderPhaseTD(sp); 
%phsFD = calcZeroOrderPhaseFDmaxmin(sp); 
%phs = mean([phsTD; phsFD], 1);

sp = phaseShift(sp, -phs); 




