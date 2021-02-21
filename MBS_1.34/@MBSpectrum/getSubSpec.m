% ****************************************************************************** 
% getSubSpec - Returns a portion of the spectrum specified by limits (ppm)
%
% [spec, freq] = getSubSpec(sp, limits) - returns two double arrays (not
% spectrum objects!) representing the spectrum and frequency axis over the
% range specicfied by limits (ppm). Note the it does not interpolate; it
% takes the closest points to the limit values

% Created: 8/7/2002 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function [spec, freq] = getSubSpec(sp, limits)

limleft = max(limits(1:2));
limright = min(limits(1:2));

jdxleft = interp1(sp.freqPPM,1:sp.N,limleft,'nearest');
jdxright = interp1(sp.freqPPM,1:sp.N,limright,'nearest');

% bounds check
jdxleft = max(1,jdxleft);
jdxright = min(sp.N, jdxright);

spec = sp.spec(jdxleft:jdxright,:);
freq = sp.freqPPM(jdxleft:jdxright);