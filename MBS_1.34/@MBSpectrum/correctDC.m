% correctDC - Removes the DC component of the FID using the mean of the
% last 1/8 points. This fails if there is large, low-frequency signal at
% the end of the fid. 
% In this case (ie lots of residual water on resonnance) you could shift
%   it off resonance so the frequency is higher.

% Created: 8/7/2002 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function [sp, correction] = correctDC(sp)

noiserange = floor(.875 * sp.N):sp.N;

% Loop version
for idx = 1:sp.M
    correction(idx) = mean(sp.fid(noiserange,idx));
	sp.fid(:,idx) = sp.fid(:,idx) - correction(idx);
end

% This is the matrix version. For some reason, this seems to be
%	leaving a small DC offset.
%sp.fid = sp.fid - mean(sp.fid(noiserange));







