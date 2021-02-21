% hsvdRemoveInside - Removes any peaks within the limits specified
%
%   [residual, peaks, spfit] = hsvdRemoveInside(sp, limitsPPM)
%       sp - input MBSpectrum
%       limitsPPM - array of 2 values specifing range to edit, in ppm
%       maxpts - max number of points to use when calculating SVD. Defaults
%       to using all, but the calculation scales with N^2
%       residual - an MBSpectrum of the residual
%       peaks - an array of lorentzian peak values. See simulateLorentzian
%           for structure
%       spfit - an NxM MBSpectrum of all M singular values

% Created: 8/6/2010 Patrick Bolan and Timo Liimatainen
% MBS - Minnesota Breast Spectroscopy package
function [residual, peaks, spfit] = hsvdRemoveInside(sp, limitsPPM, maxpts)

% Keep the original sp for final calculation, but do SVD on sp, which is
% truncated in the time domain. This works pretty well since the signal
% information is weighted to the beginning of the FID
sporig = sp;
if nargin>2
    sp = sporig.resizeTimeDomain(maxpts);
end

% First, normal HSVD
[~, peaks, ~] = hsvd(sp);

% Figure out which components are within the range
%limleft = max(limitsPPM(1:2) * sp.frq + sp.ppmref);
%limright = min(limitsPPM(1:2) * sp.frq + sp.ppmref);

% I had ppm-> freq conversion bug
limleft = max( (limitsPPM(1:2)-sp.ppmref) * sp.frq );
limright = min( (limitsPPM(1:2)-sp.ppmref) * sp.frq);

goodIdx = [];
for idx=1:size(peaks,2);
    if (peaks(idx).freq<=limleft) && (peaks(idx).freq>=limright)
        goodIdx = [goodIdx idx]; 
    end
end

peaksToRemove = peaks(goodIdx);

spfit = simulateLorentzian(peaksToRemove, sporig.N, sporig.sw, sporig.frq);
spfit.ppmref = sp.ppmref; % PJB 20190722
residual = sporig - spfit.sum;






